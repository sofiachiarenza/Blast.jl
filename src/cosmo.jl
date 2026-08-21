"""
    Background{T}

A unified container for cosmological background quantities.
Always targets the global Blast.χ grid to ensure consistent LOS integration.

# Fields:
- `cosmo::AbstractCosmology`: The underlying cosmology parameters.
- `z::Vector{T}`: Redshifts corresponding to the global χ grid.
- `χ::Vector{T}`: The global comoving distance grid (matches Blast.χ).
- `H::Vector{T}`: Hubble parameter H(z).
- `D::Vector{T}`: Raw ACE growth factor D(z).
- `D_norm::Vector{T}`: Growth factor normalized by the exact ACE value D(0).
- `f::Vector{T}`: Growth rate f(z).
- `z_of_χ`: Pre-built Akima interpolator closure for z(χ).

# Notes
The χ grid is assumed to be **uniformly spaced**. Kernel computations
(`compute_kernel!` for `CosmicShear` and `MagnificationBias`) rely on this
assumption when constructing the Simpson weight matrices.
"""
struct Background{T, F, C<:AbstractCosmology}
    cosmo::C
    z::Vector{T}
    χ::Vector{T}
    H::Vector{T}
    D::Vector{T}
    D_norm::Vector{T}
    f::Vector{T}
    z_of_χ::F
end

"""
    _make_z_of_χ(fine_z, fine_χ)

Builds the `z(χ)` akima closure inside a function barrier so that the
captured `fine_z` / `fine_χ` land as concrete field types on the closure
struct instead of the abstract `<:AbstractVector{T}` shapes that inference
leaves behind when the closure is constructed directly in `Background`.
"""
_make_z_of_χ(fine_z::AbstractVector{Tz}, fine_χ::AbstractVector{Tχ}) where {Tz, Tχ} =
    χ -> akima_interpolation(fine_z, fine_χ, χ)

"""
    Background(cosmo::AbstractCosmology; χ_grid=Blast.χ)
    Background{T}(cosmo::AbstractCosmology; χ_grid=Blast.χ) where {T}

Construct a Background snapshot by finding the redshifts corresponding to a
target χ grid.

The **typed form** `Background{T}(cosmo)` pins the element type `T` at the
call site and produces a fully-inferred `Background{T, F, C}`. ACE 0.11
parameterizes its cosmology fields concretely, but the typed form remains useful
for explicit control in inference-sensitive call sites.

The **untyped form** `Background(cosmo)` auto-detects `T` at runtime from
the computed arrays. This is the normal public path and is required when a
cosmology carries ForwardDiff `Dual` parameters. With ACE 0.11's concretely
parameterized cosmology, ordinary Float64 calls remain fully inferred.

The complete power-spectrum pipeline currently accepts only the default
production `Blast.χ` grid. A custom `χ_grid` may be used for isolated
background and kernel calculations, but `prepare_pk_workspace` and `get_Cℓ`
reject it because their plans, integration weights, and artifacts are fixed to
the production grid.
"""
function Background(cosmo::AbstractCosmology; χ_grid=Blast.χ)
    issorted(χ_grid) || throw(ArgumentError(
        "χ_grid must be monotonically increasing (got minimum at index $(argmin(χ_grid)))"))
    # Runtime T detection: evaluate the cosmology-dependent broadcasts, then
    # promote their eltypes against χ_grid. Only used by the AD path where
    # `cosmo.h` is a `Dual` and Float64-pinning would `convert`-fail.
    fine_z = collect(LinRange(0.0, Float64(Z_MAX_BACKGROUND), N_BG_FINE_GRID))
    fine_χ_raw = r_z.(fine_z, Ref(cosmo))
    H_array_raw = H_0_CONV .* cosmo.h .* E_z.(fine_z, Ref(cosmo))
    T = promote_type(eltype(χ_grid), eltype(fine_χ_raw), eltype(H_array_raw))
    return _build_background(T, cosmo, collect(χ_grid))
end

function Background{T}(cosmo::C; χ_grid=Blast.χ) where {T, C<:AbstractCosmology}
    issorted(χ_grid) || throw(ArgumentError(
        "χ_grid must be monotonically increasing (got minimum at index $(argmin(χ_grid)))"))
    return _build_background(T, cosmo, collect(χ_grid))
end

"""
    _build_background(::Type{T}, cosmo::C, χ_grid_vec::Vector{V}) where {T, C, V}

Function barrier that actually constructs the Background, specialized on the
concrete scalar, cosmology, and grid types.
"""
function _build_background(::Type{T}, cosmo::C,
                           χ_grid_vec::Vector{V}) where {T, C<:AbstractCosmology, V}
    # `fine_z` is a pure Float64 sampling grid (no cosmology dependence), used
    # for two independent interpolations below.  `collect()` materializes the
    # LinRange as a Vector so Mooncake can build a proper tangent; the
    # @from_chainrules wrapper on `akima_interpolation` does not support the
    # RData{NamedTuple{…}} tangent shape a raw LinRange produces.
    fine_z = collect(LinRange(0.0, Float64(Z_MAX_BACKGROUND), N_BG_FINE_GRID))

    # (1) Invert z ↔ χ on the fine grid. Keep the explicit conversion so the
    # Background element type follows the caller-selected T exactly.
    fine_χ_raw = r_z.(fine_z, Ref(cosmo))
    fine_χ = convert(Vector{T}, fine_χ_raw)
    # Build the closure inside a function barrier: local closures in Julia
    # capture their free vars with the types inference gave the *caller's*
    # locals. Here `fine_χ` is a `Vector{T}` so the closure's field types are
    # concrete, giving `Background{T, F, C}` a concrete `F`.
    z_of_χ_interp = _make_z_of_χ(fine_z, fine_χ)
    z_nodes_raw = z_of_χ_interp(χ_grid_vec)
    z_nodes = convert(Vector{T}, z_nodes_raw)

    # (2) H(z) has no ODE inside — safe to broadcast directly over Dual z_nodes.
    H_array_raw = H_0_CONV .* cosmo.h .* E_z.(z_nodes, Ref(cosmo))
    H_array = convert(Vector{T}, H_array_raw)

    # (3) D(z) and f(z) are computed via an ODE solve inside ACE's D_z / f_z.
    # That solver's `saveat` can't accept Dual-valued query points (saveat of
    # Vector{Dual} doesn't reconcile with the Float64 tspan/u₀ inside the
    # ODEProblem), so `D_z.(z_nodes::Vector{Dual}, Ref(cosmo_with_Duals))`
    # raises `MethodError(Float64, Dual)`.
    #
    # Workaround: evaluate D and f on the Float64 `fine_z` grid (where ACE's
    # ForwardDiff path is exercised by its own CI — z is Float64, cosmo params
    # are Dual), then interpolate at the Dual `z_nodes` via Blast's akima.
    # Akima propagates Duals through *both* the knot values and the query
    # points, so the derivative picks up the full chain rule
    #   dD_i/dp = (∂D/∂cosmo)(z)  +  (∂D/∂z)·(∂z_nodes/∂cosmo).
    # Same story for Mooncake: ACE's CI covers D_z(::Float64, ::cosmo) with
    # Mooncake, and Blast's akima has a Mooncake-registered rrule.
    #
    # Use the vectorized `D_f_z(::Vector, cosmo)` which does a SINGLE ODE
    # solve with `saveat=fine_z`, returning both D and f arrays.
    D_fine, f_fine = D_f_z(fine_z, cosmo)
    # akima_interpolation signature is (values, knots, queries); fine_z is the
    # ascending knot grid, D_fine / f_fine are the values to interpolate.
    D_array = convert(Vector{T}, akima_interpolation(D_fine, fine_z, z_nodes))
    # `fine_z` starts at exactly zero, so this is the ODE's exact D(0), not an
    # approximation at the first positive Blast χ node. Keep both conventions:
    # raw D is useful scientifically, while IA needs D(z)/D(0).
    D_norm_array = D_array ./ D_fine[1]
    f_array = convert(Vector{T}, akima_interpolation(f_fine, fine_z, z_nodes))

    return Background{T, typeof(z_of_χ_interp), typeof(cosmo)}(
        cosmo,
        z_nodes,
        convert(Vector{T}, χ_grid_vec),
        H_array,
        D_array,
        D_norm_array,
        f_array,
        z_of_χ_interp
    )
end


@doc raw"""
    compute_hubble_factor(z, cosmo)

Returns the dimensional Hubble parameter ``H(z) = H_0 E(z) = 100 h \\cdot E(z)``
in km/s/Mpc.

This is a thin wrapper around the `AbstractCosmologicalEmulators` function
`E_z(z, cosmo)`, which evaluates the dimensionless expansion history. In the
underlying `w0waCDMCosmology` implementation,

```math
E(z) = E(a=1/(1+z)),
```

with

```math
E(a) = \sqrt{\Omega_{\gamma,0} a^{-4} + \Omega_{cb,0} a^{-3} + \Omega_{k,0} a^{-2}
+ \Omega_{\Lambda,0} \, \rho_{\mathrm{DE}}(a) + \Omega_{\nu}E^2(a)},
```
"""
function compute_hubble_factor(z::Number, cosmo::AbstractCosmology)
    return H_0_CONV * cosmo.h * E_z(z, cosmo)
end

@doc raw"""
    compute_χ(z, cosmo)

Returns the comoving distance to redshift `z` in Mpc, computed via Gaussian
quadrature.

This is a thin wrapper around the `AbstractCosmologicalEmulators` function
`r_z(z, cosmo)`. In the underlying implementation,

```math
\chi(z) = \frac{c}{H_0} \int_0^z \frac{\mathrm{d}z'}{E(z')}.
```
"""
function compute_χ(z::Number, cosmo::AbstractCosmology; order=9)
    return r_z(z, cosmo, order=order)
end


# ACE 0.11 owns the cosmology and neutrino physics. It does not yet expose a
# public present-day neutrino-density accessor, so this adapter calls the exact
# implementation used internally by ACE's E(a). Keep this single dependency on
# the private function isolated here. Once ACE exports the accessor, only this
# function needs to change.
"""
    ω_ν0(cosmo::w0waCDMCosmology)

Return the present-day physical massive-neutrino density `ω_ν = Ω_ν h²` using
the exact Fermi-Dirac calculation in registered ACE 0.11.
"""
function ω_ν0(cosmo::w0waCDMCosmology)
    Ωγ0 = 2.469e-5 / cosmo.h^2
    return cosmo.h^2 * cosmo_ext._ΩνE2(one(cosmo.h), Ωγ0, cosmo.mν)
end

"""
    ω_m0(cosmo::w0waCDMCosmology)

Return the total present-day physical matter density
`ω_m = ω_b + ω_c + ω_ν`. Blast consumes total-matter `Pmm`, so this is the
density entering its Poisson-equation prefactors.
"""
ω_m0(cosmo::w0waCDMCosmology) = cosmo.ωb + cosmo.ωc + ω_ν0(cosmo)

ω_ν0(bg::Background) = ω_ν0(bg.cosmo)
ω_m0(bg::Background) = ω_m0(bg.cosmo)

# NLA convention uses dimensionless Ωm rather than the H0²Ωm combination used
# by lensing and ISW kernels. Keep this derived quantity private: Blast's public
# cosmological parameter interface is the ACE little-ω interface.
_Ω_m0(cosmo::w0waCDMCosmology) = ω_m0(cosmo) / cosmo.h^2
_Ω_m0(bg::Background) = _Ω_m0(bg.cosmo)

"""
    get_As(cosmo::AbstractCosmology)
Returns the primordial amplitude As.
"""
function get_As(cosmo::AbstractCosmology)
    return exp(cosmo.ln10Aₛ) / 1e10
end

"""
    get_ns(cosmo::AbstractCosmology)
Returns the spectral index ns.
"""
function get_ns(cosmo::AbstractCosmology)
    return cosmo.nₛ
end

# ACE 0.11 parameterizes every cosmology field concretely. The old Background
# conversion/type-assert accessors were workarounds for the pre-0.11 abstract
# `Number` fields and are no longer needed.
get_As(bg::Background{T}) where {T} = (convert(T, exp(bg.cosmo.ln10Aₛ) / 1e10))::T
get_ns(bg::Background{T}) where {T} = (convert(T, bg.cosmo.nₛ))::T

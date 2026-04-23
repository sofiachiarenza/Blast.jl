"""
    Background{T}

A unified container for cosmological background quantities.
Always targets the global Blast.χ grid to ensure consistent LOS integration.

# Fields:
- `cosmo::AbstractCosmology`: The underlying cosmology parameters.
- `z::Vector{T}`: Redshifts corresponding to the global χ grid.
- `χ::Vector{T}`: The global comoving distance grid (matches Blast.χ).
- `H::Vector{T}`: Hubble parameter H(z).
- `D::Vector{T}`: Growth factor D(z).
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
call site and produces a fully-inferred `Background{T, F, C}`. Use it when
downstream code (e.g. `f_full`) is sensitive to type stability — the
upstream `w0waCDMCosmology` declares its parameter fields as `::Number`,
so the untyped form below cannot be inferred to a concrete `T`.

The **untyped form** `Background(cosmo)` auto-detects `T` at runtime from
the computed arrays (needed for the ForwardDiff/Mooncake paths where
`cosmo` carries `Dual` parameters). It runs identically to the typed form
but its return type is `Background` UnionAll at the caller, which will
cascade runtime dispatch downstream.
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

Function barrier that actually constructs the Background, specialized on
the concrete `T`, `C`, and `V`. The explicit `C` / `V` type parameters are
load-bearing: Julia's default specialization heuristic skips specializing
on a slot whose declared type is an abstract supertype like
`::AbstractCosmology`, which would leave inference with
`typeof(cosmo)::w0waCDMCosmology` (the UnionAll — abstract!) instead of
`typeof(cosmo)::w0waCDMCosmology{Float64}`. Capturing `C` in the where
clause forces specialization on the concrete cosmology type.
"""
function _build_background(::Type{T}, cosmo::C,
                           χ_grid_vec::Vector{V}) where {T, C<:AbstractCosmology, V}
    # `fine_z` is a pure Float64 sampling grid (no cosmology dependence), used
    # for two independent interpolations below.  `collect()` materializes the
    # LinRange as a Vector so Mooncake can build a proper tangent; the
    # @from_chainrules wrapper on `akima_interpolation` does not support the
    # RData{NamedTuple{…}} tangent shape a raw LinRange produces.
    fine_z = collect(LinRange(0.0, Float64(Z_MAX_BACKGROUND), N_BG_FINE_GRID))

    # (1) Invert z ↔ χ on the fine grid. `convert(Vector{T}, …)` re-concretizes
    # what upstream `r_z(::Float64, ::w0waCDMCosmology)::Any` leaves abstract.
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
    f_array = convert(Vector{T}, akima_interpolation(f_fine, fine_z, z_nodes))

    return Background{T, typeof(z_of_χ_interp), typeof(cosmo)}(
        cosmo,
        z_nodes,
        convert(Vector{T}, χ_grid_vec),
        H_array,
        D_array,
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


"""
    get_Ωm(cosmo::AbstractCosmology)
Returns the total matter density parameter Ωm.
"""
function get_Ωm(cosmo::AbstractCosmology)
    return (cosmo.ωb + cosmo.ωc) / cosmo.h^2
end

"""
    get_Ωb(cosmo::AbstractCosmology)
Returns the baryon density parameter Ωb.
"""
function get_Ωb(cosmo::AbstractCosmology)
    return cosmo.ωb / cosmo.h^2
end

"""
    get_Ωc(cosmo::AbstractCosmology)
Returns the cold dark matter density parameter Ωc.
"""
function get_Ωc(cosmo::AbstractCosmology)
    return cosmo.ωc / cosmo.h^2
end


"""
    get_H0(cosmo::AbstractCosmology)
Returns the Hubble constant H0 in km/s/Mpc.
"""
function get_H0(cosmo::AbstractCosmology)
    return H_0_CONV * cosmo.h
end

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

# ─────────────────────────────────────────────────────────────────────────────
# Typed accessors dispatched on Background{T}
#
# `w0waCDMCosmology` (owned by the AbstractCosmologicalEmulators extension)
# is non-parametric with every parameter field declared `::Number`. Direct
# field access therefore returns `::Any`, and the `cosmo`-dispatched
# accessors above inherit that abstract return type.
#
# Inference-poisoning cascades from those accessors into every
# `compute_kernel!` that needs `H0` or `Ωm` (CosmicShear, CMBLensing,
# MagnificationBias, IntrinsicAlignment-NLA), which in turn makes
# `evaluate_components!(::WeakLensing)` and `f_full` inferred as `Any`.
#
# `Background{T}` already carries the concrete element type `T` on its
# eltype-bearing fields (`z`, `χ`, `H`, `D`, `f`) by construction. Routing
# cosmology-parameter reads through `Background{T}` therefore pins the
# return type at `T`, re-concretizing the downstream kernels.
#
# The `convert(T, …)` is runtime-free when the input is already ::T
# (identity), and correctly promotes when `cosmo.h` is a `Dual` but the
# Background's `T` is `Dual` as well (the usual ForwardDiff path).
# ─────────────────────────────────────────────────────────────────────────────
# NOTE on the `::T` typeassert: the trailing `::T` is load-bearing, not
# decorative. Without it, `convert(T, H_0_CONV * bg.cosmo.h)` comes back as
# `::Any` because `H_0_CONV * bg.cosmo.h` has already collapsed to `::Any`
# (Float64 * Number -> Any) before `convert` runs, and inference does not
# specialize the generic `convert(::Type{T}, ::Any)` back down to `T` when
# `T` is a method type-variable captured by the outer `where`. The
# post-convert typeassert forces the return type to `T` so every downstream
# consumer sees a concrete scalar.
get_H0(bg::Background{T}) where {T} = (convert(T, H_0_CONV * bg.cosmo.h))::T
get_Ωm(bg::Background{T}) where {T} = (convert(T, (bg.cosmo.ωb + bg.cosmo.ωc) / bg.cosmo.h^2))::T
get_Ωb(bg::Background{T}) where {T} = (convert(T, bg.cosmo.ωb / bg.cosmo.h^2))::T
get_Ωc(bg::Background{T}) where {T} = (convert(T, bg.cosmo.ωc / bg.cosmo.h^2))::T
get_As(bg::Background{T}) where {T} = (convert(T, exp(bg.cosmo.ln10Aₛ) / 1e10))::T
get_ns(bg::Background{T}) where {T} = (convert(T, bg.cosmo.nₛ))::T


"""
    w0waCDM(; w0, wa, H0, Ωm, Ωb, As, ns, σ8, Ωk=0.0, Ωr=0.0)

A flexible w0–wa CDM cosmological model.  The defaults `w0 = -1, wa = 0`
reduce the model to a flat ΛCDM cosmology, so this single constructor
covers both cases.

Maps standard parameters to AbstractCosmologicalEmulators format.
"""
function w0waCDM(; w0=-1.0, wa=0.0, H0=67.27, Ωm=0.3156, Ωb=0.0492, As=2.12107e-9, ns=0.9645, Ωk=0.0)
    h = H0 / H_0_CONV
    return w0waCDMCosmology(
        ωb = Ωb * h^2,
        ωc = (Ωm - Ωb) * h^2,
        ωk = Ωk,
        h = h,
        nₛ = ns,
        ln10Aₛ = log(1e10 * As),
        mν = 0.06,
        w0 = w0,
        wa = wa
    )
end

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
struct Background{T, F}
    cosmo::AbstractCosmology
    z::Vector{T}
    χ::Vector{T}
    H::Vector{T}
    D::Vector{T}
    f::Vector{T}
    z_of_χ::F
end

"""
    Background(cosmo::AbstractCosmology; χ_grid=Blast.χ)

Construct a Background snapshot by finding the redshifts corresponding to a target χ grid.
"""
function Background(cosmo::AbstractCosmology; χ_grid=Blast.χ)
    issorted(χ_grid) || throw(ArgumentError(
        "χ_grid must be monotonically increasing (got minimum at index $(argmin(χ_grid)))"))

    # `fine_z` is a pure Float64 sampling grid (no cosmology dependence), used
    # for two independent interpolations below.  `collect()` materializes the
    # LinRange as a Vector so Mooncake can build a proper tangent; the
    # @from_chainrules wrapper on `akima_interpolation` does not support the
    # RData{NamedTuple{…}} tangent shape a raw LinRange produces.
    fine_z = collect(LinRange(0.0, Float64(Z_MAX_BACKGROUND), N_BG_FINE_GRID))

    # (1) Invert z ↔ χ on the fine grid.
    fine_χ = r_z.(fine_z, Ref(cosmo))
    z_of_χ_interp(χ) = akima_interpolation(fine_z, fine_χ, χ)
    z_nodes = z_of_χ_interp(χ_grid)

    # (2) H(z) has no ODE inside — safe to broadcast directly over Dual z_nodes.
    H_array = H_0_CONV .* cosmo.h .* E_z.(z_nodes, Ref(cosmo))

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
    D_fine = D_z.(fine_z, Ref(cosmo))
    f_fine = f_z.(fine_z, Ref(cosmo))
    # akima_interpolation signature is (values, knots, queries); fine_z is the
    # ascending knot grid, D_fine / f_fine are the values to interpolate.
    D_array = akima_interpolation(D_fine, fine_z, z_nodes)
    f_array = akima_interpolation(f_fine, fine_z, z_nodes)

    # Promote the common element type across cosmology-dependent arrays and
    # the user-supplied χ grid. When `cosmo` carries Dual fields, H/D/f/z_nodes
    # are `Vector{Dual}` while χ_grid is typically `Vector{Float64}` — using a
    # single promoted `T` upcasts χ_grid so every struct field has matching
    # eltype.
    T = promote_type(eltype(χ_grid), eltype(H_array), eltype(z_nodes),
                     eltype(D_array), eltype(f_array))

    return Background{T, typeof(z_of_χ_interp)}(
        cosmo,
        convert(Vector{T}, z_nodes),
        convert(Vector{T}, collect(χ_grid)),
        convert(Vector{T}, H_array),
        convert(Vector{T}, D_array),
        convert(Vector{T}, f_array),
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

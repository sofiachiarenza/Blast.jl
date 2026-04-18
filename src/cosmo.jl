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
    T = eltype(χ_grid)
    # Dense sampling for accurate z(χ) inversion.
    # collect() materializes the LinRange as a Vector so Mooncake can build
    # a proper tangent for the captured closure field (z_of_χ_interp). The
    # @from_chainrules wrapper on _akima_interpolation does not support
    # LinRange's RData{NamedTuple{start,stop,len,lendiv}} tangent shape.
    fine_z = collect(LinRange(T(0.0), T(5.0), 1000))
    fine_χ = r_z.(fine_z, Ref(cosmo))
    
    z_of_χ_interp(χ) = Blast._akima_interpolation(fine_z, fine_χ, χ)

    z_nodes = z_of_χ_interp(χ_grid)
    
    H_array = T(100.0) .* cosmo.h .* E_z.(z_nodes, Ref(cosmo))
    D_array = D_z.(z_nodes, Ref(cosmo))
    f_array = f_z.(z_nodes, Ref(cosmo))

    return Background{T, typeof(z_of_χ_interp)}(
        cosmo, 
        z_nodes, collect(χ_grid), H_array, D_array, f_array, 
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
    return 100.0 * cosmo.h * E_z(z, cosmo)
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
    return 100.0 * cosmo.h
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
    ΛCDM(; H0, Ωm, Ωb, As, ns, σ8, Ωk=0.0, Ωr=0.0)

A flat ΛCDM cosmological model with fixed w0=-1 and wa=0.
Maps standard parameters to AbstractCosmologicalEmulators format.
"""
function ΛCDM(; H0=67.27, Ωm=0.3156, Ωb=0.0492, As=2.12107e-9, ns=0.9645, Ωk=0.0)
    h = H0 / 100.0
    return w0waCDMCosmology(
        ωb = Ωb * h^2,
        ωc = (Ωm - Ωb) * h^2,
        ωk = Ωk,
        h = h,
        nₛ = ns,
        ln10Aₛ = log(1e10 * As),
        mν = 0.06,
        w0 = -1.0,
        wa = 0.0
    )
end

"""
    w0waCDM(; w0, wa, H0, Ωm, Ωb, As, ns, σ8, Ωk=0.0, Ωr=0.0)

A flexible w0-wa CDM cosmological model.
Maps standard parameters to AbstractCosmologicalEmulators format.
"""
function w0waCDM(; w0=-1.0, wa=0.0, H0=67.27, Ωm=0.3156, Ωb=0.0492, As=2.12107e-9, ns=0.9645, Ωk=0.0)
    h = H0 / 100.0
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

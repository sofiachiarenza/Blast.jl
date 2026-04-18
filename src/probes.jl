"""
    AbstractCosmologicalProbes

Abstract supertype for cosmological probes in BLAST.
"""
abstract type AbstractCosmologicalProbes end

"""
    AbstractComponents

Abstract supertype for physical components entering cosmological probes.
"""
abstract type AbstractComponents end

"""
    prepare_nz_matrix(nz::AbstractMatrix, z::AbstractVector, z_grid::AbstractVector)

Interpolate and normalize an n(z) matrix bin-by-bin onto the calculation grid.
"""
function prepare_nz_matrix(nz::AbstractMatrix, z::AbstractVector, z_grid::AbstractVector)
    n_bins = size(nz, 1)
    nz_normed = zeros(eltype(nz), n_bins, length(z_grid))

    for b in 1:n_bins
        nz_norm_val, _ = quadgk(x -> Blast._akima_interpolation(nz[b,:], z, x), first(z_grid), last(z_grid))
        nz_norm_val > 0 || throw(ArgumentError(
            "n(z) for bin $b integrates to zero or negative ($nz_norm_val); " *
            "check that your redshift distribution is non-negative and overlaps z_grid."))
        nz_normed[b, :] = Blast._akima_interpolation(nz[b, :], z, z_grid) ./ nz_norm_val
    end
    return nz_normed
end

"""
    smooth_nz(nz::AbstractMatrix, z::AbstractVector; z_out::AbstractVector=z, span::Real=0.02)

Smooth n(z) bin-by-bin with `Loess.jl` and evaluate on `z_out`.
"""
function smooth_nz(
    nz::AbstractMatrix,
    z::AbstractVector,
    ;
    z_out::AbstractVector=z,
    span::Real=0.02
)
    size(nz, 2) == length(z) || throw(DimensionMismatch("size(nz, 2) must match length(z)."))
    0 < span <= 1 || throw(ArgumentError("span must satisfy 0 < span <= 1."))

    n_bins = size(nz, 1)
    nz_smooth = Matrix{Float64}(undef, n_bins, length(z_out))

    for b in 1:n_bins
        model = Loess.loess(z, nz[b, :], span=span)
        vs = collect(Loess.predict(model, z_out))
        vs[vs .< 0] .= 0.0
        nz_smooth[b, :] = vs
    end

    return nz_smooth
end

"""
    has_nz(::AbstractComponents) → Bool

Trait indicating whether a component carries an n(z) redshift distribution
that must be normalised before kernel computation. Defaults to `false`.
"""
has_nz(::AbstractComponents) = false
has_nz(::Nothing) = false

"""
    check_and_normalize!(Component, grid_z)

Internal helper: ensures nz_norm is populated for the current calculation grid.
Uses the `has_nz` trait to avoid `hasfield` runtime overhead.
"""
function check_and_normalize!(Component::AbstractComponents, z_grid::AbstractVector)
    if has_nz(Component)
        if size(Component.nz_norm) != (size(Component.nz, 1), length(z_grid))
            Component.nz_norm = prepare_nz_matrix(Component.nz, Component.z, z_grid)
        end
    end
end

function check_and_normalize!(Component::Nothing, z_grid::AbstractVector)
    return nothing
end

"""
    NLA_model(bg::Background; A=1.72, C1=0.0134)

Computes the Non-Linear Alignment (NLA) model for intrinsic alignments.
Returns an array evaluated on the Background grid.
"""
function NLA_model(bg::Background; A=1.72, C1=0.0134)
    Ωm = get_Ωm(bg.cosmo)
    return @. - A * C1 * Ωm / bg.D
end

@doc raw"""
    NumberCounts <: AbstractComponents

Galaxy number-counts density component.

This component stores the redshift distribution and linear bias entering

```math
K_i^{\delta}(z) = \frac{H(z)}{c}\, b_i(z)\, \hat n_i(z).
```
"""
@kwdef mutable struct NumberCounts <: AbstractComponents
    nz::Matrix{Float64} = zeros(1, 1)
    z::Vector{Float64} = zeros(1)
    nz_norm::Matrix{Float64} = zeros(1, 1)
    bias::Matrix{Float64} = zeros(1, 1)
    Kernel::Matrix{Float64} = zeros(1, 1)
    ell_prefactor::Vector{Float64} = ones(size(Blast.full_ℓ_range, 1))
    limber_factor::Vector{Float64} = ones(size(Blast.full_ℓ_range, 1))
end
has_nz(::NumberCounts) = true

@doc raw"""
    CosmicShear <: AbstractComponents

Weak-lensing shear component with kernel

```math
K_i^{\gamma}(\chi) = \frac{3 H_0^2 \Omega_m}{2 c^2}\, \chi (1+z)
\int_z^{\infty} \mathrm{d}z'\, \frac{H(z')}{c}\, \hat n_i(z')
\frac{\chi(z')-\chi}{\chi(z')}.
```
"""
@kwdef mutable struct CosmicShear <: AbstractComponents
    nz::Matrix{Float64} = zeros(1, 1)
    z::Vector{Float64} = zeros(1)
    nz_norm::Matrix{Float64} = zeros(1, 1)
    Kernel::Matrix{Float64} = zeros(1, 1)
    ell_prefactor::Vector{Float64} = @. sqrt(factorial_frac(Blast.full_ℓ_range))
    limber_factor::Vector{Float64} = (Blast.full_ℓ_range .+ 0.5) .^ (-2)
end
has_nz(::CosmicShear) = true

@doc raw"""
    CMBLensing <: AbstractComponents

CMB lensing convergence component with source plane at recombination:

```math
K^{\kappa_{\mathrm{CMB}}}(\chi) = \frac{3 H_0^2 \Omega_m}{2 c^2}\, \chi (1+z)
\left(1 - \frac{\chi}{\chi_*}\right).
```
"""
@kwdef mutable struct CMBLensing <: AbstractComponents
    Kernel::Matrix{Float64} = zeros(1, 1)
    ell_prefactor::Vector{Float64} = @. Blast.full_ℓ_range * (Blast.full_ℓ_range + 1)
    limber_factor::Vector{Float64} = (Blast.full_ℓ_range .+ 0.5) .^ (-2)
end

@doc raw"""
    RedshiftSpaceDistortions <: AbstractComponents

Redshift-space distortion component with kernel

```math
K_i^{\mathrm{RSD}}(z) = \frac{H(z)}{c}\, f(z)\, \hat n_i(z).
```
"""
@kwdef mutable struct RedshiftSpaceDistortions <: AbstractComponents
    nz::Matrix{Float64} = zeros(1, 1)
    z::Vector{Float64} = zeros(1)
    nz_norm::Matrix{Float64} = zeros(1, 1)
    growth_rate::Vector{Float64} = zeros(1)
    Kernel::Matrix{Float64} = zeros(1, 1)
    ell_prefactor::Vector{Float64} = ones(size(Blast.full_ℓ_range, 1))
    limber_factor::Vector{Float64} = ones(size(Blast.full_ℓ_range, 1))
end
has_nz(::RedshiftSpaceDistortions) = true

@doc raw"""
    MagnificationBias <: AbstractComponents

Magnification-bias component weighted by the source-slope factor $5s-2$.

Implemented as the shear-like lensing kernel times the extra slope factor:

```math
K_i^{\mu}(\chi) = \frac{3 H_0^2 \Omega_m}{2 c^2}\, \chi (1+z)
\int_z^{\infty} \mathrm{d}z'\, \frac{H(z')}{c}\, \hat n_i(z')
\frac{\chi(z')-\chi}{\chi(z')}\,\left[5 s_i(z')-2\right].
```
"""
@kwdef mutable struct MagnificationBias <: AbstractComponents
    nz::Matrix{Float64} = zeros(1, 1)
    z::Vector{Float64} = zeros(1)
    nz_norm::Matrix{Float64} = zeros(1, 1)
    s::Matrix{Float64} = zeros(1, 1)
    Kernel::Matrix{Float64} = zeros(1, 1)
    ell_prefactor::Vector{Float64} = @. Blast.full_ℓ_range * (Blast.full_ℓ_range + 1)
    limber_factor::Vector{Float64} = (Blast.full_ℓ_range .+ 0.5) .^ (-2)
end
has_nz(::MagnificationBias) = true

@doc raw"""
    IntrinsicAlignment <: AbstractComponents

Intrinsic-alignment component with NLA-style amplitude model:

Users can provide a generic amplitude table `A_IA` (same bin/redshift shape as
`nz_norm`) to be used directly in the kernel. If `A_IA` is not provided with the
required shape, BLAST falls back to the standard NLA model controlled by `A`.

```math
K_i^{\mathrm{IA}}(z) = \frac{H(z)}{c}\, A_{\mathrm{IA},i}(z)\, \hat n_i(z),
\qquad
A_{\mathrm{IA}}(z) = - A\, C_1\, \frac{\Omega_m}{D(z)}.
```
"""
@kwdef mutable struct IntrinsicAlignment <: AbstractComponents
    nz::Matrix{Float64} = zeros(1, 1)
    z::Vector{Float64} = zeros(1)
    nz_norm::Matrix{Float64} = zeros(1, 1)
    A::Float64 = 1.72  # Standard NLA amplitude
    A_IA::Matrix{Float64} = zeros(1, 1)
    Kernel::Matrix{Float64} = zeros(1, 1)
    ell_prefactor::Vector{Float64} = @. sqrt(factorial_frac(Blast.full_ℓ_range))
    limber_factor::Vector{Float64} = (Blast.full_ℓ_range .+ 0.5) .^ (-2)
end
has_nz(::IntrinsicAlignment) = true

@doc raw"""
    IntegratedSachsWolfe <: AbstractComponents

Integrated Sachs-Wolfe component sourced by the time variation of the
gravitational potential:

```math
K^{\mathrm{ISW}}(z) = \frac{3 T_{\mathrm{CMB}} H_0^2 \Omega_m}{c^3}\, H(z)\, \left[1-f(z)\right].
```
"""
@kwdef mutable struct IntegratedSachsWolfe <: AbstractComponents
    growth_rate::Vector{Float64} = zeros(1)
    Kernel::Matrix{Float64} = zeros(1, 1)
    ell_prefactor::Vector{Float64} = ones(size(Blast.full_ℓ_range, 1))
    limber_factor::Vector{Float64} = ones(size(Blast.full_ℓ_range, 1))
end

@doc raw"""
    PrimordialNonGaussianity <: AbstractComponents

Scale-dependent bias contribution induced by local-type primordial
non-Gaussianity:

```math
K_i^{\mathrm{PNG}}(z) = \frac{H(z)}{c}\, f_{\mathrm{NL}}\, b_{\Phi,i}(z)\, \hat n_i(z),
\qquad
b_{\Phi,i}(z) = 2\,\delta_c\,\left[b_i(z)-p\right].
```
"""
@kwdef mutable struct PrimordialNonGaussianity <: AbstractComponents
    nz::Matrix{Float64} = zeros(1, 1)
    z::Vector{Float64} = zeros(1)
    nz_norm::Matrix{Float64} = zeros(1, 1)
    bias::Matrix{Float64} = zeros(1, 1)
    f_NL::Float64 = 0.0
    p::Float64 = 0.0
    Kernel::Matrix{Float64} = zeros(1, 1)
    ell_prefactor::Vector{Float64} = ones(size(Blast.full_ℓ_range, 1))
    limber_factor::Vector{Float64} = ones(size(Blast.full_ℓ_range, 1))
end
has_nz(::PrimordialNonGaussianity) = true

"""
    GalaxyClustering <: AbstractCosmologicalProbes

Galaxy-clustering probe container.

It combines the density term `δ` with optional redshift-space distortions
`RSD`, magnification bias `μ`, and primordial non-Gaussianity `PNG`.
"""
@kwdef mutable struct GalaxyClustering <: AbstractCosmologicalProbes
    δ::NumberCounts
    RSD::Union{RedshiftSpaceDistortions, Nothing} = nothing
    μ::Union{MagnificationBias, Nothing} = nothing
    PNG::Union{PrimordialNonGaussianity, Nothing} = nothing
end

"""
    WeakLensing <: AbstractCosmologicalProbes

Weak-lensing probe container combining shear `γ` and optional intrinsic
alignment `IA`.
"""
@kwdef mutable struct WeakLensing <: AbstractCosmologicalProbes
    γ::CosmicShear
    IA::Union{IntrinsicAlignment, Nothing} = nothing
end

"""
    CMB <: AbstractCosmologicalProbes

CMB probe container combining lensing convergence `κ` and optional
Integrated Sachs-Wolfe contribution `ISW`.
"""
@kwdef mutable struct CMB <: AbstractCosmologicalProbes
    κ::CMBLensing
    ISW::Union{IntegratedSachsWolfe, Nothing} = nothing
end




@doc raw"""
    compute_kernel!(Component::NumberCounts, bg::Background)

Compute the galaxy number-counts kernel on the background redshift grid:

```math
K_i^{\delta}(z) = \frac{H(z)}{c}\, b_i(z)\, \hat n_i(z).
```
"""
function compute_kernel!(Component::NumberCounts, bg::Background) 
    check_and_normalize!(Component, bg.z)
    Component.Kernel = @. Component.bias * (bg.H' / C_LIGHT) * Component.nz_norm
end

@doc raw"""
    compute_kernel!(Component::CosmicShear, bg::Background)

Compute the weak-lensing shear kernel using Simpson integration on the $\chi$
grid:

```math
K_i^{\gamma}(\chi) = \frac{3 H_0^2 \Omega_m}{2 c^2}\, \chi (1+z)
\int_z^{\infty} \mathrm{d}z'\, \frac{H(z')}{c}\, \hat n_i(z')
\frac{\chi(z')-\chi}{\chi(z')}.
```
"""
function compute_kernel!(Component::CosmicShear, bg::Background)
    check_and_normalize!(Component, bg.z)
    n_bins = size(Component.nz_norm, 1)
    
    H0 = get_H0(bg.cosmo)
    Ωm = get_Ωm(bg.cosmo)
    prefac = 1.5 * H0^2 * Ωm / C_LIGHT^2
    
    simpson_matrix = simpson_weights_matrix(length(bg.z))
    Δχ = (bg.χ[end] - bg.χ[1]) / (length(bg.χ) - 1)

    @tullio kernel[idx_b, idx_zidx] := Δχ * (bg.H[idx_zp] / C_LIGHT) * simpson_matrix[idx_zidx, idx_zp] * prefac * bg.χ[idx_zidx] * (1.0 + bg.z[idx_zidx]) * Component.nz_norm[idx_b, idx_zp] * (bg.χ[idx_zp] - bg.χ[idx_zidx]) / (bg.χ[idx_zp])
    Component.Kernel = kernel
end

@doc raw"""
    compute_kernel!(Component::CMBLensing, bg::Background)

Compute the CMB lensing convergence kernel on the background grid:

```math
K^{\kappa_{\mathrm{CMB}}}(\chi) = \frac{3 H_0^2 \Omega_m}{2 c^2}\, \chi (1+z)
\left(1 - \frac{\chi}{\chi_*}\right),
\qquad \chi_* = \chi(z_*),\ z_* \simeq 1090.
```
"""
function compute_kernel!(Component::CMBLensing, bg::Background) 
    H0 = get_H0(bg.cosmo)
    Ωm = get_Ωm(bg.cosmo)
    prefac = 1.5 * H0^2 * Ωm / C_LIGHT^2
    
    χ_CMB = compute_χ(1090, bg.cosmo, order=120) 
    kernel = @. prefac * bg.χ * (1. + bg.z) * (1 - bg.χ/χ_CMB)
    Component.Kernel = reshape(kernel, 1, size(kernel,1))
end

@doc raw"""
    compute_kernel!(Component::RedshiftSpaceDistortions, bg::Background)

Compute the redshift-space-distortion kernel from the linear growth rate and
normalized redshift distribution:

```math
K_i^{\mathrm{RSD}}(z) = \frac{H(z)}{c}\, f(z)\, \hat n_i(z).
```
"""
function compute_kernel!(Component::RedshiftSpaceDistortions, bg::Background) 
    check_and_normalize!(Component, bg.z)
    Component.Kernel = @. bg.f' * (bg.H' / C_LIGHT) * Component.nz_norm
end

@doc raw"""
    compute_kernel!(Component::MagnificationBias, bg::Background)

Compute the magnification-bias kernel with source-slope factor $5s-2$:

```math
K_i^{\mu}(\chi) = \frac{3 H_0^2 \Omega_m}{2 c^2}\, \chi (1+z)
\int_z^{\infty} \mathrm{d}z'\, \frac{H(z')}{c}\, \hat n_i(z')
\frac{\chi(z')-\chi}{\chi(z')}\, \left[5 s_i(z') - 2\right].
```
"""
function compute_kernel!(Component::MagnificationBias, bg::Background)
    check_and_normalize!(Component, bg.z)
    n_bins = size(Component.nz_norm, 1)
    
    H0 = get_H0(bg.cosmo)
    Ωm = get_Ωm(bg.cosmo)
    prefac = 1.5 * H0^2 * Ωm / C_LIGHT^2
    
    simpson_matrix = simpson_weights_matrix(length(bg.z))
    Δχ = (bg.χ[end] - bg.χ[1]) / (length(bg.χ) - 1)

    @tullio kernel[idx_b, idx_zidx] := Δχ * (bg.H[idx_zp] / C_LIGHT) * simpson_matrix[idx_zidx, idx_zp] * prefac * bg.χ[idx_zidx] * (1.0 + bg.z[idx_zidx]) * Component.nz_norm[idx_b, idx_zp] * (bg.χ[idx_zp] - bg.χ[idx_zidx]) / (bg.χ[idx_zp]) * (5.0 * Component.s[idx_b, idx_zp] - 2)
    Component.Kernel = kernel
end

@doc raw"""
    compute_kernel!(Component::IntrinsicAlignment, bg::Background)

Compute the intrinsic-alignment kernel using the NLA amplitude model:

```math
K_i^{\mathrm{IA}}(z) = \frac{H(z)}{c}\, A_{\mathrm{IA},i}(z)\, \hat n_i(z),
\qquad
A_{\mathrm{IA}}(z) = - A\, C_1\, \frac{\Omega_m}{D(z)}.
```
"""
function compute_kernel!(Component::IntrinsicAlignment, bg::Background) 
    check_and_normalize!(Component, bg.z)
    
    # Use NLA model if A_IA is uninitialized
    if size(Component.A_IA) != (size(Component.nz_norm, 1), length(bg.z))
        n_bins = size(Component.nz_norm, 1)
        nla_vals = NLA_model(bg; A=Component.A)
        Component.A_IA = repeat(nla_vals', n_bins, 1)
    end
    
    Component.Kernel = @. Component.A_IA * (bg.H' / C_LIGHT) * Component.nz_norm
end

@doc raw"""
    compute_kernel!(Component::IntegratedSachsWolfe, bg::Background)

Compute the ISW kernel from the background expansion and growth history:

```math
K^{\mathrm{ISW}}(z) = \frac{3 T_{\mathrm{CMB}} H_0^2 \Omega_m}{c^3}\, H(z)\, \left[1-f(z)\right].
```
"""
function compute_kernel!(Component::IntegratedSachsWolfe, bg::Background) 
    H0 = get_H0(bg.cosmo)
    Ωm = get_Ωm(bg.cosmo)
    T_CMB = 2.726
    prefac = 3T_CMB * H0^2 * Ωm / C_LIGHT^3
    kernel = @. prefac * bg.H * (1 - bg.f)
    Component.Kernel = reshape(kernel, 1, size(kernel, 1))
end

@doc raw"""
    compute_kernel!(Component::PrimordialNonGaussianity, bg::Background)

Compute the primordial non-Gaussianity scale-dependent-bias kernel:

```math
K_i^{\mathrm{PNG}}(z) = \frac{H(z)}{c}\, f_{\mathrm{NL}}\, b_{\Phi,i}(z)\, \hat n_i(z),
\qquad
b_{\Phi,i}(z) = 2\,\delta_c\,\left[b_i(z)-p\right].
```
"""
function compute_kernel!(Component::PrimordialNonGaussianity, bg::Background) 
    check_and_normalize!(Component, bg.z)
    b_phi_vals = bΦ(Component.bias, Component.p)
    Component.Kernel = @. (bg.H' / C_LIGHT) * Component.f_NL * b_phi_vals * Component.nz_norm
end

function compute_kernel!(Component::Nothing, bg::Background) 
    return nothing
end

"""
    evaluate_components!(probe, bg)

Compute and store all kernels for the components of `probe` on the `Background` grid.
"""
function evaluate_components!(GC::GalaxyClustering, bg::Background)
    compute_kernel!(GC.δ, bg)
    compute_kernel!(GC.RSD, bg)
    compute_kernel!(GC.μ, bg)
    compute_kernel!(GC.PNG, bg)
end

function evaluate_components!(WL::WeakLensing, bg::Background)
    compute_kernel!(WL.γ, bg)
    compute_kernel!(WL.IA, bg)
end

function evaluate_components!(cmb::CMB, bg::Background)
    compute_kernel!(cmb.κ, bg)
    compute_kernel!(cmb.ISW, bg)
end



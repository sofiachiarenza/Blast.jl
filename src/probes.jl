"""
    AbstractCosmologicalProbes

Abstract supertype for cosmological probes in BLAST.

Concrete subtypes represent observables such as galaxy clustering, weak lensing,
or CMB. Each probe is composed of one or more physical components
(e.g. number counts, shear, redshift-space distortions...) that contribute to the total
angular power spectra.
"""
abstract type AbstractCosmologicalProbes end

"""
    AbstractComponents

Abstract supertype for physical components entering cosmological probes.

Concrete subtypes represent individual contributions to angular power spectra,
such as number counts, cosmic shear, magnification bias, or primordial
non-Gaussianity. Each component stores what is necessary in the kernel computation, and
ℓ-dependent prefactors required for angular power spectrum calculations.
"""
abstract type AbstractComponents end

"""
    NumberCounts <: AbstractComponents

Galaxy number-counts component.

This component describes the contribution of galaxy density fluctuations to
angular power spectra. It includes the galaxy redshift distribution, bias, and
the associated projection kernel.

Fields:
- `nz::Array{<:Any,2}`: Redshift distribution(s) of the tracer. Must have shape (nbins, z) where nbins is the number of tomographic bins.
- `z::Array{<:Any,1}`: Redshift grid corresponding to `nz`.
- `bias::Array{<:Any,2}`: Galaxy bias as a function of redshift. The bias must have size (nbins, z).
- `Kernel::Array{<:Any,2}`: Projection kernel for the number-count contribution.
- `ell_prefactor`: ℓ-dependent prefactor applied to the angular power spectrum.
- `limber_factor`: ℓ-dependent factor used in the Limber approximation.
"""
@kwdef mutable struct NumberCounts <: AbstractComponents
    nz::Array{<:Any, 2} = zeros(1, 1)
    z::Array{<:Any, 1} = zeros(1)
    bias::Array{<:Any, 2} = zeros(1, 1)
    Kernel::Array{<:Any, 2} = zeros(1, 1)
    ell_prefactor = ones(size(Blast.full_ℓ_range, 1))
    limber_factor = ones(size(Blast.full_ℓ_range, 1))
end

"""
    CosmicShear <: AbstractComponents

Cosmic shear component.

This component describes the contribution of gravitational lensing shear to
angular power spectra.

Fields:
- `nz::Array{<:Any,2}`: Source redshift distribution(s). Must have shape (nbins, z) where nbins is the number of tomographic bins.
- `z::Array{<:Any,1}`: Redshift grid corresponding to `nz`.
- `Kernel::Array{<:Any,2}`: Lensing kernel.
- `ell_prefactor`: ℓ-dependent prefactor applied to the angular power spectrum.
- `limber_factor`: ℓ-dependent factor used in the Limber approximation.
"""
@kwdef mutable struct CosmicShear <: AbstractComponents
    nz::Array{<:Any, 2} = zeros(1, 1)
    z::Array{<:Any, 1} = zeros(1)
    Kernel::Array{<:Any, 2} = zeros(1, 1)
    ell_prefactor = @. sqrt(factorial_frac(Blast.full_ℓ_range))
    limber_factor = (Blast.full_ℓ_range .+ 0.5) .^ (-2)
end

"""
    CMBLensing <: AbstractComponents

CMB lensing convergence component.

This component represents the lensing of the cosmic microwave background by
large-scale structure, contributing to angular power spectra through the
convergence field.

Fields:
- `Kernel::Array{<:Any,2}`: CMB lensing kernel.
- `ell_prefactor`: ℓ-dependent prefactor applied to the angular power spectrum.
- `limber_factor`: ℓ-dependent factor used in the Limber approximation.
"""
@kwdef mutable struct CMBLensing <: AbstractComponents
    Kernel::Array{<:Any, 2} = zeros(1, 1)
    ell_prefactor = @. Blast.full_ℓ_range * (Blast.full_ℓ_range + 1)
    limber_factor = (Blast.full_ℓ_range .+ 0.5) .^ (-2)
end

"""
    RedshiftSpaceDistortions <: AbstractComponents

Redshift-space distortion (RSD) component.

This component accounts for distortions in observed galaxy clustering caused
by peculiar velocities along the line of sight.

Fields:
- `nz::Array{<:Any,2}`: Redshift distribution(s) of the tracer.
- `z::Array{<:Any,1}`: Redshift grid.
- `growth_rate::Array{<:Any,1}`: Linear growth rate evaluated on the redshift grid.
- `Kernel::Array{<:Any,2}`: RSD projection kernel.
- `ell_prefactor`: ℓ-dependent prefactor applied to the angular power spectrum.
- `limber_factor`: ℓ-dependent factor used in the Limber approximation.
"""
@kwdef mutable struct RedshiftSpaceDistortions <: AbstractComponents
    nz::Array{<:Any, 2} = zeros(1,1)
    z::Array{<:Any, 1} = zeros(1)
    growth_rate::Array{<:Any, 1} = zeros(1)
    Kernel::Array{<:Any, 2} = zeros(1, 1)
    ell_prefactor = ones(size(Blast.full_ℓ_range, 1))
    limber_factor = ones(size(Blast.full_ℓ_range, 1))
end

"""
    MagnificationBias <: AbstractComponents

Magnification bias component.

This component describes the effect of gravitational lensing magnification on
observed galaxy number counts, parameterized by the magnification bias slope s.

Fields:
- `nz::Array{<:Any,2}`: Redshift distribution(s) of the tracer. Must have shape (nbins, z) where nbins is the number of tomographic bins.
- `z::Array{<:Any,1}`: Redshift grid.
- `s::Array{<:Any,2}`: Magnification bias slope as a function of redshift. Must have size (nbins, z).
- `Kernel::Array{<:Any,2}`: Magnification bias kernel.
- `ell_prefactor`: ℓ-dependent prefactor applied to the angular power spectrum.
- `limber_factor`: ℓ-dependent factor used in the Limber approximation.
"""
@kwdef mutable struct MagnificationBias <: AbstractComponents
    nz::Array{<:Any, 2} = zeros(1, 1)
    z::Array{<:Any, 1} = zeros(1)
    s::Array{<:Any, 2} = zeros(1,1)
    Kernel::Array{<:Any, 2} = zeros(1, 1)
    ell_prefactor = @. Blast.full_ℓ_range * (Blast.full_ℓ_range + 1)
    limber_factor = (Blast.full_ℓ_range .+ 0.5) .^ (-2)
end

"""
    IntrinsicAlignment <: AbstractComponents

Intrinsic alignment (IA) component.

This component models correlations in galaxy shapes arising from intrinsic
alignments with the large-scale tidal field, which contaminate weak lensing
measurements.

Fields:
- `nz::Array{<:Any,2}`: Source redshift distribution(s). Must have shape (nbins, z) where nbins is the number of tomographic bins.
- `z::Array{<:Any,1}`: Redshift grid.
- `A_IA::Array{<:Any,2}`: Intrinsic alignment amplitude, can be redshift dependent and different in each tomographic bin. Must have size (nbins, z).
- `Kernel::Array{<:Any,2}`: Intrinsic alignment kernel.
- `ell_prefactor`: ℓ-dependent prefactor applied to the angular power spectrum.
- `limber_factor`: ℓ-dependent factor used in the Limber approximation.
"""
@kwdef mutable struct IntrinsicAlignment <: AbstractComponents
    nz::Array{<:Any, 2} = zeros(1, 1)
    z::Array{<:Any, 1} = zeros(1)
    A_IA::Array{<:Any, 2} = zeros(1, 1)
    Kernel::Array{<:Any, 2} = zeros(1, 1)
    ell_prefactor = @. sqrt(factorial_frac(Blast.full_ℓ_range))
    limber_factor = (Blast.full_ℓ_range .+ 0.5) .^ (-2) #TODO: check that this is correct
end

"""
    IntegratedSachsWolfe <: AbstractComponents

Integrated Sachs–Wolfe (ISW) component.

This component describes temperature anisotropies in the CMB generated by the
time evolution of gravitational potentials, typically relevant on large scales.

Fields:
- `growth_rate::Array{<:Any,1}`: Growth rate of structure as a function of redshift.
- `Kernel::Array{<:Any,2}`: ISW projection kernel.
- `ell_prefactor`: ℓ-dependent prefactor applied to the angular power spectrum.
- `limber_factor`: ℓ-dependent factor used in the Limber approximation.
"""
@kwdef mutable struct IntegratedSachsWolfe <: AbstractComponents
    growth_rate::Array{<:Any, 1} = zeros(1)
    Kernel::Array{<:Any, 2} = zeros(1, 1)
    ell_prefactor = ones(size(Blast.full_ℓ_range, 1))
    limber_factor = ones(size(Blast.full_ℓ_range, 1)) #TODO: check that this is correct
end

"""
    PrimordialNonGaussianity <: AbstractComponents

Primordial non-Gaussianity (PNG) component.

This component models scale-dependent bias induced by primordial non-Gaussianity,
parameterized by the `f_NL` parameter.

Fields:
- `nz::Array{<:Any,2}`: Redshift distribution(s) of the tracer.
- `z::Array{<:Any,1}`: Redshift grid.
- `bias::Array{<:Any,2}`: Galaxy bias as a function of redshift. The bias must have size (nbins, z).
- `f_NL::Number`: Amplitude of primordial non-Gaussianity.
- `p::Number`: Universality relation parameter.
- `Kernel::Array{<:Any,2}`: PNG projection kernel.
- `ell_prefactor`: ℓ-dependent prefactor applied to the angular power spectrum.
- `limber_factor`: ℓ-dependent factor used in the Limber approximation.
"""
@kwdef mutable struct PrimordialNonGaussianity <: AbstractComponents
    nz::Array{<:Any, 2} = zeros(1,1)
    z::Array{<:Any, 1} = zeros(1)
    bias::Array{<:Any, 2} = zeros(1, 1)
    f_NL::Number = 0
    p::Number = 0
    Kernel::Array{<:Any, 2} = zeros(1, 1)
    ell_prefactor = ones(size(Blast.full_ℓ_range, 1))
    limber_factor = ones(size(Blast.full_ℓ_range, 1)) #TODO: check that this is correct
end

"""
    GalaxyClustering <: AbstractCosmologicalProbes

Galaxy clustering probe.

This probe combines multiple components contributing to observed galaxy clustering,
including density fluctuations, redshift-space distortions, magnification bias,
and primordial non-Gaussianity.

Fields:
- `δ::NumberCounts`: Number counts component.
- `RSD::Union{RedshiftSpaceDistortions,Nothing}`: Redshift-space distortion component.
- `μ::Union{MagnificationBias,Nothing}`: Magnification bias component.
- `PNG::Union{PrimordialNonGaussianity,Nothing}`: Primordial non-Gaussianity component.
"""
@kwdef mutable struct GalaxyClustering <: AbstractCosmologicalProbes
    δ::NumberCounts
    RSD::Union{RedshiftSpaceDistortions, Nothing} = nothing
    μ::Union{MagnificationBias, Nothing} = nothing
    PNG::Union{PrimordialNonGaussianity, Nothing} = nothing
end

"""
    WeakLensing <: AbstractCosmologicalProbes

Weak lensing probe.

This probe describes cosmic shear observations, optionally including intrinsic
alignment contributions.

Fields:
- `γ::CosmicShear`: Cosmic shear component.
- `IA::Union{IntrinsicAlignment,Nothing}`: Intrinsic alignment component.
"""
@kwdef mutable struct WeakLensing <: AbstractCosmologicalProbes
    γ::CosmicShear
    IA::Union{IntrinsicAlignment, Nothing} = nothing
end

"""
    CMB <: AbstractCosmologicalProbes

Cosmic microwave background (CMB) probe.

This probe combines CMB lensing and, optionally, the integrated Sachs–Wolfe effect.

Fields:
- `κ::CMBLensing`: CMB lensing convergence component.
- `ISW::Union{IntegratedSachsWolfe,Nothing}`: Integrated Sachs–Wolfe component.
"""
@kwdef mutable struct CMB <: AbstractCosmologicalProbes
    κ::CMBLensing
    ISW::Union{IntegratedSachsWolfe, Nothing} = nothing
end

"""
    compute_kernel!(Component::NumberCounts, grid::CosmologicalGrid, bg::BackgroundQuantities, cosmo::AbstractCosmology)

Compute and store the projection kernel for galaxy number counts across multiple redshift bins.
"""
function compute_kernel!(Component::NumberCounts, grid::CosmologicalGrid, bg::BackgroundQuantities, cosmo::AbstractCosmology) 

    n_bins = size(Component.nz, 1)
    kernel = zeros(n_bins, length(grid.z_range))

    for b in 1:n_bins
        nz_func = DataInterpolations.AkimaInterpolation(Component.nz[b,:], Component.z, extrapolation = ExtrapolationType.Extension)
        nz_norm, _ = quadgk(x->nz_func(x), first(grid.z_range), last(grid.z_range))

        kernel[b,:] = @. Component.bias[b,:] * (bg.Hz_array / C_LIGHT) * (nz_func.(grid.z_range) / nz_norm)
    end
    Component.Kernel = kernel
end

function compute_kernel_safe!(Component::CosmicShear, grid::CosmologicalGrid, bg::BackgroundQuantities, cosmo::AbstractCosmology) 

    n_bins = size(Component.nz, 1)
    kernel = zeros(n_bins, length(grid.z_range))

    for b in 1:n_bins
        nz_func = DataInterpolations.AkimaInterpolation(Component.nz[b,:], Component.z, extrapolation=ExtrapolationType.Extension)
        nz_norm, _ = quadgk(x->nz_func(x), first(grid.z_range), last(grid.z_range))

        prefac = 1.5 * cosmo.H0^2 * cosmo.Ωm / C_LIGHT^2
        χ_interp = DataInterpolations.AkimaInterpolation(bg.χz_array, grid.z_range, extrapolation=ExtrapolationType.Extension)

        for z_idx in 1:length(grid.z_range)
            integrand(x) = nz_func(x) * (1. - bg.χz_array[z_idx]/χ_interp(x))
            z_low = grid.z_range[z_idx]
            z_top = grid.z_range[end]
            int, _ = quadgk(x -> integrand(x), z_low, z_top) 

            kernel[b, z_idx] = prefac * bg.χz_array[z_idx] * (1. + grid.z_range[z_idx]) * int / nz_norm
        end
    end
    Component.Kernel = kernel
end

function compute_kernel!(Component::CosmicShear, grid::CosmologicalGrid, bg::BackgroundQuantities, cosmo::AbstractCosmology)

    n_bins = size(Component.nz, 1)
    kernel = zeros(n_bins, length(grid.z_range))

    n_z_array = zeros(n_bins, length(grid.z_range))

    χz_array = bg.χz_array
    z_range = grid.z_range

    prefac = 1.5 * cosmo.H0^2 * cosmo.Ωm / C_LIGHT^2
    simpson_matrix = simpson_weights_matrix(length(grid.z_range))
    Δχ = (χz_array[end] - χz_array[1]) / (length(χz_array) - 1)

    for b in 1:n_bins
        nz_func = DataInterpolations.AkimaInterpolation(Component.nz[b, :], Component.z, extrapolation=ExtrapolationType.Extension)
        nz_norm, _ = quadgk(x->nz_func(x), first(grid.z_range), last(grid.z_range))
        for (zidx, myz) in enumerate(grid.z_range)
            n_z_array[b, zidx] = nz_func(myz) / nz_norm
        end
    end

    @tullio kernel[idx_b, idx_zidx] := Δχ * (bg.Hz_array[idx_zp] / C_LIGHT) * simpson_matrix[idx_zidx, idx_zp] * prefac * χz_array[idx_zidx] * (1.0 + z_range[idx_zidx]) * n_z_array[idx_b, idx_zp] * (χz_array[idx_zp] - χz_array[idx_zidx]) / (χz_array[idx_zp] + 1e-18)

    Component.Kernel = kernel

end

function compute_kernel!(Component::CMBLensing, grid::CosmologicalGrid, bg::BackgroundQuantities, cosmo::AbstractCosmology) 

    prefac = 1.5 * cosmo.H0^2 * cosmo.Ωm / C_LIGHT^2
    χ_CMB = compute_χ(1090, cosmo) 

    kernel = @. prefac * bg.χz_array * (1. + grid.z_range) * (1 - bg.χz_array/χ_CMB)
    Component.Kernel = reshape(kernel, 1, size(kernel,1))
end

function compute_kernel!(Component::RedshiftSpaceDistortions, grid::CosmologicalGrid,  bg::BackgroundQuantities, cosmo::AbstractCosmology) 

    n_bins = size(Component.nz, 1)
    kernel = zeros(n_bins, length(grid.z_range))

    for b in 1:n_bins
        nz_func = DataInterpolations.AkimaInterpolation(Component.nz[b,:], Component.z, extrapolation=ExtrapolationType.Extension)
        nz_norm, _ = quadgk(x->nz_func(x), first(grid.z_range), last(grid.z_range))

        kernel[b,:] = @. Component.growth_rate * (bg.Hz_array / C_LIGHT) * (nz_func.(grid.z_range) / nz_norm) 
    end
    Component.Kernel = kernel
end

function compute_kernel_safe!(Component::MagnificationBias, grid::CosmologicalGrid, bg::BackgroundQuantities, cosmo::AbstractCosmology) 

    n_bins = size(Component.nz, 1)
    kernel = zeros(n_bins, length(grid.z_range))

    for b in 1:n_bins
        nz_func = DataInterpolations.AkimaInterpolation(Component.nz[b,:], Component.z, extrapolation=ExtrapolationType.Extension)
        nz_norm, _ = quadgk(x->nz_func(x), first(grid.z_range), last(grid.z_range))

        s_z = DataInterpolations.AkimaInterpolation(Component.s[b,:], grid.z_range, extrapolation=ExtrapolationType.Extension)

        prefac = 1.5 * cosmo.H0^2 * cosmo.Ωm / C_LIGHT^2
        χ_interp = DataInterpolations.AkimaInterpolation(bg.χz_array, grid.z_range, extrapolation=ExtrapolationType.Extension)

        for z_idx in 1:length(grid.z_range)
            integrand(x) = nz_func(x) * (1. - bg.χz_array[z_idx]/χ_interp(x)) * (5 .* s_z(x) .- 2)
            z_low = grid.z_range[z_idx]
            z_top =  grid.z_range[end]
            int, _ = quadgk(x -> integrand(x), z_low, z_top) 

            kernel[b, z_idx] = prefac * bg.χz_array[z_idx] * (1. + grid.z_range[z_idx]) * int / nz_norm
        end
    end
    Component.Kernel = kernel 
end

function compute_kernel!(Component::MagnificationBias, grid::CosmologicalGrid, bg::BackgroundQuantities, cosmo::AbstractCosmology)

    n_bins = size(Component.nz, 1)
    s_z_array = zeros(n_bins, length(grid.z_range))
    n_z_array = zeros(n_bins, length(grid.z_range))

    χz_array = bg.χz_array
    z_range = grid.z_range

    prefac = 1.5 * cosmo.H0^2 * cosmo.Ωm / C_LIGHT^2
    simpson_matrix = simpson_weights_matrix(length(grid.z_range))
    Δχ = (χz_array[end] - χz_array[1]) / (length(χz_array) - 1)

    for b in 1:n_bins
        nz_func = DataInterpolations.AkimaInterpolation(Component.nz[b, :], Component.z, extrapolation=ExtrapolationType.Extension)
        nz_norm, _ = quadgk(x->nz_func(x), first(grid.z_range), last(grid.z_range))
        
        s_z = DataInterpolations.AkimaInterpolation(Component.s[b, :], grid.z_range, extrapolation=ExtrapolationType.Extension)
        for (zidx, myz) in enumerate(grid.z_range)
            s_z_array[b, zidx] = s_z(myz)
            n_z_array[b, zidx] = nz_func(myz) / nz_norm
        end
    end

    @tullio kernel[idx_b, idx_zidx] := Δχ * (bg.Hz_array[idx_zp] / C_LIGHT) * simpson_matrix[idx_zidx, idx_zp] * prefac * χz_array[idx_zidx] * (1.0 + z_range[idx_zidx]) * n_z_array[idx_b, idx_zp] * (χz_array[idx_zp] - χz_array[idx_zidx]) / (χz_array[idx_zp] + 1e-18) * (5.0 * s_z_array[idx_b, idx_zp] - 2)

    Component.Kernel = kernel

end

function compute_kernel!(Component::IntrinsicAlignment, grid::CosmologicalGrid, bg::BackgroundQuantities, cosmo::AbstractCosmology) 

    n_bins = size(Component.nz, 1)
    kernel = zeros(n_bins, length(grid.z_range))

    for b in 1:n_bins
        nz_func = DataInterpolations.AkimaInterpolation(Component.nz[b,:], Component.z, extrapolation=ExtrapolationType.Extension)
        nz_norm, _ = quadgk(x->nz_func(x), first(grid.z_range), last(grid.z_range))

        kernel[b,:] = @. Component.A_IA[b,:] * (bg.Hz_array / C_LIGHT) * nz_func.(grid.z_range) / nz_norm 
    end
    Component.Kernel = kernel
end

function compute_kernel!(Component::IntegratedSachsWolfe, grid::CosmologicalGrid, bg::BackgroundQuantities, cosmo::AbstractCosmology) 
    T_CMB = 2.726
    prefac = 3T_CMB * cosmo.H0^2 * cosmo.Ωm / C_LIGHT^3

    kernel = @. prefac * bg.Hz_array * (1 - Component.growth_rate)
    Component.Kernel = reshape(kernel, 1, size(kernel, 1))
end

function compute_kernel!(Component::PrimordialNonGaussianity, grid::CosmologicalGrid, bg::BackgroundQuantities, cosmo::AbstractCosmology) 

    n_bins = size(Component.nz, 1)
    kernel = zeros(n_bins, length(grid.z_range))

    for b in 1:n_bins
        nz_func = DataInterpolations.AkimaInterpolation(Component.nz[b,:], Component.z, extrapolation = ExtrapolationType.Extension)
        nz_norm, _ = quadgk(x->nz_func(x), first(grid.z_range), last(grid.z_range))

        kernel[b,:] = (bg.Hz_array / C_LIGHT) .* Component.f_NL .* bΦ(Component.bias[b,:], Component.p) .* (nz_func.(grid.z_range) / nz_norm)
    end
    Component.Kernel = kernel
end

function compute_kernel!(Component::Nothing, grid::CosmologicalGrid, bg::BackgroundQuantities, cosmo::AbstractCosmology) 
    return nothing
end

function compute_kernel_safe!(Component::Nothing, grid::CosmologicalGrid, bg::BackgroundQuantities, cosmo::AbstractCosmology) 
    return nothing
end

function evaluate_components!(GC::GalaxyClustering, grid::CosmologicalGrid, bg::BackgroundQuantities, cosmo::AbstractCosmology) 
    compute_kernel!(GC.δ, grid, bg, cosmo)
    compute_kernel!(GC.RSD, grid, bg, cosmo)
    compute_kernel_safe!(GC.μ, grid, bg, cosmo)
    compute_kernel!(GC.PNG, grid, bg, cosmo)
end

function evaluate_components!(WL::WeakLensing, grid::CosmologicalGrid, bg::BackgroundQuantities, cosmo::AbstractCosmology)
    compute_kernel_safe!(WL.γ, grid, bg, cosmo)
    compute_kernel!(WL.IA, grid, bg, cosmo)
end

function evaluate_components!(cmb::CMB, grid::CosmologicalGrid, bg::BackgroundQuantities, cosmo::AbstractCosmology)
    compute_kernel!(cmb.κ, grid, bg, cosmo)
    compute_kernel!(cmb.ISW, grid, bg, cosmo)
end

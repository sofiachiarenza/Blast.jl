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
    compute_kernel!(component, grid, bg, cosmo)

Compute and store the projection kernel for a cosmological component.

This function mutates `component` by filling its `Kernel` field with the
redshift-dependent projection kernel evaluated on the redshift grid
`grid.z_range`.

Concrete implementations exist for each subtype of `AbstractComponents`,
corresponding to the physical equations defining that component.
"""
function compute_kernel! end

"""
    compute_kernel!(Component::NumberCounts, grid::CosmologicalGrid, bg::BackgroundQuantities, cosmo::AbstractCosmology)

Compute and store the projection kernel for galaxy number counts across multiple redshift bins.

For each bin `i`, the kernel is calculated as:
```math
W_i^{\\delta}(z) = \\frac{H(z)}{c}n_i(z)b_i(z)
```
where `n_i(z)` is the redshift distribution of the tracer, normalized to 1, and `b_i(z)` is the linear galaxy bias in the bin.
The resulting kernel has shape (n_bins, length(grid.z_range)).

# Parameters:
- `Component`: An instance of `NumberCounts`, in which the computed kernel values for each redshift bin will be stored.
- `grid`: A `CosmologicalGrid` object specifying the redshift range and grid points for kernel computation.
- `bg`: A `BackgroundQuantities` object containing arrays of Hubble parameter (`Hz_array`) and comoving distance (`χz_array`), precomputed over the grid.
- `cosmo`: An instance of `AbstractCosmology`, containing the cosmological model used to calculate the background quantities.

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

#TODO: Old safe version, figure out where to put this. I think I'll keep it for CI.
function compute_kernel_safe!(Component::CosmicShear, grid::CosmologicalGrid, bg::BackgroundQuantities, cosmo::AbstractCosmology) 

    n_bins = size(Component.nz, 1)
    kernel = zeros(n_bins, length(grid.z_range))

    for b in 1:n_bins
        nz_func = DataInterpolations.AkimaInterpolation(Component.nz[b,:], Component.z, extrapolation=ExtrapolationType.Extension)
        nz_norm, _ = quadgk(x->nz_func(x), first(grid.z_range), last(grid.z_range))

        prefac = 1.5 * cosmo.H0^2 * cosmo.Ωm / C_LIGHT^2

        for z_idx in 1:length(grid.z_range)
            integrand(x) = nz_func(x) * (1. - bg.χz_array[z_idx]/compute_χ(x, cosmo))
            z_low = grid.z_range[z_idx]
            z_top = grid.z_range[end]*1.1 #TODO: check max redshift, with n5k bins, lensing5 fallisce se uso valore diverso da 3.5
            int, err = quadgk(x -> integrand(x), z_low, z_top) #int is the lensing efficiency

            kernel[b, z_idx] = prefac * bg.χz_array[z_idx] * (1. + grid.z_range[z_idx]) * int / nz_norm
        end
    end
    Component.Kernel = kernel
end

"""
    compute_kernel!(Component::CosmicShear, grid::CosmologicalGrid, bg::BackgroundQuantities, cosmo::AbstractCosmology)

Computes the weak lensing shear kernel based on a redshift distribution `nz` and stores it in the `CosmicShear` struct. 
The kernel is defined as: 
```math
W_i^{\\gamma}(z) = \\frac{3}{2}\\frac{H_0^2}{c^2}\\Omega_m \\chi(z)(1+z)\\int_{z}^{\\infty}dz'n_i(z')(1-\\frac{\\chi(z)-}{\\chi(z')})
```
where `n_i(z)` is the redshift distribution of the tracer, normalized to 1.
The resulting kernel has shape (n_bins, length(grid.z_range)).

# Parameters:
- `Component`: An instance of `CosmicShear`, in which the computed kernel values for each redshift bin will be stored.
- `grid`: A `CosmologicalGrid` object specifying the redshift range and grid points for kernel computation.
- `bg`: A `BackgroundQuantities` object containing arrays of Hubble parameter (`Hz_array`) and comoving distance (`χz_array`), precomputed over the grid.
- `cosmo`: An instance of `AbstractCosmology`, containing the cosmological model used to calculate the background quantities.

"""
function compute_kernel!(Component::CosmicShear, grid::CosmologicalGrid, bg::BackgroundQuantities, cosmo::AbstractCosmology)

    n_bins = size(Component.nz, 1)
    kernel = zeros(n_bins, length(grid.z_range))

    n_z_array = zeros(n_bins, length(grid.z_range))

    χz_array = bg.χz_array
    z_range = grid.z_range

    prefac = 1.5 * cosmo.H0^2 * cosmo.Ωm / C_LIGHT^2
    simpson_matrix = simpson_weights_matrix(length(grid.z_range))
    Δχ = χz_array[2] - χz_array[1]

    for b in 1:n_bins
        nz_func = DataInterpolations.AkimaInterpolation(Component.nz[b, :], Component.z, extrapolation=ExtrapolationType.Extension)
        nz_norm = sum(nz_func(grid.z_range) .* bg.Hz_array / C_LIGHT) * Δχ
        for (zidx, myz) in enumerate(grid.z_range)
            n_z_array[b, zidx] = nz_func(myz) / nz_norm
        end
    end

    @tullio kernel[b, zidx] := Δχ * bg.Hz_array[zp] / C_LIGHT * simpson_matrix[zidx, zp] * prefac * χz_array[zidx] * (1.0 + z_range[zidx]) * n_z_array[b, zp] * (1.0 - χz_array[zidx] / χz_array[zp])

    Component.Kernel = kernel

end

"""
    compute_kernel!(Component::CMBLensing, grid::CosmologicalGrid, bg::BackgroundQuantities, cosmo::AbstractCosmology)

Computes the weak lensing shear kernel based on a redshift distribution `nz` and stores it in the `CosmicShear` struct. 
The kernel is defined as: 
```math
W_i^{\\kappa}(z) = \\frac{3}{2}\\frac{H_0^2}{c^2}\\Omega_m \\chi(z)(1+z)(1-\\frac{\\chi(z)-}{\\chi^*})
```
where `χ^*` is the comoving distance at the CMB redshift `χ(z=1090)`. The resulting kernel has shape (1, length(grid.z_range)).

# Parameters:
- `Component`: An instance of `CMBLensing`, in which the computed kernel values for each redshift bin will be stored.
- `grid`: A `CosmologicalGrid` object specifying the redshift range and grid points for kernel computation.
- `bg`: A `BackgroundQuantities` object containing arrays of Hubble parameter (`Hz_array`) and comoving distance (`χz_array`), precomputed over the grid.
- `cosmo`: An instance of `AbstractCosmology`, containing the cosmological model used to calculate the background quantities.

"""
function compute_kernel!(Component::CMBLensing, grid::CosmologicalGrid, bg::BackgroundQuantities, cosmo::AbstractCosmology) 

    prefac = 1.5 * cosmo.H0^2 * cosmo.Ωm / C_LIGHT^2
    χ_CMB = compute_χ(1090, cosmo) #From DESI fiducial cosmology

    kernel = @. prefac * bg.χz_array * (1. + grid.z_range) * (1 - bg.χz_array/χ_CMB)
    Component.Kernel = reshape(kernel, 1, size(kernel,1))
end

"""
    compute_kernel!(Component::RedshiftSpaceDistortions, grid::CosmologicalGrid, bg::BackgroundQuantities, cosmo::AbstractCosmology)

Compute and store the projection kernel for Redshift-Space Distortions (RSD) across multiple redshift bins.

This function implements the RSD contribution to the density contrast in the Kaiser limit. 
For each bin `i`, the kernel is calculated as:

```math
W_i^{\\mathrm{RSD}}(z) =  \\frac{H(z)}{c} n_i(z)f(z)
```
where `f(z)` is the linear growth rate and `n_i(z)` is the redshift distribution of the tracer, normalized to 1.

# Parameters:
- `Component`: An instance of `CMBLensing`, in which the computed kernel values for each redshift bin will be stored.
- `grid`: A `CosmologicalGrid` object specifying the redshift range and grid points for kernel computation.
- `bg`: A `BackgroundQuantities` object containing arrays of Hubble parameter (`Hz_array`) and comoving distance (`χz_array`), precomputed over the grid.
- `cosmo`: An instance of `AbstractCosmology`, containing the cosmological model used to calculate the background quantities.

"""
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

        for z_idx in 1:length(grid.z_range)
            integrand(x) = nz_func(x) * (1. - bg.χz_array[z_idx]/compute_χ(x, cosmo)) * (5 .* s_z(x) .- 2)
            z_low = grid.z_range[z_idx]
            z_top =  grid.z_range[end]*1.1 #TODO: this is horrible.
            int, _ = quadgk(x -> integrand(x), z_low, z_top) 

            kernel[b, z_idx] = prefac * bg.χz_array[z_idx] * (1. + grid.z_range[z_idx]) * int / nz_norm
        end
    end
    Component.Kernel = kernel 
end

"""
    compute_kernel!(Component::MagnificationBias, grid::CosmologicalGrid, bg::BackgroundQuantities, cosmo::AbstractCosmology)

Computes the magnification bias kernel based on a redshift distribution `nz` and stores it in the `MagnificationBias` struct. 
The kernel is defined as: 
```math
W_i^{\\mu}(z) = \\frac{3}{2}\\frac{H_0^2}{c^2}\\Omega_m \\chi(z)(1+z)\\int_{z}^{\\infty}dz'n_i(z')(1-\\frac{\\chi(z)-}{\\chi(z')})(5s(z')-2)
```
where `n_i(z)` is the redshift distribution of the tracer, normalized to 1. And `s_i(z)` is magnification bias slope in the bin.
The resulting kernel has shape (n_bins, length(grid.z_range)).

# Parameters:
- `Component`: An instance of `MagnificationBias`, in which the computed kernel values for each redshift bin will be stored.
- `grid`: A `CosmologicalGrid` object specifying the redshift range and grid points for kernel computation.
- `bg`: A `BackgroundQuantities` object containing arrays of Hubble parameter (`Hz_array`) and comoving distance (`χz_array`), precomputed over the grid.
- `cosmo`: An instance of `AbstractCosmology`, containing the cosmological model used to calculate the background quantities.

"""
function compute_kernel!(Component::MagnificationBias, grid::CosmologicalGrid, bg::BackgroundQuantities, cosmo::AbstractCosmology)

    n_bins = size(Component.nz, 1)
    kernel = zeros(n_bins, length(grid.z_range))

    s_z_array = zeros(n_bins, length(grid.z_range))
    n_z_array = zeros(n_bins, length(grid.z_range))

    χz_array = bg.χz_array
    z_range = grid.z_range

    prefac = 1.5 * cosmo.H0^2 * cosmo.Ωm / C_LIGHT^2
    simpson_matrix = simpson_weights_matrix(length(grid.z_range))
    Δχ = χz_array[2] - χz_array[1]

    for b in 1:n_bins
        nz_func = DataInterpolations.AkimaInterpolation(Component.nz[b, :], Component.z, extrapolation=ExtrapolationType.Extension)
        nz_norm = sum(nz_func(grid.z_range) .* bg.Hz_array / C_LIGHT) * Δχ
        #notice that, since the normalization is in z but we assume the regular grid in chi, we integrate in χ
        # and include the jacobian
        s_z = DataInterpolations.AkimaInterpolation(Component.s[b, :], grid.z_range, extrapolation=ExtrapolationType.Extension)
        for (zidx, myz) in enumerate(grid.z_range)
            s_z_array[b, zidx] = s_z(myz)
            n_z_array[b, zidx] = nz_func(myz) / nz_norm
        end
    end

    @tullio kernel[b, zidx] := Δχ * bg.Hz_array[zp] / C_LIGHT * simpson_matrix[zidx, zp] * prefac * χz_array[zidx] * (1.0 + z_range[zidx]) * n_z_array[b, zp] * (1.0 - χz_array[zidx] / χz_array[zp]) * (5.0 * s_z_array[zp] - 2)

    Component.Kernel = kernel

end

"""
    compute_kernel!(Component::IntrinsicAlignment, grid::CosmologicalGrid, bg::BackgroundQuantities, cosmo::AbstractCosmology)

Compute and store the projection kernel for Intrinsic Alignments (IA) across multiple redshift bins.

This function calculates the IA kernel, which describes how galaxy shapes are intrinsically 
aligned with the tidal field of the large-scale structure. For each bin `i`, the kernel is:
```math
W_i^{IA}(z) =  \\frac{H(z)}{c} n_i(z) A_{IA}(z)
```
where n_i(z) is the redshift distribution of the tracer, normalized to 1, and A_{IA}(z) is a general alignment amplitude, allowing the user to implement the preferred IA model.

The resulting kernel has shape (n_bins, length(grid.z_range)).

# Parameters:
- `Component`: An instance of `IntrinsicAlignment`, in which the computed kernel values for each redshift bin will be stored.
- `grid`: A `CosmologicalGrid` object specifying the redshift range and grid points for kernel computation.
- `bg`: A `BackgroundQuantities` object containing arrays of Hubble parameter (`Hz_array`) and comoving distance (`χz_array`), precomputed over the grid.
- `cosmo`: An instance of `AbstractCosmology`, containing the cosmological model used to calculate the background quantities. 
"""
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

"""
    compute_kernel!(Component::IntegratedSachsWolfe, grid::CosmologicalGrid, bg::BackgroundQuantities, cosmo::AbstractCosmology)

Computes the Integrated Sachs-Wolfe (ISW) effect kernel and stores it in the `IntegratedSachsWolfe` struct.

The ISW kernel describes the secondary anisotropy in the CMB temperature caused by time-varying gravitational potentials as photons traverse large-scale structures. It is defined as:

```math
W^{ISW}(z) = \\frac{3 T_{CMB} H_0^2 \\Omega_m}{c^3} H(z) [1 - f(z)]
```
where `T_{CMB} = 2.726K` is the CMB temperature and `f(z)` is the linear growth rate. The resulting kernel has shape (1, length(grid.z_range)).

# Parameters:
- `Component`: An instance of `IntegratedSachsWolfe`, in which the computed kernel values for each redshift bin will be stored.
- `grid`: A `CosmologicalGrid` object specifying the redshift range and grid points for kernel computation.
- `bg`: A `BackgroundQuantities` object containing arrays of Hubble parameter (`Hz_array`) and comoving distance (`χz_array`), precomputed over the grid.
- `cosmo`: An instance of `AbstractCosmology`, containing the cosmological model used to calculate the background quantities. 
"""
function compute_kernel!(Component::IntegratedSachsWolfe, grid::CosmologicalGrid, bg::BackgroundQuantities, cosmo::AbstractCosmology) 
    T_CMB = 2.726
    prefac = 3T_CMB * cosmo.H0^2 * cosmo.Ωm / C_LIGHT^3

    kernel = @. prefac * bg.Hz_array * (1 - Component.growth_rate)
    Component.Kernel = reshape(kernel, 1, size(kernel, 1))
end

"""
    compute_kernel!(Component::PrimordialNonGaussianity, grid::CosmologicalGrid, bg::BackgroundQuantities, cosmo::AbstractCosmology)

Compute and store the projection kernel for Primordial Non-Gaussianity (PNG) across multiple redshift bins.

This function calculates the kernel for the scale-dependent bias contribution arising from 
local-type primordial non-Gaussianity. For each bin `i`, the kernel is calculated as:

```math
W_i^{f_{NL}}(z) = f_{NL} b_{\\Phi}(z) \\frac{H(z)}{c} n_i(z)
```
where `f_{NL}` is the PNG parameter, `b_\\Phi` is the PNG bias, parametrized with the Universality relation `b_\\Phi = 2 \\delta_c (b_i(z) - p)`.
Finally, n_i(z) is the redshift distribution of the tracer, normalized to 1. The resulting kernel has shape (n_bins, length(grid.z_range)).
# Parameters:
- `Component`: An instance of `PrimordialNonGaussianity`, in which the computed kernel values for each redshift bin will be stored.
- `grid`: A `CosmologicalGrid` object specifying the redshift range and grid points for kernel computation.
- `bg`: A `BackgroundQuantities` object containing arrays of Hubble parameter (`Hz_array`) and comoving distance (`χz_array`), precomputed over the grid.
- `cosmo`: An instance of `AbstractCosmology`, containing the cosmological model used to calculate the background quantities. 
"""
function compute_kernel!(Component::PrimordialNonGaussianity, grid::CosmologicalGrid, bg::BackgroundQuantities, cosmo::AbstractCosmology) 

    n_bins = size(Component.nz, 1)
    kernel = zeros(n_bins, length(grid.z_range))

    for b in 1:n_bins
        nz_func = DataInterpolations.AkimaInterpolation(Component.nz[b,:], Component.z, extrapolation = ExtrapolationType.Extension)
        nz_norm, _ = quadgk(x->nz_func(x), first(grid.z_range), last(grid.z_range))

        kernel[b,:] = (bg.Hz_array / C_LIGHT) .* Component.f_NL .* bΦ(Component.bias[b,:], Component.p) .* (nz_func.(grid.z_range) ./ nz_norm)
    end
    Component.Kernel = kernel
end

function compute_kernel!(Component::Nothing, grid::CosmologicalGrid, bg::BackgroundQuantities, cosmo::AbstractCosmology) 
    return nothing
end

"""
    evaluate_components!(GC::GalaxyClustering, grid::CosmologicalGrid, bg::BackgroundQuantities, cosmo::AbstractCosmology)

Compute all individual projection kernels associated with galaxy clustering observable.

This function computes the kernels of the active components in the observable, namely:
- **Number counts** (`GC.δ`): The primary galaxy clustering signal.
- **Redshift-Space Distortions** (`GC.RSD`): The anisotropic clustering due to peculiar velocities.
- **Magnification Bias** (`GC.μ`): The change in observed density due to gravitational lensing.
- **Primordial Non-Gaussianity** (`GC.PNG`): The scale-dependent bias contribution.

# Parameters:
- `GC`: A `GalaxyClustering` object containing the sub-components to be updated.
- `grid`: The `CosmologicalGrid` defining the redshift range.
- `bg`: `BackgroundQuantities` containing precomputed distances and Hubble factors.
- `cosmo`: The `AbstractCosmology` model.
"""
function evaluate_components!(GC::GalaxyClustering, grid::CosmologicalGrid, bg::BackgroundQuantities, cosmo::AbstractCosmology) 
    compute_kernel!(GC.δ, grid, bg, cosmo)
    compute_kernel!(GC.RSD, grid, bg, cosmo)
    compute_kernel!(GC.μ, grid, bg, cosmo)
    compute_kernel!(GC.PNG, grid, bg, cosmo)
end

"""
    evaluate_components!(WL::WeakLensing, grid::CosmologicalGrid, bg::BackgroundQuantities, cosmo::AbstractCosmology)

Compute all projection kernels associated with Weak Lensing observable.

This function computes the kernels of the active components in the observable, namely:
- **Cosmic Shear** (`WL.γ`): The deflection of light by large-scale structure.
- **Intrinsic Alignment** (`WL.IA`): The intrinsic correlation of galaxy shapes with the tidal field.

# Parameters:
- `WL`: A `WeakLensing` object containing the lensing and alignment sub-components.
- `grid`: The `CosmologicalGrid` defining the redshift range.
- `bg`: `BackgroundQuantities` for distance and evolution calculations.
- `cosmo`: The `AbstractCosmology` model.
"""
function evaluate_components!(WL::WeakLensing, grid::CosmologicalGrid, bg::BackgroundQuantities, cosmo::AbstractCosmology)
    compute_kernel!(WL.γ, grid, bg, cosmo)
    compute_kernel!(WL.IA, grid, bg, cosmo)
end

"""
    evaluate_components!(cmb::CMB, grid::CosmologicalGrid, bg::BackgroundQuantities, cosmo::AbstractCosmology)

Compute all projection kernels associated with Cosmic Microwave Background (CMB) secondary anisotropies.

This function computes the kernels of the active components in the observable, namely:
- **CMB Lensing** (`cmb.κ`): The gravitational lensing of the CMB temperature and polarization.
- **Integrated Sachs-Wolfe** (`cmb.ISW`): The temperature fluctuations from time-varying potentials.

# Parameters:
- `cmb`: A `CMB` object containing the lensing and ISW sub-components.
- `grid`: The `CosmologicalGrid` defining the redshift range.
- `bg`: `BackgroundQuantities` for background evolution.
- `cosmo`: The `AbstractCosmology` model.
"""
function evaluate_components!(cmb::CMB, grid::CosmologicalGrid, bg::BackgroundQuantities, cosmo::AbstractCosmology)
    compute_kernel!(cmb.κ, grid, bg, cosmo)
    compute_kernel!(cmb.ISW, grid, bg, cosmo)
end
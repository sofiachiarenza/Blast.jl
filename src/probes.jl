"""
    AbstractCosmologicalProbes{T}
An abstract type for the shear, clustering and CMB lensing kernels.
"""
abstract type AbstractCosmologicalProbes end

abstract type AbstractComponents end

@kwdef mutable struct NullComponent <: AbstractComponents
    Kernel::AbstractArray{<:Any, 2} = zeros(1, 1)
end

function compute_kernel!(Component::NullComponent, grid::CosmologicalGrid, bg::BackgroundQuantities, cosmo::AbstractCosmology) 
    return Component.Kernel = zeros(1, length(grid.z_range))
end

"""
    GalaxyKernel{T} <: AbstractCosmologicalProbes{T}

Represents a galaxy kernel in cosmological calculations, where kernel values are provided for multiple tomographic bins.

# Parameters
- `Kernel::AbstractArray{T, 2}`: A 2D array of type `T`, with dimensions `(n_bins, nχ)`. This stores kernel values for each tomographic bin and a grid of χ values.

# Constructors
- `GalaxyKernel{T}(n_bins::Int, nχ::Int)`: Creates a `GalaxyKernel` with the specified `n_bins` and `nχ` values, initializing the kernel values to zeros of type `T`.
- `GalaxyKernel(n_bins::Int, nχ::Int)`: Creates a `GalaxyKernel` with the specified `n_bins` and `nχ` values, initializing the kernel values to zeros of type `Float64`.
"""
@kwdef mutable struct NumberCounts <: AbstractComponents
    nz::Array{<:Any, 2} = zeros(1, 1)
    z::Array{<:Any, 1} = zeros(1)
    bias::Array{<:Any, 2} = zeros(1, 1)
    Kernel::Array{<:Any, 2} = zeros(1, 1)
end


"""
    ShearKernel{T} <: AbstractCosmologicalProbes{T}

Represents a shear kernel used in cosmological lensing calculations. The kernel is defined over multiple tomographic bins.

# Parameters
- `Kernel::AbstractArray{T, 2}`: A 2D array of type `T`, with dimensions `(n_bins, nχ)`. Stores the kernel values for each tomographic bin and a grid of χ values.

# Constructors
- `ShearKernel{T}(n_bins::Int, nχ::Int)`: Creates a `ShearKernel` with the specified `n_bins` and `nχ` values, initializing the kernel values to zeros of type `T`.
- `ShearKernel(n_bins::Int, nχ::Int)`: Creates a `ShearKernel` with the specified `n_bins` and `nχ` values, initializing the kernel values to zeros of type `Float64`.
"""
@kwdef mutable struct CosmicShear <: AbstractComponents
    nz::Array{<:Any, 2} = zeros(1, 1)
    z::Array{<:Any, 1} = zeros(1)
    Kernel::Array{<:Any, 2} = zeros(1, 1)
end


"""
    CMBLensingKernel{T} <: AbstractCosmologicalProbes{T}

Represents a CMB lensing kernel.

# Parameters
- `Kernel::AbstractArray{T, 1}`: A 1D array of type `T`, with dimension `(nχ)`. Note that CMB Lensing by definition only has a single tomographic bin.

# Constructors
- `CMBLensingKernel{T}(nχ::Int)`: Creates a `CMBLensingKernel` with the specified `nχ` value, initializing the kernel values to zeros of type `T`.
- `CMBLensingKernel(nχ::Int)`: Creates a `CMBLensingKernel` with the specified `nχ` value, initializing the kernel values to zeros of type `Float64`.
"""
@kwdef mutable struct CMBLensing <: AbstractComponents
    Kernel::Array{<:Any, 2} = zeros(1, 1)
end


@kwdef mutable struct RedshiftSpaceDistortions <: AbstractComponents
    nz::Array{<:Any, 2} = zeros(1,1)
    z::Array{<:Any, 1} = zeros(1)
    growth_rate::Array{<:Any, 1} = zeros(1)
    Kernel::Array{<:Any, 2} = zeros(1, 1)
end

@kwdef mutable struct MagnificationBias <: AbstractComponents
    nz::Array{<:Any, 2} = zeros(1, 1)
    z::Array{<:Any, 1} = zeros(1)
    s::Array{<:Any, 2} = zeros(1,1)
    Kernel::Array{<:Any, 2} = zeros(1, 1)
end

@kwdef mutable struct IntrinsicAlignment <: AbstractComponents
    nz::Array{<:Any, 2} = zeros(1, 1)
    z::Array{<:Any, 1} = zeros(1)
    A_IA::Array{<:Any, 2} = zeros(1, 1)
    Kernel::Array{<:Any, 2} = zeros(1, 1)
end

@kwdef mutable struct IntegratedSachsWolfe <: AbstractComponents
    growth_rate::Array{<:Any, 1} = zeros(1)
    Kernel::Array{<:Any, 2} = zeros(1, 1)
end

@kwdef mutable struct PrimordialNonGaussianity <: AbstractComponents
    nz::Array{<:Any, 2} = zeros(1,1)
    z::Array{<:Any, 1} = zeros(1)
    bias::Array{<:Any, 2} = zeros(1, 1)
    f_NL::Number = 0
    p::Number = 0
    Kernel::Array{<:Any, 2} = zeros(1, 1)
end

@kwdef mutable struct GalaxyClustering <: AbstractCosmologicalProbes
    δ::NumberCounts
    RSD::Union{RedshiftSpaceDistortions, NullComponent} = NullComponent()
    μ::Union{MagnificationBias, NullComponent} = NullComponent()
    PNG::Union{PrimordialNonGaussianity, NullComponent} = NullComponent()
end

@kwdef mutable struct WeakLensing <: AbstractCosmologicalProbes
    γ::CosmicShear
    IA::Union{IntrinsicAlignment, NullComponent} = NullComponent()
end

@kwdef mutable struct CMB <: AbstractCosmologicalProbes
    κ::CMBLensing
    ISW::Union{IntegratedSachsWolfe, NullComponent} = NullComponent()
end


"""
    compute_kernel!(nz::AbstractArray{T, 2}, Probe::GalaxyKernel, z::AbstractArray{T, 1},
                    bias::AbstractArray{T,2}, grid::CosmologicalGrid, bg::BackgroundQuantities, 
                    cosmo::AbstractCosmology) where T

Computes the galaxy clustering kernel based on a redshift distribution `nz` and stores it in the `GalaxyKernel` struct. 
The kernel is defined as: 
```math
W_g(\\chi) = \\frac{H(z)}{c}n(z)b(z)
```

# Parameters:
- `nz`: A 2D array of type `T` where each row represents the redshift distribution of galaxies for a specific redshift bin.
- `z`: The redshift grid corresponding to the `nz` array.
- `Probe`: An instance of `GalaxyKernel`, in which the computed kernel values for each redshift bin will be stored.
- `bias`: A 2D array of type `T` where each row represents the value of the b(z) in each tomographic bin.
- `grid`: A `CosmologicalGrid` object specifying the redshift range and grid points for kernel computation.
- `bg`: A struct containing arrays of Hubble parameter (`Hz_array`) and comoving distance (`χz_array`), precomputed over the grid.
- `cosmo`: An instance of a cosmological model used to calculate the background quantities if not already provided.

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

"""
    compute_kernel!(nz::AbstractArray{T, 2}, Probe::ShearKernel, z::AbstractArray{T, 1},
                    grid::CosmologicalGrid, bg::BackgroundQuantities, 
                    cosmo::AbstractCosmology) where T

Computes the weak lensing shear kernel based on a redshift distribution `nz` and stores it in the `ShearKernel` struct. 
The kernel is defined as: 
```math
W_{\\gamma}(\\chi) = \\frac{3}{2}\\frac{H_0^2}{c^2}\\Omega_m \\frac{\\chi}{a(\\chi)}\\int_{z(\\chi)}^{\\infty}dz'n(z')\\frac{\\chi(z')-\\chi}{\\chi(z')}
```

# Parameters:
- `nz`: A 2D array of type `T` where each row corresponds to the redshift distribution for a specific shear redshift bin.
- `z`: The redshift grid corresponding to the nz array.
- `Probe`: An instance of `ShearKernel`, where computed kernel values for each redshift bin will be stored.
- `grid`: A `CosmologicalGrid` object defining the redshift range and grid points for kernel computation.
- `bg`: A struct containing precomputed Hubble parameter (`Hz_array`) and comoving distance (`χz_array`) arrays over the grid.
- `cosmo`: An instance of a cosmological model that provides background parameters needed for lensing kernel calculations.
"""
function compute_kernel!(Component::CosmicShear, grid::CosmologicalGrid, bg::BackgroundQuantities, cosmo::AbstractCosmology) 

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
    compute_kernel!(Probe::CMBLensingKernel, 
                    grid::CosmologicalGrid, bg::BackgroundQuantities, 
                    cosmo::AbstractCosmology)

Computes the CMB lensing kernel and stores it in the `CMBLensingKernel` struct.
The kernel is defined as: 
```math
W_{\\kappa}(\\chi) = \\frac{3}{2}\\frac{H_0^2}{c^2}\\Omega_m \\frac{\\chi}{a(\\chi)}\\frac{\\chi^*-\\chi}{\\chi^*},
```
where ``\\chi^* = \\chi(z_{\\mathrm{CMB}} = 1100)``


# Parameters:
- `Probe`: An instance of `CMBLensingKernel` to store the computed kernel values.
- `grid`: A grid specifying the redshift range over which the kernel is computed.
- `bg`: A struct containing precomputed Hubble parameter and comoving distance values.
- `cosmo`: A cosmological model.
"""
function compute_kernel!(Component::CMBLensing, grid::CosmologicalGrid, bg::BackgroundQuantities, cosmo::AbstractCosmology) 

    prefac = 1.5 * cosmo.H0^2 * cosmo.Ωm / C_LIGHT^2
    χ_CMB = compute_χ(1090, cosmo) #From DESI fiducial cosmology

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

function compute_kernel!(Component::MagnificationBias, grid::CosmologicalGrid, bg::BackgroundQuantities, cosmo::AbstractCosmology) 

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
            int, err = quadgk(x -> integrand(x), z_low, z_top) 

            kernel[b, z_idx] = prefac * bg.χz_array[z_idx] * (1. + grid.z_range[z_idx]) * int / nz_norm
        end
    end
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

        kernel[b,:] = Component.f_NL .* bΦ(Component.bias[b,:], Component.p) .* (nz_func.(grid.z_range) ./ nz_norm)
    end
    Component.Kernel = kernel
end

function evaluate_components!(GC::GalaxyClustering, grid::CosmologicalGrid, bg::BackgroundQuantities, cosmo::AbstractCosmology) 
    compute_kernel!(GC.δ, grid, bg, cosmo)
    compute_kernel!(GC.RSD, grid, bg, cosmo)
    compute_kernel!(GC.μ, grid, bg, cosmo)
    compute_kernel!(GC.PNG, grid, bg, cosmo)
end

function evaluate_components!(WL::WeakLensing, grid::CosmologicalGrid, bg::BackgroundQuantities, cosmo::AbstractCosmology)
    compute_kernel!(WL.γ, grid, bg, cosmo)
    compute_kernel!(WL.IA, grid, bg, cosmo)
end

function evaluate_components!(cmb::CMB, grid::CosmologicalGrid, bg::BackgroundQuantities, cosmo::AbstractCosmology)
    compute_kernel!(cmb.κ, grid, bg, cosmo)
    compute_kernel!(cmb.ISW, grid, bg, cosmo)
end
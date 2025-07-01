abstract type AbstractOLDCosmologicalProbes{T} end

@kwdef mutable struct GalaxyKernel{T} <: AbstractOLDCosmologicalProbes{T}
    Kernel::AbstractArray{T, 2} = zeros(1, 1)
end

GalaxyKernel{T}(n_bins::Int, nχ::Int) where T = GalaxyKernel{T}(Kernel = zeros(T, n_bins, nχ))

GalaxyKernel(n_bins::Int, nχ::Int) = GalaxyKernel{Float64}(n_bins, nχ)


"""
    ShearKernel{T} <: AbstractCosmologicalProbes{T}

Represents a shear kernel used in cosmological lensing calculations. The kernel is defined over multiple tomographic bins.

# Parameters
- `Kernel::AbstractArray{T, 2}`: A 2D array of type `T`, with dimensions `(n_bins, nχ)`. Stores the kernel values for each tomographic bin and a grid of χ values.

# Constructors
- `ShearKernel{T}(n_bins::Int, nχ::Int)`: Creates a `ShearKernel` with the specified `n_bins` and `nχ` values, initializing the kernel values to zeros of type `T`.
- `ShearKernel(n_bins::Int, nχ::Int)`: Creates a `ShearKernel` with the specified `n_bins` and `nχ` values, initializing the kernel values to zeros of type `Float64`.
"""
@kwdef mutable struct ShearKernel{T} <: AbstractOLDCosmologicalProbes{T}
    Kernel::AbstractArray{T, 2} = zeros(1, 1)
end

ShearKernel{T}(n_bins::Int, nχ::Int) where T = ShearKernel{T}(Kernel = zeros(T, n_bins, nχ))

ShearKernel(n_bins::Int, nχ::Int) = ShearKernel{Float64}(n_bins, nχ)

#TODO: missing docs
@kwdef mutable struct RSDKernel{T} <: AbstractOLDCosmologicalProbes{T}
    Kernel::AbstractArray{T, 2} = zeros(1, 1)
end

RSDKernel{T}(n_bins::Int, nχ::Int) where T = RSDKernel{T}(Kernel = zeros(T, n_bins, nχ))

RSDKernel(n_bins::Int, nχ::Int) = RSDKernel{Float64}(n_bins, nχ)

@kwdef mutable struct LensMagKernel{T} <: AbstractOLDCosmologicalProbes{T}
    Kernel::AbstractArray{T, 2} = zeros(1, 1)
end

LensMagKernel{T}(n_bins::Int, nχ::Int) where T = LensMagKernel{T}(Kernel = zeros(T, n_bins, nχ))

LensMagKernel(n_bins::Int, nχ::Int) = LensMagKernel{Float64}(n_bins, nχ)

"""
    CMBLensingKernel{T} <: AbstractCosmologicalProbes{T}

Represents a CMB lensing kernel.

# Parameters
- `Kernel::AbstractArray{T, 1}`: A 1D array of type `T`, with dimension `(nχ)`. Note that CMB Lensing by definition only has a single tomographic bin.

# Constructors
- `CMBLensingKernel{T}(nχ::Int)`: Creates a `CMBLensingKernel` with the specified `nχ` value, initializing the kernel values to zeros of type `T`.
- `CMBLensingKernel(nχ::Int)`: Creates a `CMBLensingKernel` with the specified `nχ` value, initializing the kernel values to zeros of type `Float64`.
"""
@kwdef mutable struct CMBLensingKernel{T} <: AbstractOLDCosmologicalProbes{T}
    Kernel::AbstractArray{T, 1} = zeros(1)
end

CMBLensingKernel{T}(nχ::Int) where T = CMBLensingKernel{T}(Kernel = zeros(T, nχ))

CMBLensingKernel(nχ::Int) = CMBLensingKernel{Float64}(nχ)

function compute_kernel!(nz::AbstractArray{T, 2}, z::AbstractArray{T, 1}, Probe::GalaxyKernel, 
                        grid::CosmologicalGrid, bg::BackgroundQuantities, 
                        cosmo::AbstractCosmology) where T

    #TODO: this test will suck for autodiff, will need fixing
    if all(iszero, bg.Hz_array) || all(iszero, bg.χz_array)
        evaluate_background_quantities!(grid, bg, cosmo)
    end

    n_bins = size(Probe.Kernel, 1)
    
    for b in 1:n_bins
        nz_func = DataInterpolations.AkimaInterpolation(nz[b,:], z, extrapolation=ExtrapolationType.Extension)
        nz_norm, _ = quadgk(x->nz_func(x), first(grid.z_range), last(grid.z_range))

        Probe.Kernel[b,:] = @. (bg.Hz_array / C_LIGHT) * (nz_func.(grid.z_range) / nz_norm)
    end
end

function compute_kernel!(nz::AbstractArray{T, 2}, z::AbstractArray{T, 1}, Probe::GalaxyKernel, 
    bias::AbstractArray{T,2}, grid::CosmologicalGrid, bg::BackgroundQuantities, 
    cosmo::AbstractCosmology) where T

    #TODO: this test will suck for autodiff, will need fixing
    if all(iszero, bg.Hz_array) || all(iszero, bg.χz_array)
    evaluate_background_quantities!(grid, bg, cosmo)
    end

    n_bins = size(Probe.Kernel, 1)

    for b in 1:n_bins
        nz_func = DataInterpolations.AkimaInterpolation(nz[b,:], z, extrapolation = ExtrapolationType.Extension)
        nz_norm, _ = quadgk(x->nz_func(x), first(grid.z_range), last(grid.z_range))

        Probe.Kernel[b,:] = @. bias[b,:] * (bg.Hz_array / C_LIGHT) * (nz_func.(grid.z_range) / nz_norm)
    end
end


function compute_kernel!(nz::AbstractArray{T, 2}, z::AbstractArray{T, 1}, Probe::ShearKernel, grid::CosmologicalGrid,
    bg::BackgroundQuantities, cosmo::AbstractCosmology) where T

    #TODO: this test will suck for autodiff, will need fixing
    if all(iszero, bg.Hz_array) || all(iszero, bg.χz_array)
        evaluate_background_quantities!(grid, bg, cosmo)
    end

    n_bins = size(Probe.Kernel, 1)

    for b in 1:n_bins
        nz_func = DataInterpolations.AkimaInterpolation(nz[b,:], z, extrapolation=ExtrapolationType.Extension)
        nz_norm, _ = quadgk(x->nz_func(x), first(grid.z_range), last(grid.z_range))

        prefac = 1.5 * cosmo.H0^2 * cosmo.Ωm / C_LIGHT^2

        for z_idx in 1:length(grid.z_range)
            integrand(x) = nz_func(x) * (1. - bg.χz_array[z_idx]/compute_χ(x, cosmo))
            z_low = grid.z_range[z_idx]
            z_top = grid.z_range[end]*1.1 #TODO: check max redshift, with n5k bins, lensing5 fallisce se uso valore diverso da 3.5
            int, err = quadgk(x -> integrand(x), z_low, z_top) 

            Probe.Kernel[b, z_idx] = prefac * bg.χz_array[z_idx] * (1. + grid.z_range[z_idx]) * int / nz_norm
        end
    end
end


function compute_kernel!(Probe::CMBLensingKernel, grid::CosmologicalGrid, 
    bg::BackgroundQuantities, cosmo::AbstractCosmology) 

    #TODO: this test will suck for autodiff, will need fixing
    if all(iszero, bg.Hz_array) || all(iszero, bg.χz_array)
        evaluate_background_quantities!(grid, bg, cosmo)
    end

    prefac = 1.5 * cosmo.H0^2 * cosmo.Ωm / C_LIGHT^2
    χ_CMB = compute_χ(1100., cosmo)

    Probe.Kernel = @. prefac * bg.χz_array * (1. + grid.z_range) * (1 - bg.χz_array/χ_CMB)
end

function compute_kernel!(nz::AbstractArray{T, 2}, z::AbstractArray{T, 1}, Probe::RSDKernel, 
    growth_factor::AbstractArray{T,1}, grid::CosmologicalGrid,  bg::BackgroundQuantities) where T

    n_bins = size(Probe.Kernel, 1)

    for b in 1:n_bins
        nz_func = DataInterpolations.AkimaInterpolation(nz[b,:], z, extrapolation=ExtrapolationType.Extension)
        nz_norm, _ = quadgk(x->nz_func(x), first(grid.z_range), last(grid.z_range))

        Probe.Kernel[b,:] = @. growth_factor * (bg.Hz_array / C_LIGHT) * (nz_func.(grid.z_range) / nz_norm) #TODO: might be missing C factors
    end
end


function compute_kernel!(nz::AbstractArray{T, 2}, z::AbstractArray{T, 1}, Probe::LensMagKernel, s::AbstractArray{T, 2}, grid::CosmologicalGrid,
    bg::BackgroundQuantities, cosmo::AbstractCosmology) where T

    #TODO: this test will suck for autodiff, will need fixing
    if all(iszero, bg.Hz_array) || all(iszero, bg.χz_array)
        evaluate_background_quantities!(grid, bg, cosmo)
    end

    n_bins = size(Probe.Kernel, 1)

    for b in 1:n_bins
        nz_func = DataInterpolations.AkimaInterpolation(nz[b,:], z, extrapolation=ExtrapolationType.Extension)
        nz_norm, _ = quadgk(x->nz_func(x), first(grid.z_range), last(grid.z_range))

        s_z = DataInterpolations.AkimaInterpolation(s[b,:], z, extrapolation=ExtrapolationType.Extension)

        prefac = 1.5 * cosmo.H0^2 * cosmo.Ωm / C_LIGHT^2

        for z_idx in 1:length(grid.z_range)
            integrand(x) = nz_func(x) * (1. - bg.χz_array[z_idx]/compute_χ(x, cosmo)) * (5 .* s_z(x) .- 2)
            z_low = grid.z_range[z_idx]
            z_top =  grid.z_range[end]*1.1 #TODO: check max redshift, with n5k bins, lensing5 fallisce se uso valore diverso da 3.5
            int, err = quadgk(x -> integrand(x), z_low, z_top) 

            Probe.Kernel[b, z_idx] = prefac * bg.χz_array[z_idx] * (1. + grid.z_range[z_idx]) * int / nz_norm
        end
    end
end

"""
    w_ell_tullio(c::AbstractArray,T::AbstractArray)
Compute the tensor contraction of the chebyshev coefficients of the power spectrum 'c' and the precomputed
integrals 'T' to obtain the projected matter densities.
"""
function w_ell_tullio(c::AbstractArray, T::AbstractArray)
    return @tullio w[i,j,k] := c[j,k,l] * T[i,j,k,l]
end

function P_phi(k::AbstractArray{<:Any,1}, cosmo::AbstractCosmology)
    return @. 9/25 * 2 * π^2 * cosmo.As / (k^3) * (k/0.05)^(cosmo.ns - 1.)
end

function extract_transfer_function(pk::AbstractArray{<:Any,2}, k::AbstractArray{<:Any, 1}, cosmo::AbstractCosmology)
    prim_pk = get_PΦ(k , cosmo)
    return sqrt.(pk ./ reshape(prim_pk, 1, :))
end
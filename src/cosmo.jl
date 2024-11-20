"""
    AbstractCosmology{T}
An abstract type representing a general cosmological model.
"""
abstract type AbstractCosmology{T} end

"""
    AbstractCosmologicalGrid{T}
An abstract type representing a grid on which cosmological quantities are evaluated.
"""
abstract type AbstractCosmologicalGrid{T} end

"""
    AbstractBackgroundQuantities{T}
An abstract type for background quantities in cosmology, such as the Hubble parameter (`H`), comoving distance (`χ`), 
and the growth factor (`D`).
"""
abstract type AbstractBackgroundQuantities{T} end

"""
    AbstractCosmologicalProbes{T}
An abstract type for the shear, clustering and CMB lensing kernels.
"""
abstract type AbstractCosmologicalProbes{T} end


# Define the flat ΛCDM cosmological model with default parameters based on the fiducial N5K cosmology.

"""
    FlatΛCDM{T}(; w0 = -1.0, wa = 0.0, H0 = 67.27, Ωm = 0.3156, Ωb = 0.0492, 
        Ωde = 0.6844, As = 2.12107e-9, σ8 = 0.816, Ωk = 0.0, Ωr = 0.0, ns = 0.9645
    )

# Parameters:
- `w0`: Dark energy equation of state parameter at present time (default: -1).
- `wa`: Time evolution of the dark energy equation of state (default: 0).
- `H0`: Hubble constant in km/s/Mpc (default: 67.27).
- `Ωm`: Matter density parameter (default: 0.3156).
- `Ωb`: Baryon density parameter (default: 0.0492).
- `Ωde`: Dark energy density parameter (default: 0.6844).
- `As`: Scalar amplitude of the primordial power spectrum (default: 2.12107e-9).
- `σ8`: Root-mean-square density fluctuation in spheres of radius 8 Mpc (default: 0.816).
- `Ωk`: Curvature density parameter (default: 0, for flat universe).
- `Ωr`: Radiation density parameter (default: 0, since radiation is negligible at low redshift).
- `ns`: Scalar spectral index (default: 0.9645).
"""
@kwdef mutable struct FlatΛCDM{T} <: AbstractCosmology{T}
    #TODO: comology will need updates, A_s and sigma8 are not independent of eachother, need more classes...
    w0::T  = -1.0
    wa::T  = 0.0
    H0::T  = 67.27
    Ωm::T  = 0.3156
    Ωb::T  = 0.0492
    Ωde::T = 0.6844
    As::T  = 2.12107e-9
    σ8::T  = 0.816
    Ωk::T  = 0.0
    Ωr::T  = 0.0
    ns::T  = 0.9645
end

"""
    CosmologicalGrid{T}(; z_range, k_range)
# Parameters:
- `z_range`: Array of redshift values where quantities like the Hubble parameter are evaluated (default: LinRange(0.001, 2.5, 300)).
- `k_range`: Array of wavenumbers for evaluating power spectra or other k-dependent quantities (default: LogSpaced(1e-5, 50, 1000)).
"""
@kwdef mutable struct CosmologicalGrid{T} <: AbstractCosmologicalGrid{T}
    z_range::AbstractArray{T} = LinRange(0.001, 2.5, 300)
    k_range::AbstractArray{T} = LinRange(1e-5, 50., 1000) # TODO: Switch to Chebyshev points for better interpolation.
end


"""
    BackgroundQuantities{T}(; Hz_array, χz_array)

# Parameters:
- `Hz_array`: Array of Hubble parameter values, evaluated on a grid of redshift values (default: zeros(500)).
- `χz_array`: Array of comoving distance values, evaluated on a grid of redshift values (default: zeros(500)).
"""
@kwdef mutable struct BackgroundQuantities{T} <: AbstractBackgroundQuantities{T}
    Hz_array::Vector{T} = zeros(500)  # TODO: How do I make it adaptable to general needs?
    χz_array::Vector{T} = zeros(500)
    # Dz_array::Vector{T} = zeros(500)  # Growth factor array (could be added when power spectrum information is available).
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
@kwdef mutable struct GalaxyKernel{T} <: AbstractCosmologicalProbes{T}
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
@kwdef mutable struct ShearKernel{T} <: AbstractCosmologicalProbes{T}
    Kernel::AbstractArray{T, 2} = zeros(1, 1)
end

ShearKernel{T}(n_bins::Int, nχ::Int) where T = ShearKernel{T}(Kernel = zeros(T, n_bins, nχ))

ShearKernel(n_bins::Int, nχ::Int) = ShearKernel{Float64}(n_bins, nχ)


"""
    CMBLensingKernel{T} <: AbstractCosmologicalProbes{T}

Represents a CMB lensing kernel.

# Parameters
- `Kernel::AbstractArray{T, 1}`: A 1D array of type `T`, with dimension `(nχ)`. Note that CMB Lensing by definition only has a single tomographic bin.

# Constructors
- `CMBLensingKernel{T}(nχ::Int)`: Creates a `CMBLensingKernel` with the specified `nχ` value, initializing the kernel values to zeros of type `T`.
- `CMBLensingKernel(nχ::Int)`: Creates a `CMBLensingKernel` with the specified `nχ` value, initializing the kernel values to zeros of type `Float64`.
"""
@kwdef mutable struct CMBLensingKernel{T} <: AbstractCosmologicalProbes{T}
    Kernel::AbstractArray{T, 1} = zeros(1)
end

CMBLensingKernel{T}(nχ::Int) where T = CMBLensingKernel{T}(Kernel = zeros(T, nχ))

CMBLensingKernel(nχ::Int) = CMBLensingKernel{Float64}(nχ)





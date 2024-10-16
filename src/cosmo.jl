# Just a starting point, mostly adapted from cosmocentral.jl. Needs a lot of expansions. 

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



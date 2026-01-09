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
    z_range::AbstractArray{T} = LinRange(0, 5, 300)
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
    # Dz_array::Vector{T} = zeros(500)  # TODO: Growth factor array (could be added when power spectrum information is available).
end


"""
    compute_adimensional_hubble_factor(z::T, cosmo::FlatΛCDM) -> T

Computes the adimensional Hubble factor `E(z)` for a given redshift `z`, using the 
cosmological parameters from a `FlatΛCDM` model.
The analitycal expression is given by:
```math
E(z)=\\sqrt{\\Omega_m(1+z)^3+\\Omega_r(1+z)^4+
\\Omega_{de}(1+z)^{3(1+w_0+w_a)}\\exp(-3w_a \\frac{z}{1+z})+\\Omega_k(1+z)^2}
```

# Parameters:
- `z`: Redshift at which to evaluate the Hubble factor.
- `cosmo`: A `FlatΛCDM` cosmological model containing parameters like Ωm, Ωr, Ωde, etc.

# Returns:
- `E_z`: The adimensional Hubble factor at redshift `z`.
"""
function compute_adimensional_hubble_factor(z::T, cosmo::FlatΛCDM) where T
    E_z = compute_adimensional_hubble_factor(z, cosmo.Ωm, cosmo.Ωr,
        cosmo.Ωde, cosmo.Ωk, cosmo.w0, cosmo.wa)
    return E_z
end

"""
    compute_adimensional_hubble_factor(z::T, Ωm::T, Ωr::T, Ωde::T, Ωk::T, w0::T, wa::T) -> T

Computes the adimensional Hubble factor `E(z)` given the redshift `z` and individual cosmological parameters.
The analitycal expression is given by:
```math
E(z)=\\sqrt{\\Omega_m(1+z)^3+\\Omega_r(1+z)^4+
\\Omega_{de}(1+z)^{3(1+w_0+w_a)}\\exp(-3w_a \\frac{z}{1+z})+\\Omega_k(1+z)^2}
```

# Parameters:
- `z`: Redshift at which to evaluate the Hubble factor.
- `Ωm`: Matter density parameter.
- `Ωr`: Radiation density parameter.
- `Ωde`: Dark energy density parameter.
- `Ωk`: Curvature density parameter.
- `w0`: Dark energy equation of state parameter at the present time.
- `wa`: Time evolution of the dark energy equation of state.

# Returns:
- `E_z`: The adimensional Hubble factor at redshift `z`.
"""
function compute_adimensional_hubble_factor(z::T, Ωm::T, Ωr::T, Ωde::T, Ωk::T, w0::T, wa::T) where T
    E_z = sqrt(Ωm*(1+z)^3 + Ωr*(1+z)^4 + Ωk*(1+z)^2 +
        Ωde*(1+z)^(3*(1+w0+wa))*exp(-3*wa*z/(1+z)))
    return E_z
end

"""
    compute_hubble_factor(z::T, cosmo::AbstractCosmology) -> T

Computes the Hubble parameter `H(z)` at a given redshift `z` using the Hubble constant `H0` and the adimensional 
Hubble factor `E(z)`.

# Parameters:
- `z`: Redshift at which to compute the Hubble parameter.
- `cosmo`: A cosmological model that contains `H0` and other necessary parameters.

# Returns:
- `H_z`: The Hubble parameter at redshift `z`.
"""
function compute_hubble_factor(z::T, cosmo::AbstractCosmology) where T
    H_z = cosmo.H0 * compute_adimensional_hubble_factor(z, cosmo)
    return H_z
end

"""
    compute_χ(z::T, cosmo::AbstractCosmology) -> T

Computes the comoving distance `χ(z)` to a given redshift `z` by numerically integrating 
the inverse of the adimensional Hubble factor `E(z)`:
```math
\\chi(z)=\\frac{c}{H_0}\\int_0^z \\frac{dz'}{E(z')}
```

# Parameters:
- `z`: Redshift up to which the comoving distance is computed.
- `cosmo`: A cosmological model containing the necessary parameters (e.g., Ωm, H0).

# Returns:
- `χ_z`: The comoving distance at redshift `z` in units of Mpc.
"""
function compute_χ(z::T, cosmo::AbstractCosmology) where T
    integral, err = quadgk(x -> 1. / compute_adimensional_hubble_factor(x, cosmo), 0., z, rtol=1e-12)
    return integral * C_LIGHT / cosmo.H0
end

"""
    evaluate_background_quantities!(grid::CosmologicalGrid, bg::BackgroundQuantities, cosmo::AbstractCosmology)

Populates the `BackgroundQuantities` struct with values for the Hubble parameter `H(z)` and comoving distance `χ(z)` 
over the redshift range specified by the `CosmologicalGrid`.

# Parameters:
- `grid`: A grid specifying the redshift range over which to evaluate the background quantities.
- `bg`: A mutable struct where the computed `H(z)` and `χ(z)` values will be stored.
- `cosmo`: A cosmological model containing the necessary parameters (e.g., H0, Ωm).

# Notes:
This function modifies the `BackgroundQuantities` struct in place by filling the arrays with the computed values.
"""
function evaluate_background_quantities!(grid::CosmologicalGrid,
    #TODO: works for now, will need vectorization and rethinking in the future.
    bg::BackgroundQuantities, cosmo::AbstractCosmology)
    for z_idx in 1:length(grid.z_range)
        # Compute the Hubble parameter H(z)
        bg.Hz_array[z_idx] = compute_hubble_factor(grid.z_range[z_idx], cosmo)
        
        # Compute the comoving distance χ(z)
        bg.χz_array[z_idx] = compute_χ(grid.z_range[z_idx], cosmo)
    end
end

#TODO I think this is not used anywhere anynore, might get rid of it.
"""
    resample_redshifts(bg::BackgroundQuantities, grid::AbstractCosmologicalGrid, new_χ::AbstractArray{T,1}) where T

Resamples redshift values corresponding to a new set of comoving distances using Akima interpolation.

# Arguments
- `bg::BackgroundQuantities`: An object containing background cosmological quantities, 
  including the mapping between redshift (`z`) and comoving distance (`χ`).
- `grid::AbstractCosmologicalGrid`: An object defining the range of redshifts (`z_range`) 
  and associated cosmological grid quantities.
- `new_χ::AbstractArray{T,1}`: A 1D array of comoving distances for which corresponding redshift values are desired.

# Returns
- `resampled_z::AbstractArray{T,1}`: A 1D array of resampled redshift values corresponding to the input comoving distances `new_χ`.
"""
function resample_redshifts(bg::BackgroundQuantities, grid::AbstractCosmologicalGrid, new_χ::AbstractArray{T,1}) where T
    z_of_χ = AkimaInterpolation(grid.z_range, bg.χz_array, extrapolation=ExtrapolationType.Extension)
    return z_of_χ.(new_χ)
end
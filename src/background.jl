"""
    compute_adimensional_hubble_factor(z::T, cosmo::FlatΛCDM) -> T

Computes the adimensional Hubble factor `E(z)` for a given redshift `z`, using the 
cosmological parameters from a `FlatΛCDM` model.

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
    compute_hubble_factor(z::T, AbstractCosmology::AbstractCosmology) -> T

Computes the Hubble parameter `H(z)` at a given redshift `z` using the Hubble constant `H0` and the adimensional 
Hubble factor `E(z)`.

# Parameters:
- `z`: Redshift at which to compute the Hubble parameter.
- `AbstractCosmology`: A cosmological model that contains `H0` and other necessary parameters.

# Returns:
- `H_z`: The Hubble parameter at redshift `z`.
"""
function compute_hubble_factor(z::T, AbstractCosmology::AbstractCosmology) where T
    H_z = AbstractCosmology.H0 * compute_adimensional_hubble_factor(z, AbstractCosmology)
    return H_z
end

"""
    compute_χ(z::T, AbstractCosmology::AbstractCosmology) -> T

Computes the comoving distance `χ(z)` to a given redshift `z` by numerically integrating 
the inverse of the adimensional Hubble factor `E(z)`.

# Parameters:
- `z`: Redshift up to which the comoving distance is computed.
- `AbstractCosmology`: A cosmological model containing the necessary parameters (e.g., Ωm, H0).

# Returns:
- `χ_z`: The comoving distance at redshift `z` in units of Mpc.
"""
function compute_χ(z::T, AbstractCosmology::AbstractCosmology) where T
    integral, err = quadgk(x -> 1. / compute_adimensional_hubble_factor(x, AbstractCosmology), 0., z, rtol=1e-12)
    return integral * C_LIGHT / AbstractCosmology.H0
end

"""
    evaluate_background_quantities!(CosmologicalGrid::CosmologicalGrid, BackgroundQuantities::BackgroundQuantities, AbstractCosmology::AbstractCosmology)

Populates the `BackgroundQuantities` struct with values for the Hubble parameter `H(z)` and comoving distance `χ(z)` 
over the redshift range specified by the `CosmologicalGrid`.

# Parameters:
- `CosmologicalGrid`: A grid specifying the redshift range over which to evaluate the background quantities.
- `BackgroundQuantities`: A mutable struct where the computed `H(z)` and `χ(z)` values will be stored.
- `AbstractCosmology`: A cosmological model containing the necessary parameters (e.g., H0, Ωm).

# Notes:
This function modifies the `BackgroundQuantities` struct in place by filling the arrays with the computed values.
"""
function evaluate_background_quantities!(CosmologicalGrid::CosmologicalGrid,
    BackgroundQuantities::BackgroundQuantities,
    AbstractCosmology::AbstractCosmology)
    for z_idx in 1:length(CosmologicalGrid.z_range)
        # Compute the Hubble parameter H(z)
        BackgroundQuantities.Hz_array[z_idx] = compute_hubble_factor(
            CosmologicalGrid.z_range[z_idx], AbstractCosmology)
        
        # Compute the comoving distance χ(z)
        BackgroundQuantities.χz_array[z_idx] = compute_χ(
            CosmologicalGrid.z_range[z_idx], AbstractCosmology)
    end
end


function compute_kernel!(nz::Vector{T}, AbstractCosmologicalProbes::GalaxyKernel, CosmologicalGrid::CosmologicalGrid,
    BackgroundQuantities::BackgroundQuantities,
    AbstractCosmology::AbstractCosmology) where T

    nz_func = DataInterpolations.AkimaInterpolation(nz, CosmologicalGrid.z_range, extrapolate=true)
    nz_norm, _ = quadgk(x->nz_func(x), first(CosmologicalGrid.z_range), last(CosmologicalGrid.z_range))

    #TODO: check if the background quantities have already been computed or not

    AbstractCosmologicalProbes.Kernel = @. (BackgroundQuantities.Hz_array / C_LIGHT) * (nz/ nz_norm)
end


function compute_kernel!(nz::Vector{T}, AbstractCosmologicalProbes::ShearKernel, CosmologicalGrid::CosmologicalGrid,
    BackgroundQuantities::BackgroundQuantities,
    AbstractCosmology::AbstractCosmology) where T

    #TODO: check if the background quantities have already been computed or not

    nz_func = DataInterpolations.AkimaInterpolation(nz, CosmologicalGrid.z_range, extrapolate=true)
    nz_norm, _ = quadgk(x->nz_func(x), first(CosmologicalGrid.z_range), last(CosmologicalGrid.z_range))

    prefac = 1.5 * AbstractCosmology.H0^2 * AbstractCosmology.Ωm / C_LIGHT^2

    for z_idx in 1:length(CosmologicalGrid.z_range)
        integrand(x) = nz_func(x) / nz_norm * (1. - BackgroundQuantities.χz_array[z_idx]/compute_χ(x, AbstractCosmology))
        z_low = CosmologicalGrid.z_range[z_idx]
        z_top = 5 #TODO: check max redshift
        int, err = quadgk(x -> integrand(x), z_low, z_top) 

        AbstractCosmologicalProbes.Kernel[z_idx] = prefac * BackgroundQuantities.χz_array[z_idx] * (1. + CosmologicalGrid.z_range[z_idx]) * int 
    end
end

function compute_kernel!(AbstractCosmologicalProbes::CMBLensingKernel, CosmologicalGrid::CosmologicalGrid,
    BackgroundQuantities::BackgroundQuantities,
    AbstractCosmology::AbstractCosmology)

    #TODO: check if the background quantities have already been computed or not

    prefac = 1.5 * AbstractCosmology.H0^2 * AbstractCosmology.Ωm / C_LIGHT^2
    χ_CMB = compute_χ(1100., AbstractCosmology)

    AbstractCosmologicalProbes.Kernel = @. prefac * BackgroundQuantities.χz_array * (1. + CosmologicalGrid.z_range) * (1 - BackgroundQuantities.χz_array/χ_CMB)
end









#Il growth factor D(z) va estratto da un power spectrum

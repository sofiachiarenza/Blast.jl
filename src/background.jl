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
function evaluate_background_quantities!(grid::CosmologicalGrid,
    bkgq::BackgroundQuantities, cosmo::AbstractCosmology)
    for z_idx in 1:length(grid.z_range)
        # Compute the Hubble parameter H(z)
        bkgq.Hz_array[z_idx] = compute_hubble_factor(grid.z_range[z_idx], cosmo)
        
        # Compute the comoving distance χ(z)
        bkgq.χz_array[z_idx] = compute_χ(grid.z_range[z_idx], cosmo)
    end
end


"""
    compute_kernel!(nz::AbstractArray{T, 2}, AbstractCosmologicalProbes::GalaxyKernel, 
                    CosmologicalGrid::CosmologicalGrid, BackgroundQuantities::BackgroundQuantities, 
                    AbstractCosmology::AbstractCosmology) where T

Computes the galaxy clustering kernel based on a redshift distribution `nz` and stores it in the `GalaxyKernel` struct. 

# Parameters:
- `nz`: A 2D array of type `T` where each row represents the redshift distribution of galaxies for a specific redshift bin.
- `AbstractCosmologicalProbes`: An instance of `GalaxyKernel`, in which the computed kernel values for each redshift bin will be stored.
- `CosmologicalGrid`: A `CosmologicalGrid` object specifying the redshift range and grid points for kernel computation.
- `BackgroundQuantities`: A struct containing arrays of Hubble parameter (`Hz_array`) and comoving distance (`χz_array`), precomputed over the grid.
- `AbstractCosmology`: An instance of a cosmological model used to calculate the background quantities if not already provided.

"""
function compute_kernel!(nz::AbstractArray{T, 2}, Probe::GalaxyKernel, 
                        grid::CosmologicalGrid, bg::BackgroundQuantities, 
                        cosmo::AbstractCosmology) where T

    #TODO: this test will suck for autodiff, will need fixing
    if all(iszero, bg.Hz_array) || all(iszero, bg.χz_array)
        evaluate_background_quantities!(grid, bg, cosmo)
    end

    n_bins = size(Probe.Kernel, 1)
    
    for b in 1:n_bins
        nz_func = DataInterpolations.AkimaInterpolation(nz[b,:], grid.z_range, extrapolate=true)
        nz_norm, _ = quadgk(x->nz_func(x), first(grid.z_range), last(grid.z_range))

        Probe.Kernel[b,:] = @. (bg.Hz_array / C_LIGHT) * (nz[b,:] / nz_norm)
    end
end

"""
    compute_kernel!(nz::AbstractArray{T, 2}, AbstractCosmologicalProbes::ShearKernel, 
                    CosmologicalGrid::CosmologicalGrid, BackgroundQuantities::BackgroundQuantities, 
                    AbstractCosmology::AbstractCosmology) where T

Computes the weak lensing shear kernel based on a redshift distribution `nz` and stores it in the `ShearKernel` struct. 

# Parameters:
- `nz`: A 2D array of type `T` where each row corresponds to the redshift distribution for a specific shear redshift bin.
- `AbstractCosmologicalProbes`: An instance of `ShearKernel`, where computed kernel values for each redshift bin will be stored.
- `CosmologicalGrid`: A `CosmologicalGrid` object defining the redshift range and grid points for kernel computation.
- `BackgroundQuantities`: A struct containing precomputed Hubble parameter (`Hz_array`) and comoving distance (`χz_array`) arrays over the grid.
- `AbstractCosmology`: An instance of a cosmological model that provides background parameters needed for lensing kernel calculations.
"""
function compute_kernel!(nz::AbstractArray{T, 2}, Probe::ShearKernel, grid::CosmologicalGrid,
    bg::BackgroundQuantities, cosmo::AbstractCosmology) where T

    #TODO: this test will suck for autodiff, will need fixing
    if all(iszero, bg.Hz_array) || all(iszero, bg.χz_array)
        evaluate_background_quantities!(grid, bg, cosmo)
    end

    n_bins = size(Probe.Kernel, 1)

    for b in 1:n_bins
        nz_func = DataInterpolations.AkimaInterpolation(nz[b,:], grid.z_range, extrapolate=true)
        nz_norm, _ = quadgk(x->nz_func(x), first(grid.z_range), last(grid.z_range))

        prefac = 1.5 * cosmo.H0^2 * cosmo.Ωm / C_LIGHT^2

        for z_idx in 1:length(grid.z_range)
            integrand(x) = nz_func(x) * (1. - bg.χz_array[z_idx]/compute_χ(x, cosmo))
            z_low = grid.z_range[z_idx]
            z_top = 5 #TODO: check max redshift, with n5k bins, lensing5 fallisce se uso valore diverso da 3.5
            int, err = quadgk(x -> integrand(x), z_low, z_top) 

            Probe.Kernel[b, z_idx] = prefac * bg.χz_array[z_idx] * (1. + grid.z_range[z_idx]) * int / nz_norm
        end
    end
end

"""
    compute_kernel!(AbstractCosmologicalProbes::CMBLensingKernel, 
                    CosmologicalGrid::CosmologicalGrid, BackgroundQuantities::BackgroundQuantities, 
                    AbstractCosmology::AbstractCosmology)

Computes the CMB lensing kernel and stores it in the `CMBLensingKernel` struct.

# Parameters:
- `AbstractCosmologicalProbes`: An instance of `CMBLensingKernel` to store the computed kernel values.
- `CosmologicalGrid`: A grid specifying the redshift range over which the kernel is computed.
- `BackgroundQuantities`: A struct containing precomputed Hubble parameter and comoving distance values.
- `AbstractCosmology`: A cosmological model.
"""
function compute_kernel!(Probe::CMBLensingKernel, grid::CosmologicalGrid,
    bg::BackgroundQuantities, cosmo::AbstractCosmology) 

    #TODO: this test will suck for autodiff, will need fixing
    if all(iszero, bg.Hz_array) || all(iszero, bg.χz_array)
        evaluate_background_quantities!(grid, bg, cosmo)
    end

    n_bins = size(Probe.Kernel, 1)

    if n_bins > 1
        throw(DomainError("CMB Lensing must have a single tomographic bin!"))
    end

    prefac = 1.5 * cosmo.H0^2 * cosmo.Ωm / C_LIGHT^2
    χ_CMB = compute_χ(1100., cosmo)

    Probe.Kernel[1,:] = @. prefac * bg.χz_array * (1. + grid.z_range) * (1 - bg.χz_array/χ_CMB)
end
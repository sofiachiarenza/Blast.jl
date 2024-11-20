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


"""
    compute_kernel!(nz::AbstractArray{T, 2}, Probe::GalaxyKernel, z::AbstractArray{T, 1},
                    grid::CosmologicalGrid, bg::BackgroundQuantities, 
                    cosmo::AbstractCosmology) where T

Computes the galaxy clustering kernel based on a redshift distribution `nz` and stores it in the `GalaxyKernel` struct. 
The kernel is defined as: 
```math
W_g(\\chi) = \\frac{H(z)}{c}n(z)
```

# Parameters:
- `nz`: A 2D array of type `T` where each row represents the redshift distribution of galaxies for a specific redshift bin.
- `z`: The redshift grid corresponding to the `nz` array.
- `Probe`: An instance of `GalaxyKernel`, in which the computed kernel values for each redshift bin will be stored.
- `grid`: A `CosmologicalGrid` object specifying the redshift range and grid points for kernel computation.
- `bg`: A struct containing arrays of Hubble parameter (`Hz_array`) and comoving distance (`χz_array`), precomputed over the grid.
- `cosmo`: An instance of a cosmological model used to calculate the background quantities if not already provided.

"""
function compute_kernel!(nz::AbstractArray{T, 2}, z::AbstractArray{T, 1}, Probe::GalaxyKernel, 
                        grid::CosmologicalGrid, bg::BackgroundQuantities, 
                        cosmo::AbstractCosmology) where T

    #TODO: this test will suck for autodiff, will need fixing
    if all(iszero, bg.Hz_array) || all(iszero, bg.χz_array)
        evaluate_background_quantities!(grid, bg, cosmo)
    end

    n_bins = size(Probe.Kernel, 1)
    
    for b in 1:n_bins
        nz_func = DataInterpolations.AkimaInterpolation(nz[b,:], z, extrapolate=true)
        nz_norm, _ = quadgk(x->nz_func(x), first(grid.z_range), last(grid.z_range))

        Probe.Kernel[b,:] = @. (bg.Hz_array / C_LIGHT) * (nz_func.(grid.z_range) / nz_norm)
    end
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
function compute_kernel!(nz::AbstractArray{T, 2}, z::AbstractArray{T, 1}, Probe::ShearKernel, grid::CosmologicalGrid,
    bg::BackgroundQuantities, cosmo::AbstractCosmology) where T

    #TODO: this test will suck for autodiff, will need fixing
    if all(iszero, bg.Hz_array) || all(iszero, bg.χz_array)
        evaluate_background_quantities!(grid, bg, cosmo)
    end

    n_bins = size(Probe.Kernel, 1)

    for b in 1:n_bins
        nz_func = DataInterpolations.AkimaInterpolation(nz[b,:], z, extrapolate=true)
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
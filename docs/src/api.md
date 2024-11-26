# API reference

```@docs
# Types
Blast.FlatΛCDM
Blast.BackgroundQuantities
Blast.AbstractCosmology
Blast.AbstractBackgroundQuantities
Blast.AbstractCosmologicalGrid
Blast.AbstractCosmologicalProbes
Blast.CosmologicalGrid
Blast.GalaxyKernel
Blast.ShearKernel
Blast.CMBLensingKernel

# Functions
Blast.compute_T̃
Blast.bessel_cheb_eval
Blast.get_clencurt_weights
Blast.get_clencurt_grid
Blast.w_ell_tullio
Blast.plan_fft
Blast.fast_chebcoefs
Blast.make_grid
Blast.grid_interpolator(::Union{Blast.GalaxyKernel, Blast.ShearKernel},::Blast.BackgroundQuantities, ::Vector{T}) where T
Blast.grid_interpolator(::Blast.CMBLensingKernel, ::Blast.BackgroundQuantities, ::Vector{T}) where T
Blast.get_kernel_array(::Blast.GalaxyKernel, ::Blast.BackgroundQuantities, ::Vector{T}) where T
Blast.get_kernel_array(::Blast.ShearKernel, ::Blast.BackgroundQuantities, ::Vector{T}) where T
Blast.get_kernel_array(::Blast.CMBLensingKernel, ::Blast.BackgroundQuantities, ::Vector{T}) where T
Blast.combine_kernels
Blast.factorial_frac
Blast.get_ell_prefactor(::Blast.GalaxyKernel, ::Blast.GalaxyKernel, ::Vector)
Blast.get_ell_prefactor(::Blast.GalaxyKernel, ::Blast.ShearKernel, ::Vector)
Blast.get_ell_prefactor(::Blast.ShearKernel, ::Blast.ShearKernel, ::Vector)
Blast.get_ell_prefactor(::Blast.CMBLensingKernel, ::Blast.ShearKernel, ::Vector)
Blast.get_ell_prefactor(::Blast.CMBLensingKernel, ::Blast.CMBLensingKernel, ::Vector)
Blast.get_ell_prefactor(::Blast.CMBLensingKernel, ::Blast.GalaxyKernel, ::Vector)
Blast.simpson_weight_array
Blast.get_clencurt_weights_R_integration
Blast.compute_Cℓ(::AbstractArray{T, 3}, 
               ::Union{Blast.GalaxyKernel, Blast.ShearKernel, Blast.CMBLensingKernel}, 
               ::Union{Blast.GalaxyKernel, Blast.ShearKernel, Blast.CMBLensingKernel}, 
               ::Blast.BackgroundQuantities, 
               R::AbstractVector, 
               ℓ_list::AbstractArray{T,1} = Blast.ℓ) where T
Blast.compute_Cℓ(::AbstractArray{T, 3}, 
               ::AbstractArray{T, 4}, 
               ::Blast.BackgroundQuantities, 
               ::AbstractArray{T, 1}, 
               ::AbstractArray{T, 1}, 
               ::AbstractArray{T,1}) where T
Blast.compute_adimensional_hubble_factor(::T, ::Blast.FlatΛCDM) where T
Blast.compute_adimensional_hubble_factor(::T, ::T, ::T, ::T, ::T, ::T, ::T) where T
Blast.compute_hubble_factor
Blast.compute_χ
Blast.evaluate_background_quantities!
Blast.resample_redshifts
Blast.compute_kernel!(::AbstractArray{T, 2}, ::AbstractArray{T, 1}, ::Blast.GalaxyKernel, ::Blast.CosmologicalGrid, ::Blast.BackgroundQuantities, ::Blast.AbstractCosmology) where T
Blast.compute_kernel!(::AbstractArray{T, 2}, ::AbstractArray{T, 1}, ::Blast.ShearKernel, ::Blast.CosmologicalGrid, ::Blast.BackgroundQuantities, ::Blast.AbstractCosmology) where T
Blast.compute_kernel!(::Blast.CMBLensingKernel, ::Blast.CosmologicalGrid, ::Blast.BackgroundQuantities, ::Blast.AbstractCosmology)
Blast.chebyshev_polynomials
Blast.interpolate_power_spectrum
Blast.unequal_time_power_spectrum
```
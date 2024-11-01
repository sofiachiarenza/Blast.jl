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
Blast.grid_interpolator
Blast.get_kernel_array(::Blast.GalaxyKernel, ::Blast.BackgroundQuantities, ::Vector{T}) where T
Blast.get_kernel_array(::Union{Blast.ShearKernel, Blast.CMBLensingKernel}, ::Blast.BackgroundQuantities, ::Vector{T}) where T
Blast.combine_kernels
Blast.factorial_frac
Blast.get_ell_prefactor
Blast.simpson_weight_array
Blast.compute_Cℓ
Blast.compute_adimensional_hubble_factor(z::T, cosmo::Blast.FlatΛCDM) where T
Blast.compute_adimensional_hubble_factor(z::T, Ωm::T, Ωr::T, Ωde::T, Ωk::T, w0::T, wa::T) where T
Blast.compute_hubble_factor
Blast.compute_χ
Blast.evaluate_background_quantities!
Blast.compute_kernel!(nz::Vector{T}, ::Blast.GalaxyKernel, ::Blast.CosmologicalGrid, ::Blast.BackgroundQuantities, ::Blast.AbstractCosmology) where T
Blast.compute_kernel!(nz::Vector{T}, ::Blast.ShearKernel, ::Blast.CosmologicalGrid, ::Blast.BackgroundQuantities, ::Blast.AbstractCosmology) where T
Blast.compute_kernel!(::Blast.CMBLensingKernel, ::Blast.CosmologicalGrid, ::Blast.BackgroundQuantities, ::Blast.AbstractCosmology)
```
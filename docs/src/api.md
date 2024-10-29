# API reference

```@docs
Blast.compute_T̃
Blast.bessel_cheb_eval
Blast.get_clencurt_weights
Blast.get_clencurt_grid
Blast.w_ell_tullio
Blast.plan_fft
Blast.fast_chebcoefs
Blast.simpson_weight_array
Blast.compute_Cℓ
Blast.FlatΛCDM
Blast.BackgroundQuantities
Blast.AbstractCosmology
Blast.AbstractBackgroundQuantities
Blast.AbstractCosmologicalGrid
Blast.ShearKernel
Blast.AbstractCosmologicalProbes
Blast.CosmologicalGrid
Blast.CMBLensingKernel
Blast.GalaxyKernel
Blast.compute_adimensional_hubble_factor(z::T, cosmo::Blast.FlatΛCDM) where T
Blast.compute_adimensional_hubble_factor(z::T, Ωm::T, Ωr::T, Ωde::T, Ωk::T, w0::T, wa::T) where T
Blast.compute_hubble_factor
Blast.compute_χ
Blast.evaluate_background_quantities!
Blast.compute_kernel!(nz::Vector{T}, Blast.AbstractCosmologicalProbes::Blast.GalaxyKernel, Blast.CosmologicalGrid::Blast.CosmologicalGrid, Blast.BackgroundQuantities::Blast.BackgroundQuantities, Blast.AbstractCosmology::Blast.AbstractCosmology) where T
Blast.compute_kernel!(nz::Vector{T}, Blast.AbstractCosmologicalProbes::Blast.ShearKernel, Blast.CosmologicalGrid::Blast.CosmologicalGrid, Blast.BackgroundQuantities::Blast.BackgroundQuantities, Blast.AbstractCosmology::Blast.AbstractCosmology) where T
Blast.compute_kernel!(Blast.AbstractCosmologicalProbes::Blast.CMBLensingKernel, Blast.CosmologicalGrid::Blast.CosmologicalGrid, Blast.BackgroundQuantities::Blast.BackgroundQuantities, Blast.AbstractCosmology::Blast.AbstractCosmology)
```
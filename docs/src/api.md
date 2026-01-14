# API Reference

This page documents the public API of `Blast.jl`, organized by functionality
and level of abstraction.

Most users should only interact with **high-level probe constructors** and
the `get_Cℓ` interface.

## Cosmology

Types and functions defining the background cosmology.

```@docs
Blast.AbstractCosmology
Blast.FlatΛCDM

Blast.compute_adimensional_hubble_factor
Blast.compute_hubble_factor
Blast.compute_χ
```

## Background Quantities and Grids

Background quantities derived from the cosmology and used to build kernels.
These objects store interpolants of background functions evaluated on the
global comoving distance grid.

```@docs
Blast.AbstractBackgroundQuantities
Blast.BackgroundQuantities
Blast.AbstractCosmologicalGrid
Blast.CosmologicalGrid
Blast.evaluate_background_quantities!
```

## Probes and Components

Each observable is decomposed into physical components (e.g. number density, RSD,
magnification, primordial non-Gaussianity). These components store kernels
as functions of comoving distance.
`AbstractCosmologicalProbes` objects group physical components (`AbstractComponents` objects) defining the Observables.

```@docs
Blast.AbstractCosmologicalProbes
Blast.AbstractComponents
Blast.evaluate_components!
```

### Galaxy Clustering Probe
```@docs
Blast.GalaxyClustering

Blast.NumberCounts
Blast.compute_kernel!(Component::Blast.NumberCounts, grid::Blast.CosmologicalGrid, bg::Blast.BackgroundQuantities, cosmo::Blast.AbstractCosmology)
Blast.RedshiftSpaceDistortions
Blast.compute_kernel!(Component::Blast.RedshiftSpaceDistortions, grid::Blast.CosmologicalGrid, bg::Blast.BackgroundQuantities, cosmo::Blast.AbstractCosmology)
Blast.MagnificationBias
Blast.compute_kernel!(Component::Blast.MagnificationBias, grid::Blast.CosmologicalGrid, bg::Blast.BackgroundQuantities, cosmo::Blast.AbstractCosmology)
Blast.PrimordialNonGaussianity
Blast.compute_kernel!(Component::Blast.PrimordialNonGaussianity, grid::Blast.CosmologicalGrid, bg::Blast.BackgroundQuantities, cosmo::Blast.AbstractCosmology)
```

### Weak Lensing Probe
```@docs
Blast.WeakLensing

Blast.CosmicShear
Blast.compute_kernel!(Component::Blast.CosmicShear, grid::Blast.CosmologicalGrid, bg::Blast.BackgroundQuantities, cosmo::Blast.AbstractCosmology)
Blast.IntrinsicAlignment
Blast.compute_kernel!(Component::Blast.IntrinsicAlignment, grid::Blast.CosmologicalGrid, bg::Blast.BackgroundQuantities, cosmo::Blast.AbstractCosmology)
```

### CMB Probe
```@docs
Blast.CMB

Blast.CMBLensing
Blast.compute_kernel!(Component::Blast.CMBLensing, grid::Blast.CosmologicalGrid, bg::Blast.BackgroundQuantities, cosmo::Blast.AbstractCosmology)
Blast.IntegratedSachsWolfe
Blast.compute_kernel!(Component::Blast.IntegratedSachsWolfe, grid::Blast.CosmologicalGrid, bg::Blast.BackgroundQuantities, cosmo::Blast.AbstractCosmology)
```

## Initialization and Setup

```@docs
Blast.FFTPlans
Blast.SetUp
Blast.plan_fft
Blast.fast_chebcoefs
Blast.AbstractCoeff
Blast.cϕTT
Blast.cϕT
Blast.cϕ
Blast.make_coeff
Blast.build_coeff
```

## Power Spectrum Handling 

Container holding all necessary matter power spectrum information.
This object has to be rebuilt when the cosmology or input power spectra change (i.e. in a MCMC chain).

```@docs
Blast.PowerSpectrum
Blast.prepare_pk_workspace
Blast.get_PΦ
Blast.get_Tm
Blast.transform_to_R_frame
```

## Projected Matter Density
The projected matter density stores the inner k integral, whose efficient computation is at the core of the Blast algorithm.

```@docs
Blast.ProjectedMatterDensityComponent
Blast.ProjectedMatterDensity
Blast.compute_w!
Blast.w_ell_tullio
```

## Angular Power Spectrum Computation

This is the main numerical engine of BLAST. Angular power spectra are
computed using a hybrid approach:
- non-Limber computation at low multipoles
- Limber approximation at high multipoles
- smooth interpolation between regimes

Users should generally call get_Cℓ.
```@docs
Blast.get_Cℓ
```

#### Internal routines
```@docs
Blast.make_grid
Blast.grid_interpolator
Blast.get_kernel_array
Blast.combine_kernels
Blast.compute_Cℓ(Component1::Blast.AbstractComponents, Component2::Blast.AbstractComponents, w::Blast.ProjectedMatterDensityComponent) 
Blast.compute_Cℓ(Component1::Blast.AbstractComponents, Component2::Blast.AbstractComponents, w02::Blast.ProjectedMatterDensityComponent, 
                    w20::Blast.ProjectedMatterDensityComponent) 
Blast.get_limber_kernel(Component::Blast.AbstractComponents)
Blast.get_limber_kernel(G::Blast.GalaxyClustering)
Blast.get_limber_kernel(G::Blast.WeakLensing)
Blast.get_limber_kernel(G::Blast.CMB)
Blast.get_limber_correction(Probe::Union{Blast.GalaxyClustering, Blast.WeakLensing, Blast.CMB}, pk::Blast.PowerSpectrum)
Blast.get_limber_correction(ProbeA::Union{Blast.GalaxyClustering, Blast.WeakLensing, Blast.CMB}, ProbeB::Union{Blast.GalaxyClustering, Blast.WeakLensing, Blast.CMB}, pk::Blast.PowerSpectrum)
Blast.get_limber_Cℓ(Probe::Union{Blast.GalaxyClustering, Blast.WeakLensing, Blast.CMB}, pk::Blast.PowerSpectrum)
Blast.get_limber_Cℓ(ProbeA::Union{Blast.GalaxyClustering, Blast.WeakLensing, Blast.CMB}, ProbeB::Union{Blast.GalaxyClustering, Blast.WeakLensing, Blast.CMB}, pk::Blast.PowerSpectrum)
```

## Numerical routines and Utilities
```@docs
Blast.get_clencurt_grid
Blast.get_clencurt_weights
Blast.get_clencurt_weights_R_integration
Blast.bessel_second_derivative
Blast.bessel_cheb_eval
Blast.compute_T̃
Blast.factorial_frac
Blast.bΦ
Blast._akima_interpolation(u, t, t_new)
Blast._akima_slopes
Blast._akima_coefficients
Blast._akima_eval
Blast._akima_interpolation(u::AbstractMatrix, t, t_new)
```



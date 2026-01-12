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
Blast.resample_redshifts
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
Blast.compute_kernel!
Blast.RedshiftSpaceDistortions
Blast.compute_kernel!
Blast.MagnificationBias
Blast.compute_kernel!
Blast.PrimordialNonGaussianity
Blast.compute_kernel!
```

### Weak Lensing Probe
```@docs
Blast.WeakLensing

Blast.CosmicShear
Blast.compute_kernel!
Blast.IntrinsicAlignment
Blast.compute_kernel!
```

### CMB Probe
```@docs
Blast.CMB

Blast.CMBLensing
Blast.compute_kernel!
Blast.IntegratedSachsWolfe
Blast.compute_kernel!
```

## Initialization and Setup

```@docs
Blast.SetUp
```

## Power Spectrum Handling 

Container holding all necessary matter power spectrum information.
This object has to be rebuilt when the cosmology or input power spectra change (i.e. in a MCMC chain).

```@docs
Blast.PowerSpectrum
Blast.prepare_pk_workspace
Blast.get_PΦ
Blast.get_Tm
```

## Projected Matter Density
The projected matter density stores the inner k integral, whose efficient computation is at the core of the Blast algorithm.

```@docs
Blast.ProjectedMatterDensity
Blast.compute_w!
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
Blast.compute_Cℓ
```

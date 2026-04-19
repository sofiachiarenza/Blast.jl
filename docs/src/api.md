# API Reference

This page documents the current API of `Blast.jl`.

## Cosmology

### Constructors

```@docs
Blast.w0waCDM
```

### Background and derived quantities

```@docs
Blast.Background
Blast.compute_hubble_factor
Blast.compute_χ
Blast.get_H0
Blast.get_Ωm
Blast.get_Ωb
Blast.get_Ωc
Blast.get_As
Blast.get_ns
```

## Probes

This section lists the concrete user-facing probe containers and components.
Internal abstract supertypes are intentionally omitted from the rendered API page.

```@docs
Blast.GalaxyClustering
Blast.WeakLensing
Blast.CMB
```

## Probe components

```@docs
Blast.NumberCounts
Blast.RedshiftSpaceDistortions
Blast.MagnificationBias
Blast.PrimordialNonGaussianity
Blast.CosmicShear
Blast.IntrinsicAlignment
Blast.CMBLensing
Blast.IntegratedSachsWolfe
```

### n(z) helpers

```@docs
Blast.prepare_nz_matrix
Blast.smooth_nz
Blast.NLA_model
```

### Kernel evaluation

```@docs
Blast.compute_kernel!(::Blast.NumberCounts, ::Blast.Background)
Blast.compute_kernel!(::Blast.RedshiftSpaceDistortions, ::Blast.Background)
Blast.compute_kernel!(::Blast.MagnificationBias, ::Blast.Background)
Blast.compute_kernel!(::Blast.PrimordialNonGaussianity, ::Blast.Background)
Blast.compute_kernel!(::Blast.CosmicShear, ::Blast.Background)
Blast.compute_kernel!(::Blast.IntrinsicAlignment, ::Blast.Background)
Blast.compute_kernel!(::Blast.CMBLensing, ::Blast.Background)
Blast.compute_kernel!(::Blast.IntegratedSachsWolfe, ::Blast.Background)
Blast.evaluate_components!
```

## Setup and Power Spectrum Workspace

```@docs
Blast.FFTPlans
Blast.SetUp
Blast.cϕTT
Blast.cϕT
Blast.cϕ
Blast.build_coeff
Blast.PowerSpectrum
Blast.prepare_pk_workspace
Blast.get_PΦ
Blast.get_Tm
Blast.transform_to_R_frame
```

## Projected Matter Density

```@docs
Blast.ProjectedMatterDensityComponent
Blast.ProjectedMatterDensity
Blast.compute_w!
Blast.w_ell_tullio
```

## Angular Power Spectra

```@docs
Blast.get_Cℓ
Blast.compute_Cℓ(::Blast.AbstractComponents, ::Blast.AbstractComponents, ::Blast.ProjectedMatterDensityComponent, ::Blast.Background)
Blast.compute_Cℓ(::Blast.AbstractComponents, ::Blast.AbstractComponents, ::Blast.ProjectedMatterDensityComponent, ::Blast.ProjectedMatterDensityComponent, ::Blast.Background)
```

## Limber helpers

```@docs
Blast.get_limber_kernel(::Blast.AbstractComponents)
Blast.get_limber_kernel(::Blast.GalaxyClustering)
Blast.get_limber_kernel(::Blast.WeakLensing)
Blast.get_limber_kernel(::Blast.CMB)
Blast.get_limber_correction(::Union{Blast.GalaxyClustering, Blast.WeakLensing, Blast.CMB}, ::Blast.PowerSpectrum)
Blast.get_limber_correction(::Union{Blast.GalaxyClustering, Blast.WeakLensing, Blast.CMB}, ::Union{Blast.GalaxyClustering, Blast.WeakLensing, Blast.CMB}, ::Blast.PowerSpectrum)
Blast.get_limber_Cℓ(::Union{Blast.GalaxyClustering, Blast.WeakLensing, Blast.CMB}, ::Blast.PowerSpectrum)
Blast.get_limber_Cℓ(::Union{Blast.GalaxyClustering, Blast.WeakLensing, Blast.CMB}, ::Union{Blast.GalaxyClustering, Blast.WeakLensing, Blast.CMB}, ::Blast.PowerSpectrum)
```

## Numerical utilities

```@docs
Blast.get_clencurt_grid
Blast.get_clencurt_weights
Blast.get_clencurt_weights_R_integration
Blast.bessel_second_derivative
Blast.bessel_cheb_eval
Blast.compute_T̃
```
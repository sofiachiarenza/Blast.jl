# Blast.jl

`Blast.jl` (**B**eyond **L**imber **A**ngular power **S**pectra **T**oolkit)
computes angular power spectra without the Limber approximation, via a
Chebyshev decomposition of the matter power spectrum. It covers a full
6x2pt framework — galaxy clustering (density, redshift-space distortions,
magnification bias, scale-dependent bias from primordial non-Gaussianity),
weak lensing (cosmic shear, intrinsic alignment), and CMB lensing (with the
integrated Sachs-Wolfe effect) — and is differentiable end-to-end via
`ChainRules`/`Mooncake`.

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/sofiachiarenza/Blast.jl")
```

## Quickstart

```julia
using Blast

cosmo = Blast.w0waCDM()
bg = Blast.Background(cosmo)

δ  = Blast.NumberCounts(nz=nz_g, z=bg.z, bias=bias)  # your n(z)/bias data
GC = Blast.GalaxyClustering(δ=δ)
Blast.evaluate_components!(GC, bg)

γ  = Blast.CosmicShear(nz=nz_s, z=bg.z)              # your n(z) data
WL = Blast.WeakLensing(γ=γ)
Blast.evaluate_components!(WL, bg)

W, plans = Blast.SetUp(GC, WL)
PowerSpectrum = Blast.prepare_pk_workspace(plans, pk, pk_limber_lin, pk_limber_nonlin, bg)
W = Blast.compute_w(W, PowerSpectrum)

Cℓ_GG = Blast.get_Cℓ(ell, GC, PowerSpectrum, W, bg, plans)
Cℓ_SS = Blast.get_Cℓ(ell, WL, PowerSpectrum, W, bg, plans)
Cℓ_GS = Blast.get_Cℓ(ell, GC, WL, PowerSpectrum, W, bg, plans)
```

See the [Example](example.md) page for a complete, runnable walkthrough with
every 6x2pt effect switched on, [The algorithm](alg.md) for the underlying
math, and the [API](api.md) page for the full reference.

## Authors
Sofia Chiarenza, PhD Student at Waterloo Centre for Astrophysics.
Marco Bonici, PostDoctoral Researcher at Waterloo Centre for Astrophysics.

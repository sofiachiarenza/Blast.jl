<p align="left">
<img width="300px" src="https://github.com/user-attachments/assets/dc268ab5-7ff8-40f1-bc37-9d3a1f356d99"/>
</p>

# Blast.jl

| **Documentation** | **Build Status** | **More Info** |
|:------------------:|:----------------:|:-------------:|
| [![Docs - Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://sofiachiarenza.github.io/Blast.jl/dev) [![Docs - Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://sofiachiarenza.github.io/Blast.jl/stable) | [![Build Status](https://github.com/sofiachiarenza/Blast.jl/workflows/CI/badge.svg)](https://github.com/sofiachiarenza/Blast.jl/actions) [![Code Coverage](https://codecov.io/github/sofiachiarenza/Blast.jl/graph/badge.svg?token=8QLDGERO9H)](https://codecov.io/github/sofiachiarenza/Blast.jl) | [![arXiv](https://img.shields.io/badge/arXiv-2410.03632-b31b1b.svg)](https://arxiv.org/abs/2410.03632) ![Size](https://img.shields.io/github/repo-size/sofiachiarenza/Blast.jl) |

This repo contains the Beyond Limber Angular power Spectra Toolkit, `Blast.jl`. The code is entirely written in `Julia` and provides differentiable functions to compute angular power spectra for the auto and cross correlation of galaxy clustering (with redshift-space distortions, magnification bias, and scale-dependent bias from primordial non-Gaussianity), weak lensing (cosmic shear with intrinsic alignment), and CMB lensing (with the integrated Sachs-Wolfe effect) — a full 6x2pt framework, differentiable end-to-end via `ChainRules`/`Mooncake`.

## Installation

In order to install `Blast.jl`, run from the `Julia` REPL

```julia
using Pkg
Pkg.add(url="https://github.com/sofiachiarenza/Blast.jl")
```

### Development setup (temporary)

This branch requires the unreleased `AbstractCosmologicalEmulators` 0.12
(branch `96-b-spline`), which is not in the General registry yet, so a plain
`Pkg.add`/`Pkg.instantiate` fails to resolve. From a clone of this repository,
pin it first, then instantiate/test as usual:

```bash
julia --project=. scripts/bootstrap_ace_dev.jl
julia --project=. -e 'using Pkg; Pkg.test()'
# docs environment:
julia --project=docs/ scripts/bootstrap_ace_dev.jl
```

The script (and its CI steps) will be removed once a registered ACE release
provides these APIs.

## Usage
After installing it, you can start instantiating the objects needed to compute the $C_\ell$'s. To begin, initialize the cosmological model and background quantities:

```julia
using Blast

cosmo = Blast.w0waCDMCosmology(
    ln10Aₛ=3.054505771332324,
    nₛ=0.9645,
    h=0.6727,
    ωb=0.022264244268,
    ωc=0.120552737256,
    ωk=0.0,
    mν=0.06,
    w0=-1.0,
    wa=0.0,
)
bg = Blast.Background(cosmo)
```

Next, build the galaxy clustering and weak lensing probes out of their
components, and evaluate their kernels against the background:

```julia
n_bins_g, n_bins_s = 10, 5
nz_g = rand(n_bins_g, length(bg.z))   # Example n(z), replace with actual data
nz_s = rand(n_bins_s, length(bg.z))   # Example n(z), replace with actual data
bias = ones(n_bins_g, length(bg.z))   # Example linear bias, replace with actual model

δ = Blast.NumberCounts(nz=nz_g, z=bg.z, bias=bias)
raw_GC = Blast.GalaxyClustering(δ=δ)
GC = Blast.prepare_probe(raw_GC, bg)

γ = Blast.CosmicShear(nz=nz_s, z=bg.z)
raw_WL = Blast.WeakLensing(γ=γ)
WL = Blast.prepare_probe(raw_WL, bg)
```

Set up the FFT plans/workspace (once, given the probe configuration), then
prepare the Chebyshev-decomposed power-spectrum workspace and the projected
matter density for your cosmology:

```julia
W, plans = Blast.SetUp(GC, WL)

# pk, pk_limber_lin, pk_limber_nonlin: your matter power spectrum, tabulated
# on Blast's internal (z, k) grids — see docs/src/example.md for how to build these
PowerSpectrum = Blast.prepare_pk_workspace(plans, pk, pk_limber_lin, pk_limber_nonlin, bg)
W = Blast.compute_w(W, PowerSpectrum)
```

Finally, compute the angular power spectra for clustering, shear, and their
cross-correlation:

```julia
ell = ...  # your desired multipoles

Cℓ_GG = Blast.get_Cℓ(ell, GC, PowerSpectrum, W, bg, plans)   # (n_ell, n_bins_g, n_bins_g)
Cℓ_SS = Blast.get_Cℓ(ell, WL, PowerSpectrum, W, bg, plans)   # (n_ell, n_bins_s, n_bins_s)
Cℓ_GS = Blast.get_Cℓ(ell, GC, WL, PowerSpectrum, W, bg, plans)  # (n_ell, n_bins_g, n_bins_s)
```

This example only switches on the basic density and cosmic-shear components.
See [`docs/src/example.md`](docs/src/example.md) (or the built docs) for a
complete 6x2pt example with every effect switched on (redshift-space
distortions, magnification bias, primordial non-Gaussianity, and intrinsic
alignment).

### Experimental non-Limber B-spline backends

`SetUp` still uses the legacy Chebyshev backend by default. This branch also
supports explicit direct and fiducial-ratio cubic B-spline tensor bundles.
These are experimental: passing algebra, API and gradient tests does **not**
certify physical-spectrum accuracy or establish a recommended coefficient count.
Supply a portable bundle manifest and choose its count explicitly; the package
does not generate or select a bank for you.

The following pattern assumes the prepared probes/background above, a bundle
path `manifest`, a chosen coefficient count `N`, and your own callable
`pmm(k, z)` returning **total-matter** linear power in Mpc³ with k in Mpc⁻¹:

```julia
bank = Blast.load_tensor_bank(manifest; count=N, representation=:ratio)
# Use representation=:direct for a direct-expansion bundle instead.
backend = Blast.CubicBSplineBackend(bank)
W_template, plans = Blast.SetUp(GC, WL; backend=backend)
grid = Blast.nonlimber_grid(plans)

# Non-Limber sampling is backend-specific. Array order is (z, k).
pk = [pmm(k, z) for z in grid.z, k in grid.k]
# Supply the same separate Limber-grid arrays as in the existing workflow.
PowerSpectrum = Blast.prepare_pk_workspace(
    plans, pk, pk_limber_lin, pk_limber_nonlin, bg; pk_grid=grid)
W = Blast.compute_w(W_template, PowerSpectrum)
Cℓ_GG = Blast.get_Cℓ(ell, GC, PowerSpectrum, W, bg, plans)
```

Do not reuse a Chebyshev-sampled `pk` array on a spline grid just because its
shape matches. Passing `pk_grid` checks coordinates as well as shape. Workspaces
and plans from different bank signatures cannot be mixed. No emulator is loaded
implicitly; `pmm` and the Limber inputs remain the caller's responsibility.

The ratio backend divides the **complete physical source** by fixed fiducial
normalizers before solving for coefficients: qTT = Pmm,fid,
qT = √(PΦ,fid Pmm,fid), and qPHI = PΦ,fid. It does not rescale the input `pk`
in place or recompute the fiducial model as sampled parameters change. The
bundle carries the corresponding distinct TT, T and Φ tensor families.

Banks are shared read-only. Give each concurrent chain its own `SetUp` plans,
mutable workspace and prepared differentiation cache. The functional
`compute_w` path above is suitable for generic-number checks; do not pass
ForwardDiff dual inputs into a preallocated Float64-only workspace.

## Citing 

If you use `Blast.jl` in your research, please cite:

S. Chiarenza, M. Bonici, W. Percival, M. White [_BLAST: Beyond Limber Angular power Spectra Toolkit. A fast and efficient algorithm for 3x2pt analysis_](https://arxiv.org/abs/2410.03632)

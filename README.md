<p align="left">
<img width="300px" src="https://github.com/user-attachments/assets/dc268ab5-7ff8-40f1-bc37-9d3a1f356d99"/>
</p>

# Blast.jl

| **Documentation** | **Build Status** | **More Info** |
|:------------------:|:----------------:|:-------------:|
| [![Docs - Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://sofiachiarenza.github.io/Blast.jl/dev) [![Docs - Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://sofiachiarenza.github.io/Blast.jl/stable) | [![Build Status](https://github.com/sofiachiarenza/Blast.jl/workflows/CI/badge.svg)](https://github.com/sofiachiarenza/Blast.jl/actions) [![Code Coverage](https://codecov.io/github/sofiachiarenza/Blast.jl/graph/badge.svg?token=8QLDGERO9H)](https://codecov.io/github/sofiachiarenza/Blast.jl) | [![arXiv](https://img.shields.io/badge/arXiv-2410.03632-b31b1b.svg)](https://arxiv.org/abs/2410.03632) ![Size](https://img.shields.io/github/repo-size/sofiachiarenza/Blast.jl) |

This repo contains the Beyond Limber Angular power Spectra Toolkit, `Blast.jl`. The code is entirely written in `Julia` and provides the functions to compute angular power spectra for the auto and cross correlation of three different probes: galaxy clustering, shear, and CMB lensing. 

## Installation

In order to install `Blast.jl`, run from the `Julia` REPL

```julia
using Pkg
Pkg.add(url="https://github.com/sofiachiarenza/Blast.jl")
```

## Usage
After installing it, you can start instantiating the objects needed to compute the $C_\ell$'s. To begin, initialize the cosmological model and background quantities:

```julia
cosmo = Blast.FlatΛCDM()
z_range = LinRange(0., 4.0, 1000)
grid = Blast.CosmologicalGrid(z_range=z_range)
bg = Blast.BackgroundQuantities(Hz_array=zeros(length(z_range)), χz_array=zeros(length(z_range)))
Blast.evaluate_background_quantities!(grid, bg, cosmo)
```

Next, define and compute the kernels for galaxy clustering, weak lensing, and CMB lensing:

```julia
GK = Blast.GalaxyKernel(10, length(grid.z_range))
SHK = Blast.ShearKernel(10, length(grid.z_range))
CMBK = Blast.CMBLensingKernel(length(grid.z_range))

nz = rand(n_bins, nz)  # Example n(z), replace with actual data
Blast.compute_kernel!(nz, grid.z_range, GK, grid, bg, cosmo) #compute clustering kernel, repeat for the other probes
```

Load the precomputed inner integrals $\Tilde{T}^{AB}_\ell(\chi_1,\chi_2)$ and evaluate the coefficients of the Chebyshev decomposition of the power spectrum:

```julia
T_LL = Blast.T_tilde_p2  # Lensing-Lensing
T_CL = Blast.T_tilde_0   # Clustering-Lensing
T_CC = Blast.T_tilde_m2  # Clustering-Clustering

plan = Blast.plan_fft(Pk, 1)
cheb_coeff = Blast.fast_chebcoefs(Pk, plan)
```
Using the Chebyshev coefficients, compute the projected matter densities $w_\ell$:

```julia
w_LL = Blast.w_ell_tullio(cheb_coeff, T_LL)
w_CL = Blast.w_ell_tullio(cheb_coeff, T_CL)
w_CC = Blast.w_ell_tullio(cheb_coeff, T_CC)
```

Finally, compute the angular power spectra for clustering, shear, CMB lensing and cross-correlations:

```julia
clustering_Cℓ = Blast.compute_Cℓ(w_CC, GK, GK, bg, R)
shear_Cℓ = Blast.compute_Cℓ(w_LL, SHK, SHK, bg, R)
cross_Cℓ = Blast.compute_Cℓ(w_CL, GK, SHK, bg, R);
cl_cross_cmb_Cℓ = Blast.compute_Cℓ(w_CL, GK, CMBK, bg, R)
sh_cross_cmb_Cℓ = Blast.compute_Cℓ(w_LL, SHK, CMBK, bg, R);
```

## Citing 

If you use `Blast.jl` in your research, please cite:

S. Chiarenza, M. Bonici, W. Percival, M. White [_BLAST: Beyond Limber Angular power Spectra Toolkit. A fast and efficient algorithm for 3x2pt analysis_](https://arxiv.org/abs/2410.03632)


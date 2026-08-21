# Example

This page walks through a complete 6x2pt $C_\ell$ computation with every
available effect switched on: galaxy clustering with density, redshift-space
distortions (RSD), magnification bias, and scale-dependent bias from
primordial non-Gaussianity (PNG); weak lensing with cosmic shear and
intrinsic alignment (NLA model).

## Load data

You'll need, for each probe, an $n(z)$ per tomographic bin and (for galaxy
clustering) a linear bias and magnification-bias slope $s$ per bin — plus
the total-matter power spectrum tabulated on Blast's internal grids: `pk` on
the non-Limber $(z, k)$ grid (`Blast.z_lin`, `Blast.k_cheb`), and
`pk_limber_lin`/`pk_limber_nonlin` on the Limber grid (`Blast.z_cheb`,
`Blast.k_limber`).

```julia
using Blast
using NPZ

nz_g = ...  # (n_bins_g, n_z) clustering n(z)
nz_s = ...  # (n_bins_s, n_z) lensing n(z)
zn   = ...  # shared redshift grid for both n(z) above

ell_list = ...  # multipoles to evaluate the Cℓ's at

pk               = ...  # non-Limber Pmm, shape (Blast.z_lin, Blast.k_cheb)
pk_limber_lin    = ...  # linear Pmm, shape (Blast.z_cheb, Blast.k_limber)
pk_limber_nonlin = ...  # nonlinear Pmm, same grid as above
```

## Cosmology and background

```julia
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

## Nuisance parameters

```julia
bias_vals = [1.376695, 1.451179, 1.528404, 1.607983, 1.689579,
             1.772899, 1.857700, 1.943754, 2.030887, 2.118943]
bias_matrix = ones(size(nz_g, 1), length(bg.z)) .* bias_vals

s_vals = [0.1648, 0.2496, 0.2708, 0.3300, 0.3880, 0.2960, 0.3580, 0.3960, 0.4320, 0.5680]
s_matrix = s_vals .* ones(size(nz_g, 1), length(bg.z))
```

## Galaxy clustering — all effects on

$\delta$ (density) + RSD + magnification bias $\mu$ + primordial
non-Gaussianity (PNG):

```julia
f_NL = 40
p = 1

δ   = Blast.NumberCounts(nz=nz_g, z=zn, bias=bias_matrix)
rsd = Blast.RedshiftSpaceDistortions(nz=nz_g, z=zn, growth_rate=bg.f)
μ   = Blast.MagnificationBias(nz=nz_g, z=zn, s=s_matrix)
PNG = Blast.PrimordialNonGaussianity(nz=nz_g, z=zn, bias=bias_matrix, f_NL=f_NL, p=p)

raw_GC = Blast.GalaxyClustering(δ=δ, RSD=rsd, μ=μ, PNG=PNG)
GC = Blast.prepare_probe(raw_GC, bg)
```

## Weak lensing — all effects on

Cosmic shear $\gamma$ + intrinsic alignment using the built-in NLA model,
$A_{\rm IA}(z) = -A C_1 \Omega_m/[D(z)/D(0)]$:

```julia
γ  = Blast.CosmicShear(nz=nz_s, z=zn)
IA = Blast.IntrinsicAlignment(nz=nz_s, z=zn, A=1.72, C1=0.0134)

raw_WL = Blast.WeakLensing(γ=γ, IA=IA)
WL = Blast.prepare_probe(raw_WL, bg)
```

## Compute the $C_\ell$'s

`SetUp` builds the FFT plans and an empty projected-matter-density workspace
once, given which components are active; `prepare_pk_workspace` builds the
Chebyshev-decomposed power-spectrum workspace for your cosmology;
`compute_w` fills in the projected matter density; `get_Cℓ` is then called
once per probe pair (auto- or cross-correlation):

```julia
W, plans = Blast.SetUp(GC, WL)

PowerSpectrum = Blast.prepare_pk_workspace(plans, pk, pk_limber_lin, pk_limber_nonlin, bg)
W = Blast.compute_w(W, PowerSpectrum)

Cℓ_GG = Blast.get_Cℓ(ell_list, GC, PowerSpectrum, W, bg, plans)       # (n_ell, n_bins_g, n_bins_g)
Cℓ_SS = Blast.get_Cℓ(ell_list, WL, PowerSpectrum, W, bg, plans)       # (n_ell, n_bins_s, n_bins_s)
Cℓ_GS = Blast.get_Cℓ(ell_list, GC, WL, PowerSpectrum, W, bg, plans)   # (n_ell, n_bins_g, n_bins_s)
```

Each `Cℓ_*` array is indexed as `(ℓ, bin_i, bin_j)`. To switch off any
effect, just omit the corresponding keyword when building `GalaxyClustering`
or `WeakLensing` (e.g. `Blast.GalaxyClustering(δ=δ)` for density-only
clustering, as in the [Quickstart](index.md)).

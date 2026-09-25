# ARF validation against ARFCAMB

Checks Blast's angular redshift fluctuations (`AngularRedshiftFluctuations`),
and the galaxy clustering pieces around them, against
[ARFCAMB](https://github.com/chmATiac/ARFCAMB) for a single Gaussian bin
(z₀ = 1, σ = 0.2; see `config.toml`).

| file | role |
|---|---|
| `config.toml` | cosmology, window, bias, s, ℓ range (read by both Python and Julia) |
| `environment.yml` | conda env `arfcamb` for building/running ARFCAMB |
| `export_blast_grids.jl` | dumps Blast's P(k) grids → `grids.npz` |
| `patches/arfcamb_no_gauge_term.patch` | drops CAMB's Newtonian-gauge density term (see below) |
| `generate_camb_reference.py` | runs patched + stock ARFCAMB → `test/reference_fixtures/arfcamb/gaussian_z1.npz` |
| `blast_pieces.jl` | Blast side, shared by the test and the notebook |
| `compare_blast_camb.ipynb` | plots and residual table |

The fixture is committed, so `test/test_arfcamb_validation.jl` (part of
`Pkg.test()`) needs neither CAMB nor Python. Regenerate it only when
`config.toml` or ARFCAMB changes.

## What is compared

The fixture carries CAMB's **linear** P(k, z) on Blast's grids. Blast uses it
directly, so only kernels and projections are tested. Pieces: `dd`, `dd_rsd`,
`dd_mu`, `gg`, `kd`, `kd_rsd`, `kd_mu`, `kg`, and density-only ARF `zz`, `dz`,
`kz`.

- **ARF is density-only.** Blast's kernel is (z − z̄) K^δ. ARFCAMB builds the
  ARF source with the same `SourceTerms` as counts, so with
  `counts_redshift`/`counts_lensing` on it adds ARF-weighted RSD and
  magnification. The strict comparison uses CAMB with `counts_density` only.
  `zz_full`, `dz_full` and `kz_full` are kept for information. The ARF velocity
  term dominates `zz` at ℓ ≲ 20.
- **z̄** is ARFCAMB's internal `Σ z W / Σ W` on the table grid
  (`fortran/SourceWindows.f90`, `SetTable`). It is saved in the fixture and
  passed to Blast. Both codes normalize the ARF window by ∫W.
- **P(k): CDM+baryons vs total.** CAMB sources galaxies with the CDM density
  and κ with the total, but Blast takes one P(k). The default `matter=:cb`
  matches the autos to about 0.3%. `:tot` is about 1% low with mν = 0.06 eV.
- **The reference is a patched ARFCAMB.** CAMB's counts density source is
  `b·δ_c^sync + (… − 3ℋ)σ/k`. The second term converts synchronous-gauge δ_c into
  the observed, Newtonian-gauge number counts. It is a horizon-scale
  relativistic correction that Blast's `NumberCounts` kernel does not model.
  `patches/arfcamb_no_gauge_term.patch` multiplies it by 0 (one line in
  `fortran/equations.f90`), and the strict fixture pieces come from that build.
  Stock ARFCAMB is also run and stored as `<piece>_stock` for reference (not
  plotted): the missing term is 8% on δδ, 3% on zz and 3.5% on κδ at ℓ = 2,
  falling below 0.5% by ℓ ≈ 20.
- **This is not a CAMB accuracy issue.** Against stock CAMB, the residual at
  ℓ < 10 does not move with `AccuracyBoost` 2 → 3, `lAccuracyBoost` 2 → 5,
  `Transfer.high_precision` or `accurate_massive_neutrinos`; all change C_ℓ by
  less than 0.02%. With the patch, Blast agrees down to ℓ = 2.

Current agreement against the patched build, max |Blast − CAMB| / √(C_aa C_bb):

| piece | ℓ ∈ [2, 10) | ℓ ∈ [10, 300] |
|---|---|---|
| dd, dd_rsd, dd_mu, gg | 0.17% | 0.37% |
| kd, kd_rsd, kd_mu, kg | 0.17% | 0.16% |
| zz | 0.30% | 0.33% |
| dz, kz | 0.06% | 0.06% |

The test's tolerance is 0.6% over all ℓ.

## Regenerating the fixture

```bash
# 1. conda env (call the binary directly if `conda` is not on your PATH)
conda env create -f validation/arfcamb/environment.yml

# 2. Two ARFCAMB builds: stock and patched. The code lives on the
#    version_March11th_2025 branch, not main. The repo ships an x86_64
#    camblib.so, so rebuild natively.
conda activate arfcamb
git clone -b version_March11th_2025 https://github.com/chmATiac/ARFCAMB.git ~/Desktop/PhD/ARFCAMB
git clone -b version_March11th_2025 ~/Desktop/PhD/ARFCAMB ~/Desktop/PhD/ARFCAMB_nogauge
git -C ~/Desktop/PhD/ARFCAMB_nogauge apply "$PWD/validation/arfcamb/patches/arfcamb_no_gauge_term.patch"
for d in ARFCAMB ARFCAMB_nogauge; do (cd ~/Desktop/PhD/$d/fortran && python ../setup.py make); done

# 3. fixture (from the Blast.jl root)
julia --project=validation/arfcamb -e 'using Pkg; Pkg.instantiate()'
julia --project=validation/arfcamb validation/arfcamb/export_blast_grids.jl
python validation/arfcamb/generate_camb_reference.py \
    --patched ~/Desktop/PhD/ARFCAMB_nogauge --stock ~/Desktop/PhD/ARFCAMB   # ~4 min
```

The generator imports each build in its own subprocess, by path, so neither
build needs to be `pip install`ed. It refuses a `--patched` build that lacks the
patch, or a `--stock` build that has it.
`setuptools<70` is pinned because ARFCAMB's build imports `pkg_resources`.

## Plots

Open `compare_blast_camb.ipynb` with a Julia kernel; it activates this folder's
project itself.

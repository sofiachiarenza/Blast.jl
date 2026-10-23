# fNL validation against ARFCAMB_nogauge and SwiftCℓ

Checks Blast's local-PNG galaxy-clustering pieces against the two external codes
that model scale-dependent bias, before and after the mixed-PNG transfer-leg
fix of issue #38. All three codes use CAMB's own linear P_cb(k, z), so only
kernels and projections differ.

| file | role |
|---|---|
| `config.toml` | cosmology (same as `../arfcamb`), three Gaussian bins with per-bin bias, f_NL, ℓ range |
| `generate_camb_png.py` | ARFCAMB_nogauge at f_NL = 0, ±f_NL → `camb_png.npz` (also P(k) for Blast and SwiftCℓ) |
| `blast_png.jl` | Blast at f_NL = 0, ±f_NL on CAMB's P(k); `--swift-inputs` also writes `swift_inputs.npz` |
| `swift_png.py` | SwiftCℓ ≥ 0.1.3 at f_NL = 0, ±f_NL → `swift_png.npz` |
| `compare_png.py` | splits each code into δ×PNG (A) and PNG×PNG (B), prints the table, writes `png_comparison.png` |
| `unfolded_check.jl` | standalone toy check of the folded vs unfolded (χ₁, χ₂) integral |

## Running

```bash
conda activate arfcamb
python validation/png_external/generate_camb_png.py --build ~/Desktop/PhD/ARFCAMB_nogauge
julia --project=. validation/png_external/blast_png.jl validation/png_external/blast_fixed.npz --swift-inputs
# pre-fix Blast: any checkout of the commit before the fix, e.g.
#   git archive <commit> | tar -x -C /tmp/blast_old
julia --project=/tmp/blast_old validation/png_external/blast_png.jl validation/png_external/blast_old.npz
conda activate swift_env
python validation/png_external/swift_png.py
python validation/png_external/compare_png.py
```

## Conventions that matter

- **ARFCAMB_nogauge** (CAMB 1.3.2 with `../arfcamb/patches/arfcamb_no_gauge_term.patch`)
  has a built-in local-fNL bias, `(b − 1) 3 δc f_NL Ωm H0² / (c² k² T D)`,
  written through CAMB's own φ(k, a). It fixes p = 1. Stock CAMB ≥ 1.5 can
  instead take a `k_bias` table, but it keeps the Newtonian-gauge density term
  that Blast does not model.
- **SwiftCℓ 0.1.2 is wrong for fNL**: its PNG kernel misses the growth factor
  every other SwiftCℓ kernel carries (fixed in 0.1.3).
- **SwiftCℓ's `b_fNL` is 3 δc (b − p)**, not Blast's b_Φ = 2 δc (b − p): its
  PNG kernel carries an explicit 2/3, with `T(z, k)` = δ_m / Φ (Blast's
  `get_Tm`). Passing b_Φ scales the δ×PNG piece by 2/3 and PNG×PNG by 4/9.
- **Equal b_Φ / b in every bin hides the pre-fix bug in GG.** The pre-fix
  δ×PNG and PNG×δ each put the transfer function on the wrong leg on one
  triangle, but when the PNG kernel is a common multiple of the δ kernel the
  two errors cancel in the summed spectrum. The residual on bin pair ij is
  ∝ (c_j − c_i) with c = f_NL b_Φ / b, so the bins here have different biases.
  A run with b = 1.5 everywhere gives old and fixed identical to machine
  precision. Pairs whose other kernel is not proportional to the PNG one
  (μ, RSD, γ, IA, κ, ISW) are affected for any bias.

## Results (2026-10-23, f_NL = 20, b = 1.2 / 1.8 / 3.0)

δ×PNG piece A, median |X − ref| / √(ref_ii ref_jj) over ℓ = 2–10:

| pair | Blast pre-fix vs CAMB | Blast fixed vs CAMB | SwiftCℓ vs CAMB | pre-fix vs SwiftCℓ | fixed vs SwiftCℓ |
|---|---|---|---|---|---|
| z0.6×z1.2 | 5.8e-3 | 8.2e-4 | 9.2e-4 | 6.7e-3 | 6.9e-4 |
| z0.6×z3 | 4.1e-4 | 1.5e-5 | 8.1e-5 | 4.5e-4 | 6.1e-5 |
| z1.2×z3 | 2.0e-3 | 1.7e-4 | 1.4e-3 | 3.1e-3 | 1.1e-3 |
| autos | unchanged | 4e-3 – 3e-2 | | | |

- Before the fix the i < j element carried the error: relative to CAMB the
  δ×PNG cross was 6% (z0.6×z1.2), 38% (z0.6×z3) and 10% (z1.2×z3) low at
  ℓ ≤ 3, while the j > i element was within ~1%. Blast's A_ij − A_ji was
  1.1% of max |A|; after the fix it is 6e-15.
- PNG×PNG (B) does not change and agrees with CAMB to < 1% below ℓ ≈ 100.
- SwiftCℓ is ~13–16% high for the narrow z = 3 bin in every piece, density
  only included. The paper run (different inputs) shows the same offset, so
  it is a SwiftCℓ limitation there, not an fNL one. Its PNG×PNG also departs
  at ℓ ≳ 50, where that piece is negligible.
- Fixed Blast oscillates around CAMB at ℓ ≳ 80 in both A and B (visible in
  `png_comparison.png`), unrelated to the ordering fix.

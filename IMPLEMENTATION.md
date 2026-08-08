# Reactant Non-Limber Implementation

This document serves as the engineering record for the Reactant-based non-Limber 3x2pt pipeline in `Blast.jl`.

## Scope

The Reactant forward implementation now covers the complete GC/WL calculation path: dynamic non-Limber preparation, W integration, non-Limber projections, dynamic linear/nonlinear Limber preparation, Limber corrections/high-ell contractions, and final requested-multipole interpolation. Reverse-mode validation remains complete only for the staged non-Limber path; complete-plan reverse mode is a separate follow-up.

The ordinary host `full_power_to_cls` pipeline remains the independent numerical oracle.

## Complete staged forward plan

`build_reactant_full_plan(workload)` constructs eight compiled stages:

```text
1. non-Limber preparation
2. W construction
3. non-Limber GG/GS/SS projection
4. linear/nonlinear Limber power products
5. paired Limber Chebyshev coefficients
6. paired Limber grid evaluation
7. six GG/GS/SS correction/high-ell contractions
8. paired final GG/GS/SS interpolation
```

The callable plan accepts three independently dynamic device arrays:

```julia
result = plan(pk_nonlimber_r, pk_limber_linear_r, pk_limber_nonlinear_r)
```

and returns `(; cl_gg, cl_gs, cl_ss)`. No host materialization, `to_rarray`, `@allowscalar`, or explicit GC occurs inside plan invocation.

Production-shape focused validation under local AbstractCosmologicalEmulators 0.11.0 and Reactant 0.2.278:

```text
Stage 4 power products:       7/7 passed
Stages 5-6 coefficients/grid: 12/12 passed
Stage 7 Limber terms:         25/25 passed
Stage 8 finalization:         10/10 passed
Complete plan + fixtures:     31/31 passed
```

The complete-plan gate reuses one compiled plan for baseline and three independent perturbations. Maximum host-parity errors were approximately:

```text
baseline:           GG 6.2e-16, GS 8.4e-18, SS 1.4e-19
non-Limber change:  GG 6.1e-16, GS 3.3e-18, SS 1.5e-19
linear Limber:      GG 6.1e-16, GS 5.4e-18, SS 1.7e-19
nonlinear Limber:   GG 5.9e-16, GS 3.3e-18, SS 2.5e-19
```

All three input perturbations produced nonzero changes in final spectra, proving that no public input was constant-folded. Ordinary host-generated text fixtures live under `benchmark/reference_fixtures/reactant_full/` and are regenerated only via `BLAST_REGENERATE_FULL_FIXTURES=1`.

Complete plan construction took 118-126 seconds on the CPU backend and peaked near 8.0 GB RSS. A 20+20-call memory diagnostic reproduced the known CPU PJRT lifetime behavior: without GC RSS grew by about 32 MB/call; with a full GC between calls RSS remained effectively flat. Production code does not force GC; long BenchmarkTools runs use `gcsample=true`.

The final equivalent long benchmark used `samples=1000`, `evals=1`,
`seconds=120`, `gcsample=true`, explicit synchronization, and separate
compilation accounting. It completed in 649 seconds with 8.05 GB peak RSS:

| Variant | Samples | Minimum ms | Median ms |
|---|---:|---:|---:|
| Plain Julia `full_power_to_cls` | 233 | 128.91 | 141.65 |
| Reactant device-resident | 228 | 144.74 | 154.82 |
| Reactant plus final host transfer | 228 | 124.63 | 150.90 |
| Three-input upload plus Reactant | 227 | 144.38 | 155.45 |

On this CPU backend, Reactant device-resident execution is approximately 9%
slower by median than the equivalent plain-Julia calculation. The lower
transfer minimum is timing noise and is not interpreted as a transfer
speedup. No CPU speedup is claimed.

## Artifact and Grid Decision

The calculation relies on precomputed $T$-tilde artifacts. The non-Limber numerical core uses the `6df6fd6960011063726cd28d2e72ad5ca221865c` artifact hash (`T_tildes_npz`). The global science grid dimensions for this branch have been resolved and kept unchanged from ordinary Blast behavior:
- $\chi$ nodes (`nchi`): 128
- Retained $R$ nodes (`nR`): 52
- Wavenumbers (`nk`): 161
- Exact non-Limber multipoles: 22

## Files and API Changes

- `src/reactant_api.jl`: Added `build_reactant_nonlimber_plan`, `AbstractReactantNonLimberPlan`, and `reactant_host_nonlimber_reference`.
- `ext/ReactantExt.jl`: Added `ReactantNonLimberPlan` struct for staged execution. Implemented separate compiled stages (`prep_comp`, `w_comp`, `proj_comp`) using explicit `StableHLO.einsum` representations to avoid the monolithic compilation trap.
- `src/Blast.jl`: Exported new API functions.
- (Unintended host changes to `src/constants.jl` and `src/projected_matter.jl` from earlier iterations were reverted to preserve the correct 128-node baseline).

## Static and Dynamic Data Boundaries

- **Dynamic Data**: The input power spectrum $P(k)$ is the only dynamic Reactant argument. It is converted to a `TracedRArray` once and passed through the plan.
- **Static Plan Data**: T-tilde kernels, Chebyshev transforms, survey kernels, and quadrature weights are initialized as device-resident static arrays during plan construction.

## Materialization Audit

The architecture strictly ensures zero host-materialization between stages:
- **No `@allowscalar`**: The plan does not trigger scalar iteration fallbacks or data-axis loops.
- **No Internal `Array` / `to_rarray` boundaries**: The outputs of `prep_comp` (coefficients) are fed directly into `w_comp` (W integration), whose outputs are fed directly into `proj_comp` (projection). The arrays remain `TracedRArray` throughout.
- **Explicit Materialization**: `Array(res.gg)` is only called by the user at the final test/benchmark boundary.
- **No Forced CPU Backend**: `build_reactant_nonlimber_plan` no longer forces `Reactant.set_default_backend("cpu")`, ensuring future GPU readiness.

## Fixture Evidence

We enforce a strict **no constant-folding proof** by executing the same compiled plan instance on both a baseline $P(k)$ and a deterministically perturbed $P(k)$.

The plan achieves exact parity (machine-epsilon) against a dedicated non-Limber host oracle (`reactant_host_nonlimber_reference`). The final integration fixture gate was rerun with the documented artifact override on 2026-08-04:

```bash
julia -t 8 --project=/home/marcobonici/Desktop/work/CosmologicalEmulators/blast_play/Blast.jl/benchmark \
  benchmark/test_reactant_nonlimber_fixtures.jl
```

Result: `55/55` tests passed, exit code `0`, wall time `2:48.21`, max RSS `7,893,904 kB`; see `benchmark/final_gate_fixtures.log`.

**Latest staged-forward parity errors (`benchmark/final_gate_forward_benchmark.log`):**
- GG: $1.0825335175794135 \times 10^{-16}$
- GS: $3.8116482626443515 \times 10^{-21}$
- SS: $1.257314531080602 \times 10^{-22}$

The fixtures are permanently recorded in `benchmark/reference_fixtures/reactant_nonlimber/`. They were properly regenerated to represent exact 22-ell non-Limber quantities using:
```bash
BLAST_REGENERATE_FIXTURES=1 julia -t 8 --project=benchmark benchmark/generate_reactant_nonlimber_fixtures.jl
```

## AD Evidence and Directional Sweep

The staged design supports explicit reverse-mode gradients with Enzyme, avoiding monolithic reverse graph overhead. The non-Limber gradients were validated by injecting signed, non-uniform cotangents ($U_{gg}, U_{gs}, U_{ss}$) into the non-Limber outputs and pulling back to $P(k)$.

The resulting gradient was verified via a central finite differences epsilon sweep against the non-Limber host oracle. The final directional-gradient gate was rerun on 2026-08-04:

```bash
julia -t 8 --project=/home/marcobonici/Desktop/work/CosmologicalEmulators/blast_play/Blast.jl/benchmark \
  benchmark/test_reactant_nonlimber_gradient.jl
```

Result: `3/3` tests passed, exit code `0`, wall time `3:36.32`, max RSS `8,258,212 kB`; see `benchmark/final_gate_gradient.log`.

**Directional derivative sweep:**

| $\epsilon$ | finite difference | AD | relative error |
|---:|---:|---:|---:|
| $10^{-3}$ | $-1.3658256790101062 \times 10^{-9}$ | $-1.0120045107566327 \times 10^{-9}$ | $2.590529477450654 \times 10^{-1}$ |
| $10^{-4}$ | $-9.859168060948054 \times 10^{-10}$ | $-1.0120045107566327 \times 10^{-9}$ | $2.6460350914556535 \times 10^{-2}$ |
| $10^{-5}$ | $-1.0112163927392115 \times 10^{-9}$ | $-1.0120045107566327 \times 10^{-9}$ | $7.793762275612834 \times 10^{-4}$ |
| $10^{-6}$ | $-1.0121569381238427 \times 10^{-9}$ | $-1.0120045107566327 \times 10^{-9}$ | $1.5059657397844062 \times 10^{-4}$ |

Best directional derivative relative error: $1.5059657397844062 \times 10^{-4}$ at $\epsilon = 10^{-6}$.

## Benchmarks and CPU Performance

Benchmarks were performed using `BenchmarkTools` on an 8-thread CPU configuration with `Reactant.set_default_backend("cpu")`. Unlike previous claims, the benchmarks strictly compare **equivalent work** (raw $P(k) \to \text{non-Limber } C_\ell$) and correctly exclude compilation.

| Variant | Minimum ms | Median ms | Output transfer | Compilation included | Notes |
|---|---:|---:|---|---|---|
| Host non-Limber from pk | 163.42 | 174.68 | Host output | No | equivalent reference |
| Reactant device-resident | 157.50 | 170.18 | No | No | synchronized |
| Reactant plus output transfer | 154.26 | 163.93 | Yes | No | synchronized |
| Reactant upload plus run | 151.35 | 167.43 | No | No | optional deployment cost |
| Plan construction | N/A | 109081.90 | N/A | Yes | fresh process (~109 seconds) |

Final staged-forward benchmark command:

```bash
julia -t 8 --project=/home/marcobonici/Desktop/work/CosmologicalEmulators/blast_play/Blast.jl/benchmark \
  benchmark/benchmark_staged_reactant_forward.jl
```

Result: exit code `0`, wall time `3:03.37`, max RSS `8,469,128 kB`; see `benchmark/final_gate_forward_benchmark.log`.

**Conclusion:** On the CPU backend, the Reactant implementation performs roughly identically to the highly optimized host Tullio loops. In this final run, the host/device minimum-time ratio was `1.0264124384176079` and the host/transfer minimum-time ratio was `1.0655235814379471`. The real value is that this zero-copy, fully device-resident explicit graph unlocks direct execution on GPUs/TPUs, which the current Tullio path cannot achieve natively.

## Cached Reverse Benchmark Evidence

The staged reverse path is benchmarked outside Blast's core package in `benchmark/benchmark_reactant_nonlimber_gradient.jl`, because Enzyme is an optional AD dependency. The benchmark is deliberately **non-Limber reverse only**. It uses the same real workload shapes and deterministic signed nonuniform cotangents as the directional-gradient validation:

- input $P(k)$ shape: `(50, 161)`
- coefficient shapes: `(161, 128, 52)`, `(161, 128, 52)`, `(161,)`
- lifted $\phi$ coefficient shape for the W stage: `(161, 128)`
- W tensor shapes: 17 tensors, each `(22, 128, 52)`
- non-Limber output cotangent shapes: GG `(22, 10, 10)`, GS `(22, 10, 5)`, SS `(22, 5, 5)`
- cotangents are signed and nonuniform; extrema were `((-1.0, 0.9432651586477729), (-1.0, 0.9178057468954721), (-0.602062469965807, 1.0))`

The timed reverse benchmark has explicit cache accounting:

- All static tensors and cotangents are uploaded once before cache construction/timing.
- Every reverse pullback is compiled before timing.
- One complete cached reverse pass is warmed before BenchmarkTools timing.
- Every timed result is synchronized before return.
- The device-resident complete reverse returns a synchronized Reactant gradient; the `to_host` variant additionally calls `Array(dpk)` inside the timed function.
- The compile counter was `19` before warm-up, after warm-up, before timing, and after all BenchmarkTools trials. Thus no timed call went through the benchmark's explicit `Reactant.@compile` cache path.

A reduced smoke validation was first run with `BLAST_REVERSE_SMOKE=1`, compiling one endpoint pullback plus the downstream W-to-coefficients and coefficient-to-pk pullbacks. It produced a nonzero `(50, 161)` gradient and the second call left the compile counter fixed at `3`; see `benchmark/reverse_smoke.log`.

The durable full benchmark was launched by `benchmark/run_full_nonlimber_reverse_benchmark_long.sh` with a six-hour timeout and the artifact override
`6df6fd6960011063726cd28d2e72ad5ca221865c = "/home/marcobonici/Desktop/work/CosmologicalEmulators/blast_play/for_marco/T_tildes_npz"`. It completed successfully on 2026-08-04 with exit code `0`, process wall time `3:34.84`, and max RSS `13,259,640 kB`; see `benchmark/full_nonlimber_reverse_benchmark_long.log` and `.status`.

| Reverse stage | Minimum ms | Median ms | Timed synchronization | Timed host materialization | Compile included |
|---|---:|---:|---|---|---|
| Endpoint-to-W cotangents, all 17 endpoint pullbacks | 21.80 | 24.72 | Yes | No | No |
| W-to-coefficients pullback | 144.67 | 161.16 | Yes | No | No |
| Coefficients-to-pk pullback | 20.40 | 23.31 | Yes | No | No |
| Complete cached reverse, device-resident gradient | 204.53 | 222.00 | Yes | No | No |
| Complete cached reverse plus final `Array(dpk)` | 196.45 | 209.24 | Yes | Yes | No |

Fresh reverse cache construction/compilation took `108.145488224` seconds and compiled exactly 19 pullbacks: 17 endpoint-to-W, one W-to-coefficients, and one coefficient-to-pk. The complete reverse warm-up returned a nonzero `(50, 161)` gradient with max absolute value `1.3912126914539274e-8`.

No equivalent host reverse benchmark exists in this run. Therefore this document does **not** claim a host-vs-Reactant reverse speedup or a compile-amortization break-even point. The honest conclusion is that cached staged Reactant reverse execution is now validated, compile-excluded, synchronized, and measured; host-equivalent reverse timing remains future work.

## Complete staged forward-plan reverse

The complete GC/WL forward calculation has now been differentiated with
Enzyme by reversing the same explicit stage boundaries rather than compiling
one monolithic reverse graph. The production driver is
`benchmark/run_staged_reactant_gradient.jl` and returns gradients with respect
to all three dynamic spectrum arrays:

```text
dpk_nonlimber:        (50, 161), maxabs 1.0746e-6
dpk_limber_linear:    (50, 257), maxabs 2.8630e-3
dpk_limber_nonlinear: (50, 257), maxabs 6.2934e-5
```

The scalar loss was the sum of all final GG, GS, and SS outputs. A bounded
relative random direction was used so every finite-difference point remained
inside the positive power-spectrum domain. Central finite differences showed
the expected convergence:

| epsilon | FD directional | Enzyme directional | relative error |
|---:|---:|---:|---:|
| 1e-3 | 3.6882561e-3 | 3.6844959e-3 | 1.0195e-3 |
| 1e-4 | 3.6845372e-3 | 3.6844959e-3 | 1.1223e-5 |
| 1e-5 | 3.6844963e-3 | 3.6844959e-3 | 1.0985e-7 |
| 1e-6 | 3.6844960e-3 | 3.6844959e-3 | 2.0111e-8 |

Cached reverse timings, excluding pullback compilation:

| Reverse calculation | Minimum ms | Median ms |
|---|---:|---:|
| non-Limber branch | 611.84 | 611.96 |
| Limber branch | 58.54 | 76.25 |
| complete finalization + both branches | 644.66 | 729.12 |

The durable run completed in 219 seconds with exit code 0 and 13.30 GB peak
RSS. This establishes complete forward-to-three-input reverse correctness.
The reverse driver remains benchmark-side because Enzyme is optional; a
public `ReactantFullGradientPlan` API has not yet been added.

## Known Failures and Limitations

- `Pkg.test()` completed with 907 passes, 2 failures, and 2 errors. The failures are the existing finite-difference tolerance failures for galaxy-bias derivatives; the errors are the pre-existing missing `test_type_stability_3x2.jl` file and a ForwardDiff `NumberCounts{Dual}` conversion failure. The full-plan patch does not alter those host/AD paths.
- The monolithic raw-$P(k,z) \to C_\ell$ pipeline approach was successfully rejected and replaced by explicitly staged `StableHLO.einsum` compilation blocks, preventing graph unrolling explosions.
- The cached reverse benchmark currently lives in `benchmark/benchmark_reactant_nonlimber_gradient.jl` as an external staged Enzyme/Reactant driver rather than as a public `ReactantNonLimberPlan` reverse API. Only the forward path is wrapped into `ReactantNonLimberPlan`.
- Complete-plan reverse mode is not implemented yet. The full staged forward plan is complete; reverse mode remains the next milestone.

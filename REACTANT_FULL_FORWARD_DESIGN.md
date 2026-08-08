# Complete Reactant Forward Plan Design

## Objective

Implement a production-shape staged Reactant plan equivalent to:

```julia
full_power_to_cls(workload; pass_plans=true)
```

The plan accepts three dynamic arrays:

```text
pk_nonlimber
pk_limber_linear
pk_limber_nonlinear
```

and returns complete requested-multipole `cl_gg`, `cl_gs`, and `cl_ss`.

## Approved decisions

1. Keep Akima in plain-Julia `transform_to_R_frame`; the B-spline experiment
   has been reverted.
2. All three power spectra are dynamic Reactant inputs.
3. Complete and validate forward mode before starting reverse mode.
4. Work directly on branch `reactant`; no worktree, commit, or push.
5. Use local AbstractCosmologicalEmulators 0.11.0 from
   `/home/marcobonici/Desktop/work/CosmologicalEmulators/blast_play/AbstractCosmologicalEmulators.jl`.

## Existing validated core

Preserve the production-shape non-Limber stages:

```text
pk_nonlimber -> coefficients -> 17 W tensors -> non-Limber GG/GS/SS
```

Axis contract:

```text
pk_nonlimber:  (50, 161)
coefficients:  (161, 128, 52), (161, 128, 52), (161,)
T HLO layout:  (161, 22, 128, 52)
W tensors:     17 × (22, 128, 52)
GG non-Limber: (22, 10, 10)
GS non-Limber: (22, 10, 5)
SS non-Limber: (22, 5, 5)
```

Do not redesign `reactant_chebyshev_matmul`, `reactant_w_ell_hlo`,
`_w_hlo_code`, W/endpoint StableHLO equations, science grids, or T-tilde
layout.

## Architecture

Do not compile one monolithic closure around `reactant_full_3x2pt`. The earlier
monolithic attempt produced an impractical graph. Use eight explicit compiled
stages with device-resident boundaries:

```text
1. non-Limber preparation             [existing]
2. W construction                     [existing]
3. non-Limber GG/GS/SS projection     [existing]
4. Limber power products
5. paired Limber Chebyshev coefficients
6. paired Limber grid evaluation
7. Limber correction/high-ell contractions
8. final GG/GS/SS interpolation
```

### Stage 4: Limber power products

Compile a pair-output closure around:

```julia
reactant_limber_power_products(pk_limber_linear,
                               pk_limber_nonlinear,
                               background)
```

The background is captured as static configuration. Both spectra remain
explicit dynamic arguments. Outputs are `logP_linear` and `logP_nonlinear`.

### Stage 5: paired coefficients

Compile one pair-output function:

```julia
(loglin, lognon, Mk, Mz) -> (
    reactant_limber_chebyshev_coefficients(loglin, Mk, Mz),
    reactant_limber_chebyshev_coefficients(lognon, Mk, Mz),
)
```

Outputs are `c_linear` and `c_nonlinear`.

### Stage 6: paired grid evaluation

Compile:

```julia
(clin, cnon, Tz, Tk) -> (
    reactant_limber_grid_from_coefficients(clin, Tz, Tk),
    reactant_limber_grid_from_coefficients(cnon, Tz, Tk),
)
```

Outputs are `P_linear(ell,chi)` and `P_nonlinear(ell,chi)`.

### Stage 7: Limber terms

Use device-resident static arrays:

```text
KG_low, KG_high
KL_low, KL_high
inv_chi2
quadrature weights
delta_chi
```

Compute:

```text
low  = (P_nonlinear - P_linear)[low ell, :] / chi²
high = P_nonlinear[high ell, :] / chi²
```

Return in this exact order:

```text
Cgg_correction, Cgg_limber
Cgs_correction, Cgs_limber
Css_correction, Css_limber
```

Refine `reactant_limber_terms_from_prepared` so its signature contains only
the four actual probe-kernel slices rather than duplicated aliases.

### Stage 8: finalization

Compile one function accepting non-Limber GG/GS/SS, all six Limber arrays,
and static final-interpolation arrays. Call `reactant_finalize_c_ell`
independently for GG, GS, and SS. Return:

```julia
(; cl_gg, cl_gs, cl_ss)
```

Do not call `reactant_full_3x2pt` here because that recomputes the non-Limber
endpoint already produced by Stage 3.

## API

Add to `src/reactant_api.jl`:

```julia
abstract type AbstractReactantFullPlan end
function build_reactant_full_plan end
```

Export `build_reactant_full_plan` from `src/Blast.jl`.

Add extension-owned `ReactantFullPlan` with concretely typed compiled handles
and static-data groups. Its public call contract is:

```julia
plan = build_reactant_full_plan(workload)
result = plan(pk_nonlimber_r, pk_limber_linear_r, pk_limber_nonlinear_r)
```

No `Array`, `collect`, scalar indexing, `@allowscalar`, `to_rarray`, or GC call
is allowed inside plan invocation.

## Static and dynamic data

Dynamic:

```text
pk_nonlimber, pk_limber_linear, pk_limber_nonlinear
```

Static for one plan:

```text
background object
probe kernels and prefactors
T-tilde tensors
quadrature data
Chebyshev transforms/evaluation matrices
requested ell grid
shape metadata
```

Dynamic background parameters and probe definitions are deferred.

## Plain-Julia oracle

Expected outputs must come from ordinary code, not the Reactant helpers:

```julia
full_power_to_cls(merge(workload, (
    pk=pk_nonlimber,
    pk_limber_lin=pk_limber_linear,
    pk_limber_nonlin=pk_limber_nonlinear,
)); pass_plans=true)
```

## Text fixtures

Create host-generated fixtures under:

```text
benchmark/reference_fixtures/reactant_full/
```

Required cases:

```text
baseline GG/GS/SS
non-Limber-only perturbation GG/GS/SS
linear-Limber-only perturbation GG/GS/SS
nonlinear-Limber-only perturbation GG/GS/SS
metadata.txt
```

Use `.txt`, not binary files. Fixtures are never regenerated during tests.
A separate generator requires an explicit environment variable. Perturbations
must be deterministic, nonuniform, shape-preserving, and physically safe.

## Validation sequence

### Gate 0: dependency baseline

Verify branch/HEAD, local ACE 0.11 path, restored Akima, and ordinary host
execution. Rerun the existing non-Limber fixture test under ACE 0.11 before
adding the full plan. Stop on compatibility regression.

### Gate 1: plain-Julia production-shape intermediates

Compare every helper against ordinary Blast intermediates:

```text
power products
Chebyshev coefficients
P_linear/P_nonlinear grids
GG/GS/SS correction arrays
GG/GS/SS high-ell arrays
final interpolation
```

### Gate 2: independent compiled stages

Compile Stages 4-8 separately. For each stage:

```text
assert all shapes
compare baseline output to plain Julia
reuse the same executable on perturbed dynamic input
compare perturbed output to plain Julia
assert a nonzero output change
```

Do not move to the next stage until the current stage passes.

### Gate 3: complete baseline plan

Build one plan and compare baseline output against both the ordinary host
oracle and saved text fixtures. Print maximum absolute and guarded relative
errors for GG/GS/SS.

### Gate 4: independent dynamic-input proof

Reuse the same plan for:

```text
baseline
only pk_nonlimber perturbed
only pk_limber_linear perturbed
only pk_limber_nonlinear perturbed
```

Every case must match its host fixture. Each independent perturbation must
produce a nonzero final-output change. This is the constant-folding proof for
all three public inputs.

### Gate 5: memory smoke

Run at least 20 synchronized calls outside BenchmarkTools while discarding
outputs; record RSS every five calls. Repeat with `GC.gc(true)` between calls
as a diagnostic only. Production plan code must not call GC.

### Gate 6: benchmark

After correctness passes, benchmark equivalent:

```text
ordinary host full_power_to_cls
device-resident full Reactant plan
Reactant plus final Array materialization
three-input upload plus plan
plan construction and per-stage compilation
```

Use BenchmarkTools interpolation, `evals=1`, explicit synchronization, and
`gcsample=true` on CPU PJRT. Exclude compilation from steady-state timing.

## Expensive compile procedure

Prove reduced examples first, then production shapes. Run real plan
construction detached with:

```text
setsid + nohup
timeout >= 6 hours
/usr/bin/time -v
durable log, status, and PID files
```

Record per-stage compile time, total construction time, peak RSS, CPU use,
exit status, and numerical result. Do not call a multi-minute compile a
failure. If failure occurs, preserve the first failing stage.

## Acceptance criteria

```text
[ ] Existing non-Limber fixtures pass under ACE 0.11.
[ ] Every Limber intermediate matches ordinary host code.
[ ] Stages 4-8 compile independently at production shapes.
[ ] Every stage passes a dynamic perturbed-input test.
[ ] One plan maps all three spectra to complete GG/GS/SS.
[ ] Baseline and three perturbation cases match text fixtures and host oracle.
[ ] No host materialization or scalar fallback occurs between stages.
[ ] Memory behavior is measured and bounded with documented GC protocol.
[ ] Compilation and steady-state runtime are reported separately.
[ ] No reverse-completion claim is made in this milestone.
```

## Reverse follow-up

After forward acceptance, compile and validate pullbacks in reverse stage
order. Each reduced pullback is checked against ForwardDiff before assembling
the production-shape reverse driver.

## Out of scope

```text
dynamic background/probes
CMB expansion beyond current GC/WL workload
interpolation-method changes
science-grid or artifact changes
monolithic compilation
GPU speed claims before CPU correctness
public reverse API during the forward milestone
```

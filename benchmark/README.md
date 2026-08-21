# Blast benchmarks

This directory contains the ordinary-Julia performance gates for the
`fix_marco_ad` repair. Reactant is deliberately outside this benchmark suite.

Set up the local benchmark environment once:

```bash
julia --project=benchmark -e 'using Pkg; Pkg.develop(path=pwd()); Pkg.instantiate()'
```

Run the deterministic controlled pipeline benchmark:

```bash
julia -t 8 --project=benchmark benchmark/benchmark_controlled_pipeline.jl
```

Run the realistic 10-bin GC × 5-bin WL primal/Mooncake benchmark:

```bash
BLAST_REALISTIC_DATA_ROOT=/path/to/data/root \
julia -t 8 --project=benchmark benchmark/benchmark_mooncake_realistic.jl
```

The realistic data root must contain `data/lsst_nz.npz`, `data/ell_n5k.npz`,
`data/pk_nonlimber_161.npz`, `data/pk_limber_lin.npz`, and
`data/pk_limber_nonlin.npz`.

Both scripts use BenchmarkTools, warm every measured path, record the complete
environment, and write raw timing samples plus a summary TSV under
`benchmark/results/`.

Permanent text references live in `benchmark/reference_results/`:

- `marco_ad_2cdb423/` is the frozen pre-change baseline.
- `fix_marco_ad_corrected/` is the final repaired baseline, including a direct
  stage-by-stage comparison with the frozen run.

`allocate_compute_w` returns mutable state. Benchmarks and inference code must
use one projected-matter workspace per concurrent task or chain; sharing one
workspace between concurrent evaluations races and corrupts results.

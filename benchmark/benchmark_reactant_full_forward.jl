using BenchmarkTools
using Blast
using Reactant
using NPZ
using Printf

const ROOT = "/home/marcobonici/Desktop/work/CosmologicalEmulators/blast_play/for_marco_copy"
include(joinpath(ROOT, "performance", "workload.jl"))

function run_synced(plan, args)
    result = plan(args...)
    foreach(Reactant.synchronize, values(result))
    return result
end

function run_to_host(plan, args)
    result = run_synced(plan, args)
    return map(Array, result)
end

function upload_and_run(plan, pk, lin, non)
    return run_synced(plan, map(Reactant.to_rarray, (pk, lin, non)))
end

function report(label, trial, elapsed)
    @printf("%s_samples=%d\n", label, length(trial.times))
    @printf("%s_elapsed_seconds=%.6f\n", label, elapsed)
    @printf("%s_min_ms=%.9f\n", label, minimum(trial).time / 1e6)
    @printf("%s_median_ms=%.9f\n", label, median(trial).time / 1e6)
    @printf("%s_allocs=%d\n", label, median(trial).allocs)
    @printf("%s_bytes=%d\n", label, median(trial).memory)
end

function main()
    Reactant.set_default_backend("cpu")
    workload = build_workload(ROOT, 1)
    t0 = time_ns()
    plan = Blast.build_reactant_full_plan(workload)
    println("plan_build_seconds=", (time_ns() - t0) / 1e9)
    args = map(Reactant.to_rarray,
        (workload.pk, workload.pk_limber_lin, workload.pk_limber_nonlin))
    run_synced(plan, args)
    GC.gc(true)

    t0 = time_ns()
    host = @benchmark full_power_to_cls($workload; pass_plans=true) samples=1000 evals=1 seconds=120 gcsample=true
    report("host", host, (time_ns() - t0) / 1e9)
    t0 = time_ns()
    device = @benchmark run_synced($plan, $args) samples=1000 evals=1 seconds=120 gcsample=true
    report("reactant_device", device, (time_ns() - t0) / 1e9)
    t0 = time_ns()
    transfer = @benchmark run_to_host($plan, $args) samples=1000 evals=1 seconds=120 gcsample=true
    report("reactant_transfer", transfer, (time_ns() - t0) / 1e9)
    t0 = time_ns()
    upload = @benchmark upload_and_run($plan, $(workload.pk), $(workload.pk_limber_lin), $(workload.pk_limber_nonlin)) samples=1000 evals=1 seconds=120 gcsample=true
    report("reactant_upload", upload, (time_ns() - t0) / 1e9)
end

main()

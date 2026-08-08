using AbstractCosmologicalEmulators
using BenchmarkTools
using Blast
using Reactant
using Printf

const ROOT = "/home/marcobonici/Desktop/work/CosmologicalEmulators/blast_play/for_marco_copy"
include(joinpath(ROOT, "performance", "workload.jl"))

# Equivalent host function: pk → coefficients → W → non-Limber output
# This performs the same work as the Reactant plan.
function host_nonlimber_from_pk(workload)
    spec = prepare_spectrum(workload)
    return Blast.reactant_host_nonlimber_reference(workload, spec)
end

function run_plan_synced(plan, pk_r)
    out = plan(pk_r)
    foreach(Reactant.synchronize, out)
    return out
end

function run_plan_to_host(plan, pk_r)
    out = run_plan_synced(plan, pk_r)
    return map(Array, out)
end

function upload_and_run(plan, pk)
    return run_plan_synced(plan, Reactant.to_rarray(pk))
end

function print_trial(label, trial, elapsed_seconds)
    @printf("%s_trial_count=%d\n", label, length(trial.times))
    @printf("%s_benchmark_elapsed_seconds=%.6f\n", label, elapsed_seconds)
    @printf("%s_min_ms=%.6f\n", label, minimum(trial).time / 1e6)
    @printf("%s_median_ms=%.6f\n", label, median(trial).time / 1e6)
    @printf("%s_min_allocs=%d\n", label, minimum(trial).allocs)
    @printf("%s_median_allocs=%.1f\n", label, Float64(median(trial).allocs))
    @printf("%s_min_memory_bytes=%d\n", label, minimum(trial).memory)
    @printf("%s_median_memory_bytes=%.1f\n", label, Float64(median(trial).memory))
end

function timed_benchmark(thunk, label)
    println("Benchmarking $label: samples=1000 evals=1 seconds=120")
    t0 = time_ns()
    trial = thunk()
    elapsed_seconds = (time_ns() - t0) / 1e9
    print_trial(label, trial, elapsed_seconds)
    return trial, elapsed_seconds
end

function main()
    Reactant.set_default_backend("cpu")
    Blast.FFTW.set_num_threads(8)
    workload = build_workload(ROOT, 8)
    selected_variant = get(ENV, "BLAST_FORWARD_VARIANT", "all")
    println("forward_variant_selection=", selected_variant)

    # 1. Build the plan (compilation cost)
    println("compiling_staged_forward=true")
    t0 = time_ns()
    plan = Blast.build_reactant_nonlimber_plan(workload)
    compile_ms = (time_ns() - t0) / 1e6
    println("staged_compile_ms=", compile_ms)

    # 2. Parity check
    pk_r = Reactant.to_rarray(workload.pk)
    result = run_plan_to_host(plan, pk_r)

    host = host_nonlimber_from_pk(workload)
    println("GG_parity_error=", maximum(abs, result.gg .- host.gg))
    println("GS_parity_error=", maximum(abs, result.gs .- host.gs))
    println("SS_parity_error=", maximum(abs, result.ss .- host.ss))

    # 3-6. Benchmarks.  Each branch still uses the exact corrected BenchmarkTools
    # configuration, but BLAST_FORWARD_VARIANT allows running variants in fresh
    # Julia processes to avoid retaining memory across long 120-second sections.
    host_bench = nothing
    device_bench = nothing
    transfer_bench = nothing
    upload_bench = nothing
    host_elapsed = 0.0
    device_elapsed = 0.0
    transfer_elapsed = 0.0
    upload_elapsed = 0.0

    if selected_variant in ("all", "host")
        host_bench, host_elapsed = timed_benchmark("host_nonlimber") do
            @benchmark host_nonlimber_from_pk($workload) samples=1000 evals=1 seconds=120
        end
    end

    if selected_variant in ("all", "device")
        device_bench, device_elapsed = timed_benchmark("reactant_device") do
            @benchmark run_plan_synced($plan, $pk_r) samples=1000 evals=1 seconds=120
        end
    end

    if selected_variant in ("all", "transfer")
        transfer_bench, transfer_elapsed = timed_benchmark("reactant_transfer") do
            @benchmark run_plan_to_host($plan, $pk_r) samples=1000 evals=1 seconds=120
        end
    end

    if selected_variant in ("all", "upload")
        upload_bench, upload_elapsed = timed_benchmark("reactant_upload") do
            @benchmark upload_and_run($plan, $workload.pk) samples=1000 evals=1 seconds=120
        end
    end

    selected_variant in ("all", "host", "device", "transfer", "upload") ||
        error("unknown BLAST_FORWARD_VARIANT=$selected_variant")

    # Summary
    println("---")
    println("julia_version=", VERSION)
    println("reactant_version=", pkgversion(Reactant))
    println("backend=cpu")
    println("threads=8")
    if selected_variant == "all"
        host_med = median(host_bench).time / 1e6
        device_med = median(device_bench).time / 1e6
        transfer_med = median(transfer_bench).time / 1e6
        println("host_vs_device_ratio=", host_med / device_med)
        println("host_vs_transfer_ratio=", host_med / transfer_med)
        println("forward_benchmark_total_elapsed_seconds=", host_elapsed + device_elapsed + transfer_elapsed + upload_elapsed)
    else
        println("host_vs_device_ratio=not_computed_in_single_variant_mode")
        println("host_vs_transfer_ratio=not_computed_in_single_variant_mode")
        println("forward_benchmark_total_elapsed_seconds=", host_elapsed + device_elapsed + transfer_elapsed + upload_elapsed)
    end
end

main()

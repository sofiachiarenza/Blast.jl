using BenchmarkTools
using DelimitedFiles
using DifferentiationInterface
using ForwardDiff
using Mooncake
using Statistics

include(joinpath(@__DIR__, "..", "test", "characterization_workload.jl"))
include(joinpath(@__DIR__, "realistic_workload.jl"))
include(joinpath(@__DIR__, "benchmark_instrumentation.jl"))
include(joinpath(@__DIR__, "fd_reference_fixture.jl"))

function summarize(label, trial)
    times_ms = trial.times ./ 1.0e6
    fastest = minimum(trial)
    return (; label, samples=length(times_ms), minimum_ms=minimum(times_ms),
            median_ms=median(times_ms), p95_ms=quantile(times_ms, 0.95),
            memory_bytes=fastest.memory, allocations=fastest.allocs, times_ms)
end

function main()
    root = get(ENV, "BLAST_REALISTIC_DATA_ROOT", "")
    isempty(root) && error("BLAST_REALISTIC_DATA_ROOT is required")
    samples = parse(Int, get(ENV, "BLAST_BENCH_SAMPLES", "100"))
    seconds = parse(Float64, get(ENV, "BLAST_BENCH_SECONDS", "120"))
    fftw_threads = parse(Int, get(ENV, "BLAST_FFTW_THREADS", "8"))
    output_directory = get(
        ENV,
        "BLAST_BENCH_OUTPUT",
        joinpath(@__DIR__, "results", "realistic_prechange"),
    )
    mkpath(output_directory)

    progress("build workload (cold setup) start")
    build_start = time_ns()
    workload = build_realistic_workload(root, fftw_threads)
    cold_build_seconds = (time_ns() - build_start) / 1.0e9
    progress("build workload done: $(cold_build_seconds) s, with_png=$(workload.with_png), active_components=$(length(workload.active_fields))")
    x0 = [1.0, 1.0, 1.0]
    x1 = [1.07, 0.93, 1.11]
    objectives = make_realistic_objective(workload, x0)
    # Optional saved independent FD reference (fail closed; default = compute FD).
    fixture_path = get(ENV, "BLAST_BENCH_FD_REFERENCE", "")
    saved = isempty(fixture_path) ? nothing :
        load_fd_reference(fixture_path, fd_fixture_context(root, x0, x1; with_png=workload.with_png,
                                                           active_components=length(workload.active_fields)))
    isnothing(saved) || progress("saved FD reference validated: $fixture_path sha256=$(saved.sha256)")
    backend = AutoMooncake(; config=nothing)

    progress("warmup (functional/reusable primals) start")
    functional0 = objectives.functional(x0)
    functional1 = objectives.functional(x1)
    primal0 = objectives.reusable(x0)
    primal1 = objectives.reusable(x1)
    isapprox(primal0, functional0; rtol=1.0e-10, atol=1.0e-12) ||
        error("Reusable/functional primal mismatch at x0")
    isapprox(primal1, functional1; rtol=1.0e-10, atol=1.0e-12) ||
        error("Reusable/functional primal mismatch at x1")
    progress("warmup done; DifferentiationInterface prepare_gradient start")
    if !isnothing(saved)   # fresh primals must match the saved Float64 references
        for (name, fresh) in (("x0", (functional0, primal0)), ("x1", (functional1, primal1))), v in fresh
            isapprox(v, saved.points[name].float_reference; rtol=1.0e-10, atol=1.0e-12) ||
                error("fresh primal at $name disagrees with the saved Float64 reference")
        end
    end
    preparation_start = time_ns()
    preparation = prepare_gradient(objectives.reusable, backend, x0)
    cold_preparation_seconds = (time_ns() - preparation_start) / 1.0e9
    progress("prepare_gradient done: $(cold_preparation_seconds) s (cold compile+preparation); first prepared gradient start")
    first_gradient_start = time_ns()
    mooncake0 = gradient(objectives.reusable, preparation, backend, x0)
    cold_first_gradient_seconds = (time_ns() - first_gradient_start) / 1.0e9
    progress("first prepared gradient done: $(cold_first_gradient_seconds) s (cold); ForwardDiff validation start")
    mooncake1 = gradient(objectives.reusable, preparation, backend, x1)
    # The current Float64 reusable workspace cannot hold ForwardDiff.Dual
    # values. Use the mathematically equivalent functional path as the
    # independent ForwardDiff reference; the fixture above validates that both
    # primals agree before either gradient is trusted.
    forward0 = isnothing(saved) ? ForwardDiff.gradient(objectives.functional_fd, x0) : saved.points["x0"].gradient
    forward1 = isnothing(saved) ? ForwardDiff.gradient(objectives.functional_fd, x1) : saved.points["x1"].gradient
    progress("ForwardDiff reference done; validating")
    isapprox(mooncake0, forward0; rtol=1.0e-10, atol=1.0e-12) ||
        error("Mooncake/ForwardDiff mismatch at x0")
    isapprox(mooncake1, forward1; rtol=1.0e-10, atol=1.0e-12) ||
        error("Mooncake/ForwardDiff mismatch at x1")

    objectives.reusable(x0)
    gradient(objectives.reusable, preparation, backend, x0)
    prepare_gradient(objectives.reusable, backend, x0)
    progress("validation and warmup complete; timed sections follow (hot timings)")
    GC.gc()
    thermal_gate("primal")
    primal_trial = @benchmark $(objectives.reusable)($x0) samples=samples seconds=seconds evals=1
    thermal_log("after primal")
    GC.gc()
    thermal_gate("gradient")
    gradient_trial = @benchmark gradient($(objectives.reusable), $preparation, $backend, $x0) samples=samples seconds=seconds evals=1
    thermal_log("after gradient")
    GC.gc()
    thermal_gate("gradient_preparation")
    prepare_trial = @benchmark prepare_gradient($(objectives.reusable), $backend, $x0) samples=max(5, samples ÷ 20) seconds=seconds evals=1
    thermal_log("after gradient_preparation")
    rows = (summarize("primal", primal_trial),
            summarize("gradient", gradient_trial),
            summarize("gradient_preparation", prepare_trial))

    println("julia_version=", VERSION)
    println("blast_version=", Base.pkgversion(Blast))
    println("blast_source=", pathof(Blast))
    println("mooncake_version=", Base.pkgversion(Mooncake))
    println("julia_threads=", Threads.nthreads())
    println("fftw_threads=", fftw_threads)
    println("cpu_model=", Sys.cpu_info()[1].model)
    println("gc_bins=", size(workload.galaxy.δ.nz, 1))
    println("wl_bins=", size(workload.lensing.γ.nz, 1))
    println("Pmm_contract=true")
    println("fd_reference_mode=", isnothing(saved) ? "fresh" : "saved")
    isnothing(saved) || println("fd_reference_fixture=", fixture_path, " sha256=", saved.sha256, " origin=", saved.origin)
    println("forwarddiff_reference=zero-partial Dual background (same physical cosmology) + functional pipeline")
    println("workload_with_png=", workload.with_png)
    println("active_projected_components=", length(workload.active_fields))
    println("active_projected_fields=", join(string.(workload.active_fields), ","))
    println("gradient_inputs=3 independent Pmm amplitudes (non-limber, limber-linear, limber-nonlinear); no f_NL or cosmological-parameter derivative is tested")
    println("cold_workload_build_seconds=", cold_build_seconds)
    println("cold_first_prepared_gradient_seconds=", cold_first_gradient_seconds)
    println("thermal_max_c_requested=", thermal_threshold())
    println("x0_primal=", primal0)
    println("x1_primal=", primal1)
    println("x0_mooncake_gradient=", mooncake0)
    println("x0_forwarddiff_gradient=", forward0)
    println("x1_mooncake_gradient=", mooncake1)
    println("x1_forwarddiff_gradient=", forward1)
    println("cold_gradient_preparation_seconds=", cold_preparation_seconds)
    open(joinpath(output_directory, "summary.tsv"), "w") do io
        println(io, "stage\tsamples\tminimum_ms\tmedian_ms\tp95_ms\tmemory_bytes\tallocations")
        for row in rows
            println(io, join((row.label, row.samples, row.minimum_ms, row.median_ms,
                              row.p95_ms, row.memory_bytes, row.allocations), '\t'))
            writedlm(joinpath(output_directory, row.label * "_times_ms.txt"), row.times_ms)
        end
    end
    println("results_directory=", output_directory)
end

main()

using BenchmarkTools
using DelimitedFiles
using DifferentiationInterface
using ForwardDiff
using Mooncake
using Statistics

include(joinpath(@__DIR__, "..", "test", "characterization_workload.jl"))
include(joinpath(@__DIR__, "realistic_workload.jl"))

objective(result, normalization) =
    normalization[1] * sum(abs2, result.cl_gg) +
    normalization[2] * sum(abs2, result.cl_gl) +
    normalization[3] * sum(abs2, result.cl_ll)

function make_realistic_objective(workload, x0)
    spectrum = prepare_realistic_spectrum(workload, x0)
    weights = Blast.compute_w(workload.weights_template, spectrum)
    result = project_realistic(workload, spectrum, weights)
    normalization = (inv(sum(abs2, result.cl_gg)),
                     inv(sum(abs2, result.cl_gl)),
                     inv(sum(abs2, result.cl_ll)))
    functional = let workload=workload, normalization=normalization
        x -> begin
            spectrum_x = prepare_realistic_spectrum(workload, x)
            weights_x = Blast.compute_w(workload.weights_template, spectrum_x)
            return objective(project_realistic(workload, spectrum_x, weights_x), normalization)
        end
    end
    workspace = Blast.allocate_compute_w(workload.weights_template, spectrum)
    reusable = let workload=workload, workspace=workspace, normalization=normalization
        x -> begin
            spectrum_x = prepare_realistic_spectrum(workload, x)
            Blast.compute_w!(workspace, spectrum_x)
            return objective(project_realistic(workload, spectrum_x, workspace), normalization)
        end
    end
    return (; functional, reusable, workspace, normalization)
end

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

    workload = build_realistic_workload(root, fftw_threads)
    x0 = [1.0, 1.0, 1.0]
    x1 = [1.07, 0.93, 1.11]
    objectives = make_realistic_objective(workload, x0)
    backend = AutoMooncake(; config=nothing)

    functional0 = objectives.functional(x0)
    functional1 = objectives.functional(x1)
    primal0 = objectives.reusable(x0)
    primal1 = objectives.reusable(x1)
    isapprox(primal0, functional0; rtol=1.0e-10, atol=1.0e-12) ||
        error("Reusable/functional primal mismatch at x0")
    isapprox(primal1, functional1; rtol=1.0e-10, atol=1.0e-12) ||
        error("Reusable/functional primal mismatch at x1")
    preparation_start = time_ns()
    preparation = prepare_gradient(objectives.reusable, backend, x0)
    cold_preparation_seconds = (time_ns() - preparation_start) / 1.0e9
    mooncake0 = gradient(objectives.reusable, preparation, backend, x0)
    mooncake1 = gradient(objectives.reusable, preparation, backend, x1)
    # The current Float64 reusable workspace cannot hold ForwardDiff.Dual
    # values. Use the mathematically equivalent functional path as the
    # independent ForwardDiff reference; the fixture above validates that both
    # primals agree before either gradient is trusted.
    forward0 = ForwardDiff.gradient(objectives.functional, x0)
    forward1 = ForwardDiff.gradient(objectives.functional, x1)
    isapprox(mooncake0, forward0; rtol=1.0e-10, atol=1.0e-12) ||
        error("Mooncake/ForwardDiff mismatch at x0")
    isapprox(mooncake1, forward1; rtol=1.0e-10, atol=1.0e-12) ||
        error("Mooncake/ForwardDiff mismatch at x1")

    objectives.reusable(x0)
    gradient(objectives.reusable, preparation, backend, x0)
    prepare_gradient(objectives.reusable, backend, x0)
    GC.gc()
    primal_trial = @benchmark $(objectives.reusable)($x0) samples=samples seconds=seconds evals=1
    GC.gc()
    gradient_trial = @benchmark gradient($(objectives.reusable), $preparation, $backend, $x0) samples=samples seconds=seconds evals=1
    GC.gc()
    prepare_trial = @benchmark prepare_gradient($(objectives.reusable), $backend, $x0) samples=max(5, samples ÷ 20) seconds=seconds evals=1
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

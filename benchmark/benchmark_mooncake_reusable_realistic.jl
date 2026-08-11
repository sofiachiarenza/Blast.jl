using BenchmarkTools
using DelimitedFiles
using DifferentiationInterface
using ForwardDiff
using Mooncake
using Statistics

include(joinpath(@__DIR__, "realistic_3x2pt_workload.jl"))

objective(result, normalizations) =
    normalizations[1] * sum(abs2, result.cl_gg) +
    normalizations[2] * sum(abs2, result.cl_gs) +
    normalizations[3] * sum(abs2, result.cl_ss)

function make_objectives(workload, x0)
    spectrum0 = prepare_scaled_spectrum(workload, x0)
    W0 = Blast.compute_w(workload.W_template, spectrum0)
    result0 = project_all(workload, spectrum0, W0)
    normalizations = (
        inv(sum(abs2, result0.cl_gg)),
        inv(sum(abs2, result0.cl_gs)),
        inv(sum(abs2, result0.cl_ss)),
    )

    functional = let workload=workload, normalizations=normalizations
        x -> begin
            spectrum = prepare_scaled_spectrum(workload, x)
            W = Blast.compute_w(workload.W_template, spectrum)
            return objective(project_all(workload, spectrum, W), normalizations)
        end
    end

    workspace = Blast.allocate_compute_w(workload.W_template, spectrum0)
    reusable = let workload=workload, workspace=workspace, normalizations=normalizations
        x -> begin
            spectrum = prepare_scaled_spectrum(workload, x)
            Blast.compute_w!(workspace, spectrum)
            return objective(project_all(workload, spectrum, workspace), normalizations)
        end
    end
    return (; functional, reusable, workspace, normalizations)
end

function summarize(label, trial)
    trial_times = trial.times ./ 1.0e6
    fastest = minimum(trial)
    println(
        label,
        "_samples=", length(trial_times),
        " min_ms=", round(minimum(trial_times); digits=6),
        " median_ms=", round(median(trial_times); digits=6),
        " p95_ms=", round(quantile(trial_times, 0.95); digits=6),
        " memory_mib=", round(fastest.memory / 2^20; digits=6),
        " allocs=", fastest.allocs,
    )
    return trial_times
end

function assert_close(label, actual, expected; rtol=1.0e-10, atol=1.0e-12)
    max_abs = maximum(abs, actual .- expected)
    max_rel = max_abs / max(maximum(abs, expected), eps(Float64))
    println(label, "_max_abs=", max_abs)
    println(label, "_max_rel=", max_rel)
    isapprox(actual, expected; rtol, atol) || error("$label failed")
end

function main()
    data_root = get(ENV, "EXAMPLE_ROOT", "")
    isempty(data_root) && error("EXAMPLE_ROOT must point to the realistic fixture directory")
    results_directory = get(
        ENV,
        "BENCHMARK_RESULTS_DIR",
        joinpath(@__DIR__, "results", "mooncake_reusable_realistic"),
    )
    mkpath(results_directory)

    workload = build_realistic_no_png_workload(data_root, 8)
    x0 = [1.0, 1.0, 1.0]
    x1 = [1.07, 0.93, 1.11]
    objectives = make_objectives(workload, x0)
    backend = AutoMooncake(; config=nothing)

    println("julia_version=", VERSION)
    println("blast_version=", Base.pkgversion(Blast))
    println("blast_source=", pathof(Blast))
    println("mooncake_version=", Base.pkgversion(Mooncake))
    println("differentiation_interface_version=", Base.pkgversion(DifferentiationInterface))
    println("forwarddiff_version=", Base.pkgversion(ForwardDiff))
    println("julia_threads=", Threads.nthreads())
    println("fftw_threads=8")
    println("cpu_model=", Sys.cpu_info()[1].model)
    println("gc_bins=", size(workload.galaxy.δ.nz, 1))
    println("wl_bins=", size(workload.lensing.γ.nz, 1))
    println("PNG_active=", !isnothing(workload.galaxy.PNG))
    println("nR=", length(Blast.R))

    # Validate both the reference point and a distinct point before timing.
    for (label, x) in (("x0", x0), ("x1", x1))
        functional_value = objectives.functional(x)
        reusable_value = objectives.reusable(x)
        println(label, "_functional_value=", functional_value)
        println(label, "_reusable_value=", reusable_value)
        assert_close(label * "_primal", [reusable_value], [functional_value])
    end

    preparation_start = time_ns()
    preparation = prepare_gradient(objectives.reusable, backend, x0)
    initial_preparation_seconds = (time_ns() - preparation_start) / 1.0e9
    println("initial_preparation_seconds=", initial_preparation_seconds)

    if parse(Bool, get(ENV, "SKIP_REFERENCE_VALIDATION", "false"))
        println("reference_validation_skipped=true")
    else
        for (label, x) in (("x0", x0), ("x1", x1))
            mooncake_gradient = gradient(objectives.reusable, preparation, backend, x)
            forwarddiff_gradient = ForwardDiff.gradient(objectives.functional, x)
            println(label, "_mooncake_gradient=", mooncake_gradient)
            println(label, "_forwarddiff_gradient=", forwarddiff_gradient)
            assert_close(label * "_gradient", mooncake_gradient, forwarddiff_gradient)
        end
    end

    if parse(Bool, get(ENV, "VALIDATE_ONLY", "false"))
        println("validation_only=true")
        return
    end

    # Warm every BenchmarkTools expression after all correctness checks.
    objectives.reusable(x0)
    gradient(objectives.reusable, preparation, backend, x0)
    prepare_gradient(objectives.reusable, backend, x0)
    GC.gc()

    primal_trial = @benchmark $(objectives.reusable)($x0) samples=200 seconds=120 evals=1
    GC.gc()
    gradient_trial = @benchmark gradient(
        $(objectives.reusable), $preparation, $backend, $x0,
    ) samples=200 seconds=120 evals=1
    GC.gc()
    prepare_trial = @benchmark prepare_gradient(
        $(objectives.reusable), $backend, $x0,
    ) samples=10 seconds=120 evals=1

    primal_times = summarize("reusable_primal", primal_trial)
    gradient_times = summarize("reusable_gradient", gradient_trial)
    prepare_times = summarize("reusable_prepare", prepare_trial)
    writedlm(joinpath(results_directory, "reusable_primal_times_ms.txt"), primal_times)
    writedlm(joinpath(results_directory, "reusable_gradient_times_ms.txt"), gradient_times)
    writedlm(joinpath(results_directory, "reusable_prepare_times_ms.txt"), prepare_times)
    println("results_directory=", results_directory)
end

main()

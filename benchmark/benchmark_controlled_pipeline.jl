using BenchmarkTools
using DelimitedFiles
using Statistics

include(joinpath(@__DIR__, "..", "test", "characterization_workload.jl"))

function summarize(label, trial)
    times_ms = trial.times ./ 1.0e6
    fastest = minimum(trial)
    return (
        label=label,
        samples=length(times_ms),
        minimum_ms=minimum(times_ms),
        median_ms=median(times_ms),
        p95_ms=quantile(times_ms, 0.95),
        memory_bytes=fastest.memory,
        allocations=fastest.allocs,
        times_ms=times_ms,
    )
end

function write_results(output_directory, rows)
    mkpath(output_directory)
    open(joinpath(output_directory, "summary.tsv"), "w") do io
        println(io, "stage\tsamples\tminimum_ms\tmedian_ms\tp95_ms\tmemory_bytes\tallocations")
        for row in rows
            println(io, join((row.label, row.samples, row.minimum_ms, row.median_ms,
                              row.p95_ms, row.memory_bytes, row.allocations), '\t'))
            writedlm(joinpath(output_directory, row.label * "_times_ms.txt"), row.times_ms)
        end
    end
end

function main()
    samples = parse(Int, get(ENV, "BLAST_BENCH_SAMPLES", "30"))
    seconds = parse(Float64, get(ENV, "BLAST_BENCH_SECONDS", "30"))
    fftw_threads = parse(Int, get(ENV, "BLAST_FFTW_THREADS", "8"))
    output_directory = get(
        ENV,
        "BLAST_BENCH_OUTPUT",
        joinpath(@__DIR__, "results", "controlled_prechange"),
    )
    Blast.FFTW.set_num_threads(fftw_threads)

    state = build_characterization_state()
    spectrum = prepare_characterization_spectrum(state)
    weights = Blast.compute_w(state.weights_template, spectrum)
    workspace = Blast.allocate_compute_w(state.weights_template, spectrum)

    # Correctness and compilation warmup.
    reference = project_characterization_spectra(state, spectrum, weights)
    Blast.compute_w!(workspace, spectrum)
    reusable = project_characterization_spectra(state, spectrum, workspace)
    for field in fieldnames(typeof(reference))
        isapprox(getfield(reusable, field), getfield(reference, field);
                 rtol=1.0e-10, atol=1.0e-15) || error("$field reusable mismatch")
    end
    run_characterization(state)
    build_characterization_probes(state.bg)
    Blast.Background(state.cosmo)
    Blast.SetUp(state.probes.galaxy, state.probes.lensing, state.probes.cmb)
    GC.gc()

    trials = Pair{String, Any}[
        "background" => @benchmark(Blast.Background($(state.cosmo));
                                     samples=samples, seconds=seconds, evals=1),
        "prepare_components" => @benchmark(build_characterization_probes($(state.bg));
                                           samples=samples, seconds=seconds, evals=1),
        "setup" => @benchmark(Blast.SetUp($(state.probes.galaxy),
                                           $(state.probes.lensing), $(state.probes.cmb));
                              samples=max(5, samples ÷ 3), seconds=seconds, evals=1),
        "prepare_pmm" => @benchmark(prepare_characterization_spectrum($state);
                                    samples=samples, seconds=seconds, evals=1),
        "compute_w" => @benchmark(Blast.compute_w($(state.weights_template), $spectrum);
                                   samples=samples, seconds=seconds, evals=1),
        "compute_w_reusable" => @benchmark(Blast.compute_w!($workspace, $spectrum);
                                            samples=samples, seconds=seconds, evals=1),
        "project_all" => @benchmark(project_characterization_spectra(
                                        $state, $spectrum, $workspace);
                                     samples=samples, seconds=seconds, evals=1),
        "full_cached" => @benchmark(run_characterization($state);
                                      samples=samples, seconds=seconds, evals=1),
        "full_from_cosmology" => @benchmark(run_characterization();
                                             samples=max(5, samples ÷ 3),
                                             seconds=seconds, evals=1),
    ]

    rows = [summarize(label, trial) for (label, trial) in trials]
    println("julia_version=", VERSION)
    println("blast_version=", Base.pkgversion(Blast))
    println("blast_source=", pathof(Blast))
    println("julia_threads=", Threads.nthreads())
    println("fftw_threads=", fftw_threads)
    println("cpu_model=", Sys.cpu_info()[1].model)
    for row in rows
        println(row.label, " min_ms=", row.minimum_ms,
                " median_ms=", row.median_ms,
                " p95_ms=", row.p95_ms,
                " memory_mib=", row.memory_bytes / 2.0^20,
                " allocations=", row.allocations)
    end
    write_results(output_directory, rows)
    println("results_directory=", output_directory)
end

main()

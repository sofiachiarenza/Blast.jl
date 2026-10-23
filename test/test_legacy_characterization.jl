using DelimitedFiles
using LinearAlgebra
using TOML

isdefined(@__MODULE__, :build_characterization_state) ||
    include("characterization_workload.jl")

@testset "Legacy marco_ad characterization" begin
    fixture_directory = joinpath(@__DIR__, "reference_fixtures", "legacy_marco_ad")
    metadata = TOML.parsefile(joinpath(fixture_directory, "metadata.toml"))
    result = run_characterization()
    arrays = characterization_arrays(result)

    @test metadata["pmm_contract"] ==
          "All supplied power spectra are total-matter Pmm in Mpc^3."
    @test Set(keys(arrays)) == Set(keys(metadata["shapes"]))

    # Stage 3 intentionally replaced normalization by D(first positive χ) with
    # normalization by ACE's exact D(0). Freeze the size of that scientific
    # change rather than silently blessing arbitrary drift in IA observables.
    #
    # The mixed-PNG transfer-leg fix (issue #38, test_png_transfer_ordering.jl)
    # moves every spectrum with a PNG cross term at this fixture's f_NL = 7:
    # cl_gg, cl_gc (κ×PNG, ISW×PNG) and cl_gl, whose range now covers both
    # changes (growth alone was 2.4e-4–2.7e-4).
    intentional_deltas = Dict(
        "lensing_ia_kernel" => (3.0e-3, 3.2e-3),
        "cl_gl" => (4.0e-4, 4.3e-4),
        "cl_lc" => (1.6e-4, 1.8e-4),
        "cl_ll" => (8.5e-5, 1.0e-4),
        "cl_gg" => (3.0e-4, 3.2e-4),
        "cl_gc" => (4.3e-3, 4.5e-3),
    )

    for name in sort!(collect(keys(arrays)))
        actual = vec(arrays[name])
        expected = vec(readdlm(joinpath(fixture_directory, name * ".txt"), Float64))
        @test length(actual) == length(expected)
        scale = max(norm(expected), floatmin(Float64))
        peak = max(maximum(abs, expected), floatmin(Float64))
        relative_l2 = norm(actual .- expected) / scale
        peak_relative_max = maximum(abs, actual .- expected) / peak

        if haskey(intentional_deltas, name)
            lower, upper = intentional_deltas[name]
            @test lower ≤ relative_l2 ≤ upper
            continue
        end

        # The primordial-potential-only projected weight is a cancellation-
        # dominated threaded reduction. Different thread schedules move its
        # tiny residual at roughly 1e-4 relative while final spectra remain
        # stable at much tighter tolerances. Characterize that known numerical
        # floor explicitly rather than pretending the reduction is bitwise.
        # Registered ACE 0.11 evaluates the same neutrino Fermi-Dirac integral
        # through its precomputed Akima table, whereas the legacy Blast code
        # duplicated a direct quadrature. Their 2.7e-10 density difference can
        # enter twice in lensing auto-spectra, so 1e-9 is the appropriate
        # cross-version characterization floor.
        tolerance = name == "weight_w_2_00_ϕ_sample" ? 5.0e-4 : 1.0e-9
        @test relative_l2 ≤ tolerance
        @test peak_relative_max ≤ tolerance
    end
end

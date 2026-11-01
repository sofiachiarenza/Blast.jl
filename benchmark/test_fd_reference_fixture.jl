# Tiny tests of the FD-reference fixture guards (fake context, temp files; no Blast needed):
#   julia benchmark/test_fd_reference_fixture.jl
using Test
include(joinpath(@__DIR__, "fd_reference_fixture.jl"))

const GOOD_LOG = """
Test Summary: | Pass  Total
04:06:17 objective check: workload built (with_png=true, active=17)
04:19:18 FD point x0: primal=2.9999999999999996 reference=3.0 gradient=[5.990864920486477, -5.672654584450479, 5.664722832973471]
04:32:00 FD point x1: primal=4.614815193325615 reference=4.614815193325651 gradient=[7.439665378732863, -7.06459273247574, 7.043301981173633]
EXIT=0
"""
fake_context() = Dict{String,Any}(
    "x0" => [1.0, 1.0, 1.0], "x1" => [1.07, 0.93, 1.11], "with_png" => true, "active_components" => 17, "n_R" => 64,
    "datasets" => Dict{String,Any}("a.npz" => "aa", "b.npz" => "bb"),
    "axes" => Dict{String,Any}("chi" => "c1", "R" => "r1"),
    "harness" => Dict{String,Any}("realistic_workload.jl" => "h1"))

function with_fixture(f; log=GOOD_LOG)
    dir = mktempdir()
    logpath = joinpath(dir, "log.txt"); write(logpath, log)
    path = joinpath(dir, "ref.toml")
    write_fd_reference(path, logpath, fake_context(); origin=Dict("baseline" => "e2b47d6"))
    return f(path, dir)
end
rejects(path, ctx) = try; load_fd_reference(path, ctx); false; catch e; e isa ErrorException; end
mutate(f) = (c = fake_context(); f(c); c)

@testset "FD reference fixture" begin
    with_fixture() do path, dir
        ref = load_fd_reference(path, fake_context())
        @test ref.points["x0"].gradient == [5.990864920486477, -5.672654584450479, 5.664722832973471]   # exact round trip
        @test ref.points["x1"].primal_fd == 4.614815193325615 && ref.points["x0"].float_reference == 3.0
        @test ref.origin["baseline"] == "e2b47d6" && length(ref.origin["log_sha256"]) == 64 && length(ref.sha256) == 64
        @test_throws ErrorException write_fd_reference(path, joinpath(dir, "log.txt"), fake_context(); origin=Dict())   # no overwrite
        @test rejects(joinpath(dir, "missing.toml"), fake_context())
        @test rejects(path, mutate(c -> c["x1"] = [1.07, 0.93, 1.12]))                       # wrong x
        @test rejects(path, mutate(c -> c["with_png"] = false))                              # wrong PNG
        @test rejects(path, mutate(c -> c["active_components"] = 13))
        @test rejects(path, mutate(c -> c["n_R"] = 52))
        @test rejects(path, mutate(c -> c["datasets"]["a.npz"] = "zz"))                      # stale data
        @test rejects(path, mutate(c -> c["datasets"]["extra.npz"] = "zz"))
        @test rejects(path, mutate(c -> c["axes"]["R"] = "r2"))                              # wrong axes
        @test rejects(path, mutate(c -> c["harness"]["realistic_workload.jl"] = "h2"))       # stale code
        text = read(path, String)
        write(path, replace(text, "5.990864920486477" => "nan")); @test rejects(path, fake_context())   # non-finite
        write(path, replace(text, "[point.x1]" => "[point.x9]")); @test rejects(path, fake_context())     # missing entry
        write(path, replace(text, FD_FIXTURE_FORMAT => "other")); @test rejects(path, fake_context())
        write(path, replace(text, "gradient = [5.990864920486477, -5.672654584450479, 5.664722832973471]" => "gradient = [1.0, 2.0]"))
        @test rejects(path, fake_context())
    end
    # Bad logs are not exported.
    for bad in (replace(GOOD_LOG, "EXIT=0" => "EXIT=1"), replace(GOOD_LOG, "active=17" => "active=13"),
                replace(GOOD_LOG, r"04:32:00 FD point x1.*\n" => ""), GOOD_LOG * "x FD point x0: primal=1.0 reference=1.0 gradient=[1.0, 2.0, 3.0]\n")
        dir = mktempdir(); write(joinpath(dir, "l"), bad)
        @test_throws ErrorException write_fd_reference(joinpath(dir, "o.toml"), joinpath(dir, "l"), fake_context(); origin=Dict())
        @test !isfile(joinpath(dir, "o.toml"))
    end
end

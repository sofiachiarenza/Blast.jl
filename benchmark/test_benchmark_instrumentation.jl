# Plain-Julia check of the thermal helper with fake samplers (no packages needed):
#   julia benchmark/test_benchmark_instrumentation.jl
include(joinpath(@__DIR__, "benchmark_instrumentation.jl"))
using Test
const noop = _ -> nothing
@testset "thermal gate" begin
    @test isnothing(wait_for_thermal_stability("x"; threshold=nothing))                      # disabled
    seq(v) = (i = Ref(0); () -> v[min(i[] += 1, length(v))])
    # cools then stabilizes
    t = Ref(0.0)
    @test wait_for_thermal_stability("cool"; threshold=60.0, sampler=seq([80, 70, 62, 58, 57, 57, 56.5]),
                                     sleeper=noop, clock=() -> (t[] += 1.0)) == 57.0
    # below threshold but drifting > 2 C is not stable until it settles
    @test wait_for_thermal_stability("drift"; threshold=60.0, sampler=seq([50, 55, 50, 51, 51, 50.5]),
                                     sleeper=noop, clock=() -> (t[] += 1.0)) == 51.0
    # never cools: timeout
    @test_throws ErrorException wait_for_thermal_stability("hot"; threshold=60.0, sampler=() -> 90.0,
                                     sleeper=noop, clock=() -> (t[] += 100.0), timeout=600.0)
    # no sensor requested-but-missing
    @test_throws ErrorException wait_for_thermal_stability("nosensor"; threshold=60.0, sampler=() -> nothing,
                                     sleeper=noop)
    # sensor parser on a fake hwmon tree
    root = mktempdir()
    for (d, name, temps) in (("hwmon0", "coretemp", (45000, 52000)), ("hwmon1", "nvme", (99000,)))
        mkpath(joinpath(root, d)); write(joinpath(root, d, "name"), name * "\n")
        for (i, v) in enumerate(temps); write(joinpath(root, d, "temp$(i)_input"), "$v\n"); end
    end
    @test cpu_max_temp_c(root) == 52.0
    @test isnothing(cpu_max_temp_c(joinpath(root, "missing")))
end

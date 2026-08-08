using Blast
using DelimitedFiles
using Reactant
using Test
using NPZ

const ROOT = "/home/marcobonici/Desktop/work/CosmologicalEmulators/blast_play/for_marco_copy"
const FIXTURES = joinpath(@__DIR__, "reference_fixtures", "reactant_full")
include(joinpath(ROOT, "performance", "workload.jl"))

host_result(workload, pk, lin, non) = full_power_to_cls(
    merge(workload, (; pk, pk_limber_lin=lin, pk_limber_nonlin=non));
    pass_plans=true)

function run_plan(plan, pk, lin, non)
    result = plan(Reactant.to_rarray(pk), Reactant.to_rarray(lin), Reactant.to_rarray(non))
    foreach(Reactant.synchronize, values(result))
    return (; cl_gg=Array(result.cl_gg), cl_gs=Array(result.cl_gs), cl_ss=Array(result.cl_ss))
end

function fixture(label, name, shape)
    data = vec(readdlm(joinpath(FIXTURES, "$(label)_$(name).txt"), Float64))
    return reshape(data, shape)
end

function main()
    Reactant.set_default_backend("cpu")
    workload = build_workload(ROOT, 1)
    pk0, lin0, non0 = workload.pk, workload.pk_limber_lin, workload.pk_limber_nonlin
    pk1 = pk0 .* reshape(1 .+ 1e-5 .* sin.(1:size(pk0, 1)), :, 1)
    lin1 = lin0 .* reshape(1 .+ 1e-5 .* sin.(1:size(lin0, 1)), :, 1)
    non1 = non0 .* reshape(1 .+ 1e-5 .* cos.(1:size(non0, 1)), :, 1)
    cases = (
        ("baseline", pk0, lin0, non0),
        ("nonlimber", pk1, lin0, non0),
        ("limber_linear", pk0, lin1, non0),
        ("limber_nonlinear", pk0, lin0, non1),
    )
    println("building_complete_reactant_plan=true")
    t0 = time_ns()
    plan = Blast.build_reactant_full_plan(workload)
    println("complete_plan_build_seconds=", (time_ns() - t0) / 1e9)
    baseline = nothing
    @testset "Complete Reactant forward plan" begin
        for (label, pk, lin, non) in cases
            host = host_result(workload, pk, lin, non)
            result = run_plan(plan, pk, lin, non)
            @test map(size, values(result)) == map(size, (host.cl_gg, host.cl_gs, host.cl_ss))
            @test result.cl_gg ≈ host.cl_gg rtol=1e-9 atol=1e-15
            @test result.cl_gs ≈ host.cl_gs rtol=1e-9 atol=1e-15
            @test result.cl_ss ≈ host.cl_ss rtol=1e-9 atol=1e-15
            @test result.cl_gg ≈ fixture(label, "cl_gg", size(result.cl_gg)) rtol=1e-9 atol=1e-15
            @test result.cl_gs ≈ fixture(label, "cl_gs", size(result.cl_gs)) rtol=1e-9 atol=1e-15
            @test result.cl_ss ≈ fixture(label, "cl_ss", size(result.cl_ss)) rtol=1e-9 atol=1e-15
            println(label, "_maxabs=", (
                maximum(abs, result.cl_gg .- host.cl_gg),
                maximum(abs, result.cl_gs .- host.cl_gs),
                maximum(abs, result.cl_ss .- host.cl_ss)))
            if label == "baseline"
                baseline = result
            else
                changes = (
                    maximum(abs, result.cl_gg .- baseline.cl_gg),
                    maximum(abs, result.cl_gs .- baseline.cl_gs),
                    maximum(abs, result.cl_ss .- baseline.cl_ss))
                @test maximum(changes) > 0
                println(label, "_dynamic_changes=", changes)
            end
        end
    end
end

main()

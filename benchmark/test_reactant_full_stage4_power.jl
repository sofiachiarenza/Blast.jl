using Blast
using Reactant
using Test
using NPZ

const ROOT = "/home/marcobonici/Desktop/work/CosmologicalEmulators/blast_play/for_marco_copy"
include(joinpath(ROOT, "performance", "workload.jl"))

function host_power_products(pk_lin, pk_nonlin, bg)
    k = 10 .^ Blast.k_limber
    Pphi = reshape(Blast.get_PΦ(k, bg), 1, :)
    Plin = permutedims(Pphi .* Blast.get_Tm(pk_lin, k, bg) .^ 2)
    Pnon = permutedims(Pphi .* Blast.get_Tm(pk_nonlin, k, bg) .^ 2)
    return log10.(Plin), log10.(Pnon)
end

function main()
    Reactant.set_default_backend("cpu")
    workload = build_workload(ROOT, 1)
    bg = workload.bg
    lin0, non0 = workload.pk_limber_lin, workload.pk_limber_nonlin
    lin1 = lin0 .* reshape(1 .+ 1e-4 .* sin.(1:size(lin0, 1)), :, 1)
    non1 = non0 .* reshape(1 .+ 1e-4 .* cos.(1:size(non0, 1)), :, 1)
    reference0 = host_power_products(lin0, non0, bg)
    reference1 = host_power_products(lin1, non1, bg)
    args0 = map(Reactant.to_rarray, (lin0, non0))
    args1 = map(Reactant.to_rarray, (lin1, non1))
    f = let bg=bg
        (lin, non) -> Blast.reactant_limber_power_products(lin, non, bg)
    end
    compiled = Reactant.compile(f, args0; sync=true)
    result0 = map(Array, compiled(args0...))
    result1 = map(Array, compiled(args1...))
    @testset "Reactant full Stage 4" begin
        @test map(size, result0) == map(size, reference0)
        @test result0[1] ≈ reference0[1] rtol=1e-11 atol=1e-12
        @test result0[2] ≈ reference0[2] rtol=1e-11 atol=1e-12
        @test result1[1] ≈ reference1[1] rtol=1e-11 atol=1e-12
        @test result1[2] ≈ reference1[2] rtol=1e-11 atol=1e-12
        @test maximum(abs, result1[1] .- result0[1]) > 0
        @test maximum(abs, result1[2] .- result0[2]) > 0
    end
    println("stage4_shapes=", map(size, result0))
    println("stage4_baseline_maxabs=", map((a,b) -> maximum(abs, a .- b), result0, reference0))
    println("stage4_perturbed_maxabs=", map((a,b) -> maximum(abs, a .- b), result1, reference1))
    println("stage4_dynamic_changes=", map((a,b) -> maximum(abs, a .- b), result1, result0))
end

main()

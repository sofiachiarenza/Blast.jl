using Blast
using Reactant
using Test
using NPZ

const ROOT = "/home/marcobonici/Desktop/work/CosmologicalEmulators/blast_play/for_marco_copy"
include(joinpath(ROOT, "performance", "workload.jl"))

function host_terms(Plin, Pnon, KG_low, KG_high, KL_low, KL_high, invχ2, weights, Δχ)
    nlow = size(KG_low, 1)
    low = (Pnon .- Plin)[1:nlow, :] .* reshape(invχ2, 1, :)
    high = Pnon[(nlow + 1):end, :] .* reshape(invχ2, 1, :)
    return (
        Blast._limber_contraction(low, KG_low, KG_low, weights, Δχ),
        Blast._limber_contraction(high, KG_high, KG_high, weights, Δχ),
        Blast._limber_contraction(low, KG_low, KL_low, weights, Δχ),
        Blast._limber_contraction(high, KG_high, KL_high, weights, Δχ),
        Blast._limber_contraction(low, KL_low, KL_low, weights, Δχ),
        Blast._limber_contraction(high, KL_high, KL_high, weights, Δχ),
    )
end

function main()
    Reactant.set_default_backend("cpu")
    workload = build_workload(ROOT, 1)
    spectrum = prepare_spectrum(workload)
    KG, KL = Blast.get_limber_kernel(workload.galaxy), Blast.get_limber_kernel(workload.lensing)
    KG_low, KG_high = Blast._low_ℓ_slice(KG), Blast._high_ℓ_slice(KG)
    KL_low, KL_high = Blast._low_ℓ_slice(KL), Blast._high_ℓ_slice(KL)
    invχ2, weights, Δχ = vec(Blast.LIMBER_INV_χ2_ROW), Blast.LIMBER_WEIGHTS, Blast.LIMBER_Δχ
    Plin = spectrum.Pδ_limber .- spectrum.ΔP_limber
    Pnon = spectrum.Pδ_limber
    reference0 = host_terms(Plin, Pnon, KG_low, KG_high, KL_low, KL_high, invχ2, weights, Δχ)
    ordinary0 = (
        Blast.get_limber_correction(KG, spectrum), Blast.get_limber_Cℓ(KG, spectrum),
        Blast.get_limber_correction(KG, KL, spectrum), Blast.get_limber_Cℓ(KG, KL, spectrum),
        Blast.get_limber_correction(KL, spectrum), Blast.get_limber_Cℓ(KL, spectrum),
    )
    args0 = map(Reactant.to_rarray, (Plin, Pnon, KG_low, KG_high, KL_low, KL_high, invχ2, weights, Δχ))
    compiled = Reactant.compile(Blast.reactant_limber_terms_from_prepared, args0; sync=true)
    result0 = map(Array, compiled(args0...))

    Plin1 = Plin .* reshape(1 .+ 1e-5 .* sin.(1:size(Plin, 1)), :, 1)
    Pnon1 = Pnon .* reshape(1 .+ 1e-5 .* cos.(1:size(Pnon, 1)), :, 1)
    args1 = (Reactant.to_rarray(Plin1), Reactant.to_rarray(Pnon1), args0[3:end]...)
    result1 = map(Array, compiled(args1...))
    reference1 = host_terms(Plin1, Pnon1, KG_low, KG_high, KL_low, KL_high, invχ2, weights, Δχ)

    @testset "Reactant full Stage 7" begin
        @test map(size, result0) == map(size, ordinary0)
        for i in 1:6
            @test reference0[i] ≈ ordinary0[i] rtol=1e-12 atol=1e-15
            @test result0[i] ≈ ordinary0[i] rtol=1e-10 atol=1e-15
            @test result1[i] ≈ reference1[i] rtol=1e-10 atol=1e-15
            @test maximum(abs, result1[i] .- result0[i]) > 0
        end
    end
    println("stage7_shapes=", map(size, result0))
    println("stage7_baseline_maxabs=", map((a,b)->maximum(abs,a.-b), result0, ordinary0))
    println("stage7_perturbed_maxabs=", map((a,b)->maximum(abs,a.-b), result1, reference1))
end

main()

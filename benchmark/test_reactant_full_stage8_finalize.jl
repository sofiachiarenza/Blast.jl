using AbstractCosmologicalEmulators
using Blast
using Reactant
using Test
using NPZ

const ROOT = "/home/marcobonici/Desktop/work/CosmologicalEmulators/blast_play/for_marco_copy"
include(joinpath(ROOT, "performance", "workload.jl"))

function main()
    Reactant.set_default_backend("cpu")
    workload = build_workload(ROOT, 1)
    spectrum = prepare_spectrum(workload)
    nonlimber = Blast.reactant_host_nonlimber_reference(workload, spectrum)
    KG, KL = Blast.get_limber_kernel(workload.galaxy), Blast.get_limber_kernel(workload.lensing)
    terms = (
        Blast.get_limber_correction(KG, spectrum), Blast.get_limber_Cℓ(KG, spectrum),
        Blast.get_limber_correction(KG, KL, spectrum), Blast.get_limber_Cℓ(KG, KL, spectrum),
        Blast.get_limber_correction(KL, spectrum), Blast.get_limber_Cℓ(KL, spectrum),
    )
    ordinary = full_power_to_cls(workload; pass_plans=true)
    nfull = length(Blast.full_ℓ_range)
    ell2 = Blast.FULL_ℓ2_REVERSED[1:nfull]
    M = Blast.reactant_chebyshev_matrix(workload.plans.plan_ℓ)
    Teval = chebyshev_polynomials(float.(workload.ell), 2.0, 2000.0, nfull - 1)
    invell2 = inv.(workload.ell .^ 2)
    f = (ngg, ngs, nss, cggc, cggl, cgsc, cgsl, cssc, cssl, e2, m, te, ie2) -> (
        Blast.reactant_finalize_c_ell(ngg, cggc, cggl, e2, m, te, ie2),
        Blast.reactant_finalize_c_ell(ngs, cgsc, cgsl, e2, m, te, ie2),
        Blast.reactant_finalize_c_ell(nss, cssc, cssl, e2, m, te, ie2),
    )
    host_args = (nonlimber.gg, nonlimber.gs, nonlimber.ss, terms..., ell2, M, Teval, invell2)
    args0 = map(Reactant.to_rarray, host_args)
    compiled = Reactant.compile(f, args0; sync=true)
    result0 = map(Array, compiled(args0...))
    references0 = (ordinary.cl_gg, ordinary.cl_gs, ordinary.cl_ss)

    non1 = map(x -> x .* reshape(1 .+ 1e-5 .* sin.(1:size(x, 1)), :, 1, 1),
               (nonlimber.gg, nonlimber.gs, nonlimber.ss))
    terms1 = map(x -> x .* reshape(1 .+ 1e-5 .* cos.(1:size(x, 1)), :, 1, 1), terms)
    host_args1 = (non1..., terms1..., ell2, M, Teval, invell2)
    args1 = map(Reactant.to_rarray, host_args1)
    result1 = map(Array, compiled(args1...))
    references1 = (
        Blast._finalize_Cℓ_parts(non1[1], terms1[1], terms1[2], workload.ell, 10, 10, workload.plans),
        Blast._finalize_Cℓ_parts(non1[2], terms1[3], terms1[4], workload.ell, 10, 5, workload.plans),
        Blast._finalize_Cℓ_parts(non1[3], terms1[5], terms1[6], workload.ell, 5, 5, workload.plans),
    )
    @testset "Reactant full Stage 8" begin
        @test map(size, result0) == map(size, references0)
        for i in 1:3
            @test result0[i] ≈ references0[i] rtol=1e-10 atol=1e-15
            @test result1[i] ≈ references1[i] rtol=1e-10 atol=1e-15
            @test maximum(abs, result1[i] .- result0[i]) > 0
        end
    end
    println("stage8_shapes=", map(size, result0))
    println("stage8_baseline_maxabs=", map((a,b)->maximum(abs,a.-b), result0, references0))
    println("stage8_perturbed_maxabs=", map((a,b)->maximum(abs,a.-b), result1, references1))
end

main()

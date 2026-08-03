using BenchmarkTools
using Blast
using Reactant

const ROOT = "/home/marcobonici/Desktop/work/CosmologicalEmulators/blast_play/for_marco_copy"
include(joinpath(ROOT, "performance", "workload.jl"))

function main()
    Reactant.set_default_backend("cpu")
    Blast.FFTW.set_num_threads(8)
    workload = build_workload(ROOT, 8)
    spectrum = prepare_spectrum(workload)
    W = Blast.compute_w(workload.W_template, spectrum)
    G = workload.galaxy
    bg = workload.bg
    integ = Blast._prepare_nonlimber_integration(bg)

    W_A = Blast.get_kernel_array(G.δ, bg)
    W_B = W_A
    pmj = W.w_2_00_ϕTT.w
    prefactor = Blast._pair_prefactor(G.δ, G.δ)
    reference = Blast._compute_Cℓ_fused_tullio(
        W_A, W_B, pmj, integ.w_χ, integ.w_R,
        prefactor, integ.Δχ, integ.χ_grid,
    )

    pmj1 = 1.07 .* pmj .+ 0.03
    W_A_r = Reactant.to_rarray(W_A)
    W_B_r = Reactant.to_rarray(W_B)
    pmj_r = Reactant.to_rarray(pmj)
    pmj1_r = Reactant.to_rarray(pmj1)
    prefactor_r = Reactant.to_rarray(prefactor)
    wχ_r = Reactant.to_rarray(integ.w_χ)
    wR_r = Reactant.to_rarray(integ.w_R)
    χ_r = Reactant.to_rarray(integ.χ_grid)

    compiled = Reactant.compile(
        Blast.reactant_nonlimber_c_ell,
        (W_A_r, W_B_r, pmj_r, prefactor_r, wχ_r, wR_r, χ_r, integ.Δχ);
        sync=true,
    )
    result0 = Array(compiled(W_A_r, W_B_r, pmj_r, prefactor_r, wχ_r, wR_r, χ_r, integ.Δχ))
    result1 = Array(compiled(W_A_r, W_B_r, pmj1_r, prefactor_r, wχ_r, wR_r, χ_r, integ.Δχ))
    reference1 = Blast._compute_Cℓ_fused_tullio(
        W_A, W_B, pmj1, integ.w_χ, integ.w_R,
        prefactor, integ.Δχ, integ.χ_grid,
    )

    println("Reactant=", pkgversion(Reactant))
    println("shape=", size(reference))
    println("parity_error=", maximum(abs, result0 .- reference))
    println("new_input_error=", maximum(abs, result1 .- reference1))
    println("compiled_delta=", maximum(abs, result1 .- result0))
    println("reference_delta=", maximum(abs, reference1 .- reference))

    host_call = () -> Blast._compute_Cℓ_fused_tullio(
        W_A, W_B, pmj, integ.w_χ, integ.w_R,
        prefactor, integ.Δχ, integ.χ_grid,
    )
    reactant_call = () -> compiled(W_A_r, W_B_r, pmj_r, prefactor_r, wχ_r, wR_r, χ_r, integ.Δχ)
    println("host_ms=", minimum(@benchmark $host_call() samples=5 evals=1).time / 1e6)
    println("reactant_ms=", minimum(@benchmark $reactant_call() samples=10 evals=1).time / 1e6)
end

main()

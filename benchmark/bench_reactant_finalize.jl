using BenchmarkTools
using Blast
using AbstractCosmologicalEmulators
using Reactant

const ROOT = "/home/marcobonici/Desktop/work/CosmologicalEmulators/blast_play/for_marco_copy"
include(joinpath(ROOT, "performance", "workload.jl"))

function main()
    Reactant.set_default_backend("cpu")
    workload = build_workload(ROOT, 8)
    spectrum = prepare_spectrum(workload)
    W = Blast.compute_w(workload.W_template, spectrum)
    G, bg, plans, ell = workload.galaxy, workload.bg, workload.plans, workload.ell
    integ = Blast._prepare_nonlimber_integration(bg)
    Wδ = Blast.get_kernel_array(G.δ, bg)
    nlow = Blast.LIMBER_N_NONLIMBER
    C_nonlimber = Blast._compute_Cℓ_fused_tullio(
        Wδ, Wδ, W.w_2_00_ϕTT.w,
        integ.w_χ, integ.w_R, Blast._pair_prefactor(G.δ, G.δ),
        integ.Δχ, integ.χ_grid,
    )
    K = Blast.get_limber_kernel(G)
    C_correction = Blast.get_limber_correction(K, spectrum)
    C_limber = Blast.get_limber_Cℓ(K, spectrum)
    println("correction_type=", typeof(C_correction), " correction_size=", size(C_correction))

    nfull = size(C_nonlimber, 1) + size(C_limber, 1)
    interp_plan = prepare_chebyshev_plan(2.0, 2000.0, nfull - 1; size_nd=(nfull,), dim=1)
    transform = Blast.reactant_chebyshev_matrix(interp_plan)
    T_eval = chebyshev_polynomials(float.(ell), 2.0, 2000.0, nfull - 1)
    ell2_reversed = Blast.FULL_ℓ2_REVERSED[1:nfull]
    inv_ell2 = inv.(ell .^ 2)
    reference = Blast._finalize_Cℓ_parts(
        C_nonlimber, C_correction, C_limber, ell,
        size(C_nonlimber, 2), size(C_nonlimber, 3), plans,
    )

    args = (
        Reactant.to_rarray(C_nonlimber), Reactant.to_rarray(C_correction),
        Reactant.to_rarray(C_limber), Reactant.to_rarray(ell2_reversed),
        Reactant.to_rarray(transform), Reactant.to_rarray(T_eval),
        Reactant.to_rarray(inv_ell2),
    )
    C_nonlimber1 = 1.03 .* C_nonlimber .+ 0.002
    args1 = (Reactant.to_rarray(C_nonlimber1), args[2:end]...)
    compiled = Reactant.compile(Blast.reactant_finalize_c_ell, args; sync=true)
    result = Array(compiled(args...))
    result1 = Array(compiled(args1...))
    reference1 = Blast._finalize_Cℓ_parts(
        C_nonlimber1, C_correction, C_limber, ell,
        size(C_nonlimber, 2), size(C_nonlimber, 3), plans,
    )
    println("Reactant=", pkgversion(Reactant))
    println("shape=", size(reference))
    println("parity_error=", maximum(abs, result .- reference))
    println("new_input_error=", maximum(abs, result1 .- reference1))
    println("compiled_delta=", maximum(abs, result1 .- result))
    println("reference_delta=", maximum(abs, reference1 .- reference))

    host_call = () -> Blast._finalize_Cℓ_parts(
        C_nonlimber, C_correction, C_limber, ell,
        size(C_nonlimber, 2), size(C_nonlimber, 3), plans,
    )
    reactant_call = () -> compiled(args...)
    println("host_ms=", minimum(@benchmark $host_call() samples=5 evals=1).time / 1e6)
    println("reactant_ms=", minimum(@benchmark $reactant_call() samples=10 evals=1).time / 1e6)
end

main()

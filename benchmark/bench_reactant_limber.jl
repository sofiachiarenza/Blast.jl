using BenchmarkTools
using Blast
using Reactant

const ROOT = "/home/marcobonici/Desktop/work/CosmologicalEmulators/blast_play/for_marco_copy"
include(joinpath(ROOT, "performance", "workload.jl"))

function run_case(label, Pterm, K1, K2, weights, Δχ)
    reference = Blast._limber_contraction(Pterm, K1, K2, weights, Δχ)
    P1 = 1.03 .* Pterm .+ 0.002
    args = (
        Reactant.to_rarray(Pterm),
        Reactant.to_rarray(K1),
        Reactant.to_rarray(K2),
        Reactant.to_rarray(weights),
        Δχ,
    )
    args1 = (Reactant.to_rarray(P1), args[2:end]...)
    compiled = Reactant.compile(Blast.reactant_limber_contraction, args; sync=true)
    result = Array(compiled(args...))
    result1 = Array(compiled(args1...))
    reference1 = Blast._limber_contraction(P1, K1, K2, weights, Δχ)
    host_call = () -> Blast._limber_contraction(Pterm, K1, K2, weights, Δχ)
    reactant_call = () -> compiled(args...)
    println(label, "_shape=", size(reference))
    println(label, "_parity_error=", maximum(abs, result .- reference))
    println(label, "_new_input_error=", maximum(abs, result1 .- reference1))
    println(label, "_compiled_delta=", maximum(abs, result1 .- result))
    println(label, "_reference_delta=", maximum(abs, reference1 .- reference))
    println(label, "_host_ms=", minimum(@benchmark $host_call() samples=5 evals=1).time / 1e6)
    println(label, "_reactant_ms=", minimum(@benchmark $reactant_call() samples=10 evals=1).time / 1e6)
end

function main()
    Reactant.set_default_backend("cpu")
    Blast.FFTW.set_num_threads(8)
    workload = build_workload(ROOT, 8)
    spectrum = prepare_spectrum(workload)
    G = workload.galaxy
    K = Blast.get_limber_kernel(G)
    nlow = Blast.LIMBER_N_NONLIMBER
    P_low = @views spectrum.ΔP_limber[1:nlow, :] .* Blast.LIMBER_INV_χ2_ROW
    P_high = @views spectrum.Pδ_limber[(nlow + 1):end, :] .* Blast.LIMBER_INV_χ2_ROW
    K_low = Blast._low_ℓ_slice(K)
    K_high = Blast._high_ℓ_slice(K)
    println("Reactant=", pkgversion(Reactant))
    run_case("limber_correction", P_low, K_low, K_low, Blast.LIMBER_WEIGHTS, Blast.LIMBER_Δχ)
    run_case("limber_high", P_high, K_high, K_high, Blast.LIMBER_WEIGHTS, Blast.LIMBER_Δχ)
end

main()

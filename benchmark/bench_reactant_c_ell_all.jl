using BenchmarkTools
using Blast
using Reactant

const ROOT = "/home/marcobonici/Desktop/work/CosmologicalEmulators/blast_play/for_marco_copy"
include(joinpath(ROOT, "performance", "workload.jl"))

function common_args(bg)
    integ = Blast._prepare_nonlimber_integration(bg)
    return integ, (
        Reactant.to_rarray(integ.w_χ),
        Reactant.to_rarray(integ.w_R),
        Reactant.to_rarray(integ.χ_grid),
    )
end

function run_standard(label, A, B, pmj, integ, dynamic_inputs)
    prefactor = Blast._pair_prefactor(A, B)
    W_A = Blast.get_kernel_array(A, integ.background)
    W_B = Blast.get_kernel_array(B, integ.background)
    reference = Blast._compute_Cℓ_fused_tullio(
        W_A, W_B, pmj, integ.w_χ, integ.w_R,
        prefactor, integ.Δχ, integ.χ_grid,
    )
    pmj1 = 1.07 .* pmj .+ 0.03
    args = (
        Reactant.to_rarray(W_A), Reactant.to_rarray(W_B),
        Reactant.to_rarray(pmj), Reactant.to_rarray(prefactor),
        dynamic_inputs..., integ.Δχ,
    )
    args1 = (args[1], args[2], Reactant.to_rarray(pmj1), args[4], args[5:end]...)
    compiled = Reactant.compile(Blast.reactant_nonlimber_c_ell, args; sync=true)
    result = Array(compiled(args...))
    result1 = Array(compiled(args1...))
    reference1 = Blast._compute_Cℓ_fused_tullio(
        W_A, W_B, pmj1, integ.w_χ, integ.w_R,
        prefactor, integ.Δχ, integ.χ_grid,
    )
    host_call = () -> Blast._compute_Cℓ_fused_tullio(
        W_A, W_B, pmj, integ.w_χ, integ.w_R,
        prefactor, integ.Δχ, integ.χ_grid,
    )
    reactant_call = () -> compiled(args...)
    println(label, "_shape=", size(reference))
    println(label, "_parity_error=", maximum(abs, result .- reference))
    println(label, "_new_input_error=", maximum(abs, result1 .- reference1))
    println(label, "_compiled_delta=", maximum(abs, result1 .- result))
    println(label, "_reference_delta=", maximum(abs, reference1 .- reference))
    println(label, "_host_ms=", minimum(@benchmark $host_call() samples=5 evals=1).time / 1e6)
    println(label, "_reactant_ms=", minimum(@benchmark $reactant_call() samples=10 evals=1).time / 1e6)
end

function run_rsd(label, A, B, pmj02, pmj20, integ, dynamic_inputs)
    prefactor = Blast._pair_prefactor(A, B)
    W_A = Blast.get_kernel_array(A, integ.background)
    W_B = Blast.get_kernel_array(B, integ.background)
    reference = Blast._compute_Cℓ_rsd_fused_tullio(
        W_A, W_B, pmj02, pmj20, integ.w_χ, integ.w_R,
        prefactor, integ.Δχ, integ.χ_grid,
    )
    pmj021 = 1.07 .* pmj02 .+ 0.03
    args = (
        Reactant.to_rarray(W_A), Reactant.to_rarray(W_B),
        Reactant.to_rarray(pmj02), Reactant.to_rarray(pmj20),
        Reactant.to_rarray(prefactor), dynamic_inputs..., integ.Δχ,
    )
    args1 = (args[1], args[2], Reactant.to_rarray(pmj021), args[4:end]...)
    compiled = Reactant.compile(Blast.reactant_nonlimber_rsd_c_ell, args; sync=true)
    result = Array(compiled(args...))
    result1 = Array(compiled(args1...))
    reference1 = Blast._compute_Cℓ_rsd_fused_tullio(
        W_A, W_B, pmj021, pmj20, integ.w_χ, integ.w_R,
        prefactor, integ.Δχ, integ.χ_grid,
    )
    println(label, "_shape=", size(reference))
    println(label, "_parity_error=", maximum(abs, result .- reference))
    println(label, "_new_input_error=", maximum(abs, result1 .- reference1))
    println(label, "_compiled_delta=", maximum(abs, result1 .- result))
    println(label, "_reference_delta=", maximum(abs, reference1 .- reference))
end

function main()
    Reactant.set_default_backend("cpu")
    Blast.FFTW.set_num_threads(8)
    workload = build_workload(ROOT, 8)
    spectrum = prepare_spectrum(workload)
    W = Blast.compute_w(workload.W_template, spectrum)
    G, L = workload.galaxy, workload.lensing
    integ0 = Blast._prepare_nonlimber_integration(workload.bg)
    integ = merge(integ0, (; background=workload.bg))
    dynamic_inputs = (
        Reactant.to_rarray(integ.w_χ),
        Reactant.to_rarray(integ.w_R),
        Reactant.to_rarray(integ.χ_grid),
    )

    println("Reactant=", pkgversion(Reactant))
    run_standard("GG", G.δ, G.δ, W.w_2_00_ϕTT.w, integ, dynamic_inputs)
    run_standard("SS", L.γ, L.γ, W.w_minus2_00_ϕTT.w, integ, dynamic_inputs)
    run_standard("GS", G.δ, L.γ, W.w_0_00_ϕTT.w, integ, dynamic_inputs)
    run_rsd("G_RSD", G.δ, G.RSD, W.w_2_02_ϕTT.w, W.w_2_20_ϕTT.w, integ, dynamic_inputs)
end

main()

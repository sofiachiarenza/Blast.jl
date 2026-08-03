using AbstractCosmologicalEmulators
using BenchmarkTools
using Blast
using LinearAlgebra
using Random
using Reactant

function report(label, trial)
    t = minimum(trial)
    println(label, "_ms=", t.time / 1e6, "_alloc_bytes=", t.memory)
end

function main()
    Reactant.set_default_backend("cpu")
    Blast.FFTW.set_num_threads(8)
    Random.seed!(20260802)

    nk = length(Blast.k_cheb)
    nχ = length(Blast.χ)
    nR = length(Blast.R)
    println("reactant=", pkgversion(Reactant))
    println("blast_reactant_available=", Blast.reactant_is_available())
    println("grid=(", nk, ",", nχ, ",", nR, ")")

    P0 = randn(nk)
    Tχ10 = randn(nχ, nk)
    TχR0 = randn(nχ, nR, nk)
    P1 = 0.3 .* randn(nk) .+ 0.2
    Tχ11 = 0.4 .* randn(nχ, nk) .- 0.1
    TχR1 = 0.2 .* randn(nχ, nR, nk) .+ 0.3

    ptt_ref = Blast._p_phi_TT_tullio(P0, Tχ10, TχR0)
    pt_ref = Blast._p_phi_T_tullio(P0, TχR0)
    P0r, Tχ10r, TχR0r = map(Reactant.to_rarray, (P0, Tχ10, TχR0))
    P1r, Tχ11r, TχR1r = map(Reactant.to_rarray, (P1, Tχ11, TχR1))

    ptt_compiled = Reactant.compile(Blast.reactant_p_phi_TT, (P0r, Tχ10r, TχR0r); sync=true)
    pt_compiled = Reactant.compile(Blast.reactant_p_phi_T, (P0r, TχR0r); sync=true)
    ptt0 = Array(ptt_compiled(P0r, Tχ10r, TχR0r))
    ptt1 = Array(ptt_compiled(P1r, Tχ11r, TχR1r))
    pt0 = Array(pt_compiled(P0r, TχR0r))
    pt1 = Array(pt_compiled(P1r, TχR1r))
    println("p_phi_TT_error=", maximum(abs, ptt0 .- ptt_ref))
    println("p_phi_TT_new_input_error=", maximum(abs, ptt1 .- Blast._p_phi_TT_tullio(P1, Tχ11, TχR1)))
    println("p_phi_TT_compiled_delta=", maximum(abs, ptt1 .- ptt0))
    println("p_phi_TT_reference_delta=", maximum(abs, Blast._p_phi_TT_tullio(P1, Tχ11, TχR1) .- ptt_ref))
    println("p_phi_T_error=", maximum(abs, pt0 .- pt_ref))
    println("p_phi_T_new_input_error=", maximum(abs, pt1 .- Blast._p_phi_T_tullio(P1, TχR1)))
    println("p_phi_T_compiled_delta=", maximum(abs, pt1 .- pt0))
    println("p_phi_T_reference_delta=", maximum(abs, Blast._p_phi_T_tullio(P1, TχR1) .- pt_ref))

    plan = prepare_chebyshev_plan(
        minimum(Blast.k_cheb), maximum(Blast.k_cheb), nk - 1;
        size_nd=(nk,), dim=1,
    )
    transform = Blast.reactant_chebyshev_matrix(plan)
    vals0 = randn(nk, nχ, nR)
    vals1 = 0.2 .* randn(nk, nχ, nR) .+ 0.1
    cheb_ref0 = chebyshev_decomposition(plan, vals0)
    cheb_ref1 = chebyshev_decomposition(plan, vals1)
    vals0r = Reactant.to_rarray(vals0)
    vals1r = Reactant.to_rarray(vals1)
    transformr = Reactant.to_rarray(transform)
    cheb_compiled = Reactant.compile(Blast.reactant_chebyshev_matmul, (vals0r, transformr); sync=true)
    cheb0 = Array(cheb_compiled(vals0r, transformr))
    cheb1 = Array(cheb_compiled(vals1r, transformr))
    println("chebyshev_error=", maximum(abs, cheb0 .- cheb_ref0))
    println("chebyshev_new_input_error=", maximum(abs, cheb1 .- cheb_ref1))
    println("chebyshev_compiled_delta=", maximum(abs, cheb1 .- cheb0))
    println("chebyshev_reference_delta=", maximum(abs, cheb_ref1 .- cheb_ref0))

    T = Blast.T_tildes.T_2_00
    T_lijk = permutedims(T, (4, 1, 2, 3))
    c0 = cheb_ref0
    c1 = cheb_ref1
    T1 = 0.73 .* T
    c0r = Reactant.to_rarray(c0)
    c1r = Reactant.to_rarray(c1)
    Tr = Reactant.to_rarray(T_lijk)
    T1r = Reactant.to_rarray(permutedims(T1, (4, 1, 2, 3)))
    w_compiled = Reactant.compile(Blast.reactant_w_ell_hlo, (c0r, Tr); sync=true)
    w0 = Array(w_compiled(c0r, Tr))
    w_c1 = Array(w_compiled(c1r, Tr))
    w_T1 = Array(w_compiled(c0r, T1r))
    w_ref0 = Blast.w_ell_tullio(c0, T)
    w_c1_ref = Blast.w_ell_tullio(c1, T)
    w_T1_ref = Blast.w_ell_tullio(c0, T1)
    println("w_error=", maximum(abs, w0 .- w_ref0))
    println("w_new_coefficient_error=", maximum(abs, w_c1 .- w_c1_ref))
    println("w_new_kernel_error=", maximum(abs, w_T1 .- w_T1_ref))
    println("w_coefficient_compiled_delta=", maximum(abs, w_c1 .- w0))
    println("w_coefficient_reference_delta=", maximum(abs, w_c1_ref .- w_ref0))
    println("w_kernel_compiled_delta=", maximum(abs, w_T1 .- w0))
    println("w_kernel_reference_delta=", maximum(abs, w_T1_ref .- w_ref0))

    report("p_phi_TT_host", @benchmark Blast._p_phi_TT_tullio($P0, $Tχ10, $TχR0) samples=10 evals=1)
    report("p_phi_TT_reactant", @benchmark $ptt_compiled($P0r, $Tχ10r, $TχR0r) samples=20 evals=1)
    report("chebyshev_host", @benchmark chebyshev_decomposition($plan, $vals0) samples=10 evals=1)
    report("chebyshev_reactant", @benchmark $cheb_compiled($vals0r, $transformr) samples=20 evals=1)
    report("w_host", @benchmark Blast.w_ell_tullio($c0, $T) samples=10 evals=1)
    report("w_reactant", @benchmark $w_compiled($c0r, $Tr) samples=20 evals=1)
end

main()

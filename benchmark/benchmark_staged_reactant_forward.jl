using AbstractCosmologicalEmulators
using BenchmarkTools
using Blast
using Reactant

const ROOT = "/home/marcobonici/Desktop/work/CosmologicalEmulators/blast_play/for_marco_copy"
include(joinpath(ROOT, "performance", "workload.jl"))

function w_fields(W)
    return (
        W.w_2_00_ϕTT.w, W.w_minus2_00_ϕTT.w, W.w_0_00_ϕTT.w,
        W.w_0_02_ϕTT.w, W.w_0_20_ϕTT.w, W.w_2_02_ϕTT.w,
        W.w_2_20_ϕTT.w, W.w_2_22_ϕTT.w, W.w_2_00_ϕT.w,
        W.w_2_00_ϕT_R1.w, W.w_0_00_ϕT.w, W.w_0_00_ϕT_R1.w,
        W.w_2_02_ϕT.w, W.w_2_02_ϕT_R1.w, W.w_2_20_ϕT.w,
        W.w_2_20_ϕT_R1.w, W.w_2_00_ϕ.w,
    )
end

function prefactors(G, L)
    return (
        δδ=Blast._pair_prefactor(G.δ, G.δ), δRSD=Blast._pair_prefactor(G.δ, G.RSD),
        RSDRSD=Blast._pair_prefactor(G.RSD, G.RSD), δμ=Blast._pair_prefactor(G.δ, G.μ),
        μμ=Blast._pair_prefactor(G.μ, G.μ), μRSD=Blast._pair_prefactor(G.μ, G.RSD),
        δfNL=Blast._pair_prefactor(G.δ, G.PNG), fNLδ=Blast._pair_prefactor(G.PNG, G.δ),
        fNLRSD=Blast._pair_prefactor(G.PNG, G.RSD), RSDfNL=Blast._pair_prefactor(G.RSD, G.PNG),
        μfNL=Blast._pair_prefactor(G.μ, G.PNG), fNLμ=Blast._pair_prefactor(G.PNG, G.μ),
        fNLfNL=Blast._pair_prefactor(G.PNG, G.PNG), γγ=Blast._pair_prefactor(L.γ, L.γ),
        γIA=Blast._pair_prefactor(L.γ, L.IA), IAγ=Blast._pair_prefactor(L.IA, L.γ),
        IAIA=Blast._pair_prefactor(L.IA, L.IA), δγ=Blast._pair_prefactor(G.δ, L.γ),
        δIA=Blast._pair_prefactor(G.δ, L.IA), RSDγ=Blast._pair_prefactor(G.RSD, L.γ),
        RSDIA=Blast._pair_prefactor(G.RSD, L.IA), μγ=Blast._pair_prefactor(G.μ, L.γ),
        μIA=Blast._pair_prefactor(G.μ, L.IA), fNLγ=Blast._pair_prefactor(G.PNG, L.γ),
        fNLIA=Blast._pair_prefactor(G.PNG, L.IA),
    )
end

function main()
    Reactant.set_default_backend("cpu")
    Blast.FFTW.set_num_threads(8)
    workload = build_workload(ROOT, 8)
    G, L, bg = workload.galaxy, workload.lensing, workload.bg
    ordinary = full_power_to_cls(workload; pass_plans=true)
    pref = prefactors(G, L)
    integ = Blast._prepare_nonlimber_integration(bg)
    K = (
        Blast.get_kernel_array(G.δ, bg), Blast.get_kernel_array(G.RSD, bg),
        Blast.get_kernel_array(G.μ, bg), Blast.get_kernel_array(G.PNG, bg),
        Blast.get_kernel_array(L.γ, bg), Blast.get_kernel_array(L.IA, bg),
    )
    T = map(T -> permutedims(T, (4, 1, 2, 3)), (
        Blast.T_tildes.T_2_00, Blast.T_tildes.T_minus2_00,
        Blast.T_tildes.T_0_00, Blast.T_tildes.T_0_02,
        Blast.T_tildes.T_0_20, Blast.T_tildes.T_2_02,
        Blast.T_tildes.T_2_20, Blast.T_tildes.T_2_22,
    ))
    Mnon_plan = prepare_chebyshev_plan(
        minimum(Blast.k_cheb), maximum(Blast.k_cheb), length(Blast.k_cheb) - 1;
        size_nd=(length(Blast.k_cheb),), dim=1,
    )
    Mnon = Blast.reactant_chebyshev_matrix(Mnon_plan)
    lp = workload.plans.plan_limber
    Kk, Kz = lp.K
    Mk = Blast.reactant_chebyshev_matrix(prepare_chebyshev_plan(
        minimum(lp.nodes[1]), maximum(lp.nodes[1]), Kk;
        size_nd=(Kk + 1,), dim=1,
    ))
    Mz = Blast.reactant_chebyshev_matrix(prepare_chebyshev_plan(
        minimum(lp.nodes[2]), maximum(lp.nodes[2]), Kz;
        size_nd=(Kz + 1,), dim=1,
    ))
    Tz = Blast.get_limber_coords_polynomials(lp, bg.z)
    Tk = workload.plans.T_k_limber
    KG, KL = Blast.get_limber_kernel(G), Blast.get_limber_kernel(L)
    KGlo, KGhi = Blast._low_ℓ_slice(KG), Blast._high_ℓ_slice(KG)
    KLlo, KLhi = Blast._low_ℓ_slice(KL), Blast._high_ℓ_slice(KL)
    invχ2 = vec(Blast.LIMBER_INV_χ2_ROW)
    nfull = length(Blast.ℓ_nonlimber) + length(Blast.ℓ_limber)
    final_plan = prepare_chebyshev_plan(2.0, 2000.0, nfull - 1; size_nd=(nfull,), dim=1)
    finalization = (
        Blast.FULL_ℓ2_REVERSED[1:nfull], Blast.reactant_chebyshev_matrix(final_plan),
        chebyshev_polynomials(float.(workload.ell), 2.0, 2000.0, nfull - 1),
        inv.(workload.ell .^ 2),
    )

    pk_r = Reactant.to_rarray(workload.pk)
    pkli_r = Reactant.to_rarray(workload.pk_limber_lin)
    pknl_r = Reactant.to_rarray(workload.pk_limber_nonlin)
    Mnon_r, Mk_r, Mz_r = map(Reactant.to_rarray, (Mnon, Mk, Mz))
    Tz_r, Tk_r = map(Reactant.to_rarray, (Tz, Tk))
    T_r = map(Reactant.to_rarray, T)
    K_r = map(Reactant.to_rarray, K)
    integ_r = map(Reactant.to_rarray, (integ.w_χ, integ.w_R, integ.χ_grid))
    final_r = map(Reactant.to_rarray, finalization)
    weights_r = Reactant.to_rarray(Blast.LIMBER_WEIGHTS)
    Δχ = Blast.LIMBER_Δχ

    prep_fun = let bg=bg
        (pk, M) -> Blast.reactant_prepare_nonlimber_spectrum(pk, bg, M)
    end
    limber_fun = let bg=bg
        (pkl, pknl, Mk_, Mz_, Tz_, Tk_) -> Blast.reactant_prepare_limber(pkl, pknl, bg, Mk_, Mz_, Tz_, Tk_)
    end
    w_fun = Blast.reactant_compute_w_from_spectrum
    limber_terms_fun = (ΔP, Pnon, KGGlo, KGGhi, KGlo, KLlo, KGhi, KLhi, KLLlo, KLLhi, inv2, weights) ->
        Blast.reactant_limber_terms_from_prepared(
            ΔP, Pnon, KGGlo, KGGhi, KGlo, KLlo, KGhi, KLhi, KLLlo, KLLhi, inv2, weights, Δχ,
        )
    full_fun = let pref=pref, Δχ=integ.Δχ
        (w1,w2,w3,w4,w5,w6,w7,w8,w9,w10,w11,w12,w13,w14,w15,w16,w17,
         k1,k2,k3,k4,k5,k6,wi,wr,χ,
         cgc,cgh,csc,csh,cec,ceh,ell2,M,Te,inv2) -> Blast.reactant_full_3x2pt(
            w1,w2,w3,w4,w5,w6,w7,w8,w9,w10,w11,w12,w13,w14,w15,w16,w17,
            k1,k2,k3,k4,k5,k6,wi,wr,χ,cgc,cgh,csc,csh,cec,ceh,
            ell2,M,Te,inv2,Δχ,pref,
        )
    end

    println("compiling_staged_forward=true")
    t0 = time_ns()
    prep_comp = Reactant.compile(prep_fun, (pk_r, Mnon_r); sync=true)
    limber_comp = Reactant.compile(limber_fun, (pkli_r, pknl_r, Mk_r, Mz_r, Tz_r, Tk_r); sync=true)
    c0 = prep_comp(pk_r, Mnon_r)
    w_comp = Reactant.compile(w_fun, (c0..., T_r...); sync=true)
    lim0 = limber_comp(pkli_r, pknl_r, Mk_r, Mz_r, Tz_r, Tk_r)
    limber_terms_comp = Reactant.compile(
        limber_terms_fun,
        (lim0[1], lim0[2], Reactant.to_rarray(KGlo), Reactant.to_rarray(KGhi),
         Reactant.to_rarray(KGlo), Reactant.to_rarray(KLlo), Reactant.to_rarray(KGhi),
         Reactant.to_rarray(KLhi), Reactant.to_rarray(KLlo), Reactant.to_rarray(KLhi),
         Reactant.to_rarray(invχ2), weights_r),
        sync=true,
    )
    C0 = limber_terms_comp(
        lim0[1], lim0[2], Reactant.to_rarray(KGlo), Reactant.to_rarray(KGhi),
        Reactant.to_rarray(KGlo), Reactant.to_rarray(KLlo), Reactant.to_rarray(KGhi),
        Reactant.to_rarray(KLhi), Reactant.to_rarray(KLlo), Reactant.to_rarray(KLhi),
        Reactant.to_rarray(invχ2), weights_r,
    )
    w0 = w_comp(c0..., T_r...)
    full_args = (w0..., K_r..., integ_r..., C0..., final_r...)
    full_comp = Reactant.compile(full_fun, full_args; sync=true)
    compile_ms = (time_ns() - t0) / 1e6

    function staged_forward()
        c = prep_comp(pk_r, Mnon_r)
        wv = w_comp(c..., T_r...)
        lv = limber_comp(pkli_r, pknl_r, Mk_r, Mz_r, Tz_r, Tk_r)
        cc = limber_terms_comp(
            lv[1], lv[2], Reactant.to_rarray(KGlo), Reactant.to_rarray(KGhi),
            Reactant.to_rarray(KGlo), Reactant.to_rarray(KLlo), Reactant.to_rarray(KGhi),
            Reactant.to_rarray(KLhi), Reactant.to_rarray(KLlo), Reactant.to_rarray(KLhi),
            Reactant.to_rarray(invχ2), weights_r,
        )
        return map(Array, full_comp(
            wv..., K_r..., integ_r...,
            cc..., final_r...,
        ))
    end
    result = staged_forward()
    println("staged_compile_ms=", compile_ms)
    println("GG_parity_error=", maximum(abs, result[1] .- ordinary.cl_gg))
    println("GS_parity_error=", maximum(abs, result[2] .- ordinary.cl_gs))
    println("SS_parity_error=", maximum(abs, result[3] .- ordinary.cl_ss))
    bench = @benchmark $staged_forward() samples=3 evals=1
    println("staged_forward_min_ms=", minimum(bench).time / 1e6)
    println("staged_forward_median_ms=", median(bench).time / 1e6)
end

main()

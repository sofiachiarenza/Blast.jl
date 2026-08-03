using BenchmarkTools
using AbstractCosmologicalEmulators
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

function main()
    Reactant.set_default_backend("cpu")
    Blast.FFTW.set_num_threads(8)
    workload = build_workload(ROOT, 8)
    spectrum = prepare_spectrum(workload)
    W = Blast.compute_w(workload.W_template, spectrum)
    G, L, bg = workload.galaxy, workload.lensing, workload.bg
    integ = Blast._prepare_nonlimber_integration(bg)
    kernels = (
        Blast.get_kernel_array(G.δ, bg), Blast.get_kernel_array(G.RSD, bg),
        Blast.get_kernel_array(G.μ, bg), Blast.get_kernel_array(G.PNG, bg),
        Blast.get_kernel_array(L.γ, bg), Blast.get_kernel_array(L.IA, bg),
    )
    prefactors = (
        δδ=Blast._pair_prefactor(G.δ, G.δ),
        δRSD=Blast._pair_prefactor(G.δ, G.RSD),
        RSDRSD=Blast._pair_prefactor(G.RSD, G.RSD),
        δμ=Blast._pair_prefactor(G.δ, G.μ),
        μμ=Blast._pair_prefactor(G.μ, G.μ),
        μRSD=Blast._pair_prefactor(G.μ, G.RSD),
        δfNL=Blast._pair_prefactor(G.δ, G.PNG),
        fNLδ=Blast._pair_prefactor(G.PNG, G.δ),
        fNLRSD=Blast._pair_prefactor(G.PNG, G.RSD),
        RSDfNL=Blast._pair_prefactor(G.RSD, G.PNG),
        μfNL=Blast._pair_prefactor(G.μ, G.PNG),
        fNLμ=Blast._pair_prefactor(G.PNG, G.μ),
        fNLfNL=Blast._pair_prefactor(G.PNG, G.PNG),
        γγ=Blast._pair_prefactor(L.γ, L.γ),
        γIA=Blast._pair_prefactor(L.γ, L.IA),
        IAγ=Blast._pair_prefactor(L.IA, L.γ),
        IAIA=Blast._pair_prefactor(L.IA, L.IA),
        δγ=Blast._pair_prefactor(G.δ, L.γ),
        δIA=Blast._pair_prefactor(G.δ, L.IA),
        RSDγ=Blast._pair_prefactor(G.RSD, L.γ),
        RSDIA=Blast._pair_prefactor(G.RSD, L.IA),
        μγ=Blast._pair_prefactor(G.μ, L.γ),
        μIA=Blast._pair_prefactor(G.μ, L.IA),
        fNLγ=Blast._pair_prefactor(G.PNG, L.γ),
        fNLIA=Blast._pair_prefactor(G.PNG, L.IA),
    )
    w0 = w_fields(W)
    w1 = map(w -> 1.01 .* w .+ 0.001, w0)
    w0r = map(Reactant.to_rarray, w0)
    w1r = map(Reactant.to_rarray, w1)
    kernels_r = map(Reactant.to_rarray, kernels)
    integration_r = (
        Reactant.to_rarray(integ.w_χ),
        Reactant.to_rarray(integ.w_R),
        Reactant.to_rarray(integ.χ_grid),
    )

    full_fun = let prefactors=prefactors
        (a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,
         Kδ,KRSD,Kμ,KfNL,Kγ,KIA,wχ,wR,χ_grid) -> Blast.reactant_nonlimber_3x2pt(
            a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,
            Kδ,KRSD,Kμ,KfNL,Kγ,KIA,wχ,wR,χ_grid,integ.Δχ,prefactors,
        )
    end
    host_args = (w0..., kernels..., integ.w_χ, integ.w_R, integ.χ_grid)
    reactant_args = (w0r..., kernels_r..., integration_r...)
    reactant_args1 = (w1r..., kernels_r..., integration_r...)
    reference = full_fun(host_args...)
    reference1 = full_fun((w1..., kernels..., integ.w_χ, integ.w_R, integ.χ_grid)...)
    compiled = Reactant.compile(full_fun, reactant_args; sync=true)
    result = compiled(reactant_args...)
    result1 = compiled(reactant_args1...)
    result_arrays = map(Array, result)
    result1_arrays = map(Array, result1)

    println("Reactant=", pkgversion(Reactant))
    for (label, r, r1, a, a1) in zip(("GG", "GS", "SS"), result_arrays, result1_arrays, reference, reference1)
        println(label, "_shape=", size(a))
        println(label, "_parity_error=", maximum(abs, r .- a))
        println(label, "_new_input_error=", maximum(abs, r1 .- a1))
        println(label, "_compiled_delta=", maximum(abs, r1 .- r))
        println(label, "_reference_delta=", maximum(abs, a1 .- a))
    end
    host_call = () -> full_fun(host_args...)
    reactant_call = () -> compiled(reactant_args...)
    println("host_ms=", minimum(@benchmark $host_call() samples=3 evals=1).time / 1e6)
    println("reactant_ms=", minimum(@benchmark $reactant_call() samples=5 evals=1).time / 1e6)

    K_G = Blast.get_limber_kernel(G)
    K_L = Blast.get_limber_kernel(L)
    Cgg_correction = Blast.get_limber_correction(K_G, spectrum)
    Cgg_limber = Blast.get_limber_Cℓ(K_G, spectrum)
    Cgs_correction = Blast.get_limber_correction(K_G, K_L, spectrum)
    Cgs_limber = Blast.get_limber_Cℓ(K_G, K_L, spectrum)
    Css_correction = Blast.get_limber_correction(K_L, spectrum)
    Css_limber = Blast.get_limber_Cℓ(K_L, spectrum)
    nfull = size(reference.gg, 1) + size(Cgg_limber, 1)
    interp_plan = prepare_chebyshev_plan(2.0, 2000.0, nfull - 1; size_nd=(nfull,), dim=1)
    transform = Blast.reactant_chebyshev_matrix(interp_plan)
    T_eval = chebyshev_polynomials(float.(workload.ell), 2.0, 2000.0, nfull - 1)
    ell2_reversed = Blast.FULL_ℓ2_REVERSED[1:nfull]
    inv_ell2 = inv.(workload.ell .^ 2)
    final_static = (
        Reactant.to_rarray(Cgg_correction), Reactant.to_rarray(Cgg_limber),
        Reactant.to_rarray(Cgs_correction), Reactant.to_rarray(Cgs_limber),
        Reactant.to_rarray(Css_correction), Reactant.to_rarray(Css_limber),
        Reactant.to_rarray(ell2_reversed), Reactant.to_rarray(transform),
        Reactant.to_rarray(T_eval), Reactant.to_rarray(inv_ell2),
    )
    full_fun = let prefactors=prefactors
        (a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,
         Kδ,KRSD,Kμ,KfNL,Kγ,KIA,wχ,wR,χ_grid,
         Cggc,Cggl,Cgsc,Cgsl,Cssc,Cssl,ell2,M,Tout,inv2) -> Blast.reactant_full_3x2pt(
            a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,
            Kδ,KRSD,Kμ,KfNL,Kγ,KIA,wχ,wR,χ_grid,
            Cggc,Cggl,Cgsc,Cgsl,Cssc,Cssl,ell2,M,Tout,inv2,integ.Δχ,prefactors,
        )
    end
    full_host_args = (w0..., kernels..., integ.w_χ, integ.w_R, integ.χ_grid,
                      Cgg_correction, Cgg_limber, Cgs_correction, Cgs_limber,
                      Css_correction, Css_limber, ell2_reversed, transform,
                      T_eval, inv_ell2)
    full_reactant_args = (w0r..., kernels_r..., integration_r..., final_static...)
    full_reactant_args1 = (w1r..., kernels_r..., integration_r..., final_static...)
    println("compiling_full_flat_reactant=true")
    full_compiled = Reactant.compile(full_fun, full_reactant_args; sync=true)
    full_result = map(Array, full_compiled(full_reactant_args...))
    full_result1 = map(Array, full_compiled(full_reactant_args1...))
    full_reference = full_fun(full_host_args...)
    ordinary_full = full_power_to_cls(workload; pass_plans=true)
    println("full_GG_parity_error=", maximum(abs, full_result[1] .- ordinary_full.cl_gg))
    println("full_GS_parity_error=", maximum(abs, full_result[2] .- ordinary_full.cl_gs))
    println("full_SS_parity_error=", maximum(abs, full_result[3] .- ordinary_full.cl_ss))
    println("full_GG_flat_reference_error=", maximum(abs, full_result[1] .- full_reference[1]))
    println("full_GS_flat_reference_error=", maximum(abs, full_result[2] .- full_reference[2]))
    println("full_SS_flat_reference_error=", maximum(abs, full_result[3] .- full_reference[3]))
    full_reference1 = full_fun((w1..., kernels..., integ.w_χ, integ.w_R, integ.χ_grid,
                                Cgg_correction, Cgg_limber, Cgs_correction, Cgs_limber,
                                Css_correction, Css_limber, ell2_reversed, transform,
                                T_eval, inv_ell2)...)
    println("full_GG_new_input_error=", maximum(abs, full_result1[1] .- full_reference1[1]))
    println("full_GS_new_input_error=", maximum(abs, full_result1[2] .- full_reference1[2]))
    println("full_SS_new_input_error=", maximum(abs, full_result1[3] .- full_reference1[3]))
    full_call = () -> full_compiled(full_reactant_args...)
    full_host_call = () -> full_fun(full_host_args...)
    println("full_host_ms=", minimum(@benchmark $full_host_call() samples=3 evals=1).time / 1e6)
    println("full_reactant_ms=", minimum(@benchmark $full_call() samples=3 evals=1).time / 1e6)
end

main()

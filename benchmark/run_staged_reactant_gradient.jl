using AbstractCosmologicalEmulators
using BenchmarkTools
using Blast
using Random
using Reactant

const ROOT = "/home/marcobonici/Desktop/work/CosmologicalEmulators/blast_play/for_marco_copy"
include(joinpath(ROOT, "performance", "workload.jl"))
include(joinpath(@__DIR__, "staged_reactant_gradient.jl"))

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
    workload = build_workload(ROOT, 8)
    G, L, bg = workload.galaxy, workload.lensing, workload.bg
    spectrum = prepare_spectrum(workload)
    W = Blast.compute_w(workload.W_template, spectrum)
    w = w_fields(W)
    integ = Blast._prepare_nonlimber_integration(bg)
    K = (
        Blast.get_kernel_array(G.δ, bg), Blast.get_kernel_array(G.RSD, bg),
        Blast.get_kernel_array(G.μ, bg), Blast.get_kernel_array(G.PNG, bg),
        Blast.get_kernel_array(L.γ, bg), Blast.get_kernel_array(L.IA, bg),
    )
    pref = prefactors(G, L)
    nonlimber = Blast.reactant_nonlimber_3x2pt(
        w..., K..., integ.w_χ, integ.w_R, integ.χ_grid, integ.Δχ, pref,
    )
    KG, KL = Blast.get_limber_kernel(G), Blast.get_limber_kernel(L)
    Ccorr = (
        Blast.get_limber_correction(KG, spectrum),
        Blast.get_limber_correction(KG, KL, spectrum),
        Blast.get_limber_correction(KL, spectrum),
    )
    Clim = (
        Blast.get_limber_Cℓ(KG, spectrum),
        Blast.get_limber_Cℓ(KG, KL, spectrum),
        Blast.get_limber_Cℓ(KL, spectrum),
    )
    nfull = size(nonlimber.gg, 1) + size(Clim[1], 1)
    nonlimber_plan = prepare_chebyshev_plan(
        minimum(Blast.k_cheb), maximum(Blast.k_cheb), length(Blast.k_cheb) - 1;
        size_nd=(length(Blast.k_cheb),), dim=1,
    )
    Mnonlimber = Blast.reactant_chebyshev_matrix(nonlimber_plan)
    reactant_coefficients = Blast.reactant_prepare_nonlimber_spectrum(workload.pk, bg, Mnonlimber)
    plan = prepare_chebyshev_plan(2.0, 2000.0, nfull - 1; size_nd=(nfull,), dim=1)
    finalization = (
        Blast.FULL_ℓ2_REVERSED[1:nfull],
        Blast.reactant_chebyshev_matrix(plan),
        chebyshev_polynomials(float.(workload.ell), 2.0, 2000.0, nfull - 1),
        inv.(workload.ell .^ 2),
    )
    nell = length(workload.ell)
    upstream = (ones(nell, size(nonlimber.gg, 2), size(nonlimber.gg, 3)),
                ones(nell, size(nonlimber.gs, 2), size(nonlimber.gs, 3)),
                ones(nell, size(nonlimber.ss, 2), size(nonlimber.ss, 3)))
    dCnon, dCcorr, dClim = staged_finalization_pullback(
        (nonlimber.gg, nonlimber.gs, nonlimber.ss), Ccorr, Clim, finalization, upstream,
    )
    println("dCnon_norms=", map(x -> maximum(abs, x), dCnon))
    println("dCcorr_norms=", map(x -> maximum(abs, x), dCcorr))
    println("dClim_norms=", map(x -> maximum(abs, x), dClim))
    endpoint_config = (
        kernels=K, integration=(integ.w_χ, integ.w_R, integ.χ_grid),
        C_terms=(Ccorr..., Clim...), finalization=finalization,
        Δχ=integ.Δχ, pref=pref,
    )
    println("starting 17 endpoint W pullbacks")
    cphi_2d = reshape(reactant_coefficients[3], :, 1) .* ones(1, size(w[1], 2))
    T_tuple = map(T -> permutedims(T, (4, 1, 2, 3)), (
        Blast.T_tildes.T_2_00, Blast.T_tildes.T_minus2_00,
        Blast.T_tildes.T_0_00, Blast.T_tildes.T_0_02,
        Blast.T_tildes.T_0_20, Blast.T_tildes.T_2_02,
        Blast.T_tildes.T_2_20, Blast.T_tildes.T_2_22,
    ))
    dpk = staged_nonlimber_pullback(
        workload.pk, bg, Mnonlimber,
        (reactant_coefficients[1], reactant_coefficients[2], cphi_2d),
        w, T_tuple, endpoint_config, dCnon,
    )
    println("dpk_shape=", size(dpk), " dpk_norm=", maximum(abs, dpk))

    # Limber forward intermediates are kept explicitly so the reverse driver
    # can traverse the split stages without rebuilding one monolithic graph.
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
    loglin, lognon = Blast.reactant_limber_power_products(
        workload.pk_limber_lin, workload.pk_limber_nonlin, bg,
    )
    clin = Blast.reactant_limber_chebyshev_coefficients(loglin, Mk, Mz)
    cnon = Blast.reactant_limber_chebyshev_coefficients(lognon, Mk, Mz)
    Plin = Blast.reactant_limber_grid_from_coefficients(clin, Tz, Tk)
    Pnon = Blast.reactant_limber_grid_from_coefficients(cnon, Tz, Tk)
    limber_state = LimberPullbackState(
        loglin, lognon, clin, cnon, Plin, Pnon,
        nothing, nothing, nothing, nothing,
    )
    invχ2 = vec(Blast.LIMBER_INV_χ2_ROW)
    KG_low, KG_high = Blast._low_ℓ_slice(KG), Blast._high_ℓ_slice(KG)
    KL_low, KL_high = Blast._low_ℓ_slice(KL), Blast._high_ℓ_slice(KL)
    dpk_lim_lin = zeros(size(workload.pk_limber_lin))
    dpk_lim_non = zeros(size(workload.pk_limber_nonlin))
    for (dCc, dCl, k1Low, k2Low, k1High, k2High) in (
        (dCcorr[1], dClim[1], KG_low, KG_low, KG_high, KG_high),
        (dCcorr[2], dClim[2], KG_low, KL_low, KG_high, KL_high),
        (dCcorr[3], dClim[3], KL_low, KL_low, KL_high, KL_high),
    )
        gl, gn = staged_limber_gradient(
            workload.pk_limber_lin, workload.pk_limber_nonlin, bg, limber_state,
            dCc, dCl, k1Low, k2Low, k1High, k2High, Blast.LIMBER_WEIGHTS, invχ2,
            Blast.LIMBER_Δχ, Mk, Mz, Tz, Tk,
        )
        dpk_lim_lin .+= gl
        dpk_lim_non .+= gn
    end
    println("dpk_limber_lin_shape=", size(dpk_lim_lin), " norm=", maximum(abs, dpk_lim_lin))
    println("dpk_limber_nonlin_shape=", size(dpk_lim_non), " norm=", maximum(abs, dpk_lim_non))

    Random.seed!(97531)
    rpk = randn(size(workload.pk)); rpk ./= maximum(abs, rpk)
    rlin = randn(size(workload.pk_limber_lin)); rlin ./= maximum(abs, rlin)
    rnon = randn(size(workload.pk_limber_nonlin)); rnon ./= maximum(abs, rnon)
    vpk = workload.pk .* rpk
    vlin = workload.pk_limber_lin .* rlin
    vnon = workload.pk_limber_nonlin .* rnon
    scalar_loss(x) = sum(x.cl_gg) + sum(x.cl_gs) + sum(x.cl_ss)
    ad_directional = sum(dpk .* vpk) + sum(dpk_lim_lin .* vlin) + sum(dpk_lim_non .* vnon)
    println("directional_ad=", ad_directional)
    for ϵ in (1e-3, 1e-4, 1e-5, 1e-6)
        plus = merge(workload, (; pk=workload.pk .+ ϵ .* vpk,
                                  pk_limber_lin=workload.pk_limber_lin .+ ϵ .* vlin,
                                  pk_limber_nonlin=workload.pk_limber_nonlin .+ ϵ .* vnon))
        minus = merge(workload, (; pk=workload.pk .- ϵ .* vpk,
                                   pk_limber_lin=workload.pk_limber_lin .- ϵ .* vlin,
                                   pk_limber_nonlin=workload.pk_limber_nonlin .- ϵ .* vnon))
        fd_directional = (scalar_loss(full_power_to_cls(plus; pass_plans=true)) -
                          scalar_loss(full_power_to_cls(minus; pass_plans=true))) / (2ϵ)
        relative_error = abs(ad_directional - fd_directional) / max(abs(fd_directional), eps())
        println("directional_epsilon=", ϵ,
                " fd=", fd_directional,
                " ad=", ad_directional,
                " relative_error=", relative_error)
    end

    # All pullback handles are cached by staged_reactant_gradient.jl at this
    # point. These timings therefore exclude Reactant/Enzyme compilation.
    nonlimber_gradient_call = () -> staged_nonlimber_pullback(
        workload.pk, bg, Mnonlimber,
        (reactant_coefficients[1], reactant_coefficients[2], cphi_2d),
        w, T_tuple, endpoint_config, dCnon,
    )
    limber_gradient_call = () -> begin
        gl = zeros(size(workload.pk_limber_lin)); gn = zeros(size(workload.pk_limber_nonlin))
        for (dCc, dCl, k1Low, k2Low, k1High, k2High) in (
            (dCcorr[1], dClim[1], KG_low, KG_low, KG_high, KG_high),
            (dCcorr[2], dClim[2], KG_low, KL_low, KG_high, KL_high),
            (dCcorr[3], dClim[3], KL_low, KL_low, KL_high, KL_high),
        )
            a, b = staged_limber_gradient(
                workload.pk_limber_lin, workload.pk_limber_nonlin, bg, limber_state,
                dCc, dCl, k1Low, k2Low, k1High, k2High,
                Blast.LIMBER_WEIGHTS, invχ2, Blast.LIMBER_Δχ, Mk, Mz, Tz, Tk,
            )
            gl .+= a; gn .+= b
        end
        return gl, gn
    end
    complete_gradient_call = () -> begin
        dN, dC, dL = staged_finalization_pullback(
            (nonlimber.gg, nonlimber.gs, nonlimber.ss), Ccorr, Clim,
            finalization, upstream,
        )
        gpk = staged_nonlimber_pullback(
            workload.pk, bg, Mnonlimber,
            (reactant_coefficients[1], reactant_coefficients[2], cphi_2d),
            w, T_tuple, endpoint_config, dN,
        )
        gl = zeros(size(workload.pk_limber_lin))
        gn = zeros(size(workload.pk_limber_nonlin))
        for (dCc, dCl, k1Low, k2Low, k1High, k2High) in (
            (dC[1], dL[1], KG_low, KG_low, KG_high, KG_high),
            (dC[2], dL[2], KG_low, KL_low, KG_high, KL_high),
            (dC[3], dL[3], KL_low, KL_low, KL_high, KL_high),
        )
            a, b = staged_limber_gradient(
                workload.pk_limber_lin, workload.pk_limber_nonlin, bg, limber_state,
                dCc, dCl, k1Low, k2Low, k1High, k2High,
                Blast.LIMBER_WEIGHTS, invχ2, Blast.LIMBER_Δχ, Mk, Mz, Tz, Tk,
            )
            gl .+= a
            gn .+= b
        end
        return gpk, gl, gn
    end
    nonlimber_bench = @benchmark $nonlimber_gradient_call() samples=3 evals=1
    limber_bench = @benchmark $limber_gradient_call() samples=3 evals=1
    complete_bench = @benchmark $complete_gradient_call() samples=3 evals=1
    println("nonlimber_gradient_steady_min_ms=", minimum(nonlimber_bench).time / 1e6)
    println("nonlimber_gradient_steady_median_ms=", median(nonlimber_bench).time / 1e6)
    println("limber_gradient_steady_min_ms=", minimum(limber_bench).time / 1e6)
    println("limber_gradient_steady_median_ms=", median(limber_bench).time / 1e6)
    println("complete_gradient_steady_min_ms=", minimum(complete_bench).time / 1e6)
    println("complete_gradient_steady_median_ms=", median(complete_bench).time / 1e6)
end

main()

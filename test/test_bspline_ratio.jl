# Phase 6: public FiducialRatio spline path (manufactured production-geometry banks, fixed analytic q).
using Test, Blast, LinearAlgebra, ForwardDiff
using AbstractCosmologicalEmulators, FastChebInterp
isdefined(@__MODULE__, :get_test_cosmo) || include("test_fixtures.jl")
isdefined(@__MODULE__, :build_characterization_state) || include("characterization_workload.jl")
isdefined(@__MODULE__, :manufactured_ratio_bank) || include("bspline_backend_fixture.jl")

_family(src) = src === :ϕTT ? :cϕTT : src === :ϕT ? :cϕT : :cϕ
ratio_array(bank, tfield, src) = src === :ϕTT ? getfield(bank.TT, tfield) : src === :ϕT ? bank.T[tfield] : bank.PHI

@testset "Ratio backend: policies, validation and source-aware components" begin
    rbank = manufactured_ratio_bank(); dbank = manufactured_direct_bank()
    backend = Blast.CubicBSplineBackend(rbank)
    @test backend.policy isa Blast.FiducialRatio && Blast.CubicBSplineBackend(dbank).policy isa Blast.DirectExpansion
    @test Blast.nonlimber_k(backend) === rbank.grid.sites && Blast.bank_signature(backend) == rbank.signature
    @test_throws ArgumentError Blast.CubicBSplineBackend(dbank; policy=Blast.FiducialRatio())
    @test_throws ArgumentError Blast.CubicBSplineBackend(rbank; policy=Blast.DirectExpansion())
    badq = Blast.FiducialNormalizers(copy(rbank.q.qTT), [rbank.q.qT[1:end-1]; -1.0], copy(rbank.q.qPHI))
    @test_throws ArgumentError Blast.CubicBSplineBackend(Blast.RatioSplineBank(rbank.TT, rbank.T, rbank.PHI, rbank.grid, badq, rbank.signature))
    short_phi = rbank.PHI[:, :, :, 1:end-1]
    @test_throws ArgumentError Blast.CubicBSplineBackend(Blast.RatioSplineBank(rbank.TT, rbank.T, short_phi, rbank.grid, rbank.q, rbank.signature))

    gc = Blast.GalaxyClustering(δ=Blast.NumberCounts(), RSD=Blast.RedshiftSpaceDistortions(), μ=Blast.MagnificationBias(), PNG=Blast.PrimordialNonGaussianity())
    wl = Blast.WeakLensing(γ=Blast.CosmicShear(), IA=Blast.IntrinsicAlignment()); cmb = Blast.CMB(κ=Blast.CMBLensing(), ISW=Blast.IntegratedSachsWolfe())
    W, plans = Blast.SetUp(gc, wl, cmb; backend)
    @test length(fieldnames(typeof(W))) == 17 && all(!isnothing, getfield.(Ref(W), fieldnames(typeof(W))))
    for (name, tfield, src, r1) in Blast._PROJECTED_MATTER_COMPONENTS
        @test getfield(W, name).T̃ === ratio_array(rbank, tfield, src)                  # actual family arrays, no dummies
    end
    @test W.w_2_02_ϕT_R1.T̃ === rbank.T.T_2_02 && W.w_2_02_ϕTT.T̃ === rbank.TT.T_2_02 && W.w_2_00_ϕ.T̃ === rbank.PHI
    @test Blast.bank_signature(W) == Blast.bank_signature(plans) == rbank.signature != Blast.bank_signature(dbank)
    # Every wrapper of SetUp works with the ratio backend and activates the same components as Chebyshev.
    for probes in ((gc,), (wl,), (gc, wl), (wl, gc), (gc, cmb), (cmb, gc), (wl, cmb), (cmb, wl), (gc, wl, cmb), (cmb, wl, gc))
        Wd, _ = Blast.SetUp(probes...); Wr, Pr = Blast.SetUp(probes...; backend)
        @test Tuple(n for n in fieldnames(typeof(Wr)) if !isnothing(getfield(Wr, n))) == Tuple(n for n in fieldnames(typeof(Wd)) if !isnothing(getfield(Wd, n)))
        @test Pr isa Blast.SpectralPlans{<:Blast.CubicBSplineBackend}
    end
end

@testset "Ratio backend: normalization of COMPLETE sources, dense oracle, 17 blocks, spectra, mixing" begin
    state = build_characterization_state()
    (; bg, probes, pmm, ell) = state
    gc, wl, cmb = probes.galaxy, probes.lensing, probes.cmb
    rbank = manufactured_ratio_bank(); backend = Blast.CubicBSplineBackend(rbank)
    qcopy = (copy(rbank.q.qTT), copy(rbank.q.qT), copy(rbank.q.qPHI))
    W, plans = Blast.SetUp(gc, wl, cmb; backend)
    grid = Blast.nonlimber_grid(plans); pk = fixture_pk(grid); pk0 = copy(pk)
    pkl, pkn = pmm.pk_limber_linear, pmm.pk_limber_nonlinear
    Pk = Blast.prepare_pk_workspace(plans, pk, pkl, pkn, bg; pk_grid=grid)
    @test pk == pk0                                                        # Pmm inputs are never normalized

    # Independent oracle: complete sources / q along k, then a dense collocation solve (not build_coeff).
    n = length(grid.k)
    A = basis_matrix(bspline_basis(CubicBSplinePlan(rbank.grid.logsites, rbank.grid.logsites; extrapolation=:throw)), rbank.grid.logsites)
    Pphi = Blast.get_PΦ(grid.k, bg); Tm = Blast.get_Tm(pk, grid.k, bg)
    TR = Blast.transform_to_R_frame(Tm, bg); T1 = TR[:, end, :]
    P_TT = Blast._p_phi_TT_tullio(Pphi, T1, TR); P_T = Blast._p_phi_T_tullio(Pphi, TR)
    divk(V, q) = [V[i, j, k] / q[i] for i in axes(V, 1), j in axes(V, 2), k in axes(V, 3)]
    dense(V) = reshape(A \ reshape(V, n, :), size(V))
    c_TT = dense(divk(P_TT, rbank.q.qTT)); c_T = dense(divk(P_T, rbank.q.qT)); c_ϕ = A \ (Pphi ./ rbank.q.qPHI)
    @test norm(Pk.cϕTT.coefs - c_TT) ≤ 1e-9 * norm(c_TT)
    @test norm(Pk.cϕT.coefs - c_T) ≤ 1e-9 * norm(c_T)
    @test norm(Pk.cϕ.coefs - c_ϕ) ≤ 1e-9 * norm(c_ϕ)
    @test norm(Pk.cϕTT.coefs - dense(P_TT)) > 1e-3 * norm(dense(P_TT))      # genuinely normalized, not the direct expansion
    # R1 commutes with the k-only normalization.
    @test divk(P_T, rbank.q.qT)[:, :, end] == P_T[:, :, end] ./ rbank.q.qT
    @test Pk.cϕT.coefs[:, :, end] ≈ dense(P_T[:, :, end:end] ./ reshape(rbank.q.qT, :, 1, 1))[:, :, 1] rtol=1e-9
    # q is fixed saved data: unchanged bytes/identity across different cosmologies.
    bg2 = Blast.Background(get_test_cosmo(h=0.70, nₛ=0.97))
    Blast.prepare_pk_workspace(plans, pk, pkl, pkn, bg2)
    @test (rbank.q.qTT, rbank.q.qT, rbank.q.qPHI) == qcopy && backend.bank.q === rbank.q

    # All 17 blocks vs literal contraction with the ACTUAL family arrays (R1: T family, sliced source coefficients).
    Wf = Blast.compute_w(W, Pk)
    refs = Dict(:cϕTT => c_TT, :cϕT => c_T, :cϕ => c_ϕ)
    for (name, tfield, src, r1) in Blast._PROJECTED_MATTER_COMPONENTS
        c = refs[_family(src)]; c = r1 ? c[:, :, end] : c
        expected = literal_w(c, ratio_array(rbank, tfield, src))
        @test norm(getfield(Wf, name).w - expected) ≤ 1e-9 * norm(expected)
    end
    @test !Blast._shares_tensor(W.w_2_00_ϕTT, W.w_2_00_ϕT, W.w_2_00_ϕT_R1, W.w_2_00_ϕ)
    Pk_B = Blast.prepare_pk_workspace(plans, 1.07 .* pk, 1.02 .* pkl, 0.96 .* pkn, bg)
    ws = Blast.allocate_compute_w(W, Pk)
    for P in (Pk, Pk_B, Pk, Pk_B)
        for f in fieldnames(typeof(ws)); fill!(getfield(ws, f).w, NaN); end
        Blast.compute_w!(ws, P); Wref = Blast.compute_w(W, P)
        for f in fieldnames(typeof(ws)); @test getfield(ws, f).w ≈ getfield(Wref, f).w rtol=1e-12; end
    end

    # Public spectra, all wrappers, plans and P = nothing.
    function spectra(P)
        (Blast.get_Cℓ(ell, gc, Pk, Wf, bg, P), Blast.get_Cℓ(ell, gc, wl, Pk, Wf, bg, P), Blast.get_Cℓ(ell, wl, gc, Pk, Wf, bg, P),
         Blast.get_Cℓ(ell, gc, cmb, Pk, Wf, bg, P), Blast.get_Cℓ(ell, cmb, gc, Pk, Wf, bg, P), Blast.get_Cℓ(ell, wl, Pk, Wf, bg, P),
         Blast.get_Cℓ(ell, wl, cmb, Pk, Wf, bg, P), Blast.get_Cℓ(ell, cmb, wl, Pk, Wf, bg, P))
    end
    for (a, b) in zip(spectra(plans), spectra(nothing)); @test all(isfinite, a) && isapprox(a, b; rtol=1e-8); end

    # Signature mixing is rejected before any contraction.
    dbank = manufactured_direct_bank(); dback = Blast.CubicBSplineBackend(dbank)
    Wd, plans_d = Blast.SetUp(gc, wl, cmb; backend=dback)
    Pk_d = Blast.prepare_pk_workspace(plans_d, fixture_pk(Blast.nonlimber_grid(plans_d)), pkl, pkn, bg)
    @test_throws ArgumentError Blast.compute_w(W, Pk_d)
    @test_throws ArgumentError Blast.compute_w(Wd, Pk)
    @test_throws ArgumentError Blast.compute_w!(ws, Pk_d)
    @test_throws ArgumentError Blast.get_Cℓ(ell, gc, wl, Pk_d, Wf, bg, plans)
    @test_throws ArgumentError Blast.get_Cℓ(ell, gc, Pk, Wf, bg, plans_d)
    @test_throws ArgumentError Blast.compute_w(W, prepare_characterization_spectrum(state))   # Chebyshev Pk
    rbank2 = manufactured_ratio_bank(n=8, offset=0.4)
    W2, plans2 = Blast.SetUp(gc, wl, cmb; backend=Blast.CubicBSplineBackend(rbank2))
    @test_throws ArgumentError Blast.compute_w(W2, Pk)
    @test_throws ArgumentError Blast.prepare_pk_workspace(plans, pk, pkl, pkn, bg; pk_grid=Blast.nonlimber_grid(plans2))
    @test_throws DimensionMismatch Blast.prepare_pk_workspace(plans, fixture_pk(Blast.nonlimber_grid(plans_d))[:, 1:end-1], pkl, pkn, bg)

    # Mixed scalar (fixed Float64 background, Dual Pmm) works for the ratio path, with zero-partial constant Phi.
    D = ForwardDiff.Dual{Nothing,Float64,1}
    pkd = [D(v, ForwardDiff.Partials{1,Float64}((i == 1 ? v : 0.0,))) for (i, v) in enumerate(vec(pk))] |> x -> reshape(x, size(pk))
    Pd = Blast.prepare_pk_workspace(plans, pkd, pkl, pkn, bg)
    @test Pd isa Blast.PowerSpectrum{D} && ForwardDiff.value.(Pd.cϕ.coefs) == Pk.cϕ.coefs
    @test all(x -> all(iszero, ForwardDiff.partials(x)), Pd.cϕ.coefs)
end

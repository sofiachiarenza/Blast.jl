# Phase 5: direct public cubic B-spline path (manufactured production-geometry banks).
using Test, Blast, LinearAlgebra
using AbstractCosmologicalEmulators, FastChebInterp      # same imports as the full suite (name collisions)
isdefined(@__MODULE__, :get_test_cosmo) || include("test_fixtures.jl")
isdefined(@__MODULE__, :build_characterization_state) || include("characterization_workload.jl")
isdefined(@__MODULE__, :manufactured_direct_bank) || include("bspline_backend_fixture.jl")

_active(W) = Tuple(n for n in fieldnames(typeof(W)) if !isnothing(getfield(W, n)))
_family(src) = src === :ϕTT ? :cϕTT : src === :ϕT ? :cϕT : :cϕ

@testset "B-spline backend: construction, grid contract and plans" begin
    bank = manufactured_direct_bank()
    backend = Blast.CubicBSplineBackend(bank)
    @test Blast.nonlimber_k(backend) === bank.grid.sites                  # exact stored physical sites
    @test Blast.nonlimber_sites(backend) === bank.grid.logsites
    @test Blast.bank_signature(backend) == Blast.bank_signature(bank) == bank.tensors.signature
    g = Blast.GalaxyClustering(δ=Blast.NumberCounts())
    W, plans = Blast.SetUp(g; backend)
    @test plans isa Blast.SpectralPlans{<:Blast.CubicBSplineBackend}
    @test Blast.nonlimber_grid(plans) == (z = Blast.z_lin, k = bank.grid.sites)
    @test W.w_2_00_ϕTT.T̃ === bank.tensors.T_2_00 && Blast.bank_signature(W) == Blast.bank_signature(backend)
    # Source plans are private to each setup; the bank is shared read-only.
    W2, plans2 = Blast.SetUp(g; backend)
    @test plans2.plan_ϕTT !== plans.plan_ϕTT && W2.w_2_00_ϕTT.T̃ === W.w_2_00_ϕTT.T̃
    # Legacy six-argument plan constructor never accepts spline plans (no method for non-Chebyshev plans).
    @test_throws MethodError Blast.FFTPlans(plans.plan_ϕTT, plans.plan_ϕT, plans.plan_ϕ, plans.plan_limber, plans.T_k_limber, plans.plan_ℓ)
    # Backend validation: wrong basis, inconsistent sites, perturbed knots.
    cheb_bank = Blast.legacy_tensor_bank()
    @test_throws MethodError Blast.CubicBSplineBackend(cheb_bank)
    tens = NamedTuple{Blast._BANK_GEOMETRIES}(Tuple(getfield(bank.tensors, g) for g in Blast._BANK_GEOMETRIES))
    wrong_basis = Blast.DirectSplineBank(Blast.TensorBank(tens, bank.grid.logsites; basis=:chebyshev, dataset_id="private:x"), bank.grid)
    @test_throws ArgumentError Blast.CubicBSplineBackend(wrong_basis)
    bad_knots = copy(bank.grid.knots); bad_knots[6] = nextfloat(bad_knots[6])
    @test_throws ArgumentError Blast.CubicBSplineBackend(Blast.DirectSplineBank(bank.tensors, Blast.SplineGrid(bank.grid.sites, bank.grid.logsites, bad_knots)))
    @test_throws ArgumentError Blast.CubicBSplineBackend(Blast.DirectSplineBank(bank.tensors,
        Blast.SplineGrid(bank.grid.sites .* (1 + 1e-9), bank.grid.logsites, bank.grid.knots)))
end

@testset "B-spline backend: public source preparation and all 17 blocks against dense/literal references" begin
    state = build_characterization_state()
    (; bg, probes, pmm, ell) = state
    bank = manufactured_direct_bank(); backend = Blast.CubicBSplineBackend(bank)
    gc, wl, cmb = probes.galaxy, probes.lensing, probes.cmb
    W, plans = Blast.SetUp(gc, wl, cmb; backend)
    grid = Blast.nonlimber_grid(plans)
    pk = fixture_pk(grid); @test all(>(0), pk)
    pkl, pkn = pmm.pk_limber_linear, pmm.pk_limber_nonlinear
    Pk = Blast.prepare_pk_workspace(plans, pk, pkl, pkn, bg; pk_grid=grid)
    @test Pk.signature == Blast.bank_signature(backend)

    # Independent coefficients: dense collocation solve of the SAME sampled sources (existing physics helpers).
    n = length(bank.grid.sites); A = basis_matrix(bspline_basis(CubicBSplinePlan(bank.grid.logsites, bank.grid.logsites; extrapolation=:throw)), bank.grid.logsites)
    Pphi = Blast.get_PΦ(grid.k, bg); Tm = Blast.get_Tm(pk, grid.k, bg)
    TR = Blast.transform_to_R_frame(Tm, bg); T1 = TR[:, end, :]
    P_TT = Blast._p_phi_TT_tullio(Pphi, T1, TR); P_T = Blast._p_phi_T_tullio(Pphi, TR)
    dense(V) = reshape(A \ reshape(V, n, :), size(V))
    c_TT, c_T, c_ϕ = dense(P_TT), dense(P_T), A \ Pphi
    @test size(Pk.cϕTT.coefs) == (n, length(Blast.χ), length(Blast.R))                 # (basis, χ, R)
    @test norm(Pk.cϕTT.coefs - c_TT) ≤ 1e-9 * norm(c_TT)
    @test norm(Pk.cϕT.coefs - c_T) ≤ 1e-9 * norm(c_T)
    @test norm(Pk.cϕ.coefs - c_ϕ) ≤ 1e-9 * norm(c_ϕ)
    @test norm(Pk.cϕTT.coefs - P_TT) > 1e-3 * norm(P_TT)                                # coefficients are NOT sample values

    # All 17 projected blocks vs literal contraction of the dense coefficients with the manufactured tensors.
    Wf = Blast.compute_w(W, Pk)
    @test length(_active(Wf)) == 17
    ref_coef = Dict(:cϕTT => c_TT, :cϕT => c_T, :cϕ => c_ϕ)
    for (name, tfield, src, r1) in Blast._PROJECTED_MATTER_COMPONENTS
        comp = getfield(Wf, name); T = getfield(bank.tensors, tfield)
        @test comp.T̃ === T                                                              # shared bank tensor, aliases preserved
        c = ref_coef[_family(src)]; c = r1 ? c[:, :, end] : c
        expected = literal_w(c, T)
        @test norm(comp.w - expected) ≤ 1e-9 * norm(expected)
    end
    @test getfield(Wf, :w_2_00_ϕTT).T̃ === getfield(Wf, :w_2_00_ϕT).T̃ === getfield(Wf, :w_2_00_ϕT_R1).T̃ === getfield(Wf, :w_2_00_ϕ).T̃

    # Functional vs mutating, repeated alternating inputs, independent workspaces.
    Pk_B = Blast.prepare_pk_workspace(plans, 1.07 .* pk, 1.02 .* pkl, 0.96 .* pkn, bg)
    ws = Blast.allocate_compute_w(W, Pk); ws2 = Blast.allocate_compute_w(W, Pk)
    @test ws.w_2_00_ϕTT.w !== ws2.w_2_00_ϕTT.w
    for P in (Pk, Pk_B, Pk, Pk_B, Pk)
        Blast.compute_w!(ws, P); Wref = Blast.compute_w(W, P)
        for f in _active(Wref); @test getfield(ws, f).w ≈ getfield(Wref, f).w rtol=1e-12; end
    end

    # Public spectra, every wrapper, P = nothing too; Chebyshev and two spline setups interleaved.
    function spectra(Pk_, W_, plans_)
        P = plans_
        (Blast.get_Cℓ(ell, gc, Pk_, W_, bg, P), Blast.get_Cℓ(ell, gc, wl, Pk_, W_, bg, P), Blast.get_Cℓ(ell, wl, gc, Pk_, W_, bg, P),
         Blast.get_Cℓ(ell, gc, cmb, Pk_, W_, bg, P), Blast.get_Cℓ(ell, cmb, gc, Pk_, W_, bg, P), Blast.get_Cℓ(ell, wl, Pk_, W_, bg, P),
         Blast.get_Cℓ(ell, wl, cmb, Pk_, W_, bg, P), Blast.get_Cℓ(ell, cmb, wl, Pk_, W_, bg, P))
    end
    cl_plans = spectra(Pk, Wf, plans); cl_none = spectra(Pk, Wf, nothing)
    for (a, b) in zip(cl_plans, cl_none); @test all(isfinite, a) && isapprox(a, b; rtol=1e-8); end
    cheb = run_characterization(state)                                                   # Chebyshev default
    bank2 = manufactured_direct_bank(n=8, offset=0.5); backend2 = Blast.CubicBSplineBackend(bank2)
    W2, plans2 = Blast.SetUp(gc, wl, cmb; backend=backend2)
    Pk2 = Blast.prepare_pk_workspace(plans2, fixture_pk(Blast.nonlimber_grid(plans2)), pkl, pkn, bg)
    Wf2 = Blast.compute_w(W2, Pk2)
    cl_plans_again = spectra(Pk, Wf, plans); cl2 = spectra(Pk2, Wf2, plans2)
    cheb_again = run_characterization(state)
    @test isapprox(cheb_again.spectra.cl_gg, cheb.spectra.cl_gg; rtol=1e-12)
    @test all(isapprox(a, b; rtol=1e-12) for (a, b) in zip(cl_plans, cl_plans_again))
    @test all(c -> all(isfinite, c), cl2)

    # Mismatched Pk / W / plans fail before any contraction.
    @test_throws ArgumentError Blast.compute_w(W, Pk2)                                   # other bank
    @test_throws ArgumentError Blast.compute_w(state.weights_template, Pk)               # Chebyshev W, spline Pk
    @test_throws ArgumentError Blast.compute_w(W, cheb.spectrum)                         # spline W, Chebyshev Pk
    @test_throws ArgumentError Blast.allocate_compute_w(W, Pk2)
    @test_throws ArgumentError Blast.compute_w!(ws, Pk2)
    @test_throws ArgumentError Blast.get_Cℓ(ell, gc, wl, Pk2, Wf, bg, plans)
    @test_throws ArgumentError Blast.get_Cℓ(ell, gc, Pk, Wf, bg, plans2)
    @test_throws ArgumentError Blast.get_Cℓ(ell, gc, Pk, Wf, bg, state.plans)

    # Grid contract: a plain array can only be shape-checked; supplied coordinates must match bit for bit.
    @test_throws DimensionMismatch Blast.prepare_pk_workspace(plans, pmm.pk_nonlimber, pkl, pkn, bg)    # Chebyshev-length Pmm
    @test_throws ArgumentError Blast.prepare_pk_workspace(plans, pk, pkl, pkn, bg; pk_grid=Blast.nonlimber_grid(state.plans))
    @test_throws ArgumentError Blast.prepare_pk_workspace(plans, pk, pkl, pkn, bg; pk_grid=(z=grid.z, k=Blast.nonlimber_grid(plans2).k))
    moved = copy(grid.k); moved[3] = nextfloat(moved[3])
    @test_throws ArgumentError Blast.prepare_pk_workspace(plans, pk, pkl, pkn, bg; pk_grid=(z=grid.z, k=moved))
    @test_throws DomainError Blast.prepare_pk_workspace(plans, -pk, pkl, pkn, bg)       # physical support enforced, no clamping
end

@testset "B-spline backend: every SetUp wrapper, ARF + PNG" begin
    state = build_characterization_state()
    gc, wl, cmb = state.probes.galaxy, state.probes.lensing, state.probes.cmb
    bank = manufactured_direct_bank(); backend = Blast.CubicBSplineBackend(bank)
    combos = [(gc,), (wl,), (gc, wl), (wl, gc), (gc, cmb), (cmb, gc), (wl, cmb), (cmb, wl), (gc, wl, cmb), (gc, cmb, wl),
              (wl, gc, cmb), (wl, cmb, gc), (cmb, gc, wl), (cmb, wl, gc)]
    for probes in combos
        Wd, Pd = Blast.SetUp(probes...)                                     # Chebyshev default
        Ws, Ps = Blast.SetUp(probes...; backend)
        @test _active(Wd) == _active(Ws)
        @test Ps isa Blast.SpectralPlans{<:Blast.CubicBSplineBackend}
        for f in _active(Ws)
            tfield = Blast._PROJECTED_MATTER_COMPONENTS[findfirst(c -> c[1] == f, Blast._PROJECTED_MATTER_COMPONENTS)][2]
            @test getfield(Ws, f).T̃ === getfield(bank.tensors, tfield)
        end
        @test Blast.bank_signature(Ws) == Blast.bank_signature(Ps) == bank.tensors.signature
    end
    # ARF + PNG gg through the spline path.
    g = gc
    arf = Blast.AngularRedshiftFluctuations(nz=g.δ.nz, z=g.δ.z, bias=copy(g.δ.bias), zbar=[0.8, 1.2])
    with_arf = Blast.prepare_probe(Blast.GalaxyClustering(δ=g.δ, RSD=g.RSD, μ=g.μ, PNG=g.PNG, ARF=arf), state.bg)
    W, plans = Blast.SetUp(with_arf; backend)
    pk = fixture_pk(Blast.nonlimber_grid(plans))
    Pk = Blast.prepare_pk_workspace(plans, pk, state.pmm.pk_limber_linear, state.pmm.pk_limber_nonlinear, state.bg)
    Wf = Blast.compute_w(W, Pk)
    cl_arf = Blast.get_Cℓ(state.ell, with_arf, Pk, Wf, state.bg, plans)
    W0, plans0 = Blast.SetUp(gc; backend)
    Wf0 = Blast.compute_w(W0, Blast.prepare_pk_workspace(plans0, pk, state.pmm.pk_limber_linear, state.pmm.pk_limber_nonlinear, state.bg))
    cl_no_arf = Blast.get_Cℓ(state.ell, gc, Blast.prepare_pk_workspace(plans0, pk, state.pmm.pk_limber_linear, state.pmm.pk_limber_nonlinear, state.bg), Wf0, state.bg, plans0)
    @test all(isfinite, cl_arf) && all(isfinite, cl_no_arf) && cl_arf != cl_no_arf
end

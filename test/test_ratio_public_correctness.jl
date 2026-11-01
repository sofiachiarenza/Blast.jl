# Phase 8-style public correctness of the ratio spline path on MANUFACTURED production-geometry banks
# (algebra/API evidence only; no accuracy, default-count or performance claim). Transferred from the external
# real-N033 checks. Selector: `ratio-public`. Helpers are prefixed `rpc_` to avoid collisions.
using Test, Blast, LinearAlgebra
using AbstractCosmologicalEmulators, FastChebInterp      # same imports as the full suite
isdefined(@__MODULE__, :build_characterization_state) || include("characterization_workload.jl")
isdefined(@__MODULE__, :manufactured_ratio_bank) || include("bspline_backend_fixture.jl")

rpc_relnorm(a, b) = norm(a .- b) / max(norm(b), floatmin(Float64))
rpc_maxrel(a, b) = maximum(abs.(a .- b)) / max(maximum(abs.(b)), floatmin(Float64))
rpc_nonlinear_boost(k, z) = 1 + 0.45 * (k / (0.35 + k))^2 / (1 + z)
function rpc_limber()
    kl = 10 .^ Blast.k_limber
    return [fixture_pmm(k, z) for z in Blast.z_cheb, k in kl], [fixture_pmm(k, z) * rpc_nonlinear_boost(k, z) for z in Blast.z_cheb, k in kl]
end

@testset "Ratio public: fiducial-separable factorization (manufactured bank, exact fiducial primordial)" begin
    bank = manufactured_ratio_bank(); backend = Blast.CubicBSplineBackend(bank)
    c0 = characterization_cosmology()
    # The manufactured q uses ln10_10_As = 3.044, n_s = 0.9649; the characterization primordial (n_s = 0.9645) would not match.
    cosmo = Blast.w0waCDMCosmology(; ln10Aₛ = 3.044, nₛ = 0.9649, h = c0.h, ωb = c0.ωb, ωc = c0.ωc, ωk = c0.ωk, mν = c0.mν, w0 = c0.w0, wa = c0.wa)
    bg = Blast.Background(cosmo)
    probes = build_characterization_probes(bg)
    W, plans = Blast.SetUp(probes.galaxy, probes.lensing, probes.cmb; backend)
    grid = Blast.nonlimber_grid(plans); nk = length(grid.k); q = bank.q
    @test isapprox(Blast.get_PΦ(grid.k, bg), q.qPHI; rtol = 1e-13)
    d(z) = 1 / (1 + z)
    pk = [q.qTT[k] * d(z)^2 for z in Blast.z_lin, k in 1:nk]
    pkl, pkn = rpc_limber()
    q0 = (copy(q.qTT), copy(q.qT), copy(q.qPHI))
    Pk = Blast.prepare_pk_workspace(plans, pk, pkl, pkn, bg; pk_grid = grid)
    @test (q.qTT, q.qT, q.qPHI) == q0                                   # fixed q: unchanged bytes
    # Independent growth through the existing R-frame helper: every k column equals d(z).
    dTR = Blast.transform_to_R_frame([d(z) for z in Blast.z_lin, _ in 1:nk], bg)      # (chi, R, k)
    dR, d1 = dTR[:, :, 1], dTR[:, end, 1]
    TT_const = d1 .* dR
    A = basis_matrix(bspline_basis(CubicBSplinePlan(bank.grid.logsites, bank.grid.logsites; extrapolation = :throw)), bank.grid.logsites)
    @test isfinite(cond(A))                                             # conditioning is recorded, not asserted small
    TOL = 1e-10
    for (fam, cref) in ((:cϕTT, TT_const), (:cϕT, dR))
        coefs = getfield(Pk, fam).coefs
        @test maximum(maximum(abs.(coefs[l, :, :] .- cref)) for l in axes(coefs, 1)) ≤ TOL * maximum(abs.(cref))   # constant along the basis axis
    end
    @test maximum(abs.(Pk.cϕ.coefs .- 1.0)) ≤ TOL
    # All 17 raw weights, entire arrays: growth scalar x sum over the basis axis of the ACTUAL family tensor.
    Wf = Blast.compute_w(W, Pk)
    for (name, tfield, src, r1) in Blast._PROJECTED_MATTER_COMPONENTS
        T = src === :ϕTT ? getfield(bank.TT, tfield) : src === :ϕT ? bank.T[tfield] : bank.PHI
        S = dropdims(sum(T; dims = 4); dims = 4)
        scale = src === :ϕTT ? reshape(TT_const, 1, size(TT_const)...) :
                src === :ϕT && !r1 ? reshape(dR, 1, size(dR)...) :
                src === :ϕT ? reshape(d1, 1, :, 1) : ones(1, 1, 1)
        expected = scale .* S
        @test maximum(abs, expected) > 0
        @test rpc_maxrel(getfield(Wf, name).w, expected) ≤ TOL
    end
end

@testset "Ratio public: ARF + PNG, GG and GC (b != p, zbar fixed)" begin
    backend = Blast.CubicBSplineBackend(manufactured_ratio_bank())
    bg = Blast.Background(characterization_cosmology())
    probes = build_characterization_probes(bg)
    g, wl, cmb = probes.galaxy, probes.lensing, probes.cmb                 # NC bias >= 1.2, PNG b-p = 0.2, magnification nonzero
    arf = Blast.AngularRedshiftFluctuations(nz = g.δ.nz, z = g.δ.z, bias = copy(g.δ.bias), zbar = [0.8, 1.2])
    g_arf = Blast.prepare_probe(Blast.GalaxyClustering(δ = g.δ, RSD = g.RSD, μ = g.μ, PNG = g.PNG, ARF = arf), bg)
    W, plans = Blast.SetUp(g, wl, cmb; backend); Wa, plans_a = Blast.SetUp(g_arf, wl, cmb; backend)
    active(Wx) = Tuple(n for n in fieldnames(typeof(Wx)) if !isnothing(getfield(Wx, n)))
    @test active(W) == active(Wa) && length(active(W)) == 17
    @test Blast.bank_signature(W) == Blast.bank_signature(Wa) == Blast.bank_signature(plans) == Blast.bank_signature(plans_a)
    grid = Blast.nonlimber_grid(plans); pk = fixture_pk(grid); pkl, pkn = rpc_limber()
    Pk = Blast.prepare_pk_workspace(plans, pk, pkl, pkn, bg; pk_grid = grid)
    Wf = Blast.compute_w(W, Pk); ell = Blast.full_ℓ_range
    gg0, gg1 = Blast.get_Cℓ(ell, g, Pk, Wf, bg, plans), Blast.get_Cℓ(ell, g_arf, Pk, Wf, bg, plans)
    gc0, gc1 = Blast.get_Cℓ(ell, g, cmb, Pk, Wf, bg, plans), Blast.get_Cℓ(ell, g_arf, cmb, Pk, Wf, bg, plans)
    gl0, gl1 = Blast.get_Cℓ(ell, g, wl, Pk, Wf, bg, plans), Blast.get_Cℓ(ell, g_arf, wl, Pk, Wf, bg, plans)
    @test all(isfinite, gg1) && all(isfinite, gc1) && all(isfinite, gl1)
    @test rpc_relnorm(gg1, gg0) > 1e-8 && rpc_relnorm(gc1, gc0) > 1e-8       # ARF has a nonzero effect on GG and GC
    @test isapprox(gl1, gl0; rtol = 1e-12)                                    # ARF does not enter GL
    cmb_off = Blast.prepare_probe(Blast.CMB(κ = Blast.CMBLensing()), bg)      # kappa kept (required); ISW toggled
    gc0_off, gc1_off = Blast.get_Cℓ(ell, g, cmb_off, Pk, Wf, bg, plans), Blast.get_Cℓ(ell, g_arf, cmb_off, Pk, Wf, bg, plans)
    @test rpc_relnorm(gc0, gc0_off) > 1e-10                                   # ISW matters for the baseline ...
    @test isapprox(gc1 .- gc0, gc1_off .- gc0_off; rtol = 1e-9, atol = 1e-12 * maximum(abs, gc1 .- gc0))   # ... not for the ARF increment
    g_np = Blast.prepare_probe(Blast.GalaxyClustering(δ = g.δ, RSD = g.RSD, μ = g.μ), bg)
    g_np_arf = Blast.prepare_probe(Blast.GalaxyClustering(δ = g.δ, RSD = g.RSD, μ = g.μ, ARF = arf), bg)
    Wn, plans_n = Blast.SetUp(g_np_arf, wl, cmb; backend); Wnf = Blast.compute_w(Wn, Pk)
    gg0_np = Blast.get_Cℓ(ell, g_np, Pk, Wnf, bg, plans_n)
    d_png, d_nopng = gg1 .- gg0, Blast.get_Cℓ(ell, g_np_arf, Pk, Wnf, bg, plans_n) .- gg0_np
    @test rpc_relnorm(gg0, gg0_np) > 1e-10 && rpc_relnorm(d_png, d_nopng) > 1e-8      # PNG matters; ARF x PNG term nonzero
    # Reversed wrappers return the same axes and values (not a transpose).
    @test Blast.get_Cℓ(ell, wl, g, Pk, Wf, bg, plans) == gl0
    @test Blast.get_Cℓ(ell, cmb, g, Pk, Wf, bg, plans) == gc0
    @test Blast.get_Cℓ(ell, cmb, g_arf, Pk, Wf, bg, plans) == gc1
end

@testset "Ratio public: interleaved private chains, coexistence with a manufactured Chebyshev bank, spawned tasks" begin
    bg = Blast.Background(characterization_cosmology())
    probes = build_characterization_probes(bg); gc, wl, cmb = probes.galaxy, probes.lensing, probes.cmb
    pkl, pkn = rpc_limber(); ell = Blast.full_ℓ_range
    function make_chain(b)
        W, plans = Blast.SetUp(gc, wl, cmb; backend = b)
        grid = Blast.nonlimber_grid(plans); base = fixture_pk(grid)
        ws = Blast.allocate_compute_w(W, Blast.prepare_pk_workspace(plans, base, pkl, pkn, bg))
        return (; W, plans, grid, base, ws)
    end
    function run_chain(c, amp)
        Pk = Blast.prepare_pk_workspace(c.plans, amp .* c.base, pkl, pkn, bg; pk_grid = c.grid)
        Blast.compute_w!(c.ws, Pk)
        return Blast.get_Cℓ(ell, gc, Pk, c.ws, bg, c.plans)
    end
    rm_bank, rm2_bank, dm_bank = manufactured_ratio_bank(), manufactured_ratio_bank(n = 8, offset = 0.4), manufactured_direct_bank()
    # A small explicit manufactured Chebyshev bank (private TensorBank, canonical Chebyshev sites), not the legacy artifact.
    cheb_sites = Blast.chebpoints(length(dm_bank.grid.logsites) - 1, minimum(dm_bank.grid.logsites), maximum(dm_bank.grid.logsites))
    cheb_tables = NamedTuple{Blast._BANK_GEOMETRIES}(Tuple(getfield(dm_bank.tensors, g) for g in Blast._BANK_GEOMETRIES))
    cheb_bank = Blast.TensorBank(cheb_tables, cheb_sites; basis = :chebyshev, dataset_id = "private:ratio-public-chebyshev9")
    rb = Blast.CubicBSplineBackend(rm_bank)
    A1, A2 = make_chain(rb), make_chain(rb)
    R2, D, C = make_chain(Blast.CubicBSplineBackend(rm2_bank)), make_chain(Blast.CubicBSplineBackend(dm_bank)), make_chain(Blast.ChebyshevBackend(cheb_bank))
    @test A1.plans !== A2.plans && A1.plans.plan_ϕTT !== A2.plans.plan_ϕTT && A1.ws.w_2_00_ϕTT.w !== A2.ws.w_2_00_ϕTT.w
    @test A1.W.w_2_00_ϕTT.T̃ === A2.W.w_2_00_ϕTT.T̃ === rm_bank.TT.T_2_00       # one read-only bank
    amps = (1.0, 1.3)
    chains = Dict("A1" => A1, "A2" => A2, "R2" => R2, "D" => D, "C" => C)
    serial = Dict((n, a) => run_chain(c, a) for (n, c) in chains for a in amps)
    @test isapprox(serial[("A1", 1.0)], serial[("A2", 1.0)]; rtol = 1e-12)
    @test rpc_relnorm(serial[("A1", 1.3)], serial[("A1", 1.0)]) > 1e-3
    order = [("A1", 1.0), ("D", 1.0), ("C", 1.0), ("A2", 1.3), ("R2", 1.0), ("A1", 1.3), ("C", 1.3), ("D", 1.3), ("R2", 1.3), ("A2", 1.0), ("C", 1.0), ("A1", 1.0)]
    for (n, a) in order
        @test isapprox(run_chain(chains[n], a), serial[(n, a)]; rtol = 1e-12)       # each equals its isolated reference; returns to x0
    end
    Pk_d = Blast.prepare_pk_workspace(D.plans, D.base, pkl, pkn, bg)
    @test_throws ArgumentError Blast.compute_w(A1.W, Pk_d)                      # cross-mix rejected; no global backend state
    @test_throws ArgumentError Blast.compute_w(C.W, Blast.prepare_pk_workspace(A1.plans, A1.base, pkl, pkn, bg))
    # Four spawned tasks with their own setups/plans/workspaces; ownership is checked at any thread count
    # (real multithreaded execution is only exercised when Threads.nthreads() > 1; nothing is claimed otherwise).
    task_amps = [(1.0, 1.3, 1.0), (1.3, 1.0, 1.3), (1.0, 1.0, 1.3), (1.3, 1.3, 1.0)]
    tasks = map(task_amps) do amps_i
        Threads.@spawn begin
            c = make_chain(rb)
            [run_chain(c, a) for a in amps_i]
        end
    end
    results = fetch.(tasks)
    for t in eachindex(task_amps), i in 1:3
        @test isapprox(results[t][i], serial[("A1", task_amps[t][i])]; rtol = 1e-12)
    end
end

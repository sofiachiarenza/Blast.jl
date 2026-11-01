# Phase 5 AD through the PUBLIC direct spline path (prepared Mooncake reverse mode vs ForwardDiff).
# Selector: `bspline-backend-ad`. Heavy: the parent launches it outside any session timeout.
using Test, Blast, ForwardDiff, Mooncake, LinearAlgebra, Random
using AbstractCosmologicalEmulators, FastChebInterp
using DifferentiationInterface: prepare_gradient, gradient
using ADTypes: AutoMooncake, AutoForwardDiff
isdefined(@__MODULE__, :get_test_cosmo) || include("test_fixtures.jl")
isdefined(@__MODULE__, :build_characterization_state) || include("characterization_workload.jl")
isdefined(@__MODULE__, :manufactured_direct_bank) || include("bspline_backend_fixture.jl")

ad_progress(msg) = (println(msg); flush(stdout))

@testset "B-spline backend AD: fixed background, Pmm amplitudes and Limber-only input (all 17 blocks)" begin
    state = build_characterization_state()
    (; bg, probes, pmm, ell) = state
    gc, wl, cmb = probes.galaxy, probes.lensing, probes.cmb
    backend = Blast.CubicBSplineBackend(manufactured_direct_bank())
    W, plans = Blast.SetUp(gc, wl, cmb; backend)
    @test count(!isnothing, getfield.(Ref(W), fieldnames(typeof(W)))) == 17
    pk = fixture_pk(Blast.nonlimber_grid(plans)); pkl, pkn = pmm.pk_limber_linear, pmm.pk_limber_nonlinear
    spectrum(a, b, c) = Blast.prepare_pk_workspace(plans, a .* pk, b .* pkl, c .* pkn, bg)
    function scalar(cls, w)           # several projections with distinct weights, normalized
        cl_gg, cl_gl, cl_ll, cl_gc, cl_lc = cls
        return w[1] * sum(abs2, cl_gg) + w[2] * sum(abs2, cl_gl) + w[3] * sum(abs2, cl_ll) + w[4] * sum(abs2, cl_gc) + w[5] * sum(abs2, cl_lc)
    end
    project(Pk, Wc) = (Blast.get_Cℓ(ell, gc, Pk, Wc, bg, plans), Blast.get_Cℓ(ell, gc, wl, Pk, Wc, bg, plans), Blast.get_Cℓ(ell, wl, Pk, Wc, bg, plans),
                       Blast.get_Cℓ(ell, gc, cmb, Pk, Wc, bg, plans), Blast.get_Cℓ(ell, wl, cmb, Pk, Wc, bg, plans))
    base = project(spectrum(1.0, 1.0, 1.0), Blast.compute_w(W, spectrum(1.0, 1.0, 1.0)))
    weights = Tuple(inv(sum(abs2, c)) for c in base)
    ws = Blast.allocate_compute_w(W, spectrum(1.0, 1.0, 1.0))

    # (a) three Pmm amplitudes, fixed Float64 background: previously a MethodError in the public path.
    functional_a = x -> (P = spectrum(x[1], x[2], x[3]); scalar(project(P, Blast.compute_w(W, P)), weights))
    reusable_a = x -> (P = spectrum(x[1], x[2], x[3]); Blast.compute_w!(ws, P); scalar(project(P, ws), weights))
    # (b) Limber-only differentiated input with constant non-Limber Pmm.
    functional_b = x -> (P = spectrum(1.0, x[1], x[2]); scalar(project(P, Blast.compute_w(W, P)), weights))
    reusable_b = x -> (P = spectrum(1.0, x[1], x[2]); Blast.compute_w!(ws, P); scalar(project(P, ws), weights))
    ad_progress("AD (a,b): setup and primal done; starting prepared Mooncake/ForwardDiff comparisons")
    for (functional, reusable, x0, x1) in ((functional_a, reusable_a, [1.0, 1.0, 1.0], [1.07, 0.93, 1.11]),
                                           (functional_b, reusable_b, [1.0, 1.0], [0.9, 1.08]))
        @test isapprox(reusable(x0), functional(x0); rtol=1e-10)
        ad_progress("AD: prepare_gradient (cold) for objective with $(length(x0)) inputs")
        mc = AutoMooncake(; config=nothing); prep = prepare_gradient(reusable, mc, x0)       # prepared once, reused
        rng = MersenneTwister(20261020)
        x_random = 0.9 .+ 0.2 .* rand(rng, length(x0))
        for x in (x0, x1, x_random, x0)
            ad_progress("AD: point $x")
            g_mc = gradient(reusable, prep, mc, x); g_fd = ForwardDiff.gradient(functional, x)
            @test all(isfinite, g_fd) && any(!iszero, g_fd)
            @test isapprox(g_mc, g_fd; rtol=1e-10, atol=1e-12)         # the unchanged strict baseline criterion
        end
    end
end

@testset "B-spline backend AD: differentiated background / cosmology through the direct source preparation (PNG active)" begin
    backend = Blast.CubicBSplineBackend(manufactured_direct_bank())
    nz_grid = ones(1, 50); z_grid = collect(LinRange(0.0, 3.6, 50))
    # b_phi ∝ bias - p: unit bias with p=1 would silently null PNG.
    png_template = Blast.PrimordialNonGaussianity(nz=nz_grid, z=z_grid, bias=fill(1.5, 1, 200), f_NL=7.0, p=1.0)
    gc_template = Blast.GalaxyClustering(δ=Blast.NumberCounts(nz=nz_grid, z=z_grid, bias=fill(1.5, 1, 200)), PNG=png_template)
    W, plans = Blast.SetUp(gc_template; backend)
    # With PNG active the T and Φ source families genuinely contribute (n_s cancels in the TT family alone).
    @test !isnothing(plans.plan_ϕT) && !isnothing(plans.plan_ϕ) && !isnothing(W.w_2_00_ϕT) && !isnothing(W.w_2_00_ϕ)
    pk = fixture_pk(Blast.nonlimber_grid(plans))
    pkl = ones(length(Blast.z_cheb), length(Blast.k_limber)); pkn = copy(pkl)
    bias_fixed = fill(1.5, 1, length(get_test_bg(get_test_cosmo()).z))
    function f_cosmo(p)               # p = (h, n_s): distances, kernels AND primordial coefficients (bg-derived)
        bg = Blast.Background(get_test_cosmo(h=p[1], nₛ=p[2]))
        raw = Blast.GalaxyClustering(δ=Blast.NumberCounts(nz=nz_grid, z=z_grid, bias=bias_fixed),
                                     PNG=Blast.PrimordialNonGaussianity(nz=nz_grid, z=z_grid, bias=bias_fixed, f_NL=7.0, p=1.0))
        gc = Blast.prepare_probe(raw, bg)
        Pk = Blast.prepare_pk_workspace(plans, pk, pkl, pkn, bg)
        Wc = Blast.compute_w(W, Pk)
        return sum(Blast.get_Cℓ(Blast.full_ℓ_range, gc, Pk, Wc, bg, plans))
    end
    p0 = [0.6727, 0.9645]
    @test isfinite(f_cosmo(p0))
    # Physical n_s dependence by finite variation of the Float64 objective (not a derivative artifact).
    @test abs(f_cosmo([p0[1], p0[2] + 1e-2]) - f_cosmo(p0)) > 1e-8 * abs(f_cosmo(p0))
    mc, fd = AutoMooncake(; config=nothing), AutoForwardDiff()
    ad_progress("AD (cosmology): prepare_gradient Mooncake and ForwardDiff at p0 (cold)")
    prep_mc = prepare_gradient(f_cosmo, mc, p0); prep_fd = prepare_gradient(f_cosmo, fd, p0)
    rng = MersenneTwister(20261021)
    p_random = [0.64 + 0.06rand(rng), 0.96 + 0.02rand(rng)]
    for p in (p0, [0.70, 0.96], [0.64, 0.98], p_random, p0)     # seeded fixture-domain draw, not a Turing prior claim
        ad_progress("AD (cosmology): point $p")
        g_mc = gradient(f_cosmo, prep_mc, mc, p); g_fd = gradient(f_cosmo, prep_fd, fd, p)
        @test all(isfinite, g_mc) && all(isfinite, g_fd) && all(!iszero, g_fd)
        @test isapprox(g_mc, g_fd; rtol=1e-9, atol=1e-12)
    end
end

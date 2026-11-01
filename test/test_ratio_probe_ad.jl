# Prepared Mooncake vs independent ForwardDiff for probe nuisance parameters through the ratio path.
# Selector: `ratio-probe-ad` (self-contained imports). Five-vector, but only FOUR are active: the stored
# RSD growth_rate constructor argument is IGNORED by the existing preparation (src/probes.jl uses bg.f), so
# its derivative is exactly zero and it serves as a negative control, not as a new nuisance parameter.
using Test, Blast, ForwardDiff, Mooncake, LinearAlgebra
using AbstractCosmologicalEmulators, FastChebInterp
using DifferentiationInterface: prepare_gradient, gradient
using ADTypes: AutoMooncake
isdefined(@__MODULE__, :build_characterization_state) || include("characterization_workload.jl")
isdefined(@__MODULE__, :manufactured_ratio_bank) || include("bspline_backend_fixture.jl")

rpa_stage(msg) = (println(msg); flush(stdout))
rpa_nonlinear_boost(k, z) = 1 + 0.45 * (k / (0.35 + k))^2 / (1 + z)

function rpa_make_probes(bg, sb, sf, sr, sm, sia)
    s = characterization_source_distributions(); n_bg = length(bg.z)
    T = promote_type(typeof(sb), typeof(sf), typeof(sr), typeof(sm), typeof(sia), eltype(bg.z))
    bias = Matrix{T}(undef, 2, n_bg); slope = similar(bias)
    @views begin
        bias[1, :] .= sb .* (1.20 .+ 0.10 .* bg.z); bias[2, :] .= sb .* (1.55 .+ 0.08 .* bg.z)
        slope[1, :] .= sm .* (0.24 .+ 0.015 .* bg.z); slope[2, :] .= sm .* (0.34 .+ 0.010 .* bg.z)
    end
    galaxy = Blast.GalaxyClustering(δ = Blast.NumberCounts(nz = s.nz_g, z = s.z, bias = bias),
        RSD = Blast.RedshiftSpaceDistortions(nz = s.nz_g, z = s.z, growth_rate = sr .* bg.f),
        μ = Blast.MagnificationBias(nz = s.nz_g, z = s.z, s = slope),
        PNG = Blast.PrimordialNonGaussianity(nz = s.nz_g, z = s.z, bias = bias, f_NL = sf, p = 1.0))
    lensing = Blast.WeakLensing(γ = Blast.CosmicShear(nz = s.nz_s, z = s.z), IA = Blast.IntrinsicAlignment(nz = s.nz_s, z = s.z, A = sia * 1.72))
    cmb = Blast.CMB(κ = Blast.CMBLensing(), ISW = Blast.IntegratedSachsWolfe(growth_rate = bg.f))
    return Blast.prepare_probe(galaxy, bg), Blast.prepare_probe(lensing, bg), Blast.prepare_probe(cmb, bg)
end

@testset "Ratio probe AD: prepared Mooncake vs ForwardDiff, four active nuisances + RSD negative control" begin
    backend = Blast.CubicBSplineBackend(manufactured_ratio_bank())
    bg = Blast.Background(characterization_cosmology())
    gc0, wl0, cmb0 = rpa_make_probes(bg, 1.0, 7.0, 1.0, 1.0, 1.0)
    W, plans = Blast.SetUp(gc0, wl0, cmb0; backend)
    grid = Blast.nonlimber_grid(plans); pk = fixture_pk(grid)
    kl = 10 .^ Blast.k_limber
    pkl = [fixture_pmm(k, z) for z in Blast.z_cheb, k in kl]; pkn = [fixture_pmm(k, z) * rpa_nonlinear_boost(k, z) for z in Blast.z_cheb, k in kl]
    Pk = Blast.prepare_pk_workspace(plans, pk, pkl, pkn, bg; pk_grid = grid)       # Pmm fixed, coefficients Float64
    Wf = Blast.compute_w(W, Pk); ell = Blast.full_ℓ_range
    spectra(x) = ((gc, wl, cmb) = rpa_make_probes(bg, x[1], x[2], x[3], x[4], x[5]);
        (Blast.get_Cℓ(ell, gc, Pk, Wf, bg, plans), Blast.get_Cℓ(ell, gc, wl, Pk, Wf, bg, plans), Blast.get_Cℓ(ell, wl, Pk, Wf, bg, plans),
         Blast.get_Cℓ(ell, gc, cmb, Pk, Wf, bg, plans), Blast.get_Cℓ(ell, wl, cmb, Pk, Wf, bg, plans)))
    x0 = [1.0, 7.0, 1.0, 1.0, 1.0]; x1 = [1.05, 5.0, 0.95, 1.1, 0.9]
    base = spectra(x0)
    @test all(c -> all(isfinite, c), base)
    w = Tuple(inv(sum(abs2, c)) for c in base)
    f = x -> (cls = spectra(x); sum(w[i] * sum(abs2, cls[i]) for i in 1:5))
    @test isfinite(f(x0)) && isfinite(f(x1))
    # Negative control: the stored RSD growth argument changes, the kernel (built from bg.f) and the primal do not.
    x_rsd = copy(x0); x_rsd[3] = 1.37
    gc_rsd, _, _ = rpa_make_probes(bg, x_rsd...)
    @test gc_rsd.RSD.growth_rate != gc0.RSD.growth_rate
    @test gc_rsd.RSD.Kernel == gc0.RSD.Kernel
    @test f(x_rsd) == f(x0)
    mc = AutoMooncake(; config = nothing)
    rpa_stage("ratio probe AD: prepare_gradient (cold)")
    prep = prepare_gradient(f, mc, x0)
    for (label, x) in (("x0", x0), ("x1", x1), ("x0_again", x0))
        rpa_stage("ratio probe AD: point $label")
        g_mc = copy(gradient(f, prep, mc, x)); g_fd = ForwardDiff.gradient(f, x)
        rpa_stage("  mooncake    = $g_mc\n  forwarddiff = $g_fd")
        @test all(isfinite, g_mc) && all(isfinite, g_fd)
        @test all(i -> !iszero(g_fd[i]), (1, 2, 4, 5))                       # the four active parameters
        @test iszero(g_mc[3]) && iszero(g_fd[3])                              # RSD argument: exactly zero (not an amplitude)
        @test isapprox(g_mc, g_fd; rtol = 1e-10, atol = 1e-12)                # strict standard
    end
end

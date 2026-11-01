# Phase 3: explicit Chebyshev bank/backend/plan ownership.
using Test, Blast, LinearAlgebra
# Load both packages that export `chebpoints` (as the full suite does) so the
# helpers below must qualify it.
using AbstractCosmologicalEmulators, FastChebInterp
isdefined(@__MODULE__, :get_test_cosmo) || include("test_fixtures.jl")
isdefined(@__MODULE__, :check_gradient) || include("test_diff_helpers.jl")
isdefined(@__MODULE__, :build_characterization_state) || include("characterization_workload.jl")

# Frozen value of the versioned digest; a change means the identity contract changed.
const PINNED_LEGACY_DIGEST = (0xe646e33677f02588, 0x3bf9ff7fa419f634, 0x9b16abac438209a7, 0x0aff2843f75bbf72)

_active_fields(W) = Tuple(name for name in fieldnames(typeof(W)) if !isnothing(getfield(W, name)))

function _spectral_probes()
    gc = Blast.GalaxyClustering(δ=Blast.NumberCounts(), RSD=Blast.RedshiftSpaceDistortions(),
                                μ=Blast.MagnificationBias(), PNG=Blast.PrimordialNonGaussianity())
    gc_arf = Blast.GalaxyClustering(δ=Blast.NumberCounts(), RSD=Blast.RedshiftSpaceDistortions(),
                                    μ=Blast.MagnificationBias(), PNG=Blast.PrimordialNonGaussianity(),
                                    ARF=Blast.AngularRedshiftFluctuations())
    gc_min = Blast.GalaxyClustering(δ=Blast.NumberCounts())
    wl = Blast.WeakLensing(γ=Blast.CosmicShear(), IA=Blast.IntrinsicAlignment())
    cmb = Blast.CMB(κ=Blast.CMBLensing(), ISW=Blast.IntegratedSachsWolfe())
    return (; gc, gc_arf, gc_min, wl, cmb)
end

# A tiny foreign Chebyshev bank (5 sites, zero tensors) with a different signature.
function _foreign_backend()
    s = Blast.chebpoints(4, log10(1e-3), log10(1.0))
    Z() = zeros(length(Blast.ℓ_nonlimber), length(Blast.χ), length(Blast.R), length(s))
    tensors = (T_2_00=Z(), T_0_00=Z(), T_minus2_00=Z(), T_0_02=Z(), T_0_20=Z(),
               T_2_02=Z(), T_2_20=Z(), T_2_22=Z())
    return Blast.ChebyshevBackend(Blast.TensorBank(tensors, s; basis=:chebyshev,
                                                   dataset_id="private:test-foreign-bank"))
end

@testset "Spectral setup: lazy fixed legacy bank" begin
    code = "using Blast; exit(isassigned(Blast._LEGACY_BANK) ? 1 : 0)"
    @test success(`$(Base.julia_cmd()) --project=$(Base.active_project()) -e $code`)

    bank = Blast.legacy_tensor_bank()
    @test bank === Blast.legacy_tensor_bank()
    @test bank.sites == Blast.k_cheb
    @test bank.signature == Blast.LEGACY_BANK_SIGNATURE
    backend = Blast.ChebyshevBackend()
    @test Blast.tensor_bank(backend) === bank
    grid = Blast.nonlimber_grid(backend)
    @test grid.k == 10 .^ Blast.k_cheb
    @test grid.z == Blast.z_lin
    @test @inferred(Blast.nonlimber_grid(backend)) isa NamedTuple
    @test Blast.FFTPlans === Blast.SpectralPlans

    p = _spectral_probes()
    W1, _ = Blast.SetUp(p.gc)
    W2, _ = Blast.SetUp(p.gc_min)
    @test W1.w_2_00_ϕTT.T̃ === W2.w_2_00_ϕTT.T̃ === bank.T_2_00
    # Source families sharing a geometry alias one bank tensor (fusion relies on it).
    @test W1.w_2_00_ϕTT.T̃ === W1.w_2_00_ϕT.T̃ === W1.w_2_00_ϕT_R1.T̃ === W1.w_2_00_ϕ.T̃
    @test W1.w_2_02_ϕTT.T̃ === W1.w_2_02_ϕT.T̃ === W1.w_2_02_ϕT_R1.T̃
    @test @inferred(Blast.bank_signature(W1)) == Blast.LEGACY_BANK_SIGNATURE
end

@testset "Spectral setup: canonical sites and deterministic signatures" begin
    bank = Blast.legacy_tensor_bank()
    tensors = NamedTuple{Blast._BANK_GEOMETRIES}(Tuple(getfield(bank, g) for g in Blast._BANK_GEOMETRIES))
    private(sites) = Blast.TensorBank(tensors, sites; basis=:chebyshev, dataset_id="private:perturbation-test")
    # An interior site moved by 1e-8 is a different grid from the canonical plan nodes: reject it,
    # even though the bank itself is internally consistent (its signature matches its own sites).
    moved = copy(Blast.k_cheb); moved[80] += 1e-8
    @test_throws ArgumentError Blast.ChebyshevBackend(private(moved))
    @test_throws ArgumentError Blast.ChebyshevBackend(private(fill(-2.0, length(Blast.k_cheb))))
    nan_sites = copy(Blast.k_cheb); nan_sites[10] = NaN
    @test_throws ArgumentError Blast.ChebyshevBackend(private(nan_sites))
    # The artifact namespace is reserved for the trusted loader.
    @test_throws ArgumentError Blast.TensorBank(tensors, Blast.k_cheb; basis=:chebyshev,
                                                dataset_id=Blast._LEGACY_DATASET_ID)
    # Sorting or reversing the declared order is not silently accepted either.
    @test_throws ArgumentError Blast.ChebyshevBackend(private(reverse(Blast.k_cheb)))
    # A few ulp of roundoff between native constructions is tolerated.
    jitter = copy(Blast.k_cheb); jitter[80] = nextfloat(jitter[80], 2)
    @test Blast.ChebyshevBackend(private(jitter)) isa Blast.ChebyshevBackend
    # The same bank reproduces the legacy-tagged contraction data only under the legacy signature.
    @test private(Blast.k_cheb).signature != Blast.LEGACY_BANK_SIGNATURE

    # Digest is deterministic, versioned, and sensitive to every identity ingredient.
    sig = Blast.LEGACY_BANK_SIGNATURE
    @test sig == Blast.bank_signature_digest(:chebyshev, Blast._LEGACY_DATASET_ID, Blast.k_cheb)
    @test startswith(Blast._LEGACY_DATASET_ID, "artifact:T_tildes:574a5b109567d1977df2376ed80379a0a2b2bd6f")
    @test sig != Blast.bank_signature_digest(:chebyshev, Blast._LEGACY_DATASET_ID, moved)
    @test sig != Blast.bank_signature_digest(:chebyshev, "artifact:T_tildes:other", Blast.k_cheb)
    @test sig != Blast.bank_signature_digest(:bspline, Blast._LEGACY_DATASET_ID, Blast.k_cheb)
    @test sig.digest == PINNED_LEGACY_DIGEST
end

@testset "Spectral setup: default equals explicit Chebyshev backend, all probe combinations" begin
    p = _spectral_probes()
    backend = Blast.ChebyshevBackend()
    combos = [(p.gc,), (p.gc_arf,), (p.gc_min,), (p.wl,), (p.gc, p.wl), (p.wl, p.gc),
              (p.gc, p.cmb), (p.cmb, p.gc), (p.gc_arf, p.cmb), (p.wl, p.cmb), (p.cmb, p.wl),
              (p.gc, p.wl, p.cmb), (p.gc, p.cmb, p.wl), (p.wl, p.gc, p.cmb),
              (p.wl, p.cmb, p.gc), (p.cmb, p.gc, p.wl), (p.cmb, p.wl, p.gc),
              (p.gc_arf, p.wl, p.cmb)]
    for probes in combos
        Wd, Pd = Blast.SetUp(probes...)
        We, Pe = Blast.SetUp(probes...; backend=backend)
        @test typeof(Wd) == typeof(We)
        @test typeof(Pd) == typeof(Pe)
        @test _active_fields(Wd) == _active_fields(We)
        for f in _active_fields(Wd)
            @test getfield(Wd, f).T̃ === getfield(We, f).T̃
            @test getfield(Wd, f).signature == getfield(We, f).signature
        end
        @test Pe.backend === backend || Pe.backend.bank === backend.bank
        @test Blast.bank_signature(Wd) == Blast.bank_signature(Pd) == Blast.LEGACY_BANK_SIGNATURE
        for plan in (:plan_ϕTT, :plan_ϕT, :plan_ϕ)
            a, b = getfield(Pd, plan), getfield(Pe, plan)
            @test isnothing(a) == isnothing(b)
            isnothing(a) || (@test a.K == b.K && a.nodes == b.nodes)
        end
        @test Pd.plan_ϕTT.nodes[1] ≈ Blast.k_cheb || Pd.plan_ϕTT.nodes[1] ≈ reverse(Blast.k_cheb)
    end
end

@testset "Spectral setup: signatures, functional vs reusable, rejection of mismatches" begin
    state = build_characterization_state()
    (; bg, probes, plans, ell) = state
    gc, wl = probes.galaxy, probes.lensing
    W, Plans = state.weights_template, plans
    pkl, pkn = state.pmm.pk_limber_linear, state.pmm.pk_limber_nonlinear
    PS = @inferred Blast.prepare_pk_workspace(Plans, state.pmm.pk_nonlimber, pkl, pkn, bg)
    @test PS.signature == Blast.bank_signature(Plans) == Blast.bank_signature(W)

    Wf = Blast.compute_w(W, PS)
    Wa = @inferred Blast.allocate_compute_w(W, PS)
    Wb = Blast.allocate_compute_w(W, PS)
    @test Blast.bank_signature(Wf) == Blast.bank_signature(Wa) == PS.signature
    @test @inferred(Blast.compute_w!(Wa, PS)) === Wa
    for f in _active_fields(W)
        @test getfield(Wa, f).w ≈ getfield(Wf, f).w rtol=1e-10 atol=1e-9
        @test getfield(Wa, f).w !== getfield(Wb, f).w          # independent workspaces
        @test getfield(Wa, f).T̃ === getfield(Wb, f).T̃          # shared read-only bank
    end
    # Independent workspaces stay independent when only one is overwritten.
    Blast.compute_w!(Wb, PS)
    fill!(Wb.w_2_00_ϕTT.w, 0.0)
    @test any(!iszero, Wa.w_2_00_ϕTT.w)

    ref = Blast.get_Cℓ(ell, gc, wl, PS, Wf, bg, Plans)
    @test Blast.get_Cℓ(ell, gc, wl, PS, Wf, bg) ≈ ref rtol=1e-8   # P = nothing is checked too
    @test Blast.get_Cℓ(ell, gc, wl, PS, Wa, bg, Plans) ≈ ref rtol=1e-10

    # Foreign configuration: different sites/bank, same everything else.
    foreign = _foreign_backend()
    Wx, Px = Blast.SetUp(gc, wl, probes.cmb; backend=foreign)
    @test Blast.bank_signature(Px) != Blast.bank_signature(Plans)
    @test size(Blast.nonlimber_grid(Px).k) == (5,)
    @test_throws DimensionMismatch Blast.prepare_pk_workspace(Px, state.pmm.pk_nonlimber, pkl, pkn, bg)
    PSx = Blast.prepare_pk_workspace(Px, rand(length(Blast.z_lin), 5) .+ 1.0, pkl, pkn, bg)
    @test PSx.signature == Blast.bank_signature(Px)
    @test_throws ArgumentError Blast.compute_w(W, PSx)
    @test_throws ArgumentError Blast.allocate_compute_w(W, PSx)
    @test_throws ArgumentError Blast.compute_w!(Wa, PSx)
    @test_throws ArgumentError Blast.get_Cℓ(ell, gc, wl, PSx, Wf, bg)
    @test_throws ArgumentError Blast.get_Cℓ(ell, gc, wl, PSx, Wf, bg, Plans)
    @test_throws ArgumentError Blast.get_Cℓ(ell, gc, wl, PS, Wf, bg, Px)   # plans from another bank
    # Interleaved setups remain valid.
    @test Blast.get_Cℓ(ell, gc, wl, PS, Wf, bg, Plans) ≈ ref rtol=1e-12

    # Legacy six-argument plan constructor: accepted only for the legacy sites.
    legacy = Blast.FFTPlans(Plans.plan_ϕTT, Plans.plan_ϕT, Plans.plan_ϕ, Plans.plan_limber,
                            Plans.T_k_limber, Plans.plan_ℓ)
    @test Blast.bank_signature(legacy) == Blast.LEGACY_BANK_SIGNATURE
    @test_throws ArgumentError Blast.FFTPlans(Px.plan_ϕTT, Px.plan_ϕT, Px.plan_ϕ, Px.plan_limber,
                                              Px.T_k_limber, Px.plan_ℓ)
    # Original @kwdef-style keyword construction, with and without the optional plans.
    kw = Blast.FFTPlans(plan_ϕTT=Plans.plan_ϕTT, plan_ϕT=Plans.plan_ϕT, plan_ϕ=Plans.plan_ϕ,
                        plan_limber=Plans.plan_limber, T_k_limber=Plans.T_k_limber, plan_ℓ=Plans.plan_ℓ)
    @test typeof(kw) == typeof(Plans)
    nophi = Blast.SpectralPlans(plan_ϕTT=Plans.plan_ϕTT, plan_limber=Plans.plan_limber,
                                T_k_limber=Plans.T_k_limber, plan_ℓ=Plans.plan_ℓ)
    @test isnothing(nophi.plan_ϕT) && isnothing(nophi.plan_ϕ)
    @test Blast.bank_signature(nophi) == Blast.LEGACY_BANK_SIGNATURE
    # A mixed-bank workspace is rejected.
    m = Wf.w_minus2_00_ϕTT
    alien = typeof(m)(m.T̃, m.w, Blast.BankSignature(:chebyshev, 5, (UInt64(0), UInt64(0), UInt64(0), UInt64(0))))
    Wmixed = Blast.ProjectedMatterDensity(w_2_00_ϕTT=Wf.w_2_00_ϕTT, w_minus2_00_ϕTT=alien)
    @test_throws ArgumentError Blast.bank_signature(Wmixed)
    # The bank constructor refuses the legacy tag for other site counts.
    c5 = Blast.cϕTT(zeros(5, 2, 2))
    @test_throws ArgumentError Blast.PowerSpectrum(c5, nothing, nothing, ones(2, 2), ones(2, 2))
end

@testset "Spectral setup: raw legacy constructors inside an AD objective" begin
    cosmo = get_test_cosmo()
    bg_ref = get_test_bg(cosmo)
    nz_grid = ones(1, 50)
    z_grid = collect(LinRange(0.0, 3.6, 50))
    pk_grid = ones(length(Blast.z_lin), length(Blast.k_cheb))
    pkl = ones(length(Blast.z_cheb), length(Blast.k_limber))
    _, Plans = Blast.SetUp(Blast.GalaxyClustering(δ=Blast.NumberCounts(nz=nz_grid, z=z_grid,
                                                                       bias=ones(1, 200))))
    bias_fixed = ones(1, length(bg_ref.z))
    function f_cosmo(q)
        bg = Background(get_test_cosmo(ωc=q[1], h=q[2], w0=q[3]))
        raw = Blast.GalaxyClustering(δ=Blast.NumberCounts(nz=nz_grid, z=z_grid, bias=bias_fixed))
        gc = Blast.prepare_probe(raw, bg)
        W_local = Blast.ProjectedMatterDensity(w_2_00_ϕTT=Blast.w_2_00_ϕTT())   # raw legacy constructor
        PS = Blast.prepare_pk_workspace(Plans, pk_grid, pkl, pkl, bg)
        W_local = Blast.compute_w(W_local, PS)
        return sum(Blast.get_Cℓ(Blast.full_ℓ_range, gc, PS, W_local, bg, Plans))
    end
    p0 = [0.12055273725599998, 0.6727, -1.0]
    g = check_gradient(f_cosmo, p0; skip_zygote=true, skip_fd=true, rtol_ad=1e-7)
    @test all(isfinite, g.ForwardDiff) && all(!iszero, g.ForwardDiff)
    @test g.Mooncake ≈ g.ForwardDiff rtol=1e-7
end

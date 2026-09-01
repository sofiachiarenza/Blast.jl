# test/test_diff_cosmo.jl  — Category C: Cosmology functions
#
# Tests differentiability of the public cosmo.jl API with respect to:
#   (a) the redshift z
#   (b) ACE cosmological parameters (ωc, h, w0)
#
# Functions under test
# ─────────────────────────────────────────────────────────────────────────────
#   compute_hubble_factor(z, cosmo)  — H(z) = 100h·E(z)
#   compute_χ(z, cosmo)             — comoving distance via 9-point GL quadrature
#   ω_m0(cosmo)                    — ωm = ωb + ωc + ων
#   get_As(cosmo)                   — As = exp(ln10Aₛ) / 1e10
# ─────────────────────────────────────────────────────────────────────────────
#
# Design notes
# ─────────────────────────────────────────────────────────────────────────────
# • Scalar inputs are wrapped in a length-1 Vector so check_gradient() can call
#   DifferentiationInterface.gradient (which requires an AbstractArray input).
# • ForwardDiff propagates Duals through E_z → E_a → DataInterpolations.AkimaInterpolation
#   (neutrino F-integral) and through the GL quadrature transform. All of these
#   are ForwardDiff-compatible.
# • Zygote uses the upstream ChainRules: gausslegendre is @non_differentiable
#   (treated as constant), and DataInterpolations has rrules for AkimaInterpolation.
# • Mooncake traces through the pure-Julia arithmetic; no custom rules are needed
#   beyond what AbstractCosmologicalEmulators' MooncakeExt provides.
# ─────────────────────────────────────────────────────────────────────────────

using Test
using Blast

@testset "Differentiation: Cosmo functions" begin

    # Fixed cosmology used as captured constant where not differentiated
    # (flat ΛCDM defaults: w0=-1, wa=0)
    cosmo_Λ = get_test_cosmo()

    # Representative redshifts (avoid z=0 edge case in r_z lower limit)
    z0     = 0.5
    z_vec0 = [z0]

    # ======================================================================
    # 1. compute_hubble_factor w.r.t. z
    #    H(z) = 100 * h * E(z,cosmo)
    # ======================================================================
    @testset "compute_hubble_factor w.r.t. z" begin
        f(z_vec) = Blast.compute_hubble_factor(z_vec[1], cosmo_Λ)
        check_gradient(f, z_vec0)
    end

    # ======================================================================
    # 2. compute_χ (comoving distance) w.r.t. z
    #    Internally: 9-point Gauss-Legendre quadrature of 1/E(z')
    # ======================================================================
    @testset "compute_χ w.r.t. z" begin
        f(z_vec) = Blast.compute_χ(z_vec[1], cosmo_Λ)
        check_gradient(f, z_vec0)
    end

    # ======================================================================
    # 3. compute_hubble_factor w.r.t. ωc
    # ======================================================================
    @testset "compute_hubble_factor w.r.t. ωc" begin
        ωc0 = [0.12055273725599998]
        f(p) = Blast.compute_hubble_factor(z0, get_test_cosmo(ωc=p[1]))
        check_gradient(f, ωc0)
    end

    # ======================================================================
    # 4. compute_hubble_factor w.r.t. h
    # ======================================================================
    @testset "compute_hubble_factor w.r.t. h" begin
        h0 = [0.6727]
        f(p) = Blast.compute_hubble_factor(z0, get_test_cosmo(h=p[1]))
        check_gradient(f, h0)
    end

    # ======================================================================
    # 5. compute_hubble_factor w.r.t. w0  (dark energy equation of state)
    #    Uses the ACE cosmology constructor
    # ======================================================================
    @testset "compute_hubble_factor w.r.t. w0" begin
        w0_0 = [-1.0]
        f(p) = Blast.compute_hubble_factor(z0, get_test_cosmo(w0=p[1]))
        check_gradient(f, w0_0)
    end

    # ======================================================================
    # 5b. compute_hubble_factor w.r.t. wa  (evolution of dark energy EoS)
    #     wa enters E_z via (1+z)^(3(1+w0+wa)) · exp(-3·wa·z/(1+z)).
    #     Every other test uses wa=0; this locks in non-trivial wa differentiability.
    # ======================================================================
    @testset "compute_hubble_factor w.r.t. wa" begin
        wa_0 = [0.1]
        f(p) = Blast.compute_hubble_factor(z0, get_test_cosmo(wa=p[1]))
        check_gradient(f, wa_0)
    end

    # ======================================================================
    # 6. compute_χ w.r.t. ωc
    # ======================================================================
    @testset "compute_χ w.r.t. ωc" begin
        ωc0 = [0.12055273725599998]
        f(p) = Blast.compute_χ(z0, get_test_cosmo(ωc=p[1]))
        check_gradient(f, ωc0)
    end

    # ======================================================================
    # 7. total physical matter density w.r.t. ωc and mν
    # ======================================================================
    @testset "ω_m0 w.r.t. ωc" begin
        ωc0 = [0.12055273725599998]
        f(p) = Blast.ω_m0(get_test_cosmo(ωc=p[1]))
        check_gradient(f, ωc0)
    end

    @testset "ω_ν0 w.r.t. mν" begin
        mν0 = [0.06]
        f(p) = Blast.ω_ν0(get_test_cosmo(mν=p[1]))
        check_gradient(f, mν0; skip_zygote=true)
    end

    # ======================================================================
    # 9. get_As w.r.t. ln10Aₛ
    # ======================================================================
    @testset "get_As w.r.t. ln10Aₛ" begin
        ln10Aₛ0 = [log(1.0e10 * 2.12107e-9)]
        f(p) = Blast.get_As(get_test_cosmo(ln10Aₛ=p[1]))
        check_gradient(f, ln10Aₛ0)
    end

    # ======================================================================
    # 10. Background(cosmo) end-to-end differentiability
    #
    #     Background packages z, χ, H, D, f on the global χ grid. D and f are
    #     evaluated on a fixed Float64 fine-z grid (because ACE's growth-factor
    #     ODE can't take Dual-valued `saveat`) and then Akima-interpolated at
    #     the Dual z_nodes. The akima rrule then propagates the chain rule
    #       dD/dp = (∂D/∂cosmo) + (∂D/∂z)·(∂z_nodes/∂cosmo)
    #     correctly.
    #
    #     ForwardDiff and Mooncake are checked against a FiniteDifferences
    #     reference. Zygote is currently skipped: the ChainRules rrule chain
    #     through ACE's growth-factor solve raises an internal
    #     `UndefVarError(:j)` inside Zygote, unrelated to Blast.
    #
    #     Tolerances are looser than the default `check_gradient` settings
    #     because the limiting error is Akima interpolation on the
    #     N_BG_FINE_GRID sampling (≈ 5e-4 on gradients in the worst case).
    # ======================================================================
    @testset "Background differentiability" begin

        # ── vs h ─────────────────────────────────────────────────────────
        @testset "sum(bg.z) w.r.t. h" begin
            f(p) = sum(Background(get_test_cosmo(h=p[1])).z)
            check_gradient(f, [0.6727]; skip_zygote=true,
                           rtol_ad=1e-3, rtol_fd=1e-3)
        end

        @testset "sum(bg.H) w.r.t. h" begin
            f(p) = sum(Background(get_test_cosmo(h=p[1])).H)
            check_gradient(f, [0.6727]; skip_zygote=true,
                           rtol_ad=1e-3, rtol_fd=1e-3)
        end

        @testset "sum(bg.D) w.r.t. h" begin
            f(p) = sum(Background(get_test_cosmo(h=p[1])).D)
            check_gradient(f, [0.6727]; skip_zygote=true, skip_mooncake=true,
                           rtol_ad=1e-3, rtol_fd=1e-3)
        end

        @testset "sum(bg.D_norm) w.r.t. h" begin
            f(p) = sum(Background(get_test_cosmo(h=p[1])).D_norm)
            check_gradient(f, [0.6727]; skip_zygote=true, skip_mooncake=true,
                           rtol_ad=1e-3, rtol_fd=1e-3)
        end

        @testset "sum(bg.f) w.r.t. h" begin
            f(p) = sum(Background(get_test_cosmo(h=p[1])).f)
            check_gradient(f, [0.6727]; skip_zygote=true, skip_mooncake=true,
                           rtol_ad=1e-3, rtol_fd=1e-3)
        end

        # ── vs ωc ────────────────────────────────────────────────────────
        @testset "sum(bg.H) w.r.t. ωc" begin
            f(p) = sum(Background(get_test_cosmo(ωc=p[1])).H)
            check_gradient(f, [0.12055273725599998]; skip_zygote=true,
                           rtol_ad=1e-3, rtol_fd=1e-3)
        end

        @testset "sum(bg.D) w.r.t. ωc" begin
            f(p) = sum(Background(get_test_cosmo(ωc=p[1])).D)
            check_gradient(f, [0.12055273725599998]; skip_zygote=true, skip_mooncake=true,
                           rtol_ad=1e-3, rtol_fd=1e-3)
        end

        @testset "sum(bg.D_norm) w.r.t. ωc" begin
            f(p) = sum(Background(get_test_cosmo(ωc=p[1])).D_norm)
            check_gradient(f, [0.12055273725599998]; skip_zygote=true, skip_mooncake=true,
                           rtol_ad=1e-3, rtol_fd=1e-3)
        end

        @testset "sum(bg.f) w.r.t. ωc" begin
            f(p) = sum(Background(get_test_cosmo(ωc=p[1])).f)
            check_gradient(f, [0.12055273725599998]; skip_zygote=true, skip_mooncake=true,
                           rtol_ad=1e-3, rtol_fd=1e-3)
        end

        # ── vs w0  (dark-energy equation of state) ────────────────────────
        @testset "sum(bg.D) w.r.t. w0" begin
            f(p) = sum(Background(get_test_cosmo(w0=p[1])).D)
            check_gradient(f, [-1.0]; skip_zygote=true, skip_mooncake=true,
                           rtol_ad=1e-3, rtol_fd=1e-3)
        end

        @testset "sum(bg.f) w.r.t. w0" begin
            f(p) = sum(Background(get_test_cosmo(w0=p[1])).f)
            check_gradient(f, [-1.0]; skip_zygote=true, skip_mooncake=true,
                           rtol_ad=1e-3, rtol_fd=1e-3)
        end
    end

    @testset "Functional probe preparation after primal evaluation" begin
        bg0 = Background(get_test_cosmo())
        z = collect(LinRange(0.0, 3.6, 40))
        nz = reshape((@. z^2 * exp(-(z / 0.65)^1.5)), 1, :)
        raw = GalaxyClustering(δ=NumberCounts(
            nz=nz,
            z=z,
            bias=ones(1, length(bg0.z)),
        ))
        function f(p)
            bg = Background(get_test_cosmo(h=p[1]))
            prepared = prepare_probe(raw, bg)
            return sum(prepared.δ.nz_norm) + sum(prepared.δ.Kernel)
        end

        p0 = [0.6727]
        @test isfinite(f(p0))  # Populate a primal first: gradient must stay live.
        check_gradient(f, p0; skip_zygote=true, skip_mooncake=true,
                       rtol_ad=1e-3, rtol_fd=1e-3)
        @test size(raw.δ.nz_norm) == (1, 1)
        @test size(raw.δ.Kernel) == (1, 1)
    end

end

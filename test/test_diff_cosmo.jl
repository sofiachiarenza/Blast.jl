# test/test_diff_cosmo.jl  — Category C: Cosmology functions
#
# Tests differentiability of the public cosmo.jl API with respect to:
#   (a) the redshift z
#   (b) cosmological parameters (Ωm, H0, w0) passed through the w0waCDM
#       constructors
#
# Functions under test
# ─────────────────────────────────────────────────────────────────────────────
#   compute_hubble_factor(z, cosmo)  — H(z) = 100h·E(z)
#   compute_χ(z, cosmo)             — comoving distance via 9-point GL quadrature
#   get_Ωm(cosmo)                   — Ωm = (ωb + ωc) / h²
#   get_H0(cosmo)                   — H0 = 100h
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
    cosmo_Λ = w0waCDM()

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
    # 3. compute_hubble_factor w.r.t. Ωm
    #    Tests AD through w0waCDM constructor + E_z
    # ======================================================================
    @testset "compute_hubble_factor (flat ΛCDM) w.r.t. Ωm" begin
        Ωm0 = [0.3156]
        f(p) = Blast.compute_hubble_factor(z0, w0waCDM(Ωm=p[1]))
        check_gradient(f, Ωm0)
    end

    # ======================================================================
    # 4. compute_hubble_factor w.r.t. H0
    #    Tests h = H0/100 propagation through struct and E_z
    # ======================================================================
    @testset "compute_hubble_factor (flat ΛCDM) w.r.t. H0" begin
        H0_0 = [67.27]
        f(p) = Blast.compute_hubble_factor(z0, w0waCDM(H0=p[1]))
        check_gradient(f, H0_0)
    end

    # ======================================================================
    # 5. compute_hubble_factor w.r.t. w0  (dark energy equation of state)
    #    Uses the w0waCDM constructor
    # ======================================================================
    @testset "compute_hubble_factor (w0waCDM) w.r.t. w0" begin
        w0_0 = [-1.0]
        f(p) = Blast.compute_hubble_factor(z0, w0waCDM(w0=p[1]))
        check_gradient(f, w0_0)
    end

    # ======================================================================
    # 5b. compute_hubble_factor w.r.t. wa  (evolution of dark energy EoS)
    #     wa enters E_z via (1+z)^(3(1+w0+wa)) · exp(-3·wa·z/(1+z)).
    #     Every other test uses wa=0; this locks in non-trivial wa differentiability.
    # ======================================================================
    @testset "compute_hubble_factor (w0waCDM) w.r.t. wa" begin
        wa_0 = [0.1]
        f(p) = Blast.compute_hubble_factor(z0, w0waCDM(wa=p[1]))
        check_gradient(f, wa_0)
    end

    # ======================================================================
    # 6. compute_χ w.r.t. Ωm
    #    GL quadrature of 1/E(z',cosmo(Ωm)) — gradient flows through
    #    the Ωm→ωc→Ωcb0 chain inside the integrand
    # ======================================================================
    @testset "compute_χ (flat ΛCDM) w.r.t. Ωm" begin
        Ωm0 = [0.3156]
        f(p) = Blast.compute_χ(z0, w0waCDM(Ωm=p[1]))
        check_gradient(f, Ωm0)
    end

    # ======================================================================
    # 7. get_Ωm w.r.t. Ωm  —  pure field accessor: (ωb + ωc) / h²
    # ======================================================================
    @testset "get_Ωm w.r.t. Ωm" begin
        Ωm0 = [0.3156]
        f(p) = Blast.get_Ωm(w0waCDM(Ωm=p[1]))
        check_gradient(f, Ωm0)
    end

    # ======================================================================
    # 8. get_H0 w.r.t. H0  —  pure field accessor: 100h
    # ======================================================================
    @testset "get_H0 w.r.t. H0" begin
        H0_0 = [67.27]
        f(p) = Blast.get_H0(w0waCDM(H0=p[1]))
        check_gradient(f, H0_0)
    end

    # ======================================================================
    # 9. get_As w.r.t. As  —  involves log/exp round-trip: exp(ln(1e10*As))/1e10
    #    FiniteDifferences is skipped: the adaptive step estimator perturbs As
    #    far outside its valid domain (As ≈ 2e-9 is tiny), causing log() to
    #    receive a negative argument.  All three AD backends agree.
    # ======================================================================
    @testset "get_As w.r.t. As" begin
        As0 = [2.12107e-9]
        f(p) = Blast.get_As(w0waCDM(As=p[1]))
        check_gradient(f, As0; skip_fd=true)
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

        # ── vs H0 ────────────────────────────────────────────────────────
        @testset "sum(bg.z) w.r.t. H0" begin
            f(p) = sum(Background(w0waCDM(H0=p[1])).z)
            check_gradient(f, [67.27]; skip_zygote=true,
                           rtol_ad=1e-3, rtol_fd=1e-3)
        end

        @testset "sum(bg.H) w.r.t. H0" begin
            f(p) = sum(Background(w0waCDM(H0=p[1])).H)
            check_gradient(f, [67.27]; skip_zygote=true,
                           rtol_ad=1e-3, rtol_fd=1e-3)
        end

        @testset "sum(bg.D) w.r.t. H0" begin
            f(p) = sum(Background(w0waCDM(H0=p[1])).D)
            check_gradient(f, [67.27]; skip_zygote=true,
                           rtol_ad=1e-3, rtol_fd=1e-3)
        end

        @testset "sum(bg.f) w.r.t. H0" begin
            f(p) = sum(Background(w0waCDM(H0=p[1])).f)
            check_gradient(f, [67.27]; skip_zygote=true,
                           rtol_ad=1e-3, rtol_fd=1e-3)
        end

        # ── vs Ωm ────────────────────────────────────────────────────────
        @testset "sum(bg.H) w.r.t. Ωm" begin
            f(p) = sum(Background(w0waCDM(Ωm=p[1])).H)
            check_gradient(f, [0.3156]; skip_zygote=true,
                           rtol_ad=1e-3, rtol_fd=1e-3)
        end

        @testset "sum(bg.D) w.r.t. Ωm" begin
            f(p) = sum(Background(w0waCDM(Ωm=p[1])).D)
            check_gradient(f, [0.3156]; skip_zygote=true,
                           rtol_ad=1e-3, rtol_fd=1e-3)
        end

        @testset "sum(bg.f) w.r.t. Ωm" begin
            f(p) = sum(Background(w0waCDM(Ωm=p[1])).f)
            check_gradient(f, [0.3156]; skip_zygote=true,
                           rtol_ad=1e-3, rtol_fd=1e-3)
        end

        # ── vs w0  (dark-energy equation of state) ────────────────────────
        @testset "sum(bg.D) w.r.t. w0" begin
            f(p) = sum(Background(w0waCDM(w0=p[1])).D)
            check_gradient(f, [-1.0]; skip_zygote=true,
                           rtol_ad=1e-3, rtol_fd=1e-3)
        end

        @testset "sum(bg.f) w.r.t. w0" begin
            f(p) = sum(Background(w0waCDM(w0=p[1])).f)
            check_gradient(f, [-1.0]; skip_zygote=true,
                           rtol_ad=1e-3, rtol_fd=1e-3)
        end
    end

end

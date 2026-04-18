# test/test_diff_cosmo.jl  — Category C: Cosmology functions
#
# Tests differentiability of the public cosmo.jl API with respect to:
#   (a) the redshift z
#   (b) cosmological parameters (Ωm, H0, w0) passed through the ΛCDM / w0waCDM
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

include("test_diff_helpers.jl")

@testset "Differentiation: Cosmo functions" begin

    # Fixed cosmology used as captured constant where not differentiated
    cosmo_Λ = ΛCDM()

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
    #    Tests AD through ΛCDM constructor + E_z
    # ======================================================================
    @testset "compute_hubble_factor (ΛCDM) w.r.t. Ωm" begin
        Ωm0 = [0.3156]
        f(p) = Blast.compute_hubble_factor(z0, ΛCDM(Ωm=p[1]))
        check_gradient(f, Ωm0)
    end

    # ======================================================================
    # 4. compute_hubble_factor w.r.t. H0
    #    Tests h = H0/100 propagation through struct and E_z
    # ======================================================================
    @testset "compute_hubble_factor (ΛCDM) w.r.t. H0" begin
        H0_0 = [67.27]
        f(p) = Blast.compute_hubble_factor(z0, ΛCDM(H0=p[1]))
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
    # 6. compute_χ w.r.t. Ωm
    #    GL quadrature of 1/E(z',cosmo(Ωm)) — gradient flows through
    #    the Ωm→ωc→Ωcb0 chain inside the integrand
    # ======================================================================
    @testset "compute_χ (ΛCDM) w.r.t. Ωm" begin
        Ωm0 = [0.3156]
        f(p) = Blast.compute_χ(z0, ΛCDM(Ωm=p[1]))
        check_gradient(f, Ωm0)
    end

    # ======================================================================
    # 7. get_Ωm w.r.t. Ωm  —  pure field accessor: (ωb + ωc) / h²
    # ======================================================================
    @testset "get_Ωm w.r.t. Ωm" begin
        Ωm0 = [0.3156]
        f(p) = Blast.get_Ωm(ΛCDM(Ωm=p[1]))
        check_gradient(f, Ωm0)
    end

    # ======================================================================
    # 8. get_H0 w.r.t. H0  —  pure field accessor: 100h
    # ======================================================================
    @testset "get_H0 w.r.t. H0" begin
        H0_0 = [67.27]
        f(p) = Blast.get_H0(ΛCDM(H0=p[1]))
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
        f(p) = Blast.get_As(ΛCDM(As=p[1]))
        check_gradient(f, As0; skip_fd=true)
    end

end

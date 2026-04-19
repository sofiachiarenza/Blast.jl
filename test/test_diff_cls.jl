# test/test_diff_cls.jl  — Category B: Tullio contractions
#
# Tests differentiability of:
#   w_ell_tullio          (3D / 2D / 1D coefficient variants, projected_matter.jl)
#   _combine_kernels_tullio                                       (cls.jl)
#   _compute_Cℓ_tullio                                            (cls.jl)
#   _compute_Cℓ_rsd_tullio                                        (cls.jl)
#   _limber_contraction                                           (limber.jl)
#
# Design notes
# ─────────────────────────────────────────────────────────────────────────────
# • ForwardDiff is SKIPPED for all Tullio functions.  `@tullio` calls
#   LoopVectorization internally; LoopVectorization does not support
#   ForwardDiff.Dual numbers.  Zygote and Mooncake use the hand-written
#   ChainRules rrules, bypassing the raw Tullio loops.
# • FiniteDifferences evaluates the primal twice with Float64±ε, so it is
#   always safe regardless of the internal implementation.
# • Only arguments that have an rrule pullback (non-NoTangent) are tested.
#   Arguments returning NoTangent() (e.g. T in w_ell_tullio, prefactor /
#   w_χ / w_R / Δχ / χ_grid in the Cℓ functions) are captured as constants.
# ─────────────────────────────────────────────────────────────────────────────

using Test
using Blast
using Random

Random.seed!(1234)

# Small dimensions — large enough to catch contraction bugs, small enough for FD
const nℓ    = 3   # multipoles
const nχ    = 4   # comoving-distance grid points
const nR    = 3   # R = χ₂/χ₁ grid points
const nC    = 5   # Chebyshev coefficient index
const nb    = 2   # probe bins (n_bins_A = n_bins_B = nb for simplicity)
const nχ_lim = 5  # Limber chi-grid points

# ─────────────────────────────────────────────────────────────────────────────
# Shared captured constants (Float64 — not differentiated)
# ─────────────────────────────────────────────────────────────────────────────
const T_ref3 = rand(nℓ, nχ, nR, nC)   # T for 3D/2D/1D w_ell_tullio

const w_χ_ref   = rand(nχ)
const w_R_ref   = rand(nR)
const pref_ref  = rand(nℓ)
const Δχ_ref    = 0.1
const χ_ref     = collect(LinRange(100.0, 3000.0, nχ))
const wlim_ref  = rand(nχ_lim)

# Pre-built kernel (used as ref for _compute_Cℓ_tullio input)
const K_ref  = rand(nb, nb, nχ, nR)
const pmj_ref = rand(nℓ, nχ, nR)

# Pre-built combined kernel (used as ref for _limber_contraction)
const P_ref  = rand(nℓ, nχ_lim)
const K1_ref = rand(nℓ, nχ_lim, nb)
const K2_ref = rand(nℓ, nχ_lim, nb)

# Pre-built W arrays (3D kernels used in RSD function)
const W_A_ref    = rand(nb, nχ, nR)
const W_B_ref    = rand(nb, nχ, nR)
const W_A_r1_ref = rand(nb, nχ)        # 2D (R=1 slice)
const W_B_r1_ref = rand(nb, nχ)        # 2D (R=1 slice)
const pmj02_ref  = rand(nℓ, nχ, nR)
const pmj20_ref  = rand(nℓ, nχ, nR)

# ─────────────────────────────────────────────────────────────────────────────

@testset "Differentiation: Tullio contractions" begin

    # ======================================================================
    # 1. w_ell_tullio  —  @tullio w[i,j,k] := c[l,...] * T[i,j,k,l]
    #    Only c is differentiable (T returns NoTangent()).
    # ======================================================================

    @testset "w_ell_tullio (3D coef) w.r.t. c" begin
        c0 = rand(nC, nχ, nR)
        f(c) = sum(Blast.w_ell_tullio(c, T_ref3))
        check_gradient(f, c0; rtol_ad=1e-8, rtol_fd=1e-5, skip_forward=true)
    end

    @testset "w_ell_tullio (2D coef) w.r.t. c" begin
        c0 = rand(nC, nχ)
        f(c) = sum(Blast.w_ell_tullio(c, T_ref3))
        check_gradient(f, c0; rtol_ad=1e-8, rtol_fd=1e-5, skip_forward=true)
    end

    @testset "w_ell_tullio (1D coef) w.r.t. c" begin
        c0 = rand(nC)
        f(c) = sum(Blast.w_ell_tullio(c, T_ref3))
        check_gradient(f, c0; rtol_ad=1e-8, rtol_fd=1e-5, skip_forward=true)
    end

    # ======================================================================
    # 2. _combine_kernels_tullio(W_A, W_B)
    #    K[i,j,c,r] = W_A_r1[i,c]*W_B[j,c,r] + W_A[i,c,r]*W_B_r1[j,c]
    #    where r1 means the last R-slice.
    # ======================================================================

    @testset "_combine_kernels_tullio w.r.t. W_A" begin
        W_A0 = rand(nb, nχ, nR)
        f(W_A) = sum(Blast._combine_kernels_tullio(W_A, W_B_ref))
        check_gradient(f, W_A0; rtol_ad=1e-8, rtol_fd=1e-5, skip_forward=true)
    end

    @testset "_combine_kernels_tullio w.r.t. W_B" begin
        W_B0 = rand(nb, nχ, nR)
        f(W_B) = sum(Blast._combine_kernels_tullio(W_A_ref, W_B))
        check_gradient(f, W_B0; rtol_ad=1e-8, rtol_fd=1e-5, skip_forward=true)
    end

    # ======================================================================
    # 3. _compute_Cℓ_tullio(K, pmj, w_χ, w_R, prefactor, Δχ, χ_grid)
    #    Cℓ[l,i,j] = prefactor[l]*χ[n]*K[i,j,n,m]*pmj[l,n,m]*w_χ[n]*w_R[m]*Δχ
    #    Pullbacks only for K and pmj.
    # ======================================================================

    @testset "_compute_Cℓ_tullio w.r.t. K" begin
        K0 = rand(nb, nb, nχ, nR)
        f(K) = sum(Blast._compute_Cℓ_tullio(K, pmj_ref, w_χ_ref, w_R_ref, pref_ref, Δχ_ref, χ_ref))
        check_gradient(f, K0; rtol_ad=1e-8, rtol_fd=1e-5, skip_forward=true)
    end

    @testset "_compute_Cℓ_tullio w.r.t. pmj" begin
        pmj0 = rand(nℓ, nχ, nR)
        f(pmj) = sum(Blast._compute_Cℓ_tullio(K_ref, pmj, w_χ_ref, w_R_ref, pref_ref, Δχ_ref, χ_ref))
        check_gradient(f, pmj0; rtol_ad=1e-8, rtol_fd=1e-5, skip_forward=true)
    end

    # ======================================================================
    # 4. _compute_Cℓ_rsd_tullio(W_A_r1, W_B, pmj02, W_A, W_B_r1, pmj20,
    #                            w_χ, w_R, prefactor, Δχ, χ_grid)
    #    Pullbacks for the first six arguments (W_A_r1 through pmj20).
    # ======================================================================

    @testset "_compute_Cℓ_rsd_tullio w.r.t. W_A_r1" begin
        W_A_r1_0 = rand(nb, nχ)
        f(W_A_r1) = sum(Blast._compute_Cℓ_rsd_tullio(
            W_A_r1, W_B_ref, pmj02_ref, W_A_ref, W_B_r1_ref, pmj20_ref,
            w_χ_ref, w_R_ref, pref_ref, Δχ_ref, χ_ref))
        check_gradient(f, W_A_r1_0; rtol_ad=1e-8, rtol_fd=1e-5, skip_forward=true)
    end

    @testset "_compute_Cℓ_rsd_tullio w.r.t. W_B" begin
        W_B0 = rand(nb, nχ, nR)
        f(W_B) = sum(Blast._compute_Cℓ_rsd_tullio(
            W_A_r1_ref, W_B, pmj02_ref, W_A_ref, W_B_r1_ref, pmj20_ref,
            w_χ_ref, w_R_ref, pref_ref, Δχ_ref, χ_ref))
        check_gradient(f, W_B0; rtol_ad=1e-8, rtol_fd=1e-5, skip_forward=true)
    end

    @testset "_compute_Cℓ_rsd_tullio w.r.t. pmj02" begin
        pmj02_0 = rand(nℓ, nχ, nR)
        f(pmj02) = sum(Blast._compute_Cℓ_rsd_tullio(
            W_A_r1_ref, W_B_ref, pmj02, W_A_ref, W_B_r1_ref, pmj20_ref,
            w_χ_ref, w_R_ref, pref_ref, Δχ_ref, χ_ref))
        check_gradient(f, pmj02_0; rtol_ad=1e-8, rtol_fd=1e-5, skip_forward=true)
    end

    @testset "_compute_Cℓ_rsd_tullio w.r.t. W_A" begin
        W_A0 = rand(nb, nχ, nR)
        f(W_A) = sum(Blast._compute_Cℓ_rsd_tullio(
            W_A_r1_ref, W_B_ref, pmj02_ref, W_A, W_B_r1_ref, pmj20_ref,
            w_χ_ref, w_R_ref, pref_ref, Δχ_ref, χ_ref))
        check_gradient(f, W_A0; rtol_ad=1e-8, rtol_fd=1e-5, skip_forward=true)
    end

    @testset "_compute_Cℓ_rsd_tullio w.r.t. W_B_r1" begin
        W_B_r1_0 = rand(nb, nχ)
        f(W_B_r1) = sum(Blast._compute_Cℓ_rsd_tullio(
            W_A_r1_ref, W_B_ref, pmj02_ref, W_A_ref, W_B_r1, pmj20_ref,
            w_χ_ref, w_R_ref, pref_ref, Δχ_ref, χ_ref))
        check_gradient(f, W_B_r1_0; rtol_ad=1e-8, rtol_fd=1e-5, skip_forward=true)
    end

    @testset "_compute_Cℓ_rsd_tullio w.r.t. pmj20" begin
        pmj20_0 = rand(nℓ, nχ, nR)
        f(pmj20) = sum(Blast._compute_Cℓ_rsd_tullio(
            W_A_r1_ref, W_B_ref, pmj02_ref, W_A_ref, W_B_r1_ref, pmj20,
            w_χ_ref, w_R_ref, pref_ref, Δχ_ref, χ_ref))
        check_gradient(f, pmj20_0; rtol_ad=1e-8, rtol_fd=1e-5, skip_forward=true)
    end

    # ======================================================================
    # 5. _limber_contraction(P_term, K1, K2, weights, Δχ)
    #    Cℓ[l,i,j] = P_term[l,m]*K1[l,m,i]*K2[l,m,j]*weights[m]*Δχ
    #    Pullbacks for P_term, K1, K2.
    # ======================================================================

    @testset "_limber_contraction w.r.t. P_term" begin
        P0 = rand(nℓ, nχ_lim)
        f(P) = sum(Blast._limber_contraction(P, K1_ref, K2_ref, wlim_ref, Δχ_ref))
        check_gradient(f, P0; rtol_ad=1e-8, rtol_fd=1e-5, skip_forward=true)
    end

    @testset "_limber_contraction w.r.t. K1" begin
        K1_0 = rand(nℓ, nχ_lim, nb)
        f(K1) = sum(Blast._limber_contraction(P_ref, K1, K2_ref, wlim_ref, Δχ_ref))
        check_gradient(f, K1_0; rtol_ad=1e-8, rtol_fd=1e-5, skip_forward=true)
    end

    @testset "_limber_contraction w.r.t. K2" begin
        K2_0 = rand(nℓ, nχ_lim, nb)
        f(K2) = sum(Blast._limber_contraction(P_ref, K1_ref, K2, wlim_ref, Δχ_ref))
        check_gradient(f, K2_0; rtol_ad=1e-8, rtol_fd=1e-5, skip_forward=true)
    end

end

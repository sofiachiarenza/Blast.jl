# test/test_diff_setup.jl  — Category D: setup-layer and probe helpers
#
# Tests differentiability of:
#   get_PΦ(k, cosmo)               — primordial power spectrum    (setup.jl)
#   get_Tm(pk, k, cosmo)           — matter transfer function     (setup.jl)
#   NLA_model(bg; A, C1)           — NLA intrinsic alignment      (probes.jl)
#   transform_to_R_frame(M, bg)    — (χ,R)–frame resampling       (setup.jl)
#
# None of these functions have hand-written ChainRules rrules — AD traces
# through underlying Julia code (algebraic broadcast / Akima interpolation).
#
# Skip flags used:
#   • get_PΦ w.r.t. k  →  skip_fd: adaptive FD step pushes k negative
#     making (-k)^(ns-1) a complex domain error
#   • get_PΦ w.r.t. As →  excluded (covered by test_diff_cosmo.jl;
#     log-domain issue applies there too)
#   • prepare_nz_matrix →  excluded: quadgk's adaptive step selection
#     is driven by primal values → Dual-mode normalization gradient is wrong
# ─────────────────────────────────────────────────────────────────────────────

using Test
using Blast
using Random
using DifferentiationInterface
using ADTypes

include("test_diff_helpers.jl")
include("test_fixtures.jl")

Random.seed!(4321)

# ─────────────────────────────────────────────────────────────────────────────
# Shared fixtures
# ─────────────────────────────────────────────────────────────────────────────
const cosmo_D   = ΛCDM()          # fixed cosmology
const bg_D      = get_test_bg(cosmo_D)  # full Background (emulator)

# A small k-grid centred well away from k=0 (avoids k^-3 divergence in FD steps)
const nk_D = 8
const k_D  = exp.(LinRange(log(0.01), log(1.0), nk_D))   # log-uniform [0.01, 1] Mpc^-1

# A small (n_z × n_k) power-spectrum matrix for get_Tm tests
const nz_pk = 3
const pk_D  = abs.(rand(nz_pk, nk_D)) .+ 1e-10   # positive, shape (n_z, n_k)

# A small n(z) matrix for prepare_nz_matrix tests
const n_bins_D = 2
const n_zpts   = 12
const z_nz     = LinRange(0.01, 2.0, n_zpts)        # redshift support
const z_grid_D = LinRange(0.05, 1.8,  20)            # evaluation grid

# ─────────────────────────────────────────────────────────────────────────────

@testset "Differentiation: Setup-layer and probe helpers" begin

    # ======================================================================
    # 1. get_PΦ(k, cosmo) w.r.t. k
    #    P_Φ(k) = (9/25)·2π²·As/k³·(k/k_pivot)^(ns-1)
    #    Gradient flows through k^-3 and (k/0.05)^(ns-1) broadcastings.
    # ======================================================================
    # FiniteDifferences skipped: adaptive step estimator can push k negative,
    # making (-k)^(ns-1) a complex-domain error for non-integer ns-1.
    @testset "get_PΦ w.r.t. k" begin
        k0 = copy(k_D)
        f(k) = sum(Blast.get_PΦ(k, cosmo_D))
        check_gradient(f, k0; skip_fd=true)
    end

    # ======================================================================
    # 2. get_PΦ(k, cosmo) w.r.t. ns  (spectral index)
    #    ns enters as the exponent of (k/0.05)^(ns-1); derivative is
    #    P_Φ·log(k/0.05).  Constructing ΛCDM(ns=p[1]) stores nₛ=p[1]
    #    directly (no log/exp), so all backends including FD are safe.
    # ======================================================================
    @testset "get_PΦ w.r.t. ns" begin
        ns0 = [0.9645]
        f(p) = sum(Blast.get_PΦ(k_D, ΛCDM(ns=p[1])))
        check_gradient(f, ns0)
    end

    # ======================================================================
    # 3. get_Tm(pk, k, cosmo) w.r.t. pk
    #    T_m(k,z) = sqrt(pk / P_Φ(k))  — gradient through element-wise sqrt
    # ======================================================================
    @testset "get_Tm w.r.t. pk" begin
        pk0 = copy(pk_D)
        f(pk) = sum(Blast.get_Tm(pk, k_D, cosmo_D))
        check_gradient(f, pk0)
    end

    # ======================================================================
    # 4. NLA_model(bg; A) w.r.t. A
    #    F = -A · C1 · Ωm / bg.D  (element-wise, bg fixed)
    #    Pure product rule — trivially differentiable.
    # ======================================================================
    @testset "NLA_model w.r.t. A" begin
        A0 = [1.72]
        f(p) = sum(Blast.NLA_model(bg_D; A=p[1]))
        check_gradient(f, A0)
    end

    # ======================================================================
    # 5. NLA_model(bg; C1) w.r.t. C1
    # ======================================================================
    @testset "NLA_model w.r.t. C1" begin
        C1_0 = [0.0134]
        f(p) = sum(Blast.NLA_model(bg_D; C1=p[1]))
        check_gradient(f, C1_0)
    end

    # ======================================================================
    # 6. transform_to_R_frame(matrix, bg) w.r.t. matrix
    #    Internally calls akima_interpolation(matrix, Blast.z_lin, x)
    #    (from AbstractCosmologicalEmulators) where x = bg.z_of_χ.(new_χs)
    #    are fixed Float64 query points.
    #    ACE's rrule for akima_interpolation covers ForwardDiff and Mooncake.
    #    Zygote is skipped: reverse-mode fails with `UndefVarError: j` when
    #    tracing through reshape(…, length(Blast.R), …) on generated code.
    #    Mooncake is skipped: the matrix Akima pullback returns a Vector
    #    tangent for the LinRange `t` argument; Mooncake has no
    #    `increment_and_get_rdata!` for (NoFData, RData{LinRange fields}, Vector)
    #    and throws.
    #    ForwardDiff traces through the primal Akima arithmetic using Duals
    #    (it doesn't invoke the rrule), so it works correctly.
    #    (prepare_nz_matrix is excluded: quadgk adaptive step selection is
    #    driven by primal values → Dual-mode normalization gradient is wrong.)
    # ======================================================================
    @testset "transform_to_R_frame w.r.t. matrix" begin
        # Shape: (n_z_lin, n_k) where n_z_lin = length(Blast.z_lin) = 50
    n_k = 6
    nz_lin = length(Blast.z_lin)
    M0 = rand(nz_lin, n_k)

    f(M) = sum(Blast.transform_to_R_frame(M, bg_D))
    check_gradient(f, M0; skip_zygote=true)
end

end

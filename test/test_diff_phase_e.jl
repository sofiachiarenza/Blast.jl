using Test
using Blast
using Random

Random.seed!(4321)

@testset "Differentiation: Phase E (Chebyshev / Limber helpers)" begin

    # ------------------------------------------------------------------
    # Shared Chebyshev/Limber fixtures
    # ------------------------------------------------------------------
    k_min, k_max = 1e-3, 1.0
    chi_min, chi_max = 0.2, 2.0
    K_k, K_chi = 4, 3

    plan_limber = Blast.prepare_chebyshev_plan((log10(k_min), chi_min), (log10(k_max), chi_max), (K_k, K_chi))

    chi_test = collect(LinRange(chi_min, chi_max, 5))
    ell_test = collect(10.0:10.0:40.0)

    @testset "chebinterp_native w.r.t. coefficients" begin
        c0 = rand(6)
        x_grid = collect(LinRange(0.0, 1.0, 8))

        f(c) = sum(Blast.chebinterp_native(c, x_grid, 0.0, 1.0))
        check_gradient(f, c0; skip_zygote=true)
    end

    @testset "build_coeff w.r.t. vals" begin
        plan_1d = Blast.prepare_chebyshev_plan(0.0, 1.0, 8)
        n_nodes = size(plan_1d.fft_plan, 1)
        vals0 = rand(n_nodes)

        f(vals) = sum(Blast.build_coeff(Blast.cϕ, vals, plan_1d).coefs)
        check_gradient(f, vals0)
    end

    @testset "limber_eval w.r.t. coeff matrix" begin
        T_chi = Blast.get_limber_coords_polynomials(plan_limber, chi_test)
        T_k = Blast.get_limber_k_polynomials(plan_limber, ell_test, chi_test; is_log_k=true)
        c0 = rand(size(T_k, 2), size(T_chi, 2))

        f(c) = sum(Blast.limber_eval(c, T_chi, T_k))
        check_gradient(f, c0)
    end

    # ----------------------------------------------------------------------
    # prepare_pk_workspace tensor products
    #   P_ϕTT[k,i,j] = P_ϕ[k] · T_χ1[i,k] · T_χR[i,j,k]
    #   P_ϕT[k,i,j]  = P_ϕ[k]              · T_χR[i,j,k]
    # ----------------------------------------------------------------------
    @testset "_p_phi_TT_tullio" begin
        nk, ni, nj = 5, 4, 3
        P0   = rand(nk)
        T1_0 = rand(ni, nk)
        TR_0 = rand(ni, nj, nk)

        check_gradient(P  -> sum(Blast._p_phi_TT_tullio(P,  T1_0, TR_0)), P0)
        check_gradient(T1 -> sum(Blast._p_phi_TT_tullio(P0, T1,   TR_0)), T1_0)
        check_gradient(TR -> sum(Blast._p_phi_TT_tullio(P0, T1_0, TR)),   TR_0)
    end

    @testset "_p_phi_T_tullio" begin
        nk, ni, nj = 5, 4, 3
        P0   = rand(nk)
        TR_0 = rand(ni, nj, nk)

        check_gradient(P  -> sum(Blast._p_phi_T_tullio(P,  TR_0)), P0)
        check_gradient(TR -> sum(Blast._p_phi_T_tullio(P0, TR)),   TR_0)
    end

    # ----------------------------------------------------------------------
    # Probe kernel contractions
    #   K[b,i] = α · Σ_p H[p] · S[i,p] · χ[i] · (1+z[i]) · nz[b,p] · (χ[p]-χ[i])/χ[p]
    #   (+ (5·s[b,p] - 2) factor for MagnificationBias)
    # Note: χ has dual role (observer χ[i] and source χ[p]) → the rrule adds
    # two contributions. `simpson_matrix` is NoTangent.
    # ----------------------------------------------------------------------
    @testset "_cosmic_shear_kernel_tullio" begin
        n_z = 8; n_bins = 3
        H0       = rand(n_z) .+ 60.0
        χ0       = collect(LinRange(0.1, 3.0, n_z))
        z0       = collect(LinRange(0.01, 2.0, n_z))
        nz0      = rand(n_bins, n_z)
        simpson0 = rand(n_z, n_z)
        Δχ0      = 0.3
        prefac0  = 1.5e-5

        cs_k(H, χ, z, nz, Δχ, prefac) =
            Blast._cosmic_shear_kernel_tullio(H, χ, z, nz, simpson0, Δχ, prefac)

        check_gradient(H  -> sum(cs_k(H,  χ0, z0, nz0, Δχ0, prefac0)), H0)
        check_gradient(χ  -> sum(cs_k(H0, χ,  z0, nz0, Δχ0, prefac0)), χ0)
        check_gradient(z  -> sum(cs_k(H0, χ0, z,  nz0, Δχ0, prefac0)), z0)
        check_gradient(nz -> sum(cs_k(H0, χ0, z0, nz,  Δχ0, prefac0)), nz0)
        check_gradient(d  -> sum(cs_k(H0, χ0, z0, nz0, d[1], prefac0)), [Δχ0])
        check_gradient(p  -> sum(cs_k(H0, χ0, z0, nz0, Δχ0, p[1])),     [prefac0])
    end

    @testset "_magnification_bias_kernel_tullio" begin
        n_z = 8; n_bins = 3
        H0       = rand(n_z) .+ 60.0
        χ0       = collect(LinRange(0.1, 3.0, n_z))
        z0       = collect(LinRange(0.01, 2.0, n_z))
        nz0      = rand(n_bins, n_z)
        s0       = rand(n_bins, n_z)
        simpson0 = rand(n_z, n_z)
        Δχ0      = 0.3
        prefac0  = 1.5e-5

        mb_k(H, χ, z, nz, s, Δχ, prefac) =
            Blast._magnification_bias_kernel_tullio(H, χ, z, nz, s, simpson0, Δχ, prefac)

        check_gradient(H  -> sum(mb_k(H,  χ0, z0, nz0, s0, Δχ0, prefac0)), H0)
        check_gradient(χ  -> sum(mb_k(H0, χ,  z0, nz0, s0, Δχ0, prefac0)), χ0)
        check_gradient(z  -> sum(mb_k(H0, χ0, z,  nz0, s0, Δχ0, prefac0)), z0)
        check_gradient(nz -> sum(mb_k(H0, χ0, z0, nz,  s0, Δχ0, prefac0)), nz0)
        check_gradient(s  -> sum(mb_k(H0, χ0, z0, nz0, s,  Δχ0, prefac0)), s0)
        check_gradient(d  -> sum(mb_k(H0, χ0, z0, nz0, s0, d[1], prefac0)), [Δχ0])
        check_gradient(p  -> sum(mb_k(H0, χ0, z0, nz0, s0, Δχ0, p[1])),     [prefac0])
    end

    # ----------------------------------------------------------------------
    # Pure-broadcast probe kernels (no contractions).
    #
    # Each kernel is a simple element-wise expression extracted from the
    # corresponding compute_kernel! method; AD backends handle these
    # natively without custom rrules. Each testset differentiates with
    # respect to every array / scalar input separately.
    # ----------------------------------------------------------------------

    # NumberCounts: K[b,i] = bias[b,i] · H[i]/c · nz[b,i]
    @testset "_number_counts_kernel" begin
        n_z = 8; n_bins = 3
        bias0 = rand(n_bins, n_z)
        H0    = rand(n_z) .+ 60.0
        nz0   = rand(n_bins, n_z)

        nc_k(bias, H, nz) = Blast._number_counts_kernel(bias, H, nz)

        check_gradient(b  -> sum(nc_k(b,     H0, nz0)), bias0)
        check_gradient(H  -> sum(nc_k(bias0, H,  nz0)), H0)
        check_gradient(nz -> sum(nc_k(bias0, H0, nz)),  nz0)
    end

    # CMBLensing: K[1,i] = prefac · χ[i] · (1+z[i]) · (1 - χ[i]/χ_CMB)
    @testset "_cmb_lensing_kernel" begin
        n_z = 8
        χ0      = collect(LinRange(0.1, 3.0, n_z))
        z0      = collect(LinRange(0.01, 2.0, n_z))
        χ_CMB0  = 9600.0
        prefac0 = 1.5e-5

        cmb_k(χ, z, χ_CMB, prefac) = Blast._cmb_lensing_kernel(χ, z, χ_CMB, prefac)

        check_gradient(χ  -> sum(cmb_k(χ,  z0, χ_CMB0, prefac0)), χ0)
        check_gradient(z  -> sum(cmb_k(χ0, z,  χ_CMB0, prefac0)), z0)
        check_gradient(c  -> sum(cmb_k(χ0, z0, c[1],   prefac0)), [χ_CMB0])
        check_gradient(p  -> sum(cmb_k(χ0, z0, χ_CMB0, p[1])),    [prefac0])
    end

    # RedshiftSpaceDistortions: K[b,i] = f[i] · H[i]/c · nz[b,i]
    @testset "_rsd_kernel" begin
        n_z = 8; n_bins = 3
        f0   = rand(n_z)
        H0   = rand(n_z) .+ 60.0
        nz0  = rand(n_bins, n_z)

        rsd_k(f, H, nz) = Blast._rsd_kernel(f, H, nz)

        check_gradient(f  -> sum(rsd_k(f,  H0, nz0)), f0)
        check_gradient(H  -> sum(rsd_k(f0, H,  nz0)), H0)
        check_gradient(nz -> sum(rsd_k(f0, H0, nz)),  nz0)
    end

    # IntrinsicAlignment — user-supplied A_IA path:
    #   K[b,i] = A_IA[b,i] · H[i]/c · nz[b,i]
    @testset "_ia_kernel (user A_IA)" begin
        n_z = 8; n_bins = 3
        A_IA0 = rand(n_bins, n_z)
        H0    = rand(n_z) .+ 60.0
        nz0   = rand(n_bins, n_z)

        ia_k(A_IA, H, nz) = Blast._ia_kernel(A_IA, H, nz)

        check_gradient(A  -> sum(ia_k(A,     H0, nz0)), A_IA0)
        check_gradient(H  -> sum(ia_k(A_IA0, H,  nz0)), H0)
        check_gradient(nz -> sum(ia_k(A_IA0, H0, nz)),  nz0)
    end

    # IntrinsicAlignment — NLA-model path (composes A_IA construction and
    # kernel multiplication into one pure function):
    #   A_NLA(z) = -A · C1 · Ωm / D(z)
    #   K[b,i]   = A_NLA[i] · H[i]/c · nz[b,i]
    # Gradients tested wrt every independent input (A, C1, Ωm, D, H, nz).
    @testset "_ia_kernel_nla (NLA path)" begin
        n_z = 8; n_bins = 3
        A0   = 1.72
        C10  = 0.0134
        Ωm0  = 0.315
        D0   = collect(LinRange(0.25, 1.0, n_z))       # growth factor
        H0   = rand(n_z) .+ 60.0
        nz0  = rand(n_bins, n_z)

        ia_nla(A, C1, Ωm, D, H, nz) =
            Blast._ia_kernel_nla(A, C1, Ωm, D, H, nz)

        check_gradient(a  -> sum(ia_nla(a[1], C10,   Ωm0,  D0, H0, nz0)), [A0])
        check_gradient(c  -> sum(ia_nla(A0,   c[1],  Ωm0,  D0, H0, nz0)), [C10])
        check_gradient(om -> sum(ia_nla(A0,   C10,   om[1],D0, H0, nz0)), [Ωm0])
        check_gradient(D  -> sum(ia_nla(A0,   C10,   Ωm0,  D,  H0, nz0)), D0)
        check_gradient(H  -> sum(ia_nla(A0,   C10,   Ωm0,  D0, H,  nz0)), H0)
        check_gradient(nz -> sum(ia_nla(A0,   C10,   Ωm0,  D0, H0, nz)),  nz0)
    end

    # IntegratedSachsWolfe: K[1,i] = prefac · H[i] · (1 - f[i])
    @testset "_isw_kernel" begin
        n_z = 8
        H0      = rand(n_z) .+ 60.0
        f0      = rand(n_z)
        prefac0 = 3.0e-12

        isw_k(H, f, prefac) = Blast._isw_kernel(H, f, prefac)

        check_gradient(H -> sum(isw_k(H,  f0, prefac0)), H0)
        check_gradient(f -> sum(isw_k(H0, f,  prefac0)), f0)
        check_gradient(p -> sum(isw_k(H0, f0, p[1])),    [prefac0])
    end

    # PrimordialNonGaussianity:
    #   b_phi[b,i] = 2 · 1.686 · (bias[b,i] - p)
    #   K[b,i]     = H[i]/c · f_NL · b_phi[b,i] · nz[b,i]
    @testset "_png_kernel" begin
        n_z = 8; n_bins = 3
        H0    = rand(n_z) .+ 60.0
        bias0 = rand(n_bins, n_z) .+ 1.0
        nz0   = rand(n_bins, n_z)
        f_NL0 = 0.5
        p0    = 1.0

        png_k(H, f_NL, bias, p, nz) = Blast._png_kernel(H, f_NL, bias, p, nz)

        check_gradient(H    -> sum(png_k(H,  f_NL0, bias0, p0, nz0)), H0)
        check_gradient(fNL  -> sum(png_k(H0, fNL[1], bias0, p0, nz0)), [f_NL0])
        check_gradient(bias -> sum(png_k(H0, f_NL0, bias, p0, nz0)), bias0)
        check_gradient(pp   -> sum(png_k(H0, f_NL0, bias0, pp[1], nz0)), [p0])
        check_gradient(nz   -> sum(png_k(H0, f_NL0, bias0, p0, nz)),  nz0)
    end
end

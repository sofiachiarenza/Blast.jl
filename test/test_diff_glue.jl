# test/test_diff_glue.jl
#
# End-to-end AD tests for the get_Cℓ glue code. The individual primitives
# (kernels, Tullio contractions, Chebyshev helpers) are all covered by
# test_diff_phase_e.jl / test_diff_setup.jl / test_diff_cosmo.jl /
# test_diff_cls.jl. This file verifies their COMPOSITION — that the
# mutating orchestration layer (evaluate_components!, get_Cℓ with its
# limber correction and Chebyshev interpolation loop) remains
# differentiable end-to-end.
#
# Zygote is skipped throughout (mutation of mutable struct fields is not
# supported). ForwardDiff is exercised where possible, but may strip
# derivative info at struct construction boundaries where fields are
# typed `Matrix{Float64}` — failures are reported not silently worked
# around.
using Test
using Blast
using Random

Random.seed!(9876)

@testset "Differentiation: glue code (end-to-end get_Cℓ)" begin

    # ---------------------------------------------------------------------
    # Shared fixtures (bg, pk grids, Plans, PS, W all precomputed outside
    # every `f` so only the probe kernel construction + Cℓ composition
    # are inside the gradient).
    # ---------------------------------------------------------------------
    cosmo = get_test_cosmo()
    bg    = get_test_bg(cosmo)

    nk              = length(Blast.k_cheb)
    nz_lin          = length(Blast.z_lin)
    nz_cheb         = length(Blast.z_cheb)
    nk_limber       = length(Blast.k_limber)
    pk_grid         = ones(nz_lin,  nk)
    pk_limber_lin   = ones(nz_cheb, nk_limber)
    pk_limber_nonl  = ones(nz_cheb, nk_limber)

    nz_bins = 1
    nz_grid = ones(nz_bins, 50)
    z_grid  = collect(LinRange(0.0, 3.6, 50))

    # ---------------------------------------------------------------------
    # Testset 1 — GalaxyClustering auto-spectrum, wrt NumberCounts.bias.
    # Exercises:
    #   NumberCounts(bias=...) construction (struct-field conversion)
    #   → evaluate_components! → compute_kernel!(::NumberCounts) mutation
    #   → get_Cℓ(GC, PS, W, bg, Plans) top-level orchestrator:
    #     → combine_kernels → get_kernel_array → akima_interpolation
    #     → _compute_Cℓ_tullio
    #     → get_limber_correction → _limber_contraction
    #     → get_limber_Cℓ → _limber_contraction
    #     → chebyshev_decomposition + chebinterp_native loop into Cℓ_final.
    # ---------------------------------------------------------------------
    @testset "GalaxyClustering auto-spectrum wrt bias" begin
        gc_template = Blast.GalaxyClustering(
            δ = Blast.NumberCounts(nz=nz_grid, z=z_grid,
                                    bias=ones(nz_bins, length(bg.z))))
        W, Plans = Blast.SetUp(gc_template)
        Blast.evaluate_components!(gc_template, bg)
        PS = Blast.prepare_pk_workspace(Plans, pk_grid,
                                         pk_limber_lin, pk_limber_nonl, bg)
        W = Blast.compute_w(W, PS)

        function f(bias_flat)
            bias = reshape(bias_flat, nz_bins, length(bg.z))
            gc = Blast.GalaxyClustering(
                δ = Blast.NumberCounts(nz=nz_grid, z=z_grid, bias=bias))
            Blast.evaluate_components!(gc, bg)
            cls = Blast.get_Cℓ(Blast.full_ℓ_range, gc, PS, W, bg, Plans)
            return sum(cls)
        end

        bias0 = ones(nz_bins * length(bg.z))
        # Post-Phase-A: struct is NumberCounts{T<:Real}, so ForwardDiff
        # Duals thread through bias → Kernel → Cℓ without being stripped.
        # Agreement with Mooncake is ~1e-12 (machine precision).
        check_gradient(f, bias0;
                       skip_forward=false,
                       skip_zygote=true,
                       rtol_fd=1e-4)
    end

    # ---------------------------------------------------------------------
    # Testset 2 — GalaxyClustering auto-spectrum, wrt NumberCounts.nz
    # (raw redshift distribution, before normalization).
    # Additionally exercises:
    #   check_and_normalize! → prepare_nz_matrix (Integrals.jl QuadGKJL
    #   native VJP) → struct-field assignment into nz_norm.
    # bias is held fixed; only nz varies.
    # ---------------------------------------------------------------------
    @testset "GalaxyClustering auto-spectrum wrt nz (Mooncake smoke test)" begin
        # The nz path is transitively exercised by Test 1 (nz_norm is an
        # input of _number_counts_kernel). This testset explicitly covers
        # the extra link: `check_and_normalize!` mutation of
        # Component.nz_norm via prepare_nz_matrix (Integrals.jl QuadGKJL).
        #
        # No FD comparison: d(sum Cℓ)/d(nz_i) spans 8+ orders of
        # magnitude (ell_prefactor × limber_factor both ~1/ℓ² cancel
        # across modes), and central_fdm(5, 1) cannot resolve the small
        # entries reliably. Mooncake is independently validated against
        # FD for bias in Test 1 through the same orchestration path.
        n_nz = 50
        nz0_mat = ones(nz_bins, n_nz)
        bias_fixed = ones(nz_bins, length(bg.z))

        gc_template = Blast.GalaxyClustering(
            δ = Blast.NumberCounts(nz=copy(nz0_mat), z=z_grid, bias=bias_fixed))
        W, Plans = Blast.SetUp(gc_template)
        Blast.evaluate_components!(gc_template, bg)
        PS = Blast.prepare_pk_workspace(Plans, pk_grid,
                                         pk_limber_lin, pk_limber_nonl, bg)
        W = Blast.compute_w(W, PS)

        function f_nz(nz_flat)
            nz = reshape(nz_flat, nz_bins, n_nz)
            gc = Blast.GalaxyClustering(
                δ = Blast.NumberCounts(nz=nz, z=z_grid, bias=bias_fixed))
            Blast.evaluate_components!(gc, bg)
            cls = Blast.get_Cℓ(Blast.full_ℓ_range, gc, PS, W, bg, Plans)
            return sum(cls)
        end

        nz0 = vec(nz0_mat)
        g_mc = DifferentiationInterface.gradient(f_nz, AutoMooncake(), nz0)

        # Smoke: Mooncake returns a finite gradient vector of the right
        # size, with at least one non-zero entry (proves the tape traced
        # all the way through the composed pipeline rather than bailing
        # out silently with zero cotangents).
        @test length(g_mc) == length(nz0)
        @test all(isfinite, g_mc)
        @test any(!iszero, g_mc)
    end

    # ---------------------------------------------------------------------
    # Testset 3 — WeakLensing auto-spectrum, wrt IntrinsicAlignment.A.
    # Exercises:
    #   CosmicShear kernel (_cosmic_shear_kernel_tullio, has rrule)
    #   + IA-NLA kernel path (_ia_kernel_nla, pure broadcast)
    #   + evaluate_components!(::WeakLensing)
    #   + top-level get_Cℓ(::WL) orchestrator (γγ, γI, Iγ, II sums)
    #   + get_limber_correction/get_limber_Cℓ for WL
    #   + the χ²_app loop inside get_kernel_array (shear branch)
    # ---------------------------------------------------------------------
    @testset "WeakLensing auto-spectrum wrt IntrinsicAlignment.A" begin
        nz_shear = ones(nz_bins, 50)
        nz_ia    = ones(nz_bins, 50)

        wl_template = Blast.WeakLensing(
            γ = Blast.CosmicShear(nz=nz_shear, z=z_grid),
            IA = Blast.IntrinsicAlignment(nz=nz_ia, z=z_grid, A=1.72))
        W, Plans = Blast.SetUp(wl_template)
        Blast.evaluate_components!(wl_template, bg)
        PS = Blast.prepare_pk_workspace(Plans, pk_grid,
                                         pk_limber_lin, pk_limber_nonl, bg)
        W = Blast.compute_w(W, PS)

        function f_A(A_vec)
            wl = Blast.WeakLensing(
                γ  = Blast.CosmicShear(nz=nz_shear, z=z_grid),
                IA = Blast.IntrinsicAlignment(nz=nz_ia, z=z_grid, A=A_vec[1]))
            Blast.evaluate_components!(wl, bg)
            cls = Blast.get_Cℓ(Blast.full_ℓ_range, wl, PS, W, bg, Plans)
            return sum(cls)
        end

        A0 = [1.72]
        # Post-Phase-A: IntrinsicAlignment{T<:Real} with A::T allows FW.
        check_gradient(f_A, A0;
                       skip_forward=false,
                       skip_zygote=true,
                       rtol_fd=1e-4)
    end

    # ---------------------------------------------------------------------
    # Testset 4 — GC × WL cross-spectrum, wrt NumberCounts.bias.
    # Exercises cross-probe composition:
    #   get_Cℓ(ℓ, ::GalaxyClustering, ::WeakLensing, ...) orchestrator
    #   → cross get_Cℓ dispatch (δγ, δI)
    #   → get_limber_correction(G, L, pk) (cross version)
    #   → get_limber_Cℓ(G, L, pk) (cross version)
    #   → final chebyshev interpolation over (nbins_A × nbins_B) pairs.
    # WL has both CosmicShear and IntrinsicAlignment; GC has only δ.
    # ---------------------------------------------------------------------
    @testset "GC × WL cross-spectrum wrt NumberCounts.bias" begin
        nz_shear = ones(nz_bins, 50)
        nz_ia    = ones(nz_bins, 50)

        gc_template = Blast.GalaxyClustering(
            δ = Blast.NumberCounts(nz=nz_grid, z=z_grid,
                                    bias=ones(nz_bins, length(bg.z))))
        wl_template = Blast.WeakLensing(
            γ  = Blast.CosmicShear(nz=nz_shear, z=z_grid),
            IA = Blast.IntrinsicAlignment(nz=nz_ia, z=z_grid, A=1.72))

        # SetUp takes (G, L) to build the cross workspace.
        W, Plans = Blast.SetUp(gc_template, wl_template)
        Blast.evaluate_components!(gc_template, bg)
        Blast.evaluate_components!(wl_template, bg)
        PS = Blast.prepare_pk_workspace(Plans, pk_grid,
                                         pk_limber_lin, pk_limber_nonl, bg)
        W = Blast.compute_w(W, PS)

        function f_cross(bias_flat)
            bias = reshape(bias_flat, nz_bins, length(bg.z))
            gc = Blast.GalaxyClustering(
                δ = Blast.NumberCounts(nz=nz_grid, z=z_grid, bias=bias))
            wl = Blast.WeakLensing(
                γ  = Blast.CosmicShear(nz=nz_shear, z=z_grid),
                IA = Blast.IntrinsicAlignment(nz=nz_ia, z=z_grid, A=1.72))
            Blast.evaluate_components!(gc, bg)
            Blast.evaluate_components!(wl, bg)
            cls = Blast.get_Cℓ(Blast.full_ℓ_range, gc, wl, PS, W, bg, Plans)
            return sum(cls)
        end

        bias0 = ones(nz_bins * length(bg.z))
        # Cross-spectrum with bias=1 against shear has gradients ~1e-11
        # element-wise; Mooncake and FD agree to 4-5 digits per entry but
        # the norm-based isapprox trips on near-zero entries with atol=0.
        # atol=1e-12 catches that without loosening rtol.
        # Post-Phase-A: ForwardDiff now works through the GC×WL cross path.
        check_gradient(f_cross, bias0;
                       skip_forward=false,
                       skip_zygote=true,
                       rtol_fd=1e-4,
                       atol=1e-12)
    end

    # ---------------------------------------------------------------------
    # Testset 5 — GC auto-spectrum wrt cosmology parameters (Ωm, H0, w0).
    #
    # The entire cosmology-dependent pipeline is inside the differentiated
    # function: Background ODE (H, χ, D, f from NumCosmo solver), kernel
    # evaluation on the fresh bg grid, full prepare_pk_workspace (includes
    # get_PΦ, get_Tm, transform_to_R_frame, build_coeff FFT path,
    # chebyshev_decomposition for the limber P_δ), compute_w mutation of
    # W, and the top-level get_Cℓ orchestrator. This is the maximal
    # Mooncake composition target for Blast.
    #
    # Plans are built outside f (FFT plans — cosmo-invariant). We build a
    # minimal W inline inside f so no mutation leaks across calls.
    #
    # Runtime: first-time Mooncake tape compile is ~8 min; ForwardDiff (post
    # Phase B) ~15 s; FD reference ~25 s. Agreement vs FD is ~1e-8 relative
    # for both AD backends (rtol_fd=1e-4 is very loose).
    # ---------------------------------------------------------------------
    @testset "GC auto-spectrum wrt cosmology (Ωm, H0, w0)" begin
        gc_template = Blast.GalaxyClustering(
            δ = Blast.NumberCounts(nz=nz_grid, z=z_grid,
                                    bias=ones(nz_bins, 200)))
        _, Plans = Blast.SetUp(gc_template)

        # length(bg.z) is cosmo-invariant in w0waCDM — safe to size the
        # bias matrix against any reference Background.
        bg_ref     = Background(w0waCDM(Ωm=0.3156))
        bias_fixed = ones(nz_bins, length(bg_ref.z))

        function f_cosmo(p)
            cosmo = w0waCDM(Ωm=p[1], H0=p[2], w0=p[3])
            bg    = Background(cosmo)
            gc = Blast.GalaxyClustering(
                δ = Blast.NumberCounts(nz=nz_grid, z=z_grid, bias=bias_fixed))
            Blast.evaluate_components!(gc, bg)
            # Fresh W inline — only δδ is active for GC-only-δ.
            W_local = Blast.ProjectedMatterDensity(
                w_2_00_ϕTT = Blast.w_2_00_ϕTT())
            PS = Blast.prepare_pk_workspace(Plans, pk_grid,
                                             pk_limber_lin, pk_limber_nonl, bg)
            W_local = Blast.compute_w(W_local, PS)
            cls = Blast.get_Cℓ(Blast.full_ℓ_range, gc, PS, W_local, bg, Plans)
            return sum(cls)
        end

        p0 = [0.3156, 67.27, -1.0]
        # ForwardDiff through cosmology is enabled as of Phase B:
        #   - PowerSpectrum{T} and 17 w_*{T} structs parametrized
        #   - evaluate_components! promotes component T against bg eltype
        #   - prepare_nz_matrix uses u-substitution for Dual bounds
        #     (QuadGK can't form Dual kronrod weights directly)
        # FW gradient vs FD agrees to ~1e-8 relative — well inside rtol_fd=1e-4.
        # FW runs in ~15s vs Mooncake's ~8min tape compile.
        check_gradient(f_cosmo, p0;
                       skip_zygote=true,
                       rtol_fd=1e-4)
    end

end

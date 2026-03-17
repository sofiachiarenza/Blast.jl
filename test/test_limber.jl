using Test
using Blast
using FastChebInterp
using StaticArrays

@testset "Limber: Polynomial helpers" begin
    cosmo = get_test_cosmo()
    bg = get_test_bg(cosmo)

    gc = Blast.GalaxyClustering(δ = Blast.NumberCounts())
    _, Plans = Blast.SetUp(gc)

    T_chi = Blast.get_limber_coords_polynomials(Plans.plan_limber, bg.z)
    T_k = Blast.get_limber_k_polynomials(Plans.plan_limber, Blast.full_ℓ_range, Blast.χ; is_log_k=true)

    @test size(T_chi, 1) == length(bg.z)
    @test size(T_k, 1) == length(Blast.full_ℓ_range)
    @test size(T_k, 3) == length(Blast.χ)
    @test all(isfinite, T_chi)
    @test all(isfinite, T_k)

    c = rand(size(T_k, 2), size(T_chi, 2))
    eval_grid = Blast.limber_eval(c, T_chi, T_k)
    @test size(eval_grid) == (length(Blast.full_ℓ_range), length(Blast.χ))
    @test all(isfinite, eval_grid)
end

@testset "Limber: FastChebInterp agreement" begin
    k_min, k_max = 1e-4, 80.0
    chi_min, chi_max = 26.0, 7000.0
    K_k, K_chi = 160, 49

    nodes_k = Blast.chebpoints(K_k, log10(k_min), log10(k_max))
    nodes_chi = Blast.chebpoints(K_chi, chi_min, chi_max)

    f_mock(lk, chi) = exp(-chi / 2000) * sin(lk) * (10^lk)^-2
    f_grid = [f_mock(lk, chi) for lk in nodes_k, chi in nodes_chi]

    lb = [log10(k_min), chi_min]
    ub = [log10(k_max), chi_max]
    interp_old = chebinterp(f_grid, lb, ub)

    ℓ_test = collect(range(2.0, 2000.0, length=100))
    chi_test = collect(LinRange(26.0, 7000.0, 96))

    result_old = zeros(length(ℓ_test), length(chi_test))
    for j in eachindex(chi_test)
        for i in eachindex(ℓ_test)
            k_val = (ℓ_test[i] + 0.5) / chi_test[j]
            result_old[i, j] = interp_old(SVector(log10(k_val), chi_test[j]))
        end
    end

    plan = Blast.prepare_chebyshev_plan((log10(k_min), chi_min), (log10(k_max), chi_max), (K_k, K_chi))
    c = Blast.chebyshev_decomposition(plan, f_grid)
    T_chi = Blast.get_limber_coords_polynomials(plan, chi_test)
    T_k = Blast.get_limber_k_polynomials(plan, ℓ_test, chi_test; is_log_k=true)
    result_new = Blast.limber_eval(c, T_chi, T_k)

    @test result_old ≈ result_new rtol=1e-6 atol=1e-7
    @test maximum(abs.(result_old .- result_new)) < 1e-6
end

@testset "Limber: Kernels and contractions" begin
    cosmo = get_test_cosmo()
    bg = get_test_bg(cosmo)

    nz = ones(1, 50)
    z_nz = collect(LinRange(0.0, 3.6, 50))

    gc = Blast.GalaxyClustering(
        δ = Blast.NumberCounts(nz=nz, z=z_nz, bias=ones(1, length(bg.z))),
        μ = Blast.MagnificationBias(nz=nz, z=z_nz, s=fill(0.4, 1, length(bg.z)))
    )
    wl = Blast.WeakLensing(γ = Blast.CosmicShear(nz=nz, z=z_nz))
    cmb = Blast.CMB(κ = Blast.CMBLensing(), ISW = Blast.IntegratedSachsWolfe())

    W, Plans = Blast.SetUp(gc, wl, cmb)
    Blast.evaluate_components!(gc, bg)
    Blast.evaluate_components!(wl, bg)
    Blast.evaluate_components!(cmb, bg)

    pk_grid = ones(length(Blast.z_lin), length(Blast.k_cheb))
    pk_limber = ones(length(Blast.z_cheb), length(Blast.k_limber))
    PS = Blast.prepare_pk_workspace(Plans, pk_grid, pk_limber, pk_limber, bg)
    Blast.compute_w!(W, PS)

    K_gc = Blast.get_limber_kernel(gc)
    K_wl = Blast.get_limber_kernel(wl)
    K_cmb = Blast.get_limber_kernel(cmb)

    @test size(K_gc) == (length(Blast.full_ℓ_range), length(Blast.χ), 1)
    @test size(K_wl) == (length(Blast.full_ℓ_range), length(Blast.χ), 1)
    @test size(K_cmb) == (length(Blast.full_ℓ_range), length(Blast.χ), 1)
    @test all(isfinite, K_gc)
    @test all(isfinite, K_wl)
    @test all(isfinite, K_cmb)

    corr_auto = Blast.get_limber_correction(gc, PS)
    corr_cross = Blast.get_limber_correction(gc, wl, PS)
    @test size(corr_auto) == (length(Blast.ℓ_nonlimber), 1, 1)
    @test size(corr_cross) == (length(Blast.ℓ_nonlimber), 1, 1)
    @test all(isfinite, corr_auto)
    @test all(isfinite, corr_cross)

    cl_auto = Blast.get_limber_Cℓ(gc, PS)
    cl_cross = Blast.get_limber_Cℓ(gc, wl, PS)
    @test size(cl_auto) == (length(Blast.full_ℓ_range) - length(Blast.ℓ_nonlimber), 1, 1)
    @test size(cl_cross) == (length(Blast.full_ℓ_range) - length(Blast.ℓ_nonlimber), 1, 1)
    @test all(isfinite, cl_auto)
    @test all(isfinite, cl_cross)

    @test Blast.get_limber_Cℓ(gc, wl, PS) ≈ Blast.get_limber_Cℓ(wl, gc, PS)
    @test Blast.get_limber_correction(gc, wl, PS) ≈ Blast.get_limber_correction(wl, gc, PS)
end
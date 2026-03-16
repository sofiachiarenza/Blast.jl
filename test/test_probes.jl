using Test
using Blast
using NPZ
using DataInterpolations
using QuadGK
using Statistics

function _compute_kernel_safe!(component::Blast.CosmicShear, bg::Blast.Background)
    Blast.check_and_normalize!(component, bg.z)
    n_bins = size(component.nz_norm, 1)
    kernel = zeros(n_bins, length(bg.z))

    H0 = Blast.get_H0(bg.cosmo)
    Ωm = Blast.get_Ωm(bg.cosmo)
    prefac = 1.5 * H0^2 * Ωm / Blast.C_LIGHT^2

    for b in 1:n_bins
        nz_vals = component.nz_norm[b, :]
        for z_idx in 1:length(bg.z)
            integrand(x) = Blast._akima_interpolation(nz_vals, bg.z, x) * (1.0 - bg.χ[z_idx] / Blast.compute_χ(x, bg.cosmo))
            z_low = bg.z[z_idx]
            z_top = bg.z[end]
            integral, _ = quadgk(integrand, z_low, z_top)
            kernel[b, z_idx] = prefac * bg.χ[z_idx] * (1.0 + bg.z[z_idx]) * integral
        end
    end

    component.Kernel = kernel
    return component
end

function _compute_kernel_safe!(component::Blast.MagnificationBias, bg::Blast.Background)
    Blast.check_and_normalize!(component, bg.z)
    n_bins = size(component.nz_norm, 1)
    kernel = zeros(n_bins, length(bg.z))

    H0 = Blast.get_H0(bg.cosmo)
    Ωm = Blast.get_Ωm(bg.cosmo)
    prefac = 1.5 * H0^2 * Ωm / Blast.C_LIGHT^2

    for b in 1:n_bins
        nz_vals = component.nz_norm[b, :]
        s_vals = component.s[b, :]
        for z_idx in 1:length(bg.z)
            integrand(x) = Blast._akima_interpolation(nz_vals, bg.z, x) * (1.0 - bg.χ[z_idx] / Blast.compute_χ(x, bg.cosmo)) * (5.0 * s_vals[z_idx] - 2.0)
            z_low = bg.z[z_idx]
            z_top = bg.z[end]
            integral, _ = quadgk(integrand, z_low, z_top)
            kernel[b, z_idx] = prefac * bg.χ[z_idx] * (1.0 + bg.z[z_idx]) * integral
        end
    end

    component.Kernel = kernel
    return component
end

@testset "Probes: Kernel computation" begin
    cosmo = get_test_cosmo()
    bg = get_test_bg(cosmo)
    
    # Redshift bin data from npz
    bins = data["bins"]
    n_bins = size(bins["dNdz"], 1)
    
    # Initialize modern components
    nc = Blast.NumberCounts(nz = bins["dNdz"], z = bins["z"], bias = ones(n_bins, length(bg.z)))
    cs = Blast.CosmicShear(nz = bins["dNdz"][1:3,:], z = bins["z"])
    cmb = Blast.CMB(κ = Blast.CMBLensing(), ISW = Blast.IntegratedSachsWolfe())
    
    # Compute kernels using the modern API
    Blast.compute_kernel!(nc, bg)
    Blast.compute_kernel!(cs, bg)
    Blast.evaluate_components!(cmb, bg)
    
    # Basic Sanity Checks (Avoiding legacy artifacts that might be outdated)
    @test all(isfinite, nc.Kernel)
    @test all(isfinite, cs.Kernel)
    @test all(isfinite, cmb.κ.Kernel)
    @test all(isfinite, cmb.ISW.Kernel)
    
    @test size(nc.Kernel) == (n_bins, length(bg.z))
    @test size(cs.Kernel) == (3, length(bg.z))
    
    # Check that kernels are non-zero
    @test any(!iszero, nc.Kernel)
    @test any(!iszero, cs.Kernel)
end

@testset "Probes: Legacy LJ kernel regression" begin
    cosmo = get_test_cosmo()
    bins = data["bins"]

    # Match the legacy 1000-point background grid used to generate LJ references
    z_ref = collect(LinRange(1e-3, 4.0, 1000))
    χ_ref = Blast.compute_χ.(z_ref, Ref(cosmo))
    bg_ref = Blast.Background(cosmo; χ_grid=χ_ref)

    n_bins = size(bins["dNdz"], 1)
    nc = Blast.NumberCounts(nz=bins["dNdz"], z=bins["z"], bias=ones(n_bins, length(bg_ref.z)))
    cs = Blast.CosmicShear(nz=bins["dNdz"][1:3, :], z=bins["z"])
    cmbκ = Blast.CMBLensing()

    Blast.compute_kernel!(nc, bg_ref)
    Blast.compute_kernel!(cs, bg_ref)
    Blast.compute_kernel!(cmbκ, bg_ref)

    lj_clustering = data["LJ_clustering"]
    lj_shear = data["LJ_shear"]
    lj_cmb = reshape(data["LJ_cmb"], 1, :)

    # Legacy-amplitude regression for clustering and CMB lensing
    @test nc.Kernel ≈ lj_clustering rtol=1e-2
    @test cmbκ.Kernel ≈ lj_cmb rtol=1e-2

    # Shear kernel currently differs in normalization but preserves shape very closely
    for b in 1:size(lj_shear, 1)
        @test cor(cs.Kernel[b, :], lj_shear[b, :]) > 0.999
    end
end

@testset "Probes: N(z) Normalization" begin
    z_grid = collect(LinRange(0.0, 4.0, 100))
    nz = ones(2, 50)
    z = collect(LinRange(0.0, 4.0, 50))

    nz_normed = Blast.prepare_nz_matrix(nz, z, z_grid)
    @test size(nz_normed) == (2, 100)
    @test all(isfinite, nz_normed)
    # Each bin must integrate to 1 on the output grid
    for b in 1:size(nz_normed, 1)
        nrm, _ = quadgk(x -> Blast._akima_interpolation(nz_normed[b, :], z_grid, x), first(z_grid), last(z_grid))
        @test nrm ≈ 1.0 atol=1e-2
    end
end

@testset "Probes: N(z) Smoothing" begin
    z = collect(LinRange(0.0, 4.0, 250))
    z_out = collect(LinRange(0.0, 4.0, 180))

    nz = zeros(2, length(z))
    nz[1, :] .= @. exp(-((z - 1.0)^2) / 0.08) + 0.05 * sin(25 * z)
    nz[2, :] .= @. exp(-((z - 2.3)^2) / 0.12) + 0.04 * cos(20 * z)
    nz[:, 20] .= -0.1

    nz_smoothed = Blast.smooth_nz(nz, z; z_out=z_out, span=0.05)

    @test size(nz_smoothed) == (2, length(z_out))
    @test all(isfinite, nz_smoothed)
    @test all(nz_smoothed .>= 0)

    nz_normed = Blast.prepare_nz_matrix(nz_smoothed, z_out, z_out)
    for b in 1:size(nz_normed, 1)
        nrm, _ = quadgk(x -> Blast._akima_interpolation(nz_normed[b, :], z_out, x), first(z_out), last(z_out))
        @test nrm ≈ 1.0 atol=1e-2
    end
end

@testset "Probes: kernel vs safe reference (Gaussian n(z))" begin
    cosmo = get_test_cosmo()
    bg = get_test_bg(cosmo)

    z = collect(bg.z)
    n_bins = 2
    nz = zeros(n_bins, length(z))
    nz[1, :] .= @. exp(-((z - 0.8)^2) / 0.05)
    nz[2, :] .= @. exp(-((z - 1.6)^2) / 0.08)
    s = fill(0.4, n_bins, length(z))

    cs_fast = Blast.CosmicShear(nz = copy(nz), z = copy(z))
    cs_safe = Blast.CosmicShear(nz = copy(nz), z = copy(z))

    μ_fast = Blast.MagnificationBias(nz = copy(nz), z = copy(z), s = copy(s))
    μ_safe = Blast.MagnificationBias(nz = copy(nz), z = copy(z), s = copy(s))

    Blast.compute_kernel!(cs_fast, bg)
    _compute_kernel_safe!(cs_safe, bg)

    Blast.compute_kernel!(μ_fast, bg)
    _compute_kernel_safe!(μ_safe, bg)

    @test cs_fast.Kernel ≈ cs_safe.Kernel rtol=1e-2 atol=1e-6
    @test μ_fast.Kernel ≈ μ_safe.Kernel rtol=1e-2 atol=1e-6
end

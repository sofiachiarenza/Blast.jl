using Test
using Blast
using NPZ
using DataInterpolations

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

@testset "Probes: N(z) Normalization" begin
    z_grid = LinRange(0.0, 4.0, 100)
    nz = ones(2, 50) 
    z = LinRange(0.0, 4.0, 50)
    
    nz_normed = Blast.prepare_nz_matrix(nz, z, z_grid)
    @test size(nz_normed) == (2, 100)
    @test all(isfinite, nz_normed)
end

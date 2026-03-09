using Test
using Blast
using NPZ

@testset "Cls: Pipeline Integration" begin
    cosmo = get_test_cosmo()
    bg = get_test_bg(cosmo)
    
    # Setup dummy power spectrum
    nk = length(Blast.k_cheb)
    nz_lin = length(Blast.z_lin)
    pk_grid = ones(nz_lin, nk)
    
    nz_cheb = length(Blast.z_cheb)
    nk_limber = length(Blast.k_limber)
    pk_limber_lin = ones(nz_cheb, nk_limber)
    pk_limber_nonlin = ones(nz_cheb, nk_limber)
    
    # Initialize GC probe
    # Simple 1-bin constant n(z)
    nz = ones(1, 50)
    z_nz = LinRange(0.0, 3.6, 50)
    gc = Blast.GalaxyClustering(δ = Blast.NumberCounts(nz=nz, z=z_nz, bias=ones(1, length(bg.z))))
    
    W, Plans = Blast.SetUp(gc)
    Blast.evaluate_components!(gc, bg)
    PS = Blast.prepare_pk_workspace(Plans, pk_grid, pk_limber_lin, pk_limber_nonlin, bg)
    
    # CRITICAL: Compute projected matter coefficients
    Blast.compute_w!(W, PS)
    
    # Compute Cls: get_Cℓ(ℓ, G, Pk, W, bg, P)
    cls = Blast.get_Cℓ(Blast.full_ℓ_range, gc, PS, W, bg, Plans)
    
    @test size(cls) == (length(Blast.full_ℓ_range), 1, 1)
    @test all(isfinite, cls)
    @test all(x -> x >= 0, cls) # Auto-spectra must be non-negative
    
    # Test correction logic
    corr = Blast.get_limber_correction(gc, PS)
    @test size(corr) == (length(Blast.ℓ_nonlimber), 1, 1)
end

@testset "Cls: Cross-probe consistency" begin
    cosmo = get_test_cosmo()
    bg = get_test_bg(cosmo)
    
    # Dummy grids
    pk_grid = ones(length(Blast.z_lin), length(Blast.k_cheb))
    pk_limber = ones(length(Blast.z_cheb), length(Blast.k_limber))
    
    # Probes
    gc = Blast.GalaxyClustering(δ = Blast.NumberCounts(nz=ones(1,50), z=LinRange(0,3.6,50), bias=ones(1,length(bg.z))))
    wl = Blast.WeakLensing(γ = Blast.CosmicShear(nz=ones(1,50), z=LinRange(0,3.6,50)))
    
    W, Plans = Blast.SetUp(gc, wl)
    Blast.evaluate_components!(gc, bg)
    Blast.evaluate_components!(wl, bg)
    PS = Blast.prepare_pk_workspace(Plans, pk_grid, pk_limber, pk_limber, bg)
    
    # CRITICAL: Compute projected matter coefficients
    Blast.compute_w!(W, PS)
    
    # Cross Cl: get_Cℓ(ℓ, G, L, Pk, W, bg, P)
    cl_xl = Blast.get_Cℓ(Blast.full_ℓ_range, gc, wl, PS, W, bg, Plans)
    
    @test size(cl_xl) == (length(Blast.full_ℓ_range), 1, 1)
    @test all(isfinite, cl_xl)
end

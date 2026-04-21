using Test
using Blast
using NPZ
using DataInterpolations

@testset "Setup: Power Spectrum Workspace" begin
    cosmo = get_test_cosmo()
    bg = get_test_bg(cosmo)
    
    # Create dummy power spectrum grids
    # Linear P(k,z) matches (z_lin, k_cheb)
    nk = length(Blast.k_cheb)
    nz_lin = length(Blast.z_lin)
    pk_grid = rand(nz_lin, nk)
    
    # Limber grids match (z_cheb, k_limber)
    nz_cheb = length(Blast.z_cheb)
    nk_limber = length(Blast.k_limber)
    pk_limber_lin = rand(nz_cheb, nk_limber)
    pk_limber_nonlin = pk_limber_lin .* 1.1 
    
    # Initialize Plans
    gc = Blast.GalaxyClustering(δ = Blast.NumberCounts())
    W, Plans = Blast.SetUp(gc)
    
    # Compute coefficients
    PS = Blast.prepare_pk_workspace(Plans, pk_grid, pk_limber_lin, pk_limber_nonlin, bg)
    
    @test PS.cϕTT isa Blast.cϕTT
    @test all(isfinite, PS.ΔP_limber)
    @test all(isfinite, PS.Pδ_limber)
    
    # Check dimensions
    @test size(PS.ΔP_limber) == (length(Blast.full_ℓ_range), length(Blast.χ))
    @test size(PS.Pδ_limber) == (length(Blast.full_ℓ_range), length(Blast.χ))
end

@testset "Setup: Probe SetUp logic" begin
    # Test that SetUp correctly returns ProjectedMatterDensity with non-nothing fields
    gc = Blast.GalaxyClustering(δ = Blast.NumberCounts(), RSD = Blast.RedshiftSpaceDistortions())
    W, Plans = Blast.SetUp(gc)
    
    # Modern fields use the w_L_MN_basis format
    @test !isnothing(W.w_2_00_ϕTT) # delta
    @test !isnothing(W.w_2_02_ϕTT) # RSD_A
    @test isnothing(W.w_0_00_ϕTT)  # μ_A (not requested)
    @test isnothing(W.w_minus2_00_ϕTT) # γ (no WL probe)
    
    # Test combined setup
    wl = Blast.WeakLensing(γ = Blast.CosmicShear())
    W_comb, Plans_comb = Blast.SetUp(gc, wl)
    @test !isnothing(W_comb.w_2_00_ϕTT)
    @test !isnothing(W_comb.w_minus2_00_ϕTT) # γ is now present
end

@testset "Setup: compute_w! populates active fields" begin
    # Direct unit test for compute_w!: for every non-nothing w_*_ϕTT/ϕT/ϕ
    # field in W, its coefs array must be finite and non-zero after
    # compute_w! runs. Catches dispatch/mutation leaks that would pass
    # through get_Cℓ (because the result still broadcasts to 0. on missing
    # components) but fail silently here.
    cosmo = get_test_cosmo()
    bg = get_test_bg(cosmo)

    nk = length(Blast.k_cheb)
    nz_lin = length(Blast.z_lin)
    pk_grid = rand(nz_lin, nk) .+ 1.0
    pk_limber_lin = rand(length(Blast.z_cheb), length(Blast.k_limber)) .+ 1.0
    pk_limber_nonlin = pk_limber_lin .* 1.05

    # Activate all GC-related components (δ, RSD, μ, PNG) plus WL + CMB
    # so every w_* basis gets exercised.
    gc = Blast.GalaxyClustering(
        δ = Blast.NumberCounts(),
        RSD = Blast.RedshiftSpaceDistortions(),
        μ = Blast.MagnificationBias(),
        PNG = Blast.PrimordialNonGaussianity())
    wl = Blast.WeakLensing(γ = Blast.CosmicShear(),
                            IA = Blast.IntrinsicAlignment())
    cmb = Blast.CMB(κ = Blast.CMBLensing(),
                     ISW = Blast.IntegratedSachsWolfe())

    W, Plans = Blast.SetUp(gc, wl, cmb)
    PS = Blast.prepare_pk_workspace(Plans, pk_grid, pk_limber_lin,
                                     pk_limber_nonlin, bg)
    Blast.compute_w!(W, PS)

    # Enumerate every field of W; skip the ones left as `nothing` for this
    # probe configuration. Every active field must have finite, non-zero
    # coefs after compute_w!.
    n_active = 0
    for name in fieldnames(typeof(W))
        field = getfield(W, name)
        field === nothing && continue
        n_active += 1
        @test all(isfinite, field.w)
        @test any(!iszero, field.w)
    end
    # Sanity: a fully-active probe configuration must light up several fields.
    @test n_active >= 5
end

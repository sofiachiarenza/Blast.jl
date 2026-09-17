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
    W = Blast.compute_w(W, PS)
    
    # Compute Cls: get_Cℓ(ℓ, G, Pk, W, bg, P)
    cls = Blast.get_Cℓ(Blast.full_ℓ_range, gc, PS, W, bg, Plans)
    
    @test size(cls) == (length(Blast.full_ℓ_range), 1, 1)
    @test all(isfinite, cls)
    @test all(x -> x >= 0, cls) # Auto-spectra must be non-negative
    
    # Test correction logic
    corr = Blast.get_limber_correction(gc, PS)
    @test size(corr) == (length(Blast.ℓ_nonlimber), 1, 1)
    @test all(isfinite, corr)
end

@testset "Cls: GC auto with RSD active" begin
    # End-to-end exercise of the RSD dispatch chain:
    #   SetUp(GC with δ+RSD) → allocates w_2_02_ϕTT / w_2_20_ϕTT / w_2_22_ϕTT
    #   evaluate_components! → populates δ.Kernel and RSD.Kernel (needs bg.f)
    #   compute_w → fills all RSD-related projected-matter coefficients
    #   get_Cℓ      → routes through _compute_Cℓ_rsd_tullio for the δ×RSD,
    #                 RSD×δ, and RSD×RSD sectors.
    cosmo = get_test_cosmo()
    bg = get_test_bg(cosmo)

    nk = length(Blast.k_cheb)
    nz_lin = length(Blast.z_lin)
    pk_grid = ones(nz_lin, nk)

    nz_cheb = length(Blast.z_cheb)
    nk_limber = length(Blast.k_limber)
    pk_limber = ones(nz_cheb, nk_limber)

    nz = ones(1, 50)
    z_nz = LinRange(0.0, 3.6, 50)

    # δ-only reference
    gc_δ = Blast.GalaxyClustering(
        δ = Blast.NumberCounts(nz=nz, z=z_nz, bias=ones(1, length(bg.z))))
    W_δ, Plans_δ = Blast.SetUp(gc_δ)
    Blast.evaluate_components!(gc_δ, bg)
    PS_δ = Blast.prepare_pk_workspace(Plans_δ, pk_grid, pk_limber, pk_limber, bg)
    W_δ = Blast.compute_w(W_δ, PS_δ)
    cls_δ = Blast.get_Cℓ(Blast.full_ℓ_range, gc_δ, PS_δ, W_δ, bg, Plans_δ)

    # δ+RSD
    gc_rsd = Blast.GalaxyClustering(
        δ = Blast.NumberCounts(nz=nz, z=z_nz, bias=ones(1, length(bg.z))),
        RSD = Blast.RedshiftSpaceDistortions(nz=nz, z=z_nz))
    W_rsd, Plans_rsd = Blast.SetUp(gc_rsd)

    # The RSD workspace fields must be allocated
    @test !isnothing(W_rsd.w_2_02_ϕTT)
    @test !isnothing(W_rsd.w_2_20_ϕTT)
    @test !isnothing(W_rsd.w_2_22_ϕTT)

    Blast.evaluate_components!(gc_rsd, bg)
    @test all(isfinite, gc_rsd.RSD.Kernel)
    @test any(!iszero, gc_rsd.RSD.Kernel)

    PS_rsd = Blast.prepare_pk_workspace(Plans_rsd, pk_grid, pk_limber, pk_limber, bg)
    W_rsd = Blast.compute_w(W_rsd, PS_rsd)
    cls_rsd = Blast.get_Cℓ(Blast.full_ℓ_range, gc_rsd, PS_rsd, W_rsd, bg, Plans_rsd)

    @test size(cls_rsd) == (length(Blast.full_ℓ_range), 1, 1)
    @test all(isfinite, cls_rsd)
    @test all(x -> x >= 0, cls_rsd)

    # RSD must actually contribute — δ-only and δ+RSD spectra must differ.
    # The RSD sectors enter as -δRSD - RSDδ + RSDRSD, so for an auto-spectrum
    # the RSD contribution is non-trivial and gc_rsd ≠ gc_δ.
    @test !(cls_rsd ≈ cls_δ)
end

@testset "Cls: GC auto with PNG active" begin
    # End-to-end exercise of the PNG dispatch chain:
    #   SetUp(GC with δ+PNG) → allocates w_2_00_ϕT / w_2_00_ϕT_R1 / w_2_00_ϕ
    #   evaluate_components! → populates PNG.Kernel via _png_kernel
    #                          (uses bg.H, bias, f_NL, p, nz_norm)
    #   compute_w → fills PNG-related projected-matter coefficients
    #   get_Cℓ      → routes through the δ×PNG / PNG×δ / PNG×PNG dispatchers.
    cosmo = get_test_cosmo()
    bg = get_test_bg(cosmo)

    nk = length(Blast.k_cheb)
    nz_lin = length(Blast.z_lin)
    pk_grid = ones(nz_lin, nk)

    nz_cheb = length(Blast.z_cheb)
    nk_limber = length(Blast.k_limber)
    pk_limber = ones(nz_cheb, nk_limber)

    nz = ones(1, 50)
    z_nz = LinRange(0.0, 3.6, 50)
    n_bg = length(bg.z)

    # δ-only reference
    gc_δ = Blast.GalaxyClustering(
        δ = Blast.NumberCounts(nz=nz, z=z_nz, bias=ones(1, n_bg)))
    W_δ, Plans_δ = Blast.SetUp(gc_δ)
    Blast.evaluate_components!(gc_δ, bg)
    PS_δ = Blast.prepare_pk_workspace(Plans_δ, pk_grid, pk_limber, pk_limber, bg)
    W_δ = Blast.compute_w(W_δ, PS_δ)
    cls_δ = Blast.get_Cℓ(Blast.full_ℓ_range, gc_δ, PS_δ, W_δ, bg, Plans_δ)

    # δ+PNG
    gc_png = Blast.GalaxyClustering(
        δ = Blast.NumberCounts(nz=nz, z=z_nz, bias=ones(1, n_bg)),
        PNG = Blast.PrimordialNonGaussianity(nz=nz, z=z_nz,
                                              bias=ones(1, n_bg),
                                              f_NL=1.0, p=0.0))
    W_png, Plans_png = Blast.SetUp(gc_png)

    # The PNG workspace fields must be allocated
    @test !isnothing(W_png.w_2_00_ϕT)
    @test !isnothing(W_png.w_2_00_ϕT_R1)
    @test !isnothing(W_png.w_2_00_ϕ)

    Blast.evaluate_components!(gc_png, bg)
    @test all(isfinite, gc_png.PNG.Kernel)
    @test any(!iszero, gc_png.PNG.Kernel)

    PS_png = Blast.prepare_pk_workspace(Plans_png, pk_grid, pk_limber, pk_limber, bg)
    W_png = Blast.compute_w(W_png, PS_png)
    cls_png = Blast.get_Cℓ(Blast.full_ℓ_range, gc_png, PS_png, W_png, bg, Plans_png)

    @test size(cls_png) == (length(Blast.full_ℓ_range), 1, 1)
    @test all(isfinite, cls_png)
    @test all(x -> x >= 0, cls_png)

    # PNG must actually contribute — δ-only and δ+PNG spectra must differ.
    @test !(cls_png ≈ cls_δ)
end

@testset "Cls: GC auto with ARF active (no PNG)" begin
    # Exercises the new _combine_gc_nonlimber_arf branch (16-arg: base + arf),
    # distinct from the pre-existing 16-arg PNG branch despite equal arity
    # (see the dispatch-collision note in src/cls.jl).
    cosmo = get_test_cosmo()
    bg = get_test_bg(cosmo)

    nk = length(Blast.k_cheb)
    nz_lin = length(Blast.z_lin)
    pk_grid = ones(nz_lin, nk)

    nz_cheb = length(Blast.z_cheb)
    nk_limber = length(Blast.k_limber)
    pk_limber = ones(nz_cheb, nk_limber)

    nz = ones(1, 50)
    z_nz = LinRange(0.0, 3.6, 50)
    n_bg = length(bg.z)

    gc_δ = Blast.GalaxyClustering(δ = Blast.NumberCounts(nz=nz, z=z_nz, bias=ones(1, n_bg)))
    W_δ, Plans_δ = Blast.SetUp(gc_δ)
    Blast.evaluate_components!(gc_δ, bg)
    PS_δ = Blast.prepare_pk_workspace(Plans_δ, pk_grid, pk_limber, pk_limber, bg)
    W_δ = Blast.compute_w(W_δ, PS_δ)
    cls_δ = Blast.get_Cℓ(Blast.full_ℓ_range, gc_δ, PS_δ, W_δ, bg, Plans_δ)

    gc_arf = Blast.GalaxyClustering(
        δ = Blast.NumberCounts(nz=nz, z=z_nz, bias=ones(1, n_bg)),
        ARF = Blast.AngularRedshiftFluctuations(nz=nz, z=z_nz, bias=ones(1, n_bg), zbar=[1.0]))
    W_arf, Plans_arf = Blast.SetUp(gc_arf)
    @test !isnothing(W_arf.w_2_00_ϕTT)  # shared with δ — no new tensor needed

    Blast.evaluate_components!(gc_arf, bg)
    @test all(isfinite, gc_arf.ARF.Kernel)
    @test any(!iszero, gc_arf.ARF.Kernel)

    PS_arf = Blast.prepare_pk_workspace(Plans_arf, pk_grid, pk_limber, pk_limber, bg)
    W_arf = Blast.compute_w(W_arf, PS_arf)
    cls_arf = Blast.get_Cℓ(Blast.full_ℓ_range, gc_arf, PS_arf, W_arf, bg, Plans_arf)

    @test size(cls_arf) == (length(Blast.full_ℓ_range), 1, 1)
    @test all(isfinite, cls_arf)
    @test !(cls_arf ≈ cls_δ)
end

@testset "Cls: GC auto with ARF and PNG combined" begin
    # Exercises the new 25-arg _combine_gc_nonlimber_arf_png branch.
    cosmo = get_test_cosmo()
    bg = get_test_bg(cosmo)

    pk_grid = ones(length(Blast.z_lin), length(Blast.k_cheb))
    pk_limber = ones(length(Blast.z_cheb), length(Blast.k_limber))

    nz = ones(1, 50)
    z_nz = LinRange(0.0, 3.6, 50)
    n_bg = length(bg.z)

    gc = Blast.GalaxyClustering(
        δ = Blast.NumberCounts(nz=nz, z=z_nz, bias=ones(1, n_bg)),
        PNG = Blast.PrimordialNonGaussianity(nz=nz, z=z_nz, bias=ones(1, n_bg), f_NL=1.0, p=0.0),
        ARF = Blast.AngularRedshiftFluctuations(nz=nz, z=z_nz, bias=ones(1, n_bg), zbar=[1.0]))
    W, Plans = Blast.SetUp(gc)
    Blast.evaluate_components!(gc, bg)
    PS = Blast.prepare_pk_workspace(Plans, pk_grid, pk_limber, pk_limber, bg)
    W = Blast.compute_w(W, PS)
    cls = Blast.get_Cℓ(Blast.full_ℓ_range, gc, PS, W, bg, Plans)

    @test size(cls) == (length(Blast.full_ℓ_range), 1, 1)
    @test all(isfinite, cls)
end

@testset "Cls: ARF pairwise dispatch matches a pre-weighted NumberCounts" begin
    # Strong correctness check for the new ARF get_Cℓ methods: since
    # ARF.Kernel == (z - zbar) .* δ.Kernel exactly (test_probes.jl), an ARF
    # component and a "fake" NumberCounts whose bias is pre-multiplied by
    # (z - zbar) must produce bit-identical kernels — and therefore
    # bit-identical Cℓ against any third partner, PROVIDED the new ARF
    # dispatch methods route to the same weight tensors as the corresponding
    # NumberCounts ones. This validates (a) no typo'd weight-tensor name in
    # the new ARF get_Cℓ methods and (b) the ARF kernel values, using the
    # already-trusted NumberCounts pathway as ground truth. It does not
    # independently re-derive the +/- sign of each term inside
    # _combine_gc_nonlimber_arf (that sign is asserted by analogy to δ's
    # sign, since ARF and δ share the same weight-tensor class).
    cosmo = get_test_cosmo()
    bg = get_test_bg(cosmo)
    z = collect(bg.z)
    n_bg = length(z)

    nz = ones(1, 50)
    z_nz = LinRange(0.0, 3.6, 50)
    bias = ones(1, n_bg)
    zbar = [1.0]
    s = fill(0.4, 1, n_bg)

    arf = Blast.AngularRedshiftFluctuations(nz=nz, z=z_nz, bias=copy(bias), zbar=copy(zbar))
    Blast.compute_kernel!(arf, bg)

    bias_eff = @. bias * (z' - zbar)
    fake_nc = Blast.NumberCounts(nz=nz, z=z_nz, bias=bias_eff)
    Blast.compute_kernel!(fake_nc, bg)
    @test arf.Kernel ≈ fake_nc.Kernel

    μ = Blast.MagnificationBias(nz=nz, z=z_nz, s=copy(s))
    Blast.compute_kernel!(μ, bg)

    W = Blast.ProjectedMatterDensity(w_2_00_ϕTT=Blast.w_2_00_ϕTT(), w_0_00_ϕTT=Blast.w_0_00_ϕTT())
    pk_grid = ones(length(Blast.z_lin), length(Blast.k_cheb))
    pk_limber = ones(length(Blast.z_cheb), length(Blast.k_limber))
    gc_for_plans = Blast.GalaxyClustering(δ = Blast.NumberCounts(nz=nz, z=z_nz, bias=bias))
    _, Plans = Blast.SetUp(gc_for_plans)
    PS = Blast.prepare_pk_workspace(Plans, pk_grid, pk_limber, pk_limber, bg)
    W = Blast.compute_w(W, PS)

    cl_arf_μ = Blast.get_Cℓ(arf, μ, W, bg)
    cl_fake_μ = Blast.get_Cℓ(fake_nc, μ, W, bg)
    @test cl_arf_μ ≈ cl_fake_μ

    cl_μ_arf = Blast.get_Cℓ(μ, arf, W, bg)
    cl_μ_fake = Blast.get_Cℓ(μ, fake_nc, W, bg)
    @test cl_μ_arf ≈ cl_μ_fake
end

@testset "Cls: KG (CMB×GC) — ARF crosses κ only, not ISW" begin
    cosmo = get_test_cosmo()
    bg = get_test_bg(cosmo)

    pk_grid = ones(length(Blast.z_lin), length(Blast.k_cheb))
    pk_limber = ones(length(Blast.z_cheb), length(Blast.k_limber))

    nz = ones(1, 50)
    z_nz = LinRange(0.0, 3.6, 50)
    n_bg = length(bg.z)

    function _kg_cls(; with_arf::Bool, with_isw::Bool)
        δ = Blast.NumberCounts(nz=nz, z=z_nz, bias=ones(1, n_bg))
        arf = with_arf ? Blast.AngularRedshiftFluctuations(nz=nz, z=z_nz, bias=ones(1, n_bg), zbar=[1.0]) : nothing
        gc = Blast.GalaxyClustering(; δ, ARF=arf)
        cmb = with_isw ? Blast.CMB(κ=Blast.CMBLensing(), ISW=Blast.IntegratedSachsWolfe()) : Blast.CMB(κ=Blast.CMBLensing())

        W, Plans = Blast.SetUp(gc, cmb)
        Blast.evaluate_components!(gc, bg)
        Blast.evaluate_components!(cmb, bg)
        PS = Blast.prepare_pk_workspace(Plans, pk_grid, pk_limber, pk_limber, bg)
        W = Blast.compute_w(W, PS)
        return Blast.get_Cℓ(Blast.full_ℓ_range, cmb, gc, PS, W, bg, Plans)
    end

    cls_κ_noarf = _kg_cls(with_arf=false, with_isw=false)
    cls_κ_arf = _kg_cls(with_arf=true, with_isw=false)
    @test size(cls_κ_arf) == (length(Blast.full_ℓ_range), 1, 1)
    @test all(isfinite, cls_κ_arf)
    # κ×ARF must actually contribute.
    @test !(cls_κ_arf ≈ cls_κ_noarf)

    cls_κisw_noarf = _kg_cls(with_arf=false, with_isw=true)
    cls_κisw_arf = _kg_cls(with_arf=true, with_isw=true)

    # The "turn on ISW" increment must be identical whether or not ARF is
    # active — this is only true if ARF contributes zero to the ISW cross
    # (i.e. no hidden Cℓ_TARF term), since the only ISW-driven contribution
    # in either case is Tδ.
    increment_noarf = cls_κisw_noarf .- cls_κ_noarf
    increment_arf = cls_κisw_arf .- cls_κ_arf
    @test increment_arf ≈ increment_noarf
end

@testset "Cls: ARF is excluded from the GC×WL (cosmic shear) cross-spectrum" begin
    cosmo = get_test_cosmo()
    bg = get_test_bg(cosmo)

    pk_grid = ones(length(Blast.z_lin), length(Blast.k_cheb))
    pk_limber = ones(length(Blast.z_cheb), length(Blast.k_limber))

    nz = ones(1, 50)
    z_nz = LinRange(0.0, 3.6, 50)
    n_bg = length(bg.z)

    wl = Blast.WeakLensing(γ = Blast.CosmicShear(nz=nz, z=z_nz))

    gc_noarf = Blast.GalaxyClustering(δ = Blast.NumberCounts(nz=nz, z=z_nz, bias=ones(1, n_bg)))
    gc_arf = Blast.GalaxyClustering(
        δ = Blast.NumberCounts(nz=nz, z=z_nz, bias=ones(1, n_bg)),
        ARF = Blast.AngularRedshiftFluctuations(nz=nz, z=z_nz, bias=ones(1, n_bg), zbar=[1.0]))

    W, Plans = Blast.SetUp(gc_arf, wl)
    Blast.evaluate_components!(gc_noarf, bg)
    Blast.evaluate_components!(gc_arf, bg)
    Blast.evaluate_components!(wl, bg)
    PS = Blast.prepare_pk_workspace(Plans, pk_grid, pk_limber, pk_limber, bg)
    W = Blast.compute_w(W, PS)

    cls_noarf = Blast.get_Cℓ(Blast.full_ℓ_range, gc_noarf, wl, PS, W, bg, Plans)
    cls_arf = Blast.get_Cℓ(Blast.full_ℓ_range, gc_arf, wl, PS, W, bg, Plans)

    # ARF must be completely absent from this cross-spectrum by design.
    @test cls_arf ≈ cls_noarf
end

@testset "Cls: CMB×GC with ISW active" begin
    # End-to-end exercise of the ISW dispatch chain through the non-limber
    # CMB×GC cross-spectrum. test_limber.jl covers ISW at the Limber-kernel
    # level; this test covers the `get_Cℓ(ℓ, ::CMB, ::GC, ...)` orchestrator,
    # which sums κ×δ, κ×{RSD,μ,PNG}, ISW×δ, ISW×{RSD,μ,PNG} contributions.
    # Only κ and ISW + δ are active here, so the Cℓ_Tδ (ISW×δ) term is the
    # only extra contribution over a pure κ×GC run.
    cosmo = get_test_cosmo()
    bg = get_test_bg(cosmo)

    nk = length(Blast.k_cheb)
    nz_lin = length(Blast.z_lin)
    pk_grid = ones(nz_lin, nk)

    nz_cheb = length(Blast.z_cheb)
    nk_limber = length(Blast.k_limber)
    pk_limber = ones(nz_cheb, nk_limber)

    nz = ones(1, 50)
    z_nz = LinRange(0.0, 3.6, 50)
    n_bg = length(bg.z)

    gc = Blast.GalaxyClustering(
        δ = Blast.NumberCounts(nz=nz, z=z_nz, bias=ones(1, n_bg)))

    # κ-only reference
    cmb_κ = Blast.CMB(κ = Blast.CMBLensing())
    W_κ, Plans_κ = Blast.SetUp(gc, cmb_κ)
    Blast.evaluate_components!(gc, bg)
    Blast.evaluate_components!(cmb_κ, bg)
    PS_κ = Blast.prepare_pk_workspace(Plans_κ, pk_grid, pk_limber, pk_limber, bg)
    W_κ = Blast.compute_w(W_κ, PS_κ)
    cls_κ_gc = Blast.get_Cℓ(Blast.full_ℓ_range, cmb_κ, gc, PS_κ, W_κ, bg, Plans_κ)

    # κ + ISW
    cmb_full = Blast.CMB(κ = Blast.CMBLensing(), ISW = Blast.IntegratedSachsWolfe())
    W_full, Plans_full = Blast.SetUp(gc, cmb_full)
    Blast.evaluate_components!(cmb_full, bg)

    @test all(isfinite, cmb_full.ISW.Kernel)
    @test any(!iszero, cmb_full.ISW.Kernel)

    PS_full = Blast.prepare_pk_workspace(Plans_full, pk_grid, pk_limber, pk_limber, bg)
    W_full = Blast.compute_w(W_full, PS_full)
    cls_full = Blast.get_Cℓ(Blast.full_ℓ_range, cmb_full, gc, PS_full, W_full, bg, Plans_full)

    @test size(cls_full) == (length(Blast.full_ℓ_range), 1, 1)
    @test all(isfinite, cls_full)

    # ISW must actually contribute — κ-only and κ+ISW cross-spectra must differ.
    @test !(cls_full ≈ cls_κ_gc)

    # Cross-direction dispatcher exists: get_Cℓ(ℓ, ::GC, ::CMB, ...)
    cls_flip = Blast.get_Cℓ(Blast.full_ℓ_range, gc, cmb_full, PS_full, W_full, bg, Plans_full)
    @test size(cls_flip) == (length(Blast.full_ℓ_range), 1, 1)
    @test all(isfinite, cls_flip)
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
    W = Blast.compute_w(W, PS)

    # Cross Cl: get_Cℓ(ℓ, G, L, Pk, W, bg, P)
    cl_gc_wl = Blast.get_Cℓ(Blast.full_ℓ_range, gc, wl, PS, W, bg, Plans)
    cl_wl_gc = Blast.get_Cℓ(Blast.full_ℓ_range, wl, gc, PS, W, bg, Plans)

    @test size(cl_gc_wl) == (length(Blast.full_ℓ_range), 1, 1)
    @test all(isfinite, cl_gc_wl)
    # Symmetry: C_ℓ(A,B) == C_ℓ(B,A)
    @test cl_gc_wl ≈ cl_wl_gc
end

@testset "Cls: GC×CMB and WL×CMB symmetry" begin
    # C_ℓ(A, B) = C_ℓ(B, A) for GC×CMB and WL×CMB. Previously only the
    # GC×WL direction was guarded by an equality test — this extends that
    # guard to the two remaining cross-probe dispatch pairs.
    cosmo = get_test_cosmo()
    bg = get_test_bg(cosmo)

    pk_grid = ones(length(Blast.z_lin), length(Blast.k_cheb))
    pk_limber = ones(length(Blast.z_cheb), length(Blast.k_limber))

    nz = ones(1, 50)
    z_nz = LinRange(0.0, 3.6, 50)
    n_bg = length(bg.z)

    gc = Blast.GalaxyClustering(
        δ = Blast.NumberCounts(nz=nz, z=z_nz, bias=ones(1, n_bg)))
    wl = Blast.WeakLensing(γ = Blast.CosmicShear(nz=nz, z=z_nz))
    cmb = Blast.CMB(κ = Blast.CMBLensing(),
                     ISW = Blast.IntegratedSachsWolfe())

    W, Plans = Blast.SetUp(gc, wl, cmb)
    Blast.evaluate_components!(gc, bg)
    Blast.evaluate_components!(wl, bg)
    Blast.evaluate_components!(cmb, bg)
    PS = Blast.prepare_pk_workspace(Plans, pk_grid, pk_limber, pk_limber, bg)
    W = Blast.compute_w(W, PS)

    # GC × CMB symmetry
    cl_gc_cmb = Blast.get_Cℓ(Blast.full_ℓ_range, gc, cmb, PS, W, bg, Plans)
    cl_cmb_gc = Blast.get_Cℓ(Blast.full_ℓ_range, cmb, gc, PS, W, bg, Plans)
    @test size(cl_gc_cmb) == (length(Blast.full_ℓ_range), 1, 1)
    @test all(isfinite, cl_gc_cmb)
    @test cl_gc_cmb ≈ cl_cmb_gc

    # WL × CMB symmetry
    cl_wl_cmb = Blast.get_Cℓ(Blast.full_ℓ_range, wl, cmb, PS, W, bg, Plans)
    cl_cmb_wl = Blast.get_Cℓ(Blast.full_ℓ_range, cmb, wl, PS, W, bg, Plans)
    @test size(cl_wl_cmb) == (length(Blast.full_ℓ_range), 1, 1)
    @test all(isfinite, cl_wl_cmb)
    @test cl_wl_cmb ≈ cl_cmb_wl
end

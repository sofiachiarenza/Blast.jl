using Test, Blast, LinearAlgebra
isdefined(@__MODULE__,:trapezoid_radial_plan) || include("test_radial_quadrature.jl")
isdefined(@__MODULE__,:manufactured_direct_bank) || include("bspline_backend_fixture.jl")

@testset "Fractional multipole bank: public preparation and pair projection" begin
    base=trapezoid_radial_plan([26.,60.,250.,800.,1600.,2800.,4200.,5500.,7000.])
    ell=[2.,2.95,80.01,219.]
    q=Blast.RadialQuadraturePlan(base.chi,base.weights,base.suffix;
        ell,R=[.2,.7,1.],R_weights=[.5,.5,0.])
    bg=Blast.Background(characterization_cosmology();radial_plan=q)
    probes=build_characterization_probes(bg);G,L=probes.galaxy,probes.lensing
    G.δ.ell_prefactor=ones(4);L.γ.ell_prefactor=sqrt.(Blast.factorial_frac.(ell))
    bank=manufactured_direct_bank(;axes=(;ell,chi=q.chi,R=q.R))
    backend=Blast.CubicBSplineBackend(bank)
    @test_throws ArgumentError Blast.SetUp(G,L;backend,radial_plan=base)
    W,P=Blast.SetUp(G,L;backend,radial_plan=q)
    grid=Blast.nonlimber_grid(P);pmm=characterization_pmm()
    ps=Blast.prepare_pk_workspace(P,fixture_pk(grid),pmm.pk_limber_linear,pmm.pk_limber_nonlinear,bg;pk_grid=grid)
    projected=Blast.compute_w(W,ps)
    component=projected.w_0_00_ϕTT
    @test size(component.w)==(4,9,3)
    ka,kb=Blast.get_kernel_array(G.δ,bg),Blast.get_kernel_array(L.γ,bg)
    actual=Blast.compute_Cℓ(G.δ,L.γ,component,bg)
    pref=2/π .* sqrt.(Blast.factorial_frac.(ell))
    @test size(actual)==(4,size(ka,1),size(kb,1))
    for l in 1:4,i in axes(ka,1),j in axes(kb,1)
        expected=sum(pref[l]*q.chi[c]*q.weights[c]*q.R_weights[r]*component.w[l,c,r]*
            (ka[i,c,end]*kb[j,c,r]+ka[i,c,r]*kb[j,c,end]) for c in 1:9,r in 1:3)
        @test actual[l,i,j]≈expected rtol=3e-12
    end
end

@testset "Public nonuniform radial path: banks, 17 blocks, spectra and mismatches" begin
    q=trapezoid_radial_plan([26.,60.,250.,800.,1600.,2800.,4200.,5500.,7000.])
    bg=Blast.Background(characterization_cosmology();radial_plan=q)
    probes=build_characterization_probes(bg)
    G,L,C=probes.galaxy,probes.lensing,probes.cmb
    bank=manufactured_direct_bank(;axes=(ell=Blast.ℓ_nonlimber,chi=q.chi,R=Blast.R))
    backend=Blast.CubicBSplineBackend(bank)
    @test_throws ArgumentError Blast.SetUp(G,L;backend)
    W,P=Blast.SetUp(G,L,C;backend,radial_plan=q)
    grid=Blast.nonlimber_grid(P);pk=fixture_pk(grid);pmm=characterization_pmm()
    ps=Blast.prepare_pk_workspace(P,pk,pmm.pk_limber_linear,pmm.pk_limber_nonlinear,bg;pk_grid=grid)
    @test ps.radial_plan===q && P.radial_plan===q
    projected=Blast.compute_w(W,ps)
    @test count(f->!isnothing(getfield(projected,f)),fieldnames(typeof(projected)))==17
    for (name,_,family,r1) in Blast._PROJECTED_MATTER_COMPONENTS
        component=getfield(projected,name)
        coefs=getfield(ps,Symbol(:c,family)).coefs
        coefs=r1 ? coefs[:,:,end] : coefs
        @test component.w≈literal_w(coefs,component.T̃) rtol=2e-12
    end
    ell=[2.,10.,100.,219.,221.,500.,1500.]
    for pair in ((G,),(G,L),(L,),(C,G),(C,L))
        x=Blast.get_Cℓ(ell,pair...,ps,projected,bg,P)
        y=Blast.get_Cℓ(ell,pair...,ps,projected,bg)
        @test all(isfinite,x)
        @test x≈y rtol=2e-10
    end
    WA,WB=Blast.get_kernel_array(G.δ,bg),Blast.get_kernel_array(L.γ,bg)
    pm=projected.w_0_00_ϕTT.w
    actual=Blast.compute_Cℓ(G.δ,L.γ,projected.w_0_00_ϕTT,bg)
    pref=Blast._pair_prefactor(G.δ,L.γ)
    for l in (1,11,22),i in axes(WA,1),j in axes(WB,1)
        expected=sum(pref[l]*q.chi[c]*q.weights[c]*q.R_weights[r]*pm[l,c,r]*
            (WA[i,c,end]*WB[j,c,r]+WA[i,c,r]*WB[j,c,end])
            for c in eachindex(q.chi),r in eachindex(q.R))
        @test actual[l,i,j]≈expected rtol=3e-12
    end
    alternate=Blast.RadialQuadraturePlan(q.chi,q.weights,1.01q.suffix)
    _,P2=Blast.SetUp(G,L,C;backend,radial_plan=alternate)
    @test_throws ArgumentError Blast.get_Cℓ(ell,G,L,ps,projected,bg,P2)
    bg2=Blast.Background(characterization_cosmology();radial_plan=alternate)
    @test_throws ArgumentError Blast.get_Cℓ(ell,G,L,ps,projected,bg2)
    wrong=copy(q.chi);wrong[4]=nextfloat(wrong[4])
    wrongq=trapezoid_radial_plan(wrong)
    @test_throws ArgumentError Blast.SetUp(G,L;backend,radial_plan=wrongq)
    wrongbg=Blast.Background(characterization_cosmology();radial_plan=wrongq)
    @test_throws ArgumentError Blast.compute_Cℓ(G.δ,L.γ,projected.w_0_00_ϕTT,wrongbg)
    for args in ((L,G),(C,G),(C,L),(G,C,L),(L,G,C),(L,C,G),(C,G,L),(C,L,G))
        _,p=Blast.SetUp(args...;backend,radial_plan=q)
        @test p.radial_plan===q
    end
end

@testset "Explicit R public path: zero-weight endpoint and all 17 blocks" begin
    base=trapezoid_radial_plan([26.,60.,250.,800.,1600.,2800.,4200.,5500.,7000.])
    q=Blast.RadialQuadraturePlan(base.chi,base.weights,base.suffix;R=[0.2,0.7,1.],R_weights=[0.5,0.5,0.])
    bg=Blast.Background(characterization_cosmology();radial_plan=q)
    probes=build_characterization_probes(bg);G,L,C=probes.galaxy,probes.lensing,probes.cmb
    bank=manufactured_direct_bank(;axes=(ell=Blast.ℓ_nonlimber,chi=q.chi,R=q.R))
    backend=Blast.CubicBSplineBackend(bank)
    @test_throws ArgumentError Blast.SetUp(G,L;backend,radial_plan=base)
    W,P=Blast.SetUp(G,L,C;backend,radial_plan=q)
    grid=Blast.nonlimber_grid(P);pk=fixture_pk(grid);pmm=characterization_pmm()
    ps=Blast.prepare_pk_workspace(P,pk,pmm.pk_limber_linear,pmm.pk_limber_nonlinear,bg;pk_grid=grid)
    @test size(ps.cϕTT.coefs)==(9,9,3)
    projected=Blast.compute_w(W,ps)
    for (name,_,family,r1) in Blast._PROJECTED_MATTER_COMPONENTS
        component=getfield(projected,name)
        coefs=getfield(ps,Symbol(:c,family)).coefs
        coefs=r1 ? coefs[:,:,end] : coefs
        @test size(component.w)==(22,9,3)
        @test component.w≈literal_w(coefs,component.T̃) rtol=2e-12
    end
    for probe in (G.δ,G.RSD,G.μ,G.PNG,L.γ,L.IA,C.κ)
        @test size(Blast.get_kernel_array(probe,bg))==(size(probe.Kernel,1),9,3)
    end
    @test Blast.get_kernel_array(G.δ,bg)[:,:,end]≈G.δ.Kernel rtol=1e-12 atol=1e-20
    for pair in ((G,),(G,L),(L,),(C,G),(C,L))
        result=Blast.get_Cℓ([2.,100.,500.],pair...,ps,projected,bg,P)
        @test all(isfinite,result)
    end
    # An explicit sum verifies that the supplied R weights reach the contraction.
    ka,kb=Blast.get_kernel_array(G.δ,bg),Blast.get_kernel_array(L.γ,bg)
    component=projected.w_0_00_ϕTT
    actual=Blast.compute_Cℓ(G.δ,L.γ,component,bg)
    pref=Blast._pair_prefactor(G.δ,L.γ)
    for ell in (1,11,22),i in axes(ka,1),j in axes(kb,1)
        expected=sum(pref[ell]*q.chi[c]*q.weights[c]*q.R_weights[r]*component.w[ell,c,r]*
            (ka[i,c,end]*kb[j,c,r]+ka[i,c,r]*kb[j,c,end]) for c in 1:9,r in 1:3)
        @test actual[ell,i,j]≈expected rtol=3e-12
    end
end

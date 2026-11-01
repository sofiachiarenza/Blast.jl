using Test, Blast, LinearAlgebra, Random
using DifferentiationInterface: prepare_gradient, gradient
using ADTypes: AutoMooncake, AutoForwardDiff
using Mooncake
isdefined(@__MODULE__,:trapezoid_radial_plan) || include("test_radial_quadrature.jl")
isdefined(@__MODULE__,:manufactured_direct_bank) || include("bspline_backend_fixture.jl")

@testset "Prepared gradients of physical signed suffix operators" begin
    q=trapezoid_radial_plan([26.,60.,250.,800.,1600.,2800.,4200.,5500.,7000.])
    s=copy(q.suffix);s[2,1]=-0.25 # algebra stress: partial-panel rows may sample to their left
    n=length(q.chi);H0=collect(LinRange(70.,350.,n));z0=collect(LinRange(.01,3.3,n))
    nz=[exp(-((z-c)/.3)^2) for c in (.5,1.5),z in z0]
    function objective(p)
        H=H0.*exp(p[1]);z=z0.*exp(p[2]);density=nz.*exp(p[3])
        slope=fill(.5+p[4],size(nz));pref=1e-8*exp(p[5])
        shear=Blast._cosmic_shear_kernel_tullio(H,q.chi,z,density,s,1.,pref)
        magnification=Blast._magnification_bias_kernel_tullio(H,q.chi,z,density,slope,s,1.,pref)
        return 100*(sum(shear)+sum(magnification))
    end
    mc,fd=AutoMooncake(;config=nothing),AutoForwardDiff()
    p0=zeros(5)
    pm=prepare_gradient(objective,mc,p0);pf=prepare_gradient(objective,fd,p0)
    for p in (p0,[.03,-.01,.1,-.1,.04],[-.02,.04,-.1,-.2,0.],p0)
        gm=gradient(objective,pm,mc,p);gf=gradient(objective,pf,fd,p)
        @test all(isfinite,gm) && all(isfinite,gf)
        @test gm≈gf rtol=1e-10 atol=1e-12
    end
end

@testset "Prepared cosmology/source gradients through nonuniform public spectra" begin
    q=trapezoid_radial_plan([26.,60.,250.,800.,1600.,2800.,4200.,5500.,7000.])
    c=characterization_cosmology()
    bg0=Blast.Background(c;radial_plan=q);probes=build_characterization_probes(bg0)
    backend=Blast.CubicBSplineBackend(manufactured_direct_bank(;axes=(ell=Blast.ℓ_nonlimber,chi=q.chi,R=Blast.R)))
    W,P=Blast.SetUp(probes.galaxy,probes.lensing;backend,radial_plan=q)
    pk0=fixture_pk(Blast.nonlimber_grid(P));pmm=characterization_pmm()
    ell=[2.,10.,100.,221.,500.]
    function outputs(p)
        cosmo=Blast.w0waCDMCosmology(ln10Aₛ=c.ln10Aₛ,nₛ=p[2],h=p[1],ωb=c.ωb,ωc=c.ωc,
            ωk=c.ωk,mν=c.mν,w0=c.w0,wa=c.wa)
        bg=Blast.Background(cosmo;radial_plan=q)
        pr=build_characterization_probes(bg)
        amp=exp(p[3])
        ps=Blast.prepare_pk_workspace(P,pk0.*amp,pmm.pk_limber_linear.*amp,pmm.pk_limber_nonlinear.*amp,bg)
        wc=Blast.compute_w(W,ps)
        gg=Blast.get_Cℓ(ell,pr.galaxy,ps,wc,bg,P)
        gl=Blast.get_Cℓ(ell,pr.galaxy,pr.lensing,ps,wc,bg,P)
        ll=Blast.get_Cℓ(ell,pr.lensing,ps,wc,bg,P)
        return (sum(abs2,gg),sum(abs2,gl),sum(abs2,ll))
    end
    p0=[c.h,c.nₛ,0.];normalization=outputs(p0)
    objective=p->sum(outputs(p)./normalization)
    mc,fd=AutoMooncake(;config=nothing),AutoForwardDiff()
    println("RADIAL_AD preparing nonuniform full-public Mooncake/ForwardDiff");flush(stdout)
    pm=prepare_gradient(objective,mc,p0);pf=prepare_gradient(objective,fd,p0)
    rng=MersenneTwister(20261022)
    random=[.65+.04rand(rng),.955+.02rand(rng),.02randn(rng)] # fixture-domain, not a Turing-prior claim
    for p in (p0,[.69,.96,.03],random,p0)
        gm=gradient(objective,pm,mc,p);gf=gradient(objective,pf,fd,p)
        @test all(isfinite,gm) && all(isfinite,gf)
        @test gm≈gf rtol=1e-9 atol=1e-12
        println("RADIAL_AD point=$p Mooncake=$gm ForwardDiff=$gf");flush(stdout)
    end
end

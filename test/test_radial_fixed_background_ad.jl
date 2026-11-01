# Isolates the new radial operators from the separately reproduced ACE growth
# ODE Float64/Dual primal drift. Full cosmology AD remains a separate, strict gate.
using Test, Blast, Random
using DifferentiationInterface: prepare_gradient, gradient
using ADTypes: AutoMooncake, AutoForwardDiff
using Mooncake
isdefined(@__MODULE__,:trapezoid_radial_plan) || include("test_radial_quadrature.jl")
isdefined(@__MODULE__,:manufactured_direct_bank) || include("bspline_backend_fixture.jl")

@testset "Nonuniform fixed-background public AD: functional and reusable weights" for custom_R in (false,true)
    base=trapezoid_radial_plan([26.,60.,250.,800.,1600.,2800.,4200.,5500.,7000.])
    q=custom_R ? Blast.RadialQuadraturePlan(base.chi,base.weights,base.suffix;
        R=[0.2,0.7,1.],R_weights=[0.5,0.5,0.]) : base
    bg=Blast.Background(characterization_cosmology();radial_plan=q)
    sources=characterization_source_distributions()
    probes=build_characterization_probes(bg)
    backend=Blast.CubicBSplineBackend(manufactured_direct_bank(;axes=(ell=Blast.ℓ_nonlimber,chi=q.chi,R=q.R)))
    W,P=Blast.SetUp(probes.galaxy,probes.lensing;backend,radial_plan=q)
    pk0=fixture_pk(Blast.nonlimber_grid(P));pmm=characterization_pmm()
    ps0=Blast.prepare_pk_workspace(P,pk0,pmm.pk_limber_linear,pmm.pk_limber_nonlinear,bg)
    cache=Blast.allocate_compute_w(W,ps0)
    ell=[2.,10.,100.,221.,500.]
    function outputs(p,reuse)
        bias=fill(1.5p[2],size(sources.nz_g,1),length(bg.z))
        slope=fill(.3+p[3],size(bias))
        g=Blast.prepare_probe(Blast.GalaxyClustering(
            δ=Blast.NumberCounts(nz=sources.nz_g,z=sources.z,bias=bias),
            RSD=Blast.RedshiftSpaceDistortions(nz=sources.nz_g,z=sources.z,growth_rate=bg.f),
            μ=Blast.MagnificationBias(nz=sources.nz_g,z=sources.z,s=slope),
            PNG=Blast.PrimordialNonGaussianity(nz=sources.nz_g,z=sources.z,bias=bias,f_NL=p[4],p=1.)),bg)
        l=Blast.prepare_probe(Blast.WeakLensing(γ=Blast.CosmicShear(nz=sources.nz_s,z=sources.z),
            IA=Blast.IntrinsicAlignment(nz=sources.nz_s,z=sources.z,A=p[5])),bg)
        amp=exp(p[1])
        ps=Blast.prepare_pk_workspace(P,pk0.*amp,pmm.pk_limber_linear.*amp,pmm.pk_limber_nonlinear.*amp,bg)
        wc=if reuse
            Blast.compute_w!(cache,ps);cache
        else
            Blast.compute_w(W,ps)
        end
        gg=Blast.get_Cℓ(ell,g,ps,wc,bg,P);gl=Blast.get_Cℓ(ell,g,l,ps,wc,bg,P)
        ll=Blast.get_Cℓ(ell,l,ps,wc,bg,P)
        return (sum(abs2,gg),sum(abs2,gl),sum(abs2,ll))
    end
    p0=[0.,1.,0.,7.,1.72];normalization=outputs(p0,false)
    functional=p->sum(outputs(p,false)./normalization)
    reusable=p->sum(outputs(p,true)./normalization)
    mc,fd=AutoMooncake(;config=nothing),AutoForwardDiff()
    println("RADIAL_FIXED_AD R_count=$(length(q.R)) preparing functional/reusable Mooncake and ForwardDiff");flush(stdout)
    pf=prepare_gradient(functional,mc,p0);pr=prepare_gradient(reusable,mc,p0)
    pd=prepare_gradient(functional,fd,p0)
    for p in (p0,[.02,.95,.02,-20.,1.5],[-.01,1.05,-.03,0.,1.9],[.01,1.,0.,20.,1.72],p0)
        gf=gradient(functional,pf,mc,p);gr=gradient(reusable,pr,mc,p);gd=gradient(functional,pd,fd,p)
        @test all(isfinite,gf) && all(isfinite,gr) && all(isfinite,gd)
        @test gf≈gd rtol=1e-9 atol=1e-12
        @test gr≈gd rtol=1e-9 atol=1e-12
        @test reusable(p)≈functional(p) rtol=1e-12
        println("RADIAL_FIXED_AD point=$p functional=$gf reusable=$gr ForwardDiff=$gd");flush(stdout)
    end
end

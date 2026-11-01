using Test, LinearAlgebra
using Blast
isdefined(@__MODULE__,:characterization_cosmology) || include("characterization_workload.jl")

function trapezoid_radial_plan(x)
    n=length(x);s=zeros(n,n)
    for i in 1:n-1,j in i:n-1
        h=(x[j+1]-x[j])/2
        s[i,j]+=h;s[i,j+1]+=h
    end
    return Blast.RadialQuadraturePlan(x,s[1,:],s)
end

@testset "Explicit multipole geometry without changing radial measure" begin
    base=Blast.LEGACY_RADIAL_PLAN
    ell=[2.0,2.95,80.01,219.0]
    q=Blast.RadialQuadraturePlan(base.chi,base.weights,base.suffix;legacy=true,ell)
    @test q.signature==base.signature
    @test q.geometry==Blast._geometry_signature(ell,base.chi,base.R)
    @test q.geometry!=base.geometry
    @test q.chiR_flat==base.chiR_flat && q.R_weights==base.R_weights
    ell[2]=3.01
    @test q.geometry!=Blast._geometry_signature(ell,base.chi,base.R)
    for bad in (Float64[],[0.0],[NaN],[Inf],[-2.0])
        @test_throws ArgumentError Blast.RadialQuadraturePlan(base.chi,base.weights,base.suffix;ell=bad)
    end
end

@testset "Fixed radial plan ownership and physical suffix kernels" begin
    original=[26.,60.,250.,800.,1600.,2800.,4200.,5500.,7000.]
    plan=trapezoid_radial_plan(original)
    @test plan.chi==original && plan.chi!==original
    original[2]=61.
    @test plan.chi[2]==60.
    @test plan.suffix[1,:]==plan.weights
    @test all(iszero,plan.suffix[end,:])
    @test plan.chiR_flat==vec(plan.chi*transpose(Blast.R))
    @test plan.signature==trapezoid_radial_plan(plan.chi).signature
    @test plan.signature!=trapezoid_radial_plan(original).signature
    @test_throws ArgumentError Blast.RadialQuadraturePlan([1.,1.],[1.,1.],ones(2,2))
    @test_throws DimensionMismatch Blast.RadialQuadraturePlan([1.,2.],[1.],ones(2,2))
    @test_throws ArgumentError Blast.RadialQuadraturePlan([1.,2.],[1.,NaN],ones(2,2))
    @test_throws ArgumentError Blast.RadialQuadraturePlan(plan.chi,plan.weights,plan.suffix;legacy=true)
    bg=Blast.Background(characterization_cosmology();radial_plan=plan)
    @test bg.χ==plan.chi && bg.radial_plan===plan
    @test_throws ArgumentError Blast.Background(characterization_cosmology();radial_plan=plan,χ_grid=original)
    integ=@inferred Blast._prepare_nonlimber_integration(bg)
    @test integ.w_χ===plan.weights && integ.Δχ==1.0
    raw=Blast.CosmicShear(nz=ones(1,101),z=collect(LinRange(0.,5.,101)))
    shear=Blast.prepare_probe(Blast.WeakLensing(γ=raw),bg).γ
    pref=1.5*Blast.H_0_CONV^2*Blast.ω_m0(bg)/Blast.C_LIGHT^2
    expected=[pref*bg.χ[i]*(1+bg.z[i])*sum(plan.suffix[i,j]*bg.H[j]/Blast.C_LIGHT*
        shear.nz_norm[1,j]*(bg.χ[j]-bg.χ[i])/bg.χ[j] for j in eachindex(bg.χ)) for i in eachindex(bg.χ)]
    @test vec(shear.Kernel)≈expected rtol=3e-13 atol=1e-20
    bare=Blast.Background(characterization_cosmology();χ_grid=plan.chi)
    @test_throws ArgumentError Blast.prepare_probe(Blast.WeakLensing(γ=raw),bare)
    default=Blast.Background(characterization_cosmology())
    @test default.radial_plan===Blast.LEGACY_RADIAL_PLAN
    s,h=Blast._kernel_quadrature(default.χ,default.radial_plan)
    @test s==Blast.simpson_weights_matrix(128)
    @test h==(7000-26)/127
    lifted=Blast.Background{Float64,typeof(default.z_of_χ),typeof(default.cosmo)}(
        default.cosmo,default.z,default.χ,default.H,default.D,default.D_norm,default.f,default.z_of_χ)
    @test lifted.radial_plan===Blast.LEGACY_RADIAL_PLAN
end

@testset "Explicit R quadrature ownership, endpoints and validation" begin
    base=Blast.LEGACY_RADIAL_PLAN
    r=[0.2,0.7,1.0];wr=[0.2,0.3,0.0]
    q=Blast.RadialQuadraturePlan(base.chi,base.weights,base.suffix;legacy=true,R=r,R_weights=wr)
    @test q.R==r && q.R!==r && q.R_weights==wr && q.R_weights!==wr
    @test sum(q.R_weights)==0.5 # Preserve supplied measure; never renormalize.
    @test q.chiR_flat==vec(q.chi*transpose(r))
    @test q.inv_chiR2_row==reshape(inv.(q.chiR_flat.^2),1,:)
    @test q.geometry!=base.geometry && q.signature!=base.signature
    same=Blast.RadialQuadraturePlan(base.chi,base.weights,base.suffix;legacy=true,R=copy(Blast.R))
    @test same.signature==base.signature && same.geometry==base.geometry
    altered=Blast.RadialQuadraturePlan(base.chi,base.weights,base.suffix;R=r,R_weights=2wr)
    @test altered.geometry==q.geometry && altered.signature!=q.signature
    r[1]=0.25;wr[1]=0.1
    @test q.R[1]==0.2 && q.R_weights[1]==0.2
    @test_throws ArgumentError Blast.RadialQuadraturePlan(base.chi,base.weights,base.suffix;R=r)
    @test_throws DimensionMismatch Blast.RadialQuadraturePlan(base.chi,base.weights,base.suffix;R=r,R_weights=[1.0])
    for bad in ([0.,0.7,1.], [0.2,0.7,0.9], [0.7,0.2,1.], [0.2,0.2,1.], [NaN,0.7,1.], [1.])
        @test_throws ArgumentError Blast.RadialQuadraturePlan(base.chi,base.weights,base.suffix;R=bad,R_weights=ones(length(bad)))
    end
    for bad in ([-0.1,0.7,0.], [0.1,Inf,0.], zeros(3))
        @test_throws ArgumentError Blast.RadialQuadraturePlan(base.chi,base.weights,base.suffix;R=r,R_weights=bad)
    end
end

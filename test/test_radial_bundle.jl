using Test, Blast
isdefined(@__MODULE__,:trapezoid_radial_plan) || include("test_radial_quadrature.jl")
isdefined(@__MODULE__,:manufactured_direct_bank) || include("bspline_backend_fixture.jl")
isdefined(@__MODULE__,:write_synthetic_bundle) || include("tensor_bundle_fixture.jl")

@testset "Canonical direct/ratio bundles retain nonuniform radial axes" begin
    q=trapezoid_radial_plan([26.,60.,250.,800.,1600.,2800.,4200.,5500.,7000.])
    bg=Blast.Background(characterization_cosmology();radial_plan=q)
    probes=build_characterization_probes(bg);pmm=characterization_pmm()
    gridaxes=(ell=copy(Blast.ℓ_nonlimber),chi=copy(q.chi),R=copy(Blast.R))
    mktempdir() do root
        write_synthetic_bundle(root;counts=(9,),kmin=5e-5,kmax=16.,grid_axes=gridaxes)
        @test_throws ArgumentError Blast.load_tensor_bank(joinpath(root,"bundle.toml");count=9,representation=:direct)
        for rep in (:direct,:ratio)
            bank=syn_load(root;rep,axes=gridaxes,expected_k_range=(5e-5,16.))
            axes=Blast.bank_axes(bank)
            @test axes.chi==q.chi && axes.ell==Blast.ℓ_nonlimber && axes.R==Blast.R
            @test axes.chi!==gridaxes.chi
            @test Blast.bank_signature(bank).geometry==q.geometry
            backend=Blast.CubicBSplineBackend(bank)
            G,L=probes.galaxy,probes.lensing
            W,P=Blast.SetUp(G,L;backend,radial_plan=q)
            pk=fixture_pk(Blast.nonlimber_grid(P))
            ps=Blast.prepare_pk_workspace(P,pk,pmm.pk_limber_linear,pmm.pk_limber_nonlinear,bg)
            @test size(ps.cϕTT.coefs)==(9,9,64)
            wc=Blast.compute_w(W,ps)
            @test all(isfinite,Blast.get_Cℓ([2.,100.,500.],G,L,ps,wc,bg,P))
            @test_throws ArgumentError Blast.SetUp(G,L;backend)
        end
    end
end

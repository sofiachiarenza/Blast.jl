using Test
using Blast

function _type_stability_workload()
    bg = get_test_bg()
    nz = ones(1, 50)
    z_nz = LinRange(0.0, 3.6, 50)
    n_bg = length(bg.z)

    galaxy = Blast.GalaxyClustering(
        δ=Blast.NumberCounts(nz=nz, z=z_nz, bias=ones(1, n_bg)),
        RSD=Blast.RedshiftSpaceDistortions(nz=nz, z=z_nz),
        μ=Blast.MagnificationBias(nz=nz, z=z_nz, s=ones(1, n_bg)),
        PNG=Blast.PrimordialNonGaussianity(
            nz=nz, z=z_nz, bias=ones(1, n_bg), f_NL=1.0, p=0.0,
        ),
    )
    lensing = Blast.WeakLensing(
        γ=Blast.CosmicShear(nz=nz, z=z_nz),
        IA=Blast.IntrinsicAlignment(nz=nz, z=z_nz, A_IA=ones(1, n_bg)),
    )
    W, plans = Blast.SetUp(galaxy, lensing)
    Blast.evaluate_components!(galaxy, bg)
    Blast.evaluate_components!(lensing, bg)

    pk = ones(length(Blast.z_lin), length(Blast.k_cheb))
    pk_limber = ones(length(Blast.z_cheb), length(Blast.k_limber))
    spectrum = Blast.prepare_pk_workspace(plans, pk, pk_limber, pk_limber, bg)
    W = Blast.compute_w(W, spectrum)
    return (; galaxy, lensing, spectrum, W, bg, plans)
end

function _plain_3x2pt(ℓ, workload)
    cl_gg = Blast.get_Cℓ(
        ℓ, workload.galaxy, workload.spectrum, workload.W,
        workload.bg, workload.plans,
    )
    cl_gs = Blast.get_Cℓ(
        ℓ, workload.galaxy, workload.lensing, workload.spectrum,
        workload.W, workload.bg, workload.plans,
    )
    cl_ss = Blast.get_Cℓ(
        ℓ, workload.lensing, workload.spectrum, workload.W,
        workload.bg, workload.plans,
    )
    return (; cl_gg, cl_gs, cl_ss)
end

@testset "plain Julia full 3x2pt inference" begin
    workload = _type_stability_workload()
    ℓ = Blast.full_ℓ_range

    cl_gg = @inferred Blast.get_Cℓ(
        ℓ, workload.galaxy, workload.spectrum, workload.W,
        workload.bg, workload.plans,
    )
    cl_gs = @inferred Blast.get_Cℓ(
        ℓ, workload.galaxy, workload.lensing, workload.spectrum,
        workload.W, workload.bg, workload.plans,
    )
    cl_ss = @inferred Blast.get_Cℓ(
        ℓ, workload.lensing, workload.spectrum, workload.W,
        workload.bg, workload.plans,
    )
    result = @inferred _plain_3x2pt(ℓ, workload)

    @test cl_gg isa Array{Float64, 3}
    @test cl_gs isa Array{Float64, 3}
    @test cl_ss isa Array{Float64, 3}
    @test result isa NamedTuple{
        (:cl_gg, :cl_gs, :cl_ss),
        Tuple{Array{Float64, 3}, Array{Float64, 3}, Array{Float64, 3}},
    }
end

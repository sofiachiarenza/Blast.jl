using Test

isdefined(@__MODULE__, :build_characterization_state) ||
    include("characterization_workload.jl")

@testset "source distributions" begin
    @test_throws DimensionMismatch NumberCounts(
        nz=ones(1, 4), z=[0.0, 0.5, 1.0], bias=ones(1, length(Blast.χ)))
    @test_throws ArgumentError NumberCounts(
        nz=ones(1, 3), z=[0.0, 1.0, 0.5], bias=ones(1, length(Blast.χ)))
    @test_throws ArgumentError NumberCounts(
        nz=zeros(1, 3), z=[0.0, 0.5, 1.0], bias=ones(1, length(Blast.χ)))
    @test_throws ArgumentError CosmicShear(
        nz=[1.0 NaN 1.0], z=[0.0, 0.5, 1.0])
end

@testset "component coordinate shapes" begin
    bg = Background(characterization_cosmology())
    sources = characterization_source_distributions()
    wrong_bias = ones(2, length(bg.z) - 1)
    raw_gc = GalaxyClustering(δ=NumberCounts(
        nz=sources.nz_g, z=sources.z, bias=wrong_bias))
    @test_throws DimensionMismatch prepare_probe(raw_gc, bg)

    wrong_A_IA = ones(2, length(bg.z) - 1)
    raw_wl = WeakLensing(
        γ=CosmicShear(nz=sources.nz_s, z=sources.z),
        IA=IntrinsicAlignment(nz=sources.nz_s, z=sources.z, A_IA=wrong_A_IA),
    )
    @test_throws DimensionMismatch prepare_probe(raw_wl, bg)
end

@testset "Pmm and production-grid boundaries" begin
    state = build_characterization_state()
    pmm = state.pmm

    @test_throws DimensionMismatch prepare_pk_workspace(
        state.plans, pmm.pk_nonlimber[1:end-1, :],
        pmm.pk_limber_linear, pmm.pk_limber_nonlinear, state.bg)

    bad_pmm = copy(pmm.pk_nonlimber)
    bad_pmm[1] = NaN
    @test_throws ArgumentError prepare_pk_workspace(
        state.plans, bad_pmm,
        pmm.pk_limber_linear, pmm.pk_limber_nonlinear, state.bg)

    zero_pmm = copy(pmm.pk_nonlimber)
    zero_pmm[1] = 0.0
    @test_throws DomainError prepare_pk_workspace(
        state.plans, zero_pmm,
        pmm.pk_limber_linear, pmm.pk_limber_nonlinear, state.bg)

    custom_χ = collect(LinRange(first(Blast.χ), last(Blast.χ), length(Blast.χ) - 2))
    custom_bg = Background(state.cosmo; χ_grid=custom_χ)
    @test_throws DimensionMismatch prepare_pk_workspace(
        state.plans, pmm.pk_nonlimber,
        pmm.pk_limber_linear, pmm.pk_limber_nonlinear, custom_bg)
end

@testset "multipole domain" begin
    state = build_characterization_state()
    spectrum = prepare_characterization_spectrum(state)
    weights = compute_w(state.weights_template, spectrum)
    @test_throws DomainError get_Cℓ(
        [0.0], state.probes.galaxy, spectrum, weights, state.bg, state.plans)
    @test_throws DomainError get_Cℓ(
        [2001.0], state.probes.galaxy, spectrum, weights, state.bg, state.plans)
    @test_throws ArgumentError get_Cℓ(
        [NaN], state.probes.galaxy, spectrum, weights, state.bg, state.plans)
end

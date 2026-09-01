using Test
using Blast

@testset "Cosmology: Background checks" begin
    cosmo = get_test_cosmo()
    bg = get_test_bg(cosmo)
    
    # Adimensional Hubble factor at z=0 is 1.0
    E0_test = Blast.E_z(0.0, cosmo)
    @test E0_test ≈ 1.0
    
    # H0 at z=0 matches 100 * h
    H0_test = Blast.compute_hubble_factor(0.0, cosmo)
    @test H0_test ≈ 100.0 * cosmo.h

    # Check that H_array, χ, and z in Background object are consistent.
    # bg.z contains the redshifts found for the Blast.χ grid via the akima
    # inversion built from an N_BG_FINE_GRID-point (fine_z, fine_χ) table.
    # The round-trip `r_z(z_nodes[i], cosmo) ≈ bg.χ[i]` is bounded by that
    # inversion's accuracy, worst near χ=26 Mpc (lowest knot density in the
    # high-dχ/dz region). At N_BG_FINE_GRID=100 the floor is ~5e-4 relative /
    # ~0.02 Mpc absolute — cosmologically negligible (redshift error ~3e-6).
    # See benchmark/bg_fine_grid_sweep.jl for the full convergence table.
    for i in 1:length(bg.z)
        @test bg.H[i] ≈ Blast.compute_hubble_factor(bg.z[i], cosmo)
        @test bg.χ[i] ≈ Blast.compute_χ(bg.z[i], cosmo) rtol=5e-4
    end
end

@testset "Cosmology: Parameter accessors" begin
    cosmo = get_test_cosmo()
    Ωγ0 = 2.469e-5 / cosmo.h^2
    expected_ων = cosmo.h^2 * Blast.cosmo_ext._ΩνE2(1.0, Ωγ0, cosmo.mν)
    @test Blast.ω_ν0(cosmo) ≈ expected_ων
    @test Blast.ω_m0(cosmo) ≈ cosmo.ωb + cosmo.ωc + expected_ων
    @test Blast.get_ns(cosmo) == cosmo.nₛ
    @test Blast.get_As(cosmo) == exp(cosmo.ln10Aₛ) / 1e10
end

@testset "Cosmology: neutrino physical-density adapter" begin
    for mν in (0.0, 0.06, 0.3, [0.01, 0.02, 0.03])
        cosmo = get_test_cosmo(mν=mν)
        Ωγ0 = 2.469e-5 / cosmo.h^2
        expected = cosmo.h^2 * Blast.cosmo_ext._ΩνE2(1.0, Ωγ0, mν)
        @test Blast.ω_ν0(cosmo) ≈ expected rtol=1e-13
        @test Blast.ω_m0(cosmo) ≈ cosmo.ωb + cosmo.ωc + expected rtol=1e-13
    end
end

@testset "Cosmology: E_z physical reasonableness" begin
    cosmo = get_test_cosmo()
    # E(0) = 1 by definition
    @test Blast.E_z(0.0, cosmo) ≈ 1.0
    # For a flat ΛCDM with Ωm < 1, E(z) > 1 for z > 0
    for z in [0.5, 1.0, 2.0]
        @test Blast.E_z(z, cosmo) > 1.0
    end
    # E_z should be monotonically increasing for these redshifts
    zs = [0.0, 0.5, 1.0, 2.0, 3.0]
    Ez = Blast.E_z.(zs, Ref(cosmo))
    @test issorted(Ez)
end

@testset "Cosmology: Grid inversion z(χ) vs χ(z)" begin
    cosmo = get_test_cosmo()
    bg = get_test_bg(cosmo)
    
    # Pick a random point in the grid
    target_z = 1.2
    target_χ = Blast.compute_χ(target_z, cosmo)
    
    # Test our pre-built interpolators
    @test bg.z_of_χ(target_χ) ≈ target_z rtol=1e-4
    @test Blast.compute_χ(target_z, cosmo) ≈ target_χ rtol=1e-6 # Analytical truth
end

@testset "Cosmology: exact growth normalization" begin
    cosmo = get_test_cosmo()
    bg = get_test_bg(cosmo)
    D0 = Blast.D_z(0.0, cosmo)
    @test bg.D_norm ≈ bg.D ./ D0 rtol=5e-13
    @test bg.D_norm[1] < 1
    @test bg.D[1] != bg.D_norm[1]
end

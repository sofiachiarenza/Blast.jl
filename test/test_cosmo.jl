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

    # Check that H_array, χ, and z in Background object are consistent
    # Note: bg.z contains the redshifts found for the Blast.χ grid
    for i in 1:length(bg.z)
        @test bg.H[i] ≈ Blast.compute_hubble_factor(bg.z[i], cosmo)
        @test bg.χ[i] ≈ Blast.compute_χ(bg.z[i], cosmo) rtol=1e-5
    end
end

@testset "Cosmology: Parameter accessors" begin
    cosmo = get_test_cosmo()
    @test Blast.get_H0(cosmo) ≈ 100.0 * cosmo.h
    @test Blast.get_Ωm(cosmo) ≈ (cosmo.ωb + cosmo.ωc) / cosmo.h^2
    @test Blast.get_ns(cosmo) == cosmo.nₛ
    @test Blast.get_As(cosmo) == exp(cosmo.ln10Aₛ) / 1e10
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

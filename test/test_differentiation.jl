using Test
using Blast
using DifferentiationInterface
using ADTypes
using ForwardDiff
using Zygote
using Mooncake
using LinearAlgebra

@testset "Differentiation: Akima Spline AD" begin
    # Test data representative of Blast usage
    x1 = collect(LinRange(0.0, 5.0, 20))
    y = sin.(x1)
    x2 = collect(LinRange(0.5, 4.5, 30))

    @testset "Gradient w.r.t. y (data values)" begin
        f_interp(y_vals) = sum(akima_interpolation(y_vals, x1, x2))
        grad_fwd = DifferentiationInterface.gradient(f_interp, AutoForwardDiff(), y)
        grad_zyg = DifferentiationInterface.gradient(f_interp, AutoZygote(), y)
        grad_mc = DifferentiationInterface.gradient(f_interp, AutoMooncake(), y)

        @test grad_fwd ≈ grad_zyg rtol=1e-9
        @test grad_mc ≈ grad_zyg rtol=1e-9
    end

    @testset "Gradient w.r.t. x1 (input grid)" begin
        f_x1(x1_vals) = sum(akima_interpolation(y, x1_vals, x2))
        grad_fwd = DifferentiationInterface.gradient(f_x1, AutoForwardDiff(), x1)
        grad_zyg = DifferentiationInterface.gradient(f_x1, AutoZygote(), x1)
        grad_mc = DifferentiationInterface.gradient(f_x1, AutoMooncake(), x1)

        @test grad_fwd ≈ grad_zyg rtol=1e-9
        @test grad_mc ≈ grad_zyg rtol=1e-9
    end

    @testset "Gradient w.r.t. x2 (query points)" begin
        f_x2(x2_vals) = sum(akima_interpolation(y, x1, x2_vals))
        grad_fwd = DifferentiationInterface.gradient(f_x2, AutoForwardDiff(), x2)
        grad_zyg = DifferentiationInterface.gradient(f_x2, AutoZygote(), x2)
        grad_mc = DifferentiationInterface.gradient(f_x2, AutoMooncake(), x2)

        @test grad_fwd ≈ grad_zyg rtol=1e-9
        @test grad_mc ≈ grad_zyg rtol=1e-9
    end

    @testset "Matrix Akima AD" begin
        Y_mat = hcat(sin.(x1), cos.(x1), sin.(2 .* x1))
        
        f_mat(Y) = sum(akima_interpolation(Y, x1, x2))

        grad_fwd = DifferentiationInterface.gradient(f_mat, AutoForwardDiff(), Y_mat)
        grad_zyg = DifferentiationInterface.gradient(f_mat, AutoZygote(), Y_mat)
        grad_mc = DifferentiationInterface.gradient(f_mat, AutoMooncake(), Y_mat)

        @test grad_fwd ≈ grad_zyg rtol=1e-9
        @test grad_mc ≈ grad_zyg rtol=1e-9
    end
end

@testset "Differentiation: Tullio Contractions" begin
    # Dimensions matching non-Limber integration
    # Use global nχ and nR to avoid DimensionMismatch in w_ell_tullio
    nℓ_nl = length(Blast.ℓ_nonlimber)
    nχ = length(Blast.χ)
    nR = length(Blast.R)
    
    # Matching the specific 3D rule: @tullio w[i,j,k] := c[l,j,k] * T[i,j,k,l]
    # where l=nR, j=nχ, k=nR, i=nℓ_nl
    c = rand(nR, nχ, nR) 
    T = rand(nℓ_nl, nχ, nR, nR)
    
    f_tullio(coeffs) = sum(Blast.w_ell_tullio(coeffs, T))
    
    # Gradient verification
    grad_zyg = DifferentiationInterface.gradient(f_tullio, AutoZygote(), c)
    grad_mc = DifferentiationInterface.gradient(f_tullio, AutoMooncake(), c)
    
    @test grad_zyg ≈ grad_mc rtol=1e-8
end


using Test
using Blast
using DifferentiationInterface
using ADTypes
using ForwardDiff
using Zygote
using Mooncake
using FiniteDifferences
using LinearAlgebra

@testset "Differentiation: Akima Spline AD" begin
    # Test data for Akima AD
    x1 = collect(LinRange(0.0, 5.0, 20))
    y = sin.(x1)
    x2 = collect(LinRange(0.5, 4.5, 30))

    # Wrapper for 1D interpolation
    f_interp(y_vals) = sum(Blast._akima_interpolation(y_vals, x1, x2))

    @testset "Gradient w.r.t. y (data values)" begin
        backend_fd = AutoFiniteDifferences(central_fdm(5, 1))
        
        grad_fd = DifferentiationInterface.gradient(f_interp, backend_fd, y)
        grad_fwd = DifferentiationInterface.gradient(f_interp, AutoForwardDiff(), y)
        grad_zyg = DifferentiationInterface.gradient(f_interp, AutoZygote(), y)
        grad_mc = DifferentiationInterface.gradient(f_interp, AutoMooncake(), y)

        @test grad_fwd ≈ grad_fd rtol=1e-8
        @test grad_zyg ≈ grad_fd rtol=1e-8
        @test grad_mc ≈ grad_fd rtol=1e-8
    end

    @testset "Gradient w.r.t. x1 (input grid)" begin
        f_x1(x1_vals) = sum(Blast._akima_interpolation(y, x1_vals, x2))
        backend_fd = AutoFiniteDifferences(central_fdm(5, 1))
        
        grad_fd = DifferentiationInterface.gradient(f_x1, backend_fd, x1)
        grad_fwd = DifferentiationInterface.gradient(f_x1, AutoForwardDiff(), x1)
        grad_zyg = DifferentiationInterface.gradient(f_x1, AutoZygote(), x1)
        grad_mc = DifferentiationInterface.gradient(f_x1, AutoMooncake(), x1)

        @test grad_fwd ≈ grad_fd rtol=1e-8
        @test grad_zyg ≈ grad_fd rtol=1e-8
        @test grad_mc ≈ grad_fd rtol=1e-8
    end

    @testset "Gradient w.r.t. x2 (query points)" begin
        f_x2(x2_vals) = sum(Blast._akima_interpolation(y, x1, x2_vals))
        backend_fd = AutoFiniteDifferences(central_fdm(5, 1))
        
        grad_fd = DifferentiationInterface.gradient(f_x2, backend_fd, x2)
        grad_fwd = DifferentiationInterface.gradient(f_x2, AutoForwardDiff(), x2)
        grad_zyg = DifferentiationInterface.gradient(f_x2, AutoZygote(), x2)
        grad_mc = DifferentiationInterface.gradient(f_x2, AutoMooncake(), x2)

        @test grad_fwd ≈ grad_fd rtol=1e-8
        @test grad_zyg ≈ grad_fd rtol=1e-8
        @test grad_mc ≈ grad_fd rtol=1e-8
    end

    @testset "Matrix Akima AD" begin
        n_cols = 3
        Y_mat = hcat([sin.(x1 .+ i) for i in 1:n_cols]...)
        
        f_mat(Y) = sum(Blast._akima_interpolation(Y, x1, x2))
        
        backend_fd = AutoFiniteDifferences(central_fdm(5, 1))
        grad_fd = DifferentiationInterface.gradient(f_mat, backend_fd, Y_mat)
        grad_fwd = DifferentiationInterface.gradient(f_mat, AutoForwardDiff(), Y_mat)
        grad_zyg = DifferentiationInterface.gradient(f_mat, AutoZygote(), Y_mat)
        grad_mc = DifferentiationInterface.gradient(f_mat, AutoMooncake(), Y_mat)

        @test grad_fwd ≈ grad_fd rtol=1e-7
        @test grad_zyg ≈ grad_fd rtol=1e-7
        @test grad_mc ≈ grad_fd rtol=1e-7
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

@testset "Differentiation: C_ℓ Pipeline AD" begin
    n_chi = length(Blast.χ)
    nR = length(Blast.R)
    n_bins = 1
    
    cosmo = get_test_cosmo()
    bg = get_test_bg(cosmo)
    
    # Wrapper that differentiates w.r.t NumberCounts kernel
    function loss_kernels(k_vals)
        k_mat = reshape(k_vals, n_bins, n_chi)
        nc = Blast.NumberCounts(Kernel = k_mat, bias = ones(n_bins, n_chi))
        
        # Mock non-limber coefficients
        nℓ_nl = 22 # Matching T_tilde artifact dimension
        PS = Blast.PowerSpectrum(
            cϕTT = Blast.cϕTT(coefs = ones(161, n_chi, nR)),
            cϕT = nothing, cϕ = nothing,
            ΔP_limber = zeros(length(Blast.full_ℓ_range), n_chi),
            Pδ_limber = zeros(length(Blast.full_ℓ_range), n_chi)
        )
        
        W_comp = Blast.w_2_00_ϕTT(w = ones(nℓ_nl, n_chi, nR))
        
        # Differentiable C_ℓ call
        cls = Blast.compute_Cℓ(nc, nc, W_comp, bg)
        return sum(cls)
    end

    k_in_vec = rand(n_bins * n_chi)
    
    # Compare AD backends on the real grid
    grad_zyg = DifferentiationInterface.gradient(loss_kernels, AutoZygote(), k_in_vec)
    grad_mc = DifferentiationInterface.gradient(loss_kernels, AutoMooncake(), k_in_vec)
    
    @test grad_zyg ≈ grad_mc rtol=1e-7
end

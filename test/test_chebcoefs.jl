using Test
using Blast
using FastChebInterp
using StaticArrays
using NPZ
using LinearAlgebra

@testset "Chebcoefs: FFT vs FastChebInterp" begin
    # Power of two nodes usually avoids plan mismatches in FFTW wrappers
    n_nodes_requested = 128
    plan = Blast.prepare_chebyshev_plan(0.0, 1.0, n_nodes_requested)
    
    # Introspect the plan to see what size it actually expects
    # AbstractCosmologicalEmulators uses 'fft_plan'
    n_nodes = size(plan.fft_plan, 1)
    
    n_cols = 5
    A = rand(Float64, n_nodes, n_cols)
    
    my_coefs = zeros(n_nodes, n_cols)
    for i in 1:n_cols
        vec_in = A[:, i]
        my_coefs[:, i] = Blast.chebyshev_decomposition(plan, vec_in)
    end

    true_coefs = zeros(n_nodes, n_cols)
    for i in 1:n_cols
        true_coefs[:,i] = FastChebInterp.chebcoefs(A[:,i])
    end
    
    @test true_coefs ≈ my_coefs
end

@testset "T tilde CC consistency" begin
    # Basic sanity check for the T_tilde integration logic
    ℓ = 100.0
    chi = LinRange(26, 7000, 10)
    R_grid = Blast.chebpoints(20, -1, 1)
    R_grid = reverse(R_grid[R_grid.>0])
    kmax = 200/13 
    kmin = 2.5/7000

    # Test the compute_T̃ logic produces finite, structured output
    res = Blast.compute_T̃(ℓ, chi, R_grid, kmin, kmax, 2, 0, 0, n_cheb=119, N=2^(10)+1)
    
    @test size(res) == (1, 10, length(R_grid), 120)
    @test all(isfinite, res)
    @test any(!iszero, res)
end

@testset "Projected-matter contraction" begin
    # Matching Tullio indices in src/projected_matter.jl:
    # @tullio w[i,j,k] := c[l,j,k] * T[i,j,k,l]
    i, j, k, l = 3, 7, 10, 8
    c = rand(l, j, k)
    T_mat = rand(i, j, k, l)

    w_blast = Blast.w_ell_tullio(c, T_mat)
    @test size(w_blast) == (i, j, k)
    @test !all(iszero, w_blast)
end

@testset "T tilde legacy artifact regression" begin
    chi = LinRange(26, 7000, 10)
    R_grid = Blast.chebpoints(20, -1, 1)
    R_grid = reverse(R_grid[R_grid .> 0])
    kmax = 200 / 13
    kmin = 2.5 / 7000

    # Legacy sectors from runtests_old.jl
    # sector tag, beta exponent, relative-norm tolerance
    configs = [
        ("CC", 2, 5e-4),
        ("CL", 0, 5e-4),
        ("LL", -2, 5e-2)
    ]

    for (sector, beta, tol) in configs
        ref = data["T_tilde"][sector]["2.0"]
        calc = Blast.compute_T̃(2.0, chi, R_grid, kmin, kmax, beta, 0, 0, n_cheb=119, N=2^12+1)

        rel_norm = norm(vec(calc .- ref)) / norm(vec(ref))
        @test rel_norm < tol
    end
end

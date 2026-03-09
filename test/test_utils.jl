using Test
using Blast
using FastTransforms
using Tullio

@testset "Utils: Simpson and CC weights" begin
    # Simpson integration check
    n = 101
    x = LinRange(0, 1, n)
    Δx = (last(x) - first(x)) / (n - 1)
    weights = Blast.simpson_weights_array(n)
    
    # \int_0^1 x dx = 0.5
    @tullio integral := x[i] * weights[i] * Δx
    @test integral ≈ 0.5
    
    # CC weights check
    min_v, max_v = -1.0, 1.0
    N = 100
    x_blast = Blast.get_clencurt_grid(min_v, max_v, N)
    x_check = FastTransforms.clenshawcurtisnodes(Float64, N)
    @test isapprox(x_check, x_blast)

    w_blast = Blast.get_clencurt_weights(min_v, max_v, N)
    object = FastTransforms.chebyshevmoments1(Float64, N)
    w_check = FastTransforms.clenshawcurtisweights(object)
    @test isapprox(w_check, w_blast)
end

@testset "Utils: Akima and Bessel" begin
    # Precomputation tests for Bessel/Cheb evaluation
    min_v = 10^(-1) * (1 + 1e-10)
    max_v = 10.0
    N = 100
    x_check = Blast.get_clencurt_grid(min_v, max_v, Int(N))

    ncheb = 2
    ell = 1
    test_chi = zeros(10)
    # Correct signature: (ell, min, max, chi, ncheb, N, deriv_order)
    T_blast, Bessel_blast = Blast.bessel_cheb_eval(ell, min_v, max_v, test_chi, Int(ncheb), Int(N), 0)

    T_check = zeros(ncheb + 1, N)
    T_check[1, :] = ones(N)
    T_check[2, :] = log10.(x_check)
    T_check[3, :] = 2 .* log10.(x_check) .^ 2 .- 1

    Bessel_check = zeros(10, N) 

    @test isapprox(Bessel_check, Bessel_blast)
    @test isapprox(T_check[1, :], T_blast[1, :])
    @test isapprox(T_check[2, :], T_blast[2, :])
    @test isapprox(T_check[3, :], T_blast[3, :])
end

@testset "Utils: Matrix products" begin
    # Matching Tullio indices in src/projected_matter.jl:
    # @tullio w[i,j,k] := c[l,j,k] * T[i,j,k,l]
    i, j, k, l = 3, 7, 10, 8
    c = rand(l, j, k)
    T_mat = rand(i, j, k, l)
    w_blast = Blast.w_ell_tullio(c, T_mat)
    @test size(w_blast) == (i, j, k)
    @test !all(iszero, w_blast)
end

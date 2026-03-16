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

    # CC weights on [-1, 1] must sum to 2 (∫₋₁¹ 1 dx = 2)
    w_unit = Blast.get_clencurt_weights(-1.0, 1.0, N)
    @test sum(w_unit) ≈ 2.0
end

@testset "Utils: ChebPoly and Bessel" begin
    # Precomputation tests for Bessel/Cheb evaluation
    min_v = 10^(-1) * (1 + 1e-10)
    max_v = 10.0
    N = 100
    x_check = Blast.get_clencurt_grid(min_v, max_v, Int(N))

    ncheb = 2
    ell = 1
    test_chi = zeros(10)
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

@testset "Utils: Akima interpolation" begin
    t = collect(LinRange(0.0, 5.0, 20))
    u = sin.(t)
    t_new = collect(LinRange(0.2, 4.8, 35))

    @testset "Vector interpolation" begin
        # Interpolant should reproduce data values at original nodes
        u_on_nodes = Blast._akima_interpolation(u, t, t)
        @test u_on_nodes ≈ u rtol=1e-12

        # Output size should match query grid
        u_interp = Blast._akima_interpolation(u, t, t_new)
        @test length(u_interp) == length(t_new)
        @test all(isfinite, u_interp)
    end

    @testset "Matrix interpolation" begin
        u_mat = hcat(sin.(t), cos.(t), sin.(2 .* t))

        interp_mat = Blast._akima_interpolation(u_mat, t, t_new)
        interp_cols = hcat([Blast._akima_interpolation(u_mat[:, i], t, t_new) for i in 1:size(u_mat, 2)]...)

        @test size(interp_mat) == (length(t_new), size(u_mat, 2))
        @test interp_mat ≈ interp_cols rtol=1e-12
    end
end

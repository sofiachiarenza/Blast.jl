using Test
using Blast

using FastTransforms

@testset "Matrix product test" begin
    i = 3
    j = 7
    l = 10 
    k = 8 
    c = zeros(j,k,l)
    T = rand(i,j,k,l)
    w_check = zeros(i,j,k)
    w_blast = Blast.w_ell_tullio(c,T)
    @test isapprox(w_check, w_blast)
end

@testset "Precomputation tests" begin
    min = -1 
    max = 1
    N = 100
    x_blast = Blast.get_clencurt_grid(min, max, N)
    x_check = FastTransforms.clenshawcurtisnodes(Float64, N)
    @test isapprox(x_check, x_blast)

    w_blast = Blast.get_clencurt_weights(min, max, N)
    object = FastTransforms.chebyshevmoments1(Float64, N)
    w_check = FastTransforms.clenshawcurtisweights(object)
    @test isapprox(w_check, w_blast)

    min = 10 ^ (-1) .* (1+1e-10)
    max = 10

    x_check = Blast.get_clencurt_grid(min, max, N)

    ncheb = 2
    ell = 1
    test_chi = zeros(10)
    T_blast, Bessel_blast = Blast.Bessel_Cheb_eval(ell, min, max, test_chi, ncheb, N)

    T_check = zeros(ncheb+1, N)
    T_check[1, :] = ones(N)
    T_check[2, :] = log10.(x_check)
    T_check[3, :] = 2 .* log10.(x_check).^2 .- 1

    Bessel_check = zeros(10, N)

    @test isapprox(Bessel_check, Bessel_blast)
    @test isapprox(T_check[1, :], T_blast[1,:])
    @test isapprox(T_check[2, :], T_blast[2,:])
    @test isapprox(T_check[3, :], T_blast[3,:])
end
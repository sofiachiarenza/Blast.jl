using Test
using Blast

using FastTransforms
using FastChebInterp
using NPZ

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

# downloading files LL
run(`wget https://zenodo.org/api/records/13870111/files-archive`)
run(`mv files-archive T_tilde_LL_check.tar.xz`)
run(`mkdir T_tilde_LL_check`)
run(`tar xvf T_tilde_LL_check.tar.xz -C T_tilde_LL_check`)
run(`rm T_tilde_LL_check.tar.xz`)

#CL
run(`wget https://zenodo.org/api/records/13870094/files-archive`)
run(`mv files-archive T_tilde_CL_check.tar.xz`)
run(`mkdir T_tilde_CL_check`)
run(`tar xvf T_tilde_CL_check.tar.xz -C T_tilde_CL_check`)
run(`rm T_tilde_CL_check.tar.xz`)

#CC
run(`wget https://zenodo.org/api/records/13870029/files-archive`)
run(`mv files-archive T_tilde_CC_check.tar.xz`)
run(`mkdir T_tilde_CC_check`)
run(`tar xvf T_tilde_CC_check.tar.xz -C T_tilde_CC_check`)
run(`rm T_tilde_CC_check.tar.xz`)

T_CC_check = zeros(3,96,48,120)
T_CC_check[1,:,:,:] = npzread("T_tilde_CC_check/T_tilde_l_2.0.npy")[:,:,:,1:120]
T_CC_check[2,:,:,:] = npzread("T_tilde_CC_check/T_tilde_l_97.1.npy")[:,:,:,1:120]
T_CC_check[3,:,:,:] = npzread("T_tilde_CC_check/T_tilde_l_211.6.npy")[:,:,:,1:120]

T_CL_check = zeros(3,96,48,120)
T_CL_check[1,:,:,:] = npzread("T_tilde_CL_check/T_tilde_l_2.0.npy")[:,:,:,1:120]
T_CL_check[2,:,:,:] = npzread("T_tilde_CL_check/T_tilde_l_97.1.npy")[:,:,:,1:120]
T_CL_check[3,:,:,:] = npzread("T_tilde_CL_check/T_tilde_l_211.6.npy")[:,:,:,1:120]

T_LL_check = zeros(3,96,48,120)
T_LL_check[1,:,:,:] = npzread("T_tilde_LL_check/T_tilde_l_2.0.npy")[:,:,:,1:120]
T_LL_check[2,:,:,:] = npzread("T_tilde_LL_check/T_tilde_l_97.1.npy")[:,:,:,1:120]
T_LL_check[3,:,:,:] = npzread("T_tilde_LL_check/T_tilde_l_211.6.npy")[:,:,:,1:120]

run(`rm -r T_tilde_CC_check`)
run(`rm -r T_tilde_CL_check`)
run(`rm -r T_tilde_LL_check`)


@testset "T tilde validation" begin
    ℓs = [2.0, 97.07777459, 211.63514264]
    T_CC_blast = zeros(3,96,48,120)
    T_CL_blast = zeros(3,96,48,120)
    T_LL_blast = zeros(3,96,48,120)

    chi = LinRange(26,7000,96)
    R = FastChebInterp.chebpoints(96, -1, 1)
    R = reverse(R[R.>0])
    kmax = 200/13 
    kmin = 2.5/7000

    for (il, l) in enumerate(ℓs)
        println("Performing integration for ℓ=$l...\n")
        T_CC_blast[il,:,:,:] = Blast.compute_T̃(l, chi, R, kmin, kmax, 2)
        T_CL_blast[il,:,:,:] = Blast.compute_T̃(l, chi, R, kmin, kmax, 0)
        T_LL_blast[il,:,:,:] = Blast.compute_T̃(l, chi, R, kmin, kmax, -2)
    end

    @test isapprox(T_CC_check, T_CC_blast)
    @test isapprox(T_CL_check, T_CL_blast)
    @test isapprox(T_LL_check, T_LL_blast)
end
using Test
using Blast

using FastTransforms
using FastChebInterp
using NPZ
using FFTW
using PhysicalConstants
using DataInterpolations
using Tullio

import PhysicalConstants.CODATA2018: c_0
const C_LIGHT = c_0.val * 10^(-3) #speed of light in Km/s

input_path = pwd()

run(`wget --content-disposition "https://zenodo.org/record/13984489/files/dNdzs_fullwidth.npz?download=1"`)
n5k_bins = npzread(input_path*"/dNdzs_fullwidth.npz")
run(`bash -c "rm dNdzs_fullwidth.npz"`)

run(`wget --content-disposition "https://zenodo.org/records/13984491/files/LJ_clustering_kernels.npz?download=1"`)
LJ_clustering_kernels = npzread(input_path*"/LJ_clustering_kernels.npz")
run(`bash -c "rm LJ_clustering_kernels.npz"`)

run(`wget --content-disposition "https://zenodo.org/records/13984495/files/LJ_shear_kernels.npz?download=1"`)
LJ_shear_kernels = npzread(input_path*"/LJ_shear_kernels.npz")
run(`bash -c "rm LJ_shear_kernels.npz"`)

run(`wget --content-disposition "https://zenodo.org/records/13984500/files/LJ_cmb_kernel.npz?download=1"`)
LJ_cmb_kernel = npzread(input_path*"/LJ_cmb_kernel.npz")
run(`bash -c "rm LJ_cmb_kernel.npz"`)

@testset "Background checks" begin
    cosmo = Blast.FlatΛCDM()
    z_range = Array(LinRange(0., 3.5, 1000))
    grid = Blast.CosmologicalGrid(z_range=z_range)
    bg = Blast.BackgroundQuantities(
    Hz_array = zeros(length(z_range)), χz_array=zeros(length(z_range)) )
    
    E0_test = Blast.compute_adimensional_hubble_factor(0., cosmo)
    @test E0_test == 1.
    H0_test = Blast.compute_hubble_factor(0., cosmo)
    @test H0_test == cosmo.H0

    #now check for the function evaluate_background_quantities
    test_H_array = zeros(length(z_range))
    test_χ_array = zeros(length(z_range))
    for (iz, z) in enumerate(z_range)
        test_H_array[iz] = Blast.compute_hubble_factor(z, cosmo)
        test_χ_array[iz] = Blast.compute_χ(z, cosmo)
    end
    Blast.evaluate_background_quantities!(grid, bg, cosmo)

    @test test_H_array ≈ bg.Hz_array
    @test test_χ_array ≈ bg.χz_array

    #testing the kernels - comparing to LimberJack 
    blast_cl_ker = zeros(10, length(grid.z_range))
    blast_sh_ker = zeros(5, length(grid.z_range))
    blast_cmb_ker = zeros(length(grid.z_range))
   
    print("Computing clustering kernels...\n")
    for i in 1:10
        interp = DataInterpolations.AkimaInterpolation(n5k_bins["dNdz_cl"][:,i], n5k_bins["z_cl"], extrapolate = true)
        GK = Blast.GalaxyKernel(zeros(length(z_range)))
        Blast.compute_kernel!(Array(interp.(z_range)), GK, grid, bg, cosmo)
        blast_cl_ker[i,:] = GK.Kernel
    end

    print("Computing lensing kernels...\n")
    for i in 1:5
        interp = DataInterpolations.AkimaInterpolation(n5k_bins["dNdz_sh"][:,i], n5k_bins["z_sh"], extrapolate = true)
        SHK = Blast.ShearKernel(zeros(length(z_range)))
        Blast.compute_kernel!(Array(interp.(z_range)), SHK, grid, bg, cosmo)
        blast_sh_ker[i,:] = SHK.Kernel
    end

    print("Computing CMB kernel...\n")
    CMBK = Blast.CMBLensingKernel(zeros(length(z_range)))
    Blast.compute_kernel!(CMBK, grid, bg, cosmo)
    blast_cmb_ker = CMBK.Kernel

    @test isapprox(blast_cl_ker, LJ_clustering_kernels, rtol=1e-5)
    @test isapprox(blast_sh_ker, LJ_shear_kernels, rtol=1e-2)
    @test isapprox(blast_cmb_ker, LJ_cmb_kernel, rtol=1e-5)
    

end

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
    T_blast, Bessel_blast = Blast.bessel_cheb_eval(ell, min, max, test_chi, ncheb, N)

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

#CC
run(`wget --content-disposition https://zenodo.org/api/records/13885803/files-archive`)
run(`unzip 13885803.zip`)
run(`rm -r 13885803.zip`)

input_path = pwd()

T_CC_check = zeros(3,10,10,120)
T_CC_check[1,:,:,:] = npzread(input_path*"/T_tilde_l_2.0.npy")
T_CC_check[2,:,:,:] = npzread(input_path*"/T_tilde_l_97.1.npy")
T_CC_check[3,:,:,:] = npzread(input_path*"/T_tilde_l_211.6.npy")

run(`bash -c "rm T_tilde_l_*"`)

@testset "T tilde CC validation" begin
    ℓs = [2.0, 97.07777459, 211.63514264]
    T_CC_blast = zeros(3,10,10,120)

    chi = LinRange(26,7000,10)
    R = FastChebInterp.chebpoints(20, -1, 1)
    R = reverse(R[R.>0])
    kmax = 200/13 
    kmin = 2.5/7000

    for (il, l) in enumerate(ℓs)
        println("Performing integration for ℓ=$l...")
        T_CC_blast[il,:,:,:] = Blast.compute_T̃(l, chi, R, kmin, kmax, 2, N=2^(14)+1)
    end

    @test isapprox(T_CC_check, T_CC_blast)
end


#CL
run(`wget --content-disposition https://zenodo.org/api/records/13885823/files-archive`)
run(`unzip 13885823.zip`)
run(`rm -r 13885823.zip`)

T_CL_check = zeros(3,10,10,120)
T_CL_check[1,:,:,:] = npzread(input_path*"/T_tilde_l_2.0.npy")
T_CL_check[2,:,:,:] = npzread(input_path*"/T_tilde_l_97.1.npy")
T_CL_check[3,:,:,:] = npzread(input_path*"/T_tilde_l_211.6.npy")

run(`bash -c "rm T_tilde_l_*"`)

@testset "T tilde CL validation" begin
    ℓs = [2.0, 97.07777459, 211.63514264]
    T_CL_blast = zeros(3,10,10,120)

    chi = LinRange(26,7000,10)
    R = FastChebInterp.chebpoints(20, -1, 1)
    R = reverse(R[R.>0])
    kmax = 200/13 
    kmin = 2.5/7000

    for (il, l) in enumerate(ℓs)
        println("Performing integration for ℓ=$l...")
        T_CL_blast[il,:,:,:] = Blast.compute_T̃(l, chi, R, kmin, kmax, 0,  N=2^(14)+1)
    end

    @test isapprox(T_CL_check, T_CL_blast)
end

# downloading files LL
run(`wget --content-disposition https://zenodo.org/api/records/13885822/files-archive`)
run(`unzip 13885822.zip`)
run(`rm -r 13885822.zip`)

T_LL_check = zeros(3,10,10,120)
T_LL_check[1,:,:,:] = npzread(input_path*"/T_tilde_l_2.0.npy")
T_LL_check[2,:,:,:] = npzread(input_path*"/T_tilde_l_97.1.npy")
T_LL_check[3,:,:,:] = npzread(input_path*"/T_tilde_l_211.6.npy")

run(`bash -c "rm T_tilde_l_*"`)

@testset "T tilde LL validation" begin
    ℓs = [2.0, 97.07777459, 211.63514264]
    T_LL_blast = zeros(3,10,10,120)

    chi = LinRange(26,7000,10)
    R = FastChebInterp.chebpoints(20, -1, 1)
    R = reverse(R[R.>0])
    kmax = 200/13 
    kmin = 2.5/7000

    for (il, l) in enumerate(ℓs)
        println("Performing integration for ℓ=$l...")
        T_LL_blast[il,:,:,:] = Blast.compute_T̃(l, chi, R, kmin, kmax, -2,  N=2^(14)+1)
    end

    @test isapprox(T_LL_check, T_LL_blast)
end

@testset "Chebyshev coefficients" begin
    dims = (2^10, 2^6)
    A = rand(Float64, dims)

    plan = Blast.plan_fft(A)
    my_coefs = Blast.fast_chebcoefs(A, plan)

    true_coefs = zeros(dims)
    for i in 1:dims[2]
        true_coefs[:,i] = FastChebInterp.chebcoefs(A[:,i])
    end
    true_coefs
    

    @test true_coefs ≈ my_coefs
end

@testset "Outer integrals tests" begin
    n = 100
    x = LinRange(0,1,n)
    Δx = ((last(x)-first(x))/(n-1))
    weights = Blast.simpson_weight_array(n)

    @tullio integral = x[i]*weights[i]*Δx

    @test integral ≈ 0.5

    pmd = ones(1, 200, 50)
    kernel = ones(1, 1, 200, 50)
    χ = LinRange(10, 100, 200) 
    R = chebpoints(100,-1,1)
    R = reverse(R[R.>0])

    cl_test = Blast.compute_Cℓ(pmd, kernel, χ, R)
    cl_true = 4950*(R[end]-R[1])

    @test isapprox(cl_test[1,1,1], cl_true, rtol = 1e-5) 

end
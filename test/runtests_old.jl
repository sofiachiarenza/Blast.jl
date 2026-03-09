using Test
using Blast

using FastTransforms
using FastChebInterp
using NPZ
using FFTW
using PhysicalConstants
using DataInterpolations
using Tullio
using StaticArrays

import PhysicalConstants.CODATA2018: c_0
const C_LIGHT = c_0.val * 10^(-3) #speed of light in Km/s

input_path = pwd()

run(`wget --content-disposition "https://zenodo.org/records/13997096/files/bins.npz?download=1"`)
bins = npzread(input_path*"/bins.npz")
run(`bash -c "rm bins.npz"`)

run(`wget --content-disposition "https://zenodo.org/records/13996320/files/LJ_clustering_kernels.npz?download=1"`)
LJ_clustering_kernels = npzread(input_path*"/LJ_clustering_kernels.npz")
run(`bash -c "rm LJ_clustering_kernels.npz"`)

run(`wget --content-disposition "https://zenodo.org/records/13996321/files/LJ_shear_kernels.npz?download=1"`)
LJ_shear_kernels = npzread(input_path*"/LJ_shear_kernels.npz")[1:3,:]
run(`bash -c "rm LJ_shear_kernels.npz"`)

run(`wget --content-disposition "https://zenodo.org/records/13997095/files/LJ_cmb_kernel.npz?download=1"`)
LJ_cmb_kernel = npzread(input_path*"/LJ_cmb_kernel.npz")
run(`bash -c "rm LJ_cmb_kernel.npz"`)

@testset "Background checks" begin
    cosmo = Blast.FlatΛCDM()
    z_range = Array(LinRange(0., 4, 1000))
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

    print("Computing clustering kernels...\n")

    nz_interp = zeros(10, length(z_range))
    for i in 1:10
        interp = DataInterpolations.AkimaInterpolation(bins["dNdz"][i,:], bins["z"], extrapolate = true)
        nz_interp[i,:] = interp.(z_range)
    end

    GK = Blast.GalaxyKernel(10, length(grid.z_range) )
    Blast.compute_kernel!(bins["dNdz"], bins["z"], GK, ones(10), grid, bg, cosmo)


    print("Computing shear kernels...\n")
    SHK = Blast.ShearKernel(3, length(grid.z_range))
    Blast.compute_kernel!(bins["dNdz"][1:3,:], bins["z"], SHK, grid, bg, cosmo)

    print("Computing CMB kernels...\n")
    CMBK = Blast.CMBLensingKernel(length(grid.z_range))
    Blast.compute_kernel!(CMBK, grid, bg, cosmo)

    @test isapprox(GK.Kernel, LJ_clustering_kernels, rtol=1e-5)
    @test isapprox(SHK.Kernel, LJ_shear_kernels, rtol=1e-3)
    @test isapprox(CMBK.Kernel, LJ_cmb_kernel, rtol=1e-5)
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

    plan = Blast.plan_fft(A, 1)
    my_coefs = Blast.fast_chebcoefs(A, plan)

    true_coefs = zeros(dims)
    for i in 1:dims[2]
        true_coefs[:,i] = FastChebInterp.chebcoefs(A[:,i])
    end
    true_coefs
    
    @test true_coefs ≈ my_coefs
end

@testset "Outer integrals tests" begin

    elle = 7

    a = Blast.factorial_frac(elle)
    b = factorial(elle + 2) / factorial(elle - 2)

    @test a ≈ b

    n = 100
    x = LinRange(0,1,n)
    Δx = ((last(x)-first(x))/(n-1))
    weights = Blast.simpson_weight_array(n)

    @tullio integral = x[i]*weights[i]*Δx

    @test integral ≈ 0.5

    pmd = ones(1, 200, 50)

    Probe1 = Blast.GalaxyKernel(1, 200)
    Probe1.Kernel = ones(1,200)
    Probe2 = Blast.GalaxyKernel(1,200)
    Probe2.Kernel = ones(1,200) 

    χ = Array(LinRange(10, 100, 200))
    R = chebpoints(100,-1,1)
    R = reverse(R[R.>0])

    cosmo = Blast.FlatΛCDM()
    bg = Blast.BackgroundQuantities(Hz_array = zeros(200), χz_array=χ )

    cl_test = Blast.compute_Cℓ(pmd, Probe1, Probe2, bg, R, Float64[1.0]) 
    cl_true = 4950*(R[end]-R[1]) * 2 / π * 2 #2/pi is the ell prefactor, the other 2 comes from the window combination!

    @test isapprox(cl_test[1,1,1], cl_true, rtol = 1e-5) 

    nχ = 200
    z = LinRange(0.01, 4, nχ)
    R = chebpoints(100,-1,1)
    R = reverse(R[R.>0])
    nR = length(R)

    cosmo = Blast.FlatΛCDM()
    bg = Blast.BackgroundQuantities(Hz_array = zeros(nχ), χz_array = zeros(nχ))
    grid = Blast.CosmologicalGrid(z_range = z)
    Blast.evaluate_background_quantities!(grid, bg, cosmo)

    Probe1 = Blast.GalaxyKernel(1, nχ)
    Probe1.Kernel = ones(1,nχ)
    Probe2 = Blast.ShearKernel(1,nχ)
    Probe2.Kernel = ones(1,nχ) 
    nz = rand(3,nχ)
    Blast.compute_kernel!(nz, z, Probe1, ones(1), grid, bg, cosmo)
    Blast.compute_kernel!(nz, z, Probe2, grid, bg, cosmo)

    w = rand(length(Blast.ℓ), nχ, nR)

    Cℓ_mod1 = Blast.compute_Cℓ(w, Probe1, Probe2, bg, R)

    w_χ = Blast.simpson_weight_array(nχ)
    w_R = Blast.get_clencurt_weights_R_integration(2*nR+1)
    pref= Blast.get_ell_prefactor(Probe1, Probe2, Blast.ℓ)

    pref_check = Blast.get_ell_prefactor(Probe2, Probe1, Blast.ℓ)

    @test pref ≈ pref_check

    pref_LL = Blast.get_ell_prefactor(Probe2, Probe2, Blast.ℓ)

    @test pref_LL[1] ≈ 2 / π * factorial(4)/factorial(0) 

    K = Blast.combine_kernels(Probe1, Probe2, bg, R)

    Cℓ_mod2 = Blast.compute_Cℓ(w, K, bg, w_χ, w_R, pref)

    @test Cℓ_mod1 ≈ Cℓ_mod2

end

@testset "Tomographic bins combination" begin
    a = Vector([1.,2.,3.,4.,5.])
    b = ones(8)

    cosmo = Blast.FlatΛCDM()
    bg = Blast.BackgroundQuantities(Hz_array = zeros(5), χz_array = a )

    true_res = vcat(fill(a, length(b))...)

    @test true_res ≈ Blast.make_grid(bg, b)

    GK = Blast.GalaxyKernel(1,length(a))
    GK.Kernel = ones(size(GK.Kernel,1), size(GK.Kernel,2))

    SHK = Blast.ShearKernel(1, length(a))
    SHK.Kernel = ones(size(SHK.Kernel,1), size(SHK.Kernel,2))

    CK = Blast.CMBLensingKernel(length(a))
    CK.Kernel = ones(size(CK.Kernel))

    theory_gal = ones(1, length(a), length(b))
    theory_sh = ones(1, length(a), length(b)) ./ reshape(a .^ 2, 1, 5, 1)
    theory_cmb = ones(1, length(a), length(b)) ./ reshape(a .^ 2, 1, 5, 1)


    @test theory_gal ≈ Blast.get_kernel_array(GK, bg, b)
    @test theory_sh ≈ Blast.get_kernel_array(SHK, bg, b)
    @test theory_cmb ≈ Blast.get_kernel_array(CK, bg, b)
end

run(`wget --content-disposition "https://zenodo.org/records/14192971/files/pk_n5k_cheb.npz?download=1"`)
pk_n5k = npzread(input_path*"/pk_n5k_cheb.npz")
run(`bash -c "rm pk_n5k_cheb.npz"`)

run(`wget --content-disposition "https://zenodo.org/records/14193379/files/n5k_zs.npz?download=1"`)
n5k_zs = npzread(input_path*"/n5k_zs.npz")
run(`bash -c "rm n5k_zs.npz"`)

@testset "Power spectrum interpolation tests" begin
    x = rand(1000) * 10.0 .+ 1.0  
    n_cheb = 120  
    #Blast evaluation of ChebPolys
    Tcheb_test = Blast.chebyshev_polynomials(x, n_cheb, minimum(x), maximum(x))

    #Standard evaluation of ChebPolys
    x_cheb = chebpoints(n_cheb-1, minimum(x), maximum(x)) 
    c = FastChebInterp.ChebPoly(x_cheb, SA[minimum(x)], SA[maximum(x)])
    T_reference = zeros(n_cheb, length(x))
    for i in 1:n_cheb
        copy_c = deepcopy(c)
        copy_c.coefs .= 0
        copy_c.coefs[i] = 1.0
        T_reference[i, :] = copy_c.(x)
    end

    @test Tcheb_test ≈ T_reference

    cosmo =Blast.FlatΛCDM()
    z_range = Array([0.0, 0.5, 1.0, 2.0])
    grid = Blast.CosmologicalGrid(z_range = z_range)
    bg = Blast.BackgroundQuantities(Hz_array = zeros(length(z_range)), χz_array = zeros(length(z_range)))
    Blast.evaluate_background_quantities!(grid, bg, cosmo)
    R = 0.5
    new_χ = bg.χz_array .* R

    z_of_χ = DataInterpolations.AkimaInterpolation(grid.z_range, bg.χz_array, extrapolate=true)

    new_z = Blast.resample_redshifts(bg, grid, new_χ)
    test_z = z_of_χ.(new_χ)

    @test new_z ≈ test_z

    nχ = 96
    χ = LinRange(26, 7000, nχ)
    R = chebpoints(96, -1, 1)
    R = reverse(R[R.>0])
    kmax = 200/13 
    kmin = 2.5/7000
    n_cheb = 119
    k_cheb = chebpoints(n_cheb, log10(kmin), log10(kmax))
    cosmo = Blast.FlatΛCDM()
    z_range = n5k_zs #redshifts corresponding to the χ array
    bgrid = Blast.CosmologicalGrid(z_range = z_range)
    bg = Blast.BackgroundQuantities(Hz_array = zeros(length(z_range)), χz_array = Array(χ))
    z_cheb = chebpoints(7, 0, 3.5)
    plan = Blast.plan_fft(pk_n5k, 1)
    blast_pk = Blast.interpolate_power_spectrum(pk_n5k, z_cheb, R, plan, bg, bgrid)

    fastcheb_check = zeros(n_cheb+1, length(χ), length(R))
    newzs = Blast.resample_redshifts(bg, bgrid, Blast.make_grid(bg, R))
   
    for i in 1:n_cheb+1
        rightinterp = chebinterp(pk_n5k[:,i], minimum(newzs), maximum(newzs))
        fastcheb_check[i,:,:] = reshape(rightinterp.(newzs) , 96, 48)
    end

    @test blast_pk ≈ fastcheb_check

    pk = ones(3, 3, 3)
    expected_output = ones(3, 3, 3)
    result = Blast.unequal_time_power_spectrum(pk)
    @test result == expected_output

end

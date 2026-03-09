using Test
using Blast
using NPZ
using DataInterpolations
using StaticArrays
using FastTransforms
using FastChebInterp
using Downloads

"""
    load_reference_data()

Downloads (if missing) and loads all reference data from Zenodo for validation.
"""
function load_reference_data()
    data = Dict{String, Any}()
    
    # 1. Bins and basic kernels
    if !isfile("bins.npz")
        Downloads.download("https://zenodo.org/records/13997096/files/bins.npz?download=1", "bins.npz")
    end
    data["bins"] = npzread("bins.npz")

    if !isfile("LJ_clustering_kernels.npz")
        Downloads.download("https://zenodo.org/records/13996320/files/LJ_clustering_kernels.npz?download=1", "LJ_clustering_kernels.npz")
    end
    data["LJ_clustering"] = npzread("LJ_clustering_kernels.npz")

    if !isfile("LJ_shear_kernels.npz")
        Downloads.download("https://zenodo.org/records/13996321/files/LJ_shear_kernels.npz?download=1", "LJ_shear_kernels.npz")
    end
    data["LJ_shear"] = npzread("LJ_shear_kernels.npz")[1:3,:]

    if !isfile("LJ_cmb_kernel.npz")
        Downloads.download("https://zenodo.org/records/13997095/files/LJ_cmb_kernel.npz?download=1", "LJ_cmb_kernel.npz")
    end
    data["LJ_cmb"] = npzread("LJ_cmb_kernel.npz")

    # 2. T_tilde artifacts (CC, CL, LL)
    if !isfile("T_tilde_l_2.0.npy")
        println("Downloading T_tilde artifacts...")
        # Note: Zenodo zip-archive links are slightly different
        Downloads.download("https://zenodo.org/api/records/13885803/files-archive", "CC.zip")
        run(`unzip -o -q CC.zip`)
        Downloads.download("https://zenodo.org/api/records/13885823/files-archive", "CL.zip")
        run(`unzip -o -q CL.zip`)
        Downloads.download("https://zenodo.org/api/records/13885822/files-archive", "LL.zip")
        run(`unzip -o -q LL.zip`)
        
        rm("CC.zip", force=true)
        rm("CL.zip", force=true)
        rm("LL.zip", force=true)
    end

    # 3. Power spectrum interpolation data
    if !isfile("pk_n5k_cheb.npz")
        Downloads.download("https://zenodo.org/records/14192971/files/pk_n5k_cheb.npz?download=1", "pk_n5k_cheb.npz")
    end
    data["pk_n5k"] = npzread("pk_n5k_cheb.npz")

    if !isfile("n5k_zs.npz")
        Downloads.download("https://zenodo.org/records/14193379/files/n5k_zs.npz?download=1", "n5k_zs.npz")
    end
    data["n5k_zs"] = npzread("n5k_zs.npz")

    return data
end

"""
    get_test_cosmo()
Returns a standard ΛCDM cosmology for testing.
"""
function get_test_cosmo()
    return ΛCDM(H0=67.27, Ωm=0.3156, Ωb=0.0492)
end

"""
    get_test_bg(cosmo=get_test_cosmo())
Returns a Background object built on the default Blast.χ grid.
"""
function get_test_bg(cosmo=get_test_cosmo())
    return Background(cosmo)
end

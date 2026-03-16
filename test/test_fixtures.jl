using Test
using Blast
using NPZ
using DataInterpolations
using StaticArrays
using FastTransforms
using FastChebInterp
using Downloads

function _load_t_tilde_reference!(data::Dict{String, Any}, record_id::Int, sector::String, ell_tag::String)
    local_file = "T_tilde_$(sector)_l_$(ell_tag).npy"

    if !isfile(local_file)
        tmpdir = mktempdir()
        archive = joinpath(tmpdir, "files-archive.zip")
        Downloads.download("https://zenodo.org/api/records/$(record_id)/files-archive", archive)
        run(`unzip -o -q $archive -d $tmpdir`)
        cp(joinpath(tmpdir, "T_tilde_l_$(ell_tag).npy"), local_file; force=true)
    end

    t_tilde = get!(data, "T_tilde", Dict{String, Dict{String, Any}}())
    sector_refs = get!(t_tilde, sector, Dict{String, Any}())
    sector_refs[ell_tag] = npzread(local_file)
end

"""
    load_reference_data()

Downloads (if missing) and loads the reference data used by the modular test suite.
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

    # 2. T_tilde artifacts used by the legacy regression tests
    _load_t_tilde_reference!(data, 13885803, "CC", "2.0")
    _load_t_tilde_reference!(data, 13885823, "CL", "2.0")
    _load_t_tilde_reference!(data, 13885822, "LL", "2.0")

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

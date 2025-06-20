module Blast

using LoopVectorization
using Tullio
using FastTransforms
using FastChebInterp
using SpecialFunctions
using DataInterpolations
using Interpolations
using StaticArrays
using FFTW
using NPZ
using QuadGK
using Artifacts
using PhysicalConstants

include("cosmo.jl")
include("deprecated.jl")
include("probes.jl")
include("chebcoefs.jl")
include("power_spectrum.jl")
include("projected_matter.jl")
include("integrals.jl")
include("cls.jl")
include("utils.jl")
include("limber.jl")

import PhysicalConstants.CODATA2018: c_0

const C_LIGHT = c_0.val * 10^(-3) #speed of light in Km/s

function load_precomputed_Ts(folder::String)
    #ell_vector = npzread(joinpath(folder, "ell_list.npy"))
    full_T = npzread(joinpath(folder, "T_tilde.npy"))
    #return ell_vector, full_T
    return full_T
end

struct T̃_data
    T_2_00::AbstractArray{<:Any, 4}
    T_0_00::AbstractArray{<:Any, 4}
    T_minus2_00::AbstractArray{<:Any, 4}
    T_0_02::AbstractArray{<:Any, 4}
    T_0_20::AbstractArray{<:Any, 4}
    T_2_02::AbstractArray{<:Any, 4}
    T_2_20::AbstractArray{<:Any, 4}
    T_2_22::AbstractArray{<:Any, 4}
end

function __init__()

    global T_tilde_m2 = load_precomputed_Ts(artifact"T_tilde_2")
    global T_tilde_0 = load_precomputed_Ts(artifact"T_tilde_0")
    global T_tilde_p2 = load_precomputed_Ts(artifact"T_tilde_-2")

    global full_ℓ_range = reverse(chebpoints(100, 2, 2000))

    global ℓ_nonlimber = full_ℓ_range[full_ℓ_range .< 220]
    global ℓ_limber = full_ℓ_range[full_ℓ_range .> 220]

    nχ = 200
    global χ = Array(LinRange(26, 7000, nχ))
    R = chebpoints(200, -1, 1)
    global R = reverse(R[R.>0])
    nR = length(R)
    kmax = 200/13 
    kmin = 2.5/7000
    n_cheb = 119
    global k_cheb = chebpoints(n_cheb, log10(kmin), log10(kmax))
    global k_limber = chebpoints(599, log10(1e-4), log10(80))
    global z_cheb = chebpoints(49, 0, 3.6)

    global T_tildes = T̃_data(
        Blast.load_Ts("/Users/sofiachiarenza/Desktop/PhD_Stuff/Blastoise/T_tildes/high_precision/T_tilde_2", nχ, nR, n_cheb+1),
        Blast.load_Ts("/Users/sofiachiarenza/Desktop/PhD_Stuff/Blastoise/T_tildes/high_precision/T_tilde_0", nχ, nR, n_cheb+1),
        Blast.load_Ts("/Users/sofiachiarenza/Desktop/PhD_Stuff/Blastoise/T_tildes/high_precision/T_tilde_-2", nχ, nR, n_cheb+1),
        Blast.load_Ts("/Users/sofiachiarenza/Desktop/PhD_Stuff/Blastoise/T_tildes/high_precision/T_tilde_0_02", nχ, nR, n_cheb+1),
        Blast.load_Ts("/Users/sofiachiarenza/Desktop/PhD_Stuff/Blastoise/T_tildes/high_precision/T_tilde_0_20", nχ, nR, n_cheb+1),
        Blast.load_Ts("/Users/sofiachiarenza/Desktop/PhD_Stuff/Blastoise/T_tildes/high_precision/T_tilde_2_02", nχ, nR, n_cheb+1),
        Blast.load_Ts("/Users/sofiachiarenza/Desktop/PhD_Stuff/Blastoise/T_tildes/high_precision/T_tilde_2_20", nχ, nR, n_cheb+1),
        Blast.load_Ts("/Users/sofiachiarenza/Desktop/PhD_Stuff/Blastoise/T_tildes/high_precision/T_tilde_2_22", nχ, nR, n_cheb+1)
    )

end

end # module Blast

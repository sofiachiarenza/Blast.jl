module Blast

using LoopVectorization
using Tullio
using FastTransforms
using FastChebInterp
using SpecialFunctions
using DataInterpolations
using StaticArrays
using FFTW
using NPZ
using QuadGK
using Artifacts
using PhysicalConstants

include("projected_matter.jl")
include("chebcoefs.jl")
include("integrals.jl")
include("cosmo.jl")
include("background.jl")

import PhysicalConstants.CODATA2018: c_0

const C_LIGHT = c_0.val * 10^(-3) #speed of light in Km/s

function load_precomputed_Ts(folder::String)
    ell_vector = npzread(joinpath(folder, "ell_list.npy"))
    full_T = npzread(joinpath(folder, "T_tilde.npy"))
    return ell_vector, full_T
end

function __init__()

    global ℓ, T_tilde_m2 = load_precomputed_Ts(artifact"T_tilde_2")
    global ℓ, T_tilde_0 = load_precomputed_Ts(artifact"T_tilde_0")
    global ℓ, T_tilde_p2 = load_precomputed_Ts(artifact"T_tilde_-2")

end

end # module Blast

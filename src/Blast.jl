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

#TODO: added limberjack to the project, remove it later and just keep it in the test's project.

include("projected_matter.jl")
include("chebcoefs.jl")
include("integrals.jl")
include("cosmo.jl")
include("background.jl")

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

module Blast

using LoopVectorization
using Tullio
using FastTransforms
using FastChebInterp
using SpecialFunctions
using StaticArrays
using FFTW
using NPZ
using Artifacts

include("projected_matter.jl")
include("chebcoefs.jl")
include("integrals.jl")

function load_precomputed_Ts(folder::String; nχ::Int = 96, nR::Int = 48, n_cheb::Int = 120)
    ell_vector = npzread(joinpath(folder, "ell_list.npy"))
    nℓ = length(ell_vector)
    full_T = zeros(nℓ, nχ, nR, n_cheb)
    for i in 1:nℓ
        l_string = string(round(ell_vector[i]; digits=1))
        filename = joinpath(folder, "T_tilde_l_$l_string.npy")
        if isfile(filename)
            full_T[i,:,:,:] = npzread(filename)
        else
            throw(LoadError(filename, "File doesn't exist."))
        end
    end
    return full_T
end

function __init__()

    global T_tilde_m2 = load_precomputed_Ts(artifact"T_tilde_2")
    global T_tilde_0 = load_precomputed_Ts(artifact"T_tilde_0")
    global T_tilde_p2 = load_precomputed_Ts(artifact"T_tilde_-2")

end

end # module Blast

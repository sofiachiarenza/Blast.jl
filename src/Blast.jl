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

function load_precomputed_Ts(folder; nℓ=22, nχ=96, nR=48, n_cheb=120)
    ell_vector =reverse(chebpoints(100, 2,2000))[1:nℓ]
    full_T = zeros(nℓ, nχ, nR, n_cheb)
    for i in 1:nℓ
        l_string = string(round(ell_vector[i]; digits=1))
        filename = joinpath(folder, "T_tilde_l_$l_string.npy")
        if isfile(filename)
            full_T[i,:,:,:] = npzread(filename)
        else
            println("Missing file!")
        end
    end
    return full_T
end

function __init__()

    global T_tilde_CC = load_precomputed_Ts(artifact"T_tilde_2")
    global T_tilde_CL = load_precomputed_Ts(artifact"T_tilde_0")
    global T_tilde_LL = load_precomputed_Ts(artifact"T_tilde_-2")

end

end # module Blast

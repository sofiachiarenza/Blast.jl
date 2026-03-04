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
using Memoization
using ChainRulesCore
using ChainRules
using Mooncake

include("chebcoefs.jl")

# Define global constants at top-level
const full_ℓ_range = reverse(chebpoints(100, 2, 2000))
const ℓ_nonlimber = full_ℓ_range[full_ℓ_range .< 220]
const ℓ_limber = full_ℓ_range[full_ℓ_range .> 220]

const nχ = 96
const χ = Array(LinRange(26, 7000, nχ))

const _R_nodes = chebpoints(64*2, -1, 1)
const R = reverse(_R_nodes[_R_nodes .> 0])

const k_cheb = chebpoints(160, log10(5e-5), log10(16))
const k_limber = chebpoints(256, log10(1e-4), log10(80))
const z_cheb = chebpoints(49, 0, 3.6)
const z_lin = LinRange(0, 3.6, 50)

# Functional include order
include("utils.jl")
include("cosmo.jl")
include("deprecated.jl")
include("probes.jl")
include("setup.jl")
include("projected_matter.jl")
include("cls.jl")
include("limber.jl")
include("chainrules.jl")

import PhysicalConstants.CODATA2018: c_0
const C_LIGHT = c_0.val * 10^(-3) #speed of light in Km/s

struct T̃{A<:AbstractArray{<:Any,4}}
    T_2_00::A
    T_0_00::A
    T_minus2_00::A
    T_0_02::A
    T_0_20::A
    T_2_02::A
    T_2_20::A
    T_2_22::A
end

function __init__()
    Tdir = artifact"T_tildes"
    global T_tildes = T̃(
        npzread(joinpath(Tdir, "T_tildes_artifact/T_2_00.npz")),
        npzread(joinpath(Tdir, "T_tildes_artifact/T_0_00.npz")),
        npzread(joinpath(Tdir, "T_tildes_artifact/T_minus2_00.npz")),
        npzread(joinpath(Tdir, "T_tildes_artifact/T_0_02.npz")),
        npzread(joinpath(Tdir, "T_tildes_artifact/T_0_20.npz")),
        npzread(joinpath(Tdir, "T_tildes_artifact/T_2_02.npz")),
        npzread(joinpath(Tdir, "T_tildes_artifact/T_2_20.npz")),
        npzread(joinpath(Tdir, "T_tildes_artifact/T_2_22.npz"))
    )
end

end # module Blast

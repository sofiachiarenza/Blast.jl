module Blast

using LoopVectorization
using Tullio
using FastTransforms
using FastChebInterp
using SpecialFunctions
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

# These dependencies are used by the BackgroundCosmologyExt extension
# They MUST be present in the project but do not necessarily need to be 'using'ed here
using DataInterpolations, FastGaussQuadrature, Integrals, LinearAlgebra, OrdinaryDiffEqTsit5, SciMLSensitivity
using AbstractCosmologicalEmulators

# --- 1. HANDLE EXTENSION AND CORE TYPES ---
# We look for the extension that provides differentiable background functions
const ext = Base.get_extension(AbstractCosmologicalEmulators, :BackgroundCosmologyExt)

if !isnothing(ext)
    using .ext: AbstractCosmology, w0waCDMCosmology, D_z, D_f_z, f_z, E_z, d̃A_z, dM_z, dA_z, dL_z, r_z
    # Re-export for user convenience
    export AbstractCosmology, w0waCDMCosmology, D_z, D_f_z, f_z, E_z, d̃A_z, dM_z, dA_z, dL_z, r_z
else
    # Fallback or error if essential differentiable background is missing
    @error "BackgroundCosmologyExt extension not loaded. Please ensure all dependencies are installed."
end

# --- 2. INCLUDE COMPONENT FILES (Order matters!) ---
include("chebcoefs.jl")

# Define global constants at top-level so they are available at include-time
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

include("utils.jl")
include("cosmo.jl")
#include("deprecated.jl")
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

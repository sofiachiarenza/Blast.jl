module Blast

using LoopVectorization
using Tullio
using SpecialFunctions
using StaticArrays
using FFTW
using NPZ
using QuadGK
using Artifacts
using PhysicalConstants
using Memoization
using FastTransforms
using ChainRulesCore
using ChainRules
using Mooncake

# These dependencies are required for the AbstractCosmologicalEmulators background extension
using DataInterpolations, FastGaussQuadrature, Integrals, LinearAlgebra, OrdinaryDiffEqTsit5, SciMLSensitivity
using AbstractCosmologicalEmulators

# --- 1. HANDLE EXTENSION SYMBOLS ---
# Access extension symbols via Base.get_extension to ensure precompilation stability.
# Using 'const' aliases avoids "undeclared binding" warnings.
const cosmo_ext = Base.get_extension(AbstractCosmologicalEmulators, :BackgroundCosmologyExt)

if !isnothing(cosmo_ext)
    const AbstractCosmology = cosmo_ext.AbstractCosmology
    const w0waCDMCosmology = cosmo_ext.w0waCDMCosmology
    const E_z = cosmo_ext.E_z
    const r_z = cosmo_ext.r_z
    const D_z = cosmo_ext.D_z
    const f_z = cosmo_ext.f_z
    
    export AbstractCosmology, w0waCDMCosmology
    export Background, prepare_nz_matrix, NLA_model, ΛCDM, w0waCDM
else
    @error "BackgroundCosmologyExt extension not loaded. Differentiable background will not be available."
end

# --- 2. CORE INFRASTRUCTURE ---
include("chebcoefs.jl")  # Native Chebyshev engine
include("constants.jl")  # Global grids and integration constants

# --- 3. PIPELINE COMPONENTS ---
include("utils.jl")      # Math and interpolation helpers
include("cosmo.jl")      # Unified Background and Parameter accessors
include("probes.jl")     # Redshift distributions and Probes
include("setup.jl")      # FFTPlans and PowerSpectrum workspace
include("projected_matter.jl")
include("cls.jl")
include("limber.jl")     # Limber grid and evaluation
include("chainrules.jl")

# --- 4. PHYSICAL CONSTANTS ---
import PhysicalConstants.CODATA2018: c_0
const C_LIGHT = c_0.val * 10^(-3) # speed of light in km/s

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

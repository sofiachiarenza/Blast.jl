module Blast

using LoopVectorization
using Tullio
using SpecialFunctions
using StaticArrays
using FFTW
using NPZ
using Artifacts
using PhysicalConstants
using Memoization
using FastTransforms
using Loess
using ChainRulesCore
using ChainRules

# These dependencies are required for the AbstractCosmologicalEmulators background extension
using DataInterpolations, FastGaussQuadrature, Integrals, LinearAlgebra, OrdinaryDiffEqTsit5, SciMLSensitivity
using AbstractCosmologicalEmulators

const cosmo_ext = Base.get_extension(AbstractCosmologicalEmulators, :BackgroundCosmologyExt)

if !isnothing(cosmo_ext)
    const AbstractCosmology = cosmo_ext.AbstractCosmology
    const w0waCDMCosmology = cosmo_ext.w0waCDMCosmology
    const E_z = cosmo_ext.E_z
    const r_z = cosmo_ext.r_z
    const D_z = cosmo_ext.D_z
    const f_z = cosmo_ext.f_z
    const D_f_z = cosmo_ext.D_f_z
else
    @error "BackgroundCosmologyExt extension not loaded. Differentiable background will not be available."
end

# Physical constants — defined before includes so all source files can use C_LIGHT
import PhysicalConstants.CODATA2018: c_0
const C_LIGHT = c_0.val * 10^(-3)   # speed of light in km/s

include("constants.jl")
include("chebcoefs.jl")
include("utils.jl")
include("cosmo.jl")
include("probes.jl")
include("setup.jl")
include("projected_matter.jl")
include("cls.jl")
include("limber.jl")
include("chainrules.jl")
include("reactant_api.jl")

# ── Public API ────────────────────────────────────────────────────────────────
# Cosmology
export AbstractCosmology, w0waCDMCosmology
export Background, w0waCDM

# Probes
export AbstractCosmologicalProbes, AbstractComponents
export GalaxyClustering, WeakLensing, CMB

# Components
export NumberCounts, CosmicShear, CMBLensing
export RedshiftSpaceDistortions, MagnificationBias
export IntrinsicAlignment, IntegratedSachsWolfe, PrimordialNonGaussianity

# Setup and workspace
export SetUp, FFTPlans, PowerSpectrum, ProjectedMatterDensity

# Core functions
export evaluate_components!, compute_kernel!, compute_w
export prepare_pk_workspace, get_Cℓ
export reactant_is_available, reactant_p_phi_TT, reactant_p_phi_T
export reactant_chebyshev_matrix, reactant_chebyshev_matmul
export reactant_w_ell_hlo, reactant_chebyshev_w_ell_hlo
export reactant_prepare_nonlimber_spectrum, reactant_compute_w
export reactant_nonlimber_c_ell, reactant_nonlimber_rsd_c_ell
export reactant_limber_contraction
export reactant_finalize_c_ell
export reactant_nonlimber_3x2pt
export reactant_full_3x2pt
export reactant_chebyshev_2d_matmul, reactant_limber_eval, reactant_prepare_limber
export reactant_limber_power_products, reactant_limber_chebyshev_coefficients
export reactant_limber_grid_from_coefficients, reactant_limber_grid_difference
export reactant_compute_w_from_spectrum, reactant_limber_terms_from_prepared
# ─────────────────────────────────────────────────────────────────────────────

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
    paths = joinpath.(Tdir, "T_tildes_artifact", string.(fieldnames(T̃)) .* ".npz")
    global T_tildes = T̃(npzread.(paths)...)
end

end # module Blast

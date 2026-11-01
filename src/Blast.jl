module Blast

using LoopVectorization
using Tullio
using SpecialFunctions
using StaticArrays
using FFTW
using NPZ
using Artifacts
using SHA
using TOML
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
include("bsplinecoefs.jl")
include("utils.jl")
include("radial_quadrature.jl")
include("tensor_banks.jl")
include("spectral_backends.jl")
include("tensor_bundle.jl")
include("bspline_backend.jl")
include("cosmo.jl")
include("probes.jl")
include("setup.jl")
include("projected_matter.jl")
include("cls.jl")
include("limber.jl")
include("chainrules.jl")

# ── Public API ────────────────────────────────────────────────────────────────
# Cosmology. The cosmology type itself is owned by
# AbstractCosmologicalEmulators; Blast does not maintain a competing parameter
# translation layer.
export AbstractCosmology, w0waCDMCosmology
export Background, ω_ν0, ω_m0

# Probes
export AbstractCosmologicalProbes, AbstractComponents
export GalaxyClustering, WeakLensing, CMB

# Components
export NumberCounts, CosmicShear, CMBLensing
export RedshiftSpaceDistortions, MagnificationBias
export IntrinsicAlignment, IntegratedSachsWolfe, PrimordialNonGaussianity
export AngularRedshiftFluctuations

# Setup and workspace
export SetUp, SpectralPlans, FFTPlans, PowerSpectrum, ProjectedMatterDensity
export ChebyshevBackend, CubicBSplineBackend, nonlimber_grid, load_tensor_bank
export RadialQuadraturePlan

# Core functions
export prepare_probe, evaluate_components!, compute_kernel!, compute_w, compute_w!, allocate_compute_w
export prepare_pk_workspace, get_Cℓ
# ─────────────────────────────────────────────────────────────────────────────

end # module Blast

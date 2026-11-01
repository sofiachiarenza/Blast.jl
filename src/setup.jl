# Parametrized on each plan's concrete type so every field access infers to
# a concrete plan, not the abstract `ChebyshevPlan` UnionAll. Without this,
# `prepare_pk_workspace(::SpectralPlans, …)` could not infer its return to
# `PowerSpectrum{T}` — only to the bare `PowerSpectrum` UnionAll — which
# poisoned `f_full`'s inference all the way to the top.
"""
    SpectralPlans

Fixed spectral configuration of one setup: the concrete non-Limber `backend`
(which owns the tensor bank and sampling sites), the source decomposition plans,
and the Limber/ℓ machinery. Constructed by [`SetUp`](@ref) and not mutated
afterwards; the bank is shared read-only, so changing basis means building new
plans and workspaces.

`FFTPlans` is the legacy name of this type.

# Fields
- `backend`: the non-Limber backend (bank and sampling contract).
- `plan_ϕTT`: plan for the unequal-time power spectrum `P(k, χ₁, χ₂)`.
- `plan_ϕT`: Optional plan for the power spectrum `P(k, χ₁)` built with a single transfer function.
- `plan_ϕ`: Optional plan for the primordial power spectrum.
- `plan_limber`: FFT plan for the Limber power spectrum.
- `T_k_limber`: Precomputed k-basis polynomials for the Limber grid.
- `plan_ℓ`: Plan for the final C_ℓ interpolation.
"""
struct SpectralPlans{B<:AbstractNonLimberBackend,
                     PϕTT,                      # source plans: the backend's concrete plan type (or Nothing)
                     PϕT,
                     Pϕ,
                     PL<:ChebyshevPlan,
                     Pℓ<:ChebyshevPlan}
    backend::B
    plan_ϕTT::PϕTT
    plan_ϕT::PϕT
    plan_ϕ::Pϕ
    plan_limber::PL
    T_k_limber::Array{Float64, 3}
    plan_ℓ::Pℓ
    radial_plan::RadialQuadraturePlan
end
SpectralPlans(backend::AbstractNonLimberBackend,a,b,c,d,e,f)=
    SpectralPlans(backend,a,b,c,d,e,f,LEGACY_RADIAL_PLAN)

function _check_radial_bank(backend,plan)
    bank_signature(backend).geometry==plan.geometry || throw(ArgumentError(
        "tensor bank axes differ from the selected radial quadrature plan"))
    return nothing
end

"Legacy name of [`SpectralPlans`](@ref)."
const FFTPlans = SpectralPlans

"""
    SpectralPlans(plan_ϕTT, plan_ϕT, plan_ϕ, plan_limber, T_k_limber, plan_ℓ)

Legacy six-argument form: wraps the default Chebyshev backend, after checking
that the supplied source plans really are the Chebyshev plans of that backend
(same site count and bounds, correct array shape). Plans for any other
configuration must be built through `SetUp(...; backend=...)`.
"""
function SpectralPlans(plan_ϕTT::ChebyshevPlan, plan_ϕT, plan_ϕ, plan_limber::ChebyshevPlan,
                       T_k_limber::Array{Float64,3}, plan_ℓ::ChebyshevPlan)
    backend = ChebyshevBackend()
    nk = length(nonlimber_sites(backend))
    _is_source_plan(backend, plan_ϕTT, (nk, length(Blast.χ), length(Blast.R))) ||
        throw(ArgumentError("plan_ϕTT is not the legacy Chebyshev (k, χ, R) plan; use SetUp(...; backend=...)"))
    isnothing(plan_ϕT) || _is_source_plan(backend, plan_ϕT, (nk, length(Blast.χ), length(Blast.R))) ||
        throw(ArgumentError("plan_ϕT is not the legacy Chebyshev (k, χ, R) plan; use SetUp(...; backend=...)"))
    isnothing(plan_ϕ) || _is_source_plan(backend, plan_ϕ, (nk,)) ||
        throw(ArgumentError("plan_ϕ is not the legacy Chebyshev k plan; use SetUp(...; backend=...)"))
    return SpectralPlans(backend, plan_ϕTT, plan_ϕT, plan_ϕ, plan_limber, T_k_limber, plan_ℓ)
end

# Legacy keyword form of the original `@kwdef FFTPlans`.
function SpectralPlans(; plan_ϕTT, plan_ϕT=nothing, plan_ϕ=nothing, plan_limber, T_k_limber, plan_ℓ)
    return SpectralPlans(plan_ϕTT, plan_ϕT, plan_ϕ, plan_limber, T_k_limber, plan_ℓ)
end

bank_signature(P::SpectralPlans) = bank_signature(P.backend)
nonlimber_grid(P::SpectralPlans) = nonlimber_grid(P.backend)

function _setup_limber_plan(radial_plan::RadialQuadraturePlan=LEGACY_RADIAL_PLAN)
    lk_min, lk_max = minimum(k_limber), maximum(k_limber)
    z_min, z_max = minimum(z_cheb), maximum(z_cheb)
    K_k, K_z = length(k_limber) - 1, length(z_cheb) - 1

    plan_limber = prepare_chebyshev_plan((lk_min, z_min), (lk_max, z_max), (K_k, K_z))
    
    T_k = get_limber_k_polynomials(plan_limber, Blast.full_ℓ_range, radial_plan.chi; is_log_k=true)
    
    plan_ℓ = prepare_chebyshev_plan(2, 2000, 100)
    
    return plan_limber, T_k, plan_ℓ
end

"""
    SetUp(probes...; backend = ChebyshevBackend())

Construct the projected matter density containers and spectral plans required
for the given set of cosmological probes. The keyword `backend` selects the
non-Limber backend that owns the tensor bank and sampling sites (default: the
legacy Chebyshev bank, loaded once on first use and shared).

# Returns
- `W::ProjectedMatterDensity`: Container holding all required inner k integrals.
- `P::SpectralPlans`: Pre-allocated plans (legacy name `FFTPlans`).
"""
function SetUp(G::GalaxyClustering; backend::AbstractNonLimberBackend=ChebyshevBackend(), radial_plan::RadialQuadraturePlan=LEGACY_RADIAL_PLAN)
    bank = tensor_bank(backend)
    _check_radial_bank(backend,radial_plan)
    plan_ϕTT = source_plan_ϕTT(backend)
    plan_limber, T_k, plan_ℓ = _setup_limber_plan(radial_plan)

    # Each `have_*` folds to a compile-time constant once G's concrete type
    # params are known, so every ternary below is const-propagated away and
    # the returned `ProjectedMatterDensity{...}` has fully concrete type
    # parameters (never `Union{T, Nothing}` nor UnionAll).
    have_RSD = !isnothing(G.RSD)
    have_μ   = !isnothing(G.μ)
    have_PNG = !isnothing(G.PNG)

    plan_ϕT = have_PNG ? plan_ϕTT : nothing
    plan_ϕ  = have_PNG ? source_plan_ϕ(backend) : nothing

    W = ProjectedMatterDensity(
        w_2_00_ϕTT      = w_2_00_ϕTT(bank),  # δ (always)
        w_minus2_00_ϕTT = have_μ ? w_minus2_00_ϕTT(bank) : nothing,
        w_0_00_ϕTT      = have_μ ? w_0_00_ϕTT(bank) : nothing,
        w_0_02_ϕTT      = (have_RSD && have_μ) ? w_0_02_ϕTT(bank) : nothing,
        w_0_20_ϕTT      = (have_RSD && have_μ) ? w_0_20_ϕTT(bank) : nothing,
        w_2_02_ϕTT      = have_RSD ? w_2_02_ϕTT(bank) : nothing,
        w_2_20_ϕTT      = have_RSD ? w_2_20_ϕTT(bank) : nothing,
        w_2_22_ϕTT      = have_RSD ? w_2_22_ϕTT(bank) : nothing,
        w_2_00_ϕT       = have_PNG ? w_2_00_ϕT(bank) : nothing,
        w_2_00_ϕT_R1    = have_PNG ? w_2_00_ϕT_R1(bank) : nothing,
        w_0_00_ϕT       = (have_μ && have_PNG) ? w_0_00_ϕT(bank) : nothing,
        w_0_00_ϕT_R1    = (have_μ && have_PNG) ? w_0_00_ϕT_R1(bank) : nothing,
        w_2_02_ϕT       = (have_RSD && have_PNG) ? w_2_02_ϕT(bank) : nothing,
        w_2_02_ϕT_R1    = (have_RSD && have_PNG) ? w_2_02_ϕT_R1(bank) : nothing,
        w_2_20_ϕT       = (have_RSD && have_PNG) ? w_2_20_ϕT(bank) : nothing,
        w_2_20_ϕT_R1    = (have_RSD && have_PNG) ? w_2_20_ϕT_R1(bank) : nothing,
        w_2_00_ϕ        = have_PNG ? w_2_00_ϕ(bank) : nothing,
    )

    Plans = SpectralPlans(backend, plan_ϕTT, plan_ϕT, plan_ϕ, plan_limber, T_k, plan_ℓ, radial_plan)
    return W, Plans
end

function SetUp(L::WeakLensing; backend::AbstractNonLimberBackend=ChebyshevBackend(), radial_plan::RadialQuadraturePlan=LEGACY_RADIAL_PLAN)
    bank = tensor_bank(backend)
    _check_radial_bank(backend,radial_plan)
    plan_ϕTT = source_plan_ϕTT(backend)
    plan_limber, T_k, plan_ℓ = _setup_limber_plan(radial_plan)

    W = ProjectedMatterDensity(
        w_minus2_00_ϕTT = w_minus2_00_ϕTT(bank),  # γ (always under L)
    )

    Plans = SpectralPlans(backend, plan_ϕTT, nothing, nothing, plan_limber, T_k, plan_ℓ, radial_plan)
    return W, Plans
end

function SetUp(G::GalaxyClustering, L::WeakLensing; backend::AbstractNonLimberBackend=ChebyshevBackend(), radial_plan::RadialQuadraturePlan=LEGACY_RADIAL_PLAN)
    bank = tensor_bank(backend)
    _check_radial_bank(backend,radial_plan)
    plan_ϕTT = source_plan_ϕTT(backend)
    plan_limber, T_k, plan_ℓ = _setup_limber_plan(radial_plan)

    have_RSD = !isnothing(G.RSD)
    have_μ   = !isnothing(G.μ)
    have_PNG = !isnothing(G.PNG)
    have_γ   = !isnothing(L.γ)  # always true — L.γ is required

    plan_ϕT = have_PNG ? plan_ϕTT : nothing
    plan_ϕ  = have_PNG ? source_plan_ϕ(backend) : nothing

    W = ProjectedMatterDensity(
        w_2_00_ϕTT      = w_2_00_ϕTT(bank),            # δ
        w_minus2_00_ϕTT = w_minus2_00_ϕTT(bank),       # γ
        w_0_00_ϕTT      = w_0_00_ϕTT(bank),            # δ×γ cross
        w_0_02_ϕTT      = (have_RSD && (have_μ || have_γ)) ? w_0_02_ϕTT(bank) : nothing,
        w_0_20_ϕTT      = (have_RSD && (have_μ || have_γ)) ? w_0_20_ϕTT(bank) : nothing,
        w_2_02_ϕTT      = have_RSD ? w_2_02_ϕTT(bank) : nothing,
        w_2_20_ϕTT      = have_RSD ? w_2_20_ϕTT(bank) : nothing,
        w_2_22_ϕTT      = have_RSD ? w_2_22_ϕTT(bank) : nothing,
        w_2_00_ϕT       = have_PNG ? w_2_00_ϕT(bank) : nothing,
        w_2_00_ϕT_R1    = have_PNG ? w_2_00_ϕT_R1(bank) : nothing,
        w_0_00_ϕT       = ((have_μ || have_γ) && have_PNG) ? w_0_00_ϕT(bank) : nothing,
        w_0_00_ϕT_R1    = ((have_μ || have_γ) && have_PNG) ? w_0_00_ϕT_R1(bank) : nothing,
        w_2_02_ϕT       = (have_RSD && have_PNG) ? w_2_02_ϕT(bank) : nothing,
        w_2_02_ϕT_R1    = (have_RSD && have_PNG) ? w_2_02_ϕT_R1(bank) : nothing,
        w_2_20_ϕT       = (have_RSD && have_PNG) ? w_2_20_ϕT(bank) : nothing,
        w_2_20_ϕT_R1    = (have_RSD && have_PNG) ? w_2_20_ϕT_R1(bank) : nothing,
        w_2_00_ϕ        = have_PNG ? w_2_00_ϕ(bank) : nothing,
    )

    Plans = SpectralPlans(backend, plan_ϕTT, plan_ϕT, plan_ϕ, plan_limber, T_k, plan_ℓ, radial_plan)
    return W, Plans
end

function SetUp(L::WeakLensing, G::GalaxyClustering; backend::AbstractNonLimberBackend=ChebyshevBackend(), radial_plan::RadialQuadraturePlan=LEGACY_RADIAL_PLAN)
    return SetUp(G, L; backend=backend,radial_plan)
end

function SetUp(G::GalaxyClustering, C::CMB; backend::AbstractNonLimberBackend=ChebyshevBackend(), radial_plan::RadialQuadraturePlan=LEGACY_RADIAL_PLAN)
    bank = tensor_bank(backend)
    _check_radial_bank(backend,radial_plan)
    plan_ϕTT = source_plan_ϕTT(backend)
    plan_limber, T_k, plan_ℓ = _setup_limber_plan(radial_plan)

    have_RSD = !isnothing(G.RSD)
    have_μ   = !isnothing(G.μ)
    have_PNG = !isnothing(G.PNG)
    have_κ   = !isnothing(C.κ)  # always true — C.κ is required

    plan_ϕT = have_PNG ? plan_ϕTT : nothing
    plan_ϕ  = have_PNG ? source_plan_ϕ(backend) : nothing

    W = ProjectedMatterDensity(
        w_2_00_ϕTT      = w_2_00_ϕTT(bank),     # δ
        w_minus2_00_ϕTT = have_μ ? w_minus2_00_ϕTT(bank) : nothing,
        w_0_00_ϕTT      = w_0_00_ϕTT(bank),     # κ×δ cross (always under G×C)
        w_0_02_ϕTT      = (have_RSD && (have_μ || have_κ)) ? w_0_02_ϕTT(bank) : nothing,
        w_0_20_ϕTT      = (have_RSD && (have_μ || have_κ)) ? w_0_20_ϕTT(bank) : nothing,
        w_2_02_ϕTT      = have_RSD ? w_2_02_ϕTT(bank) : nothing,
        w_2_20_ϕTT      = have_RSD ? w_2_20_ϕTT(bank) : nothing,
        w_2_22_ϕTT      = have_RSD ? w_2_22_ϕTT(bank) : nothing,
        w_2_00_ϕT       = have_PNG ? w_2_00_ϕT(bank) : nothing,
        w_2_00_ϕT_R1    = have_PNG ? w_2_00_ϕT_R1(bank) : nothing,
        w_0_00_ϕT       = ((have_μ || have_κ) && have_PNG) ? w_0_00_ϕT(bank) : nothing,
        w_0_00_ϕT_R1    = ((have_μ || have_κ) && have_PNG) ? w_0_00_ϕT_R1(bank) : nothing,
        w_2_02_ϕT       = (have_RSD && have_PNG) ? w_2_02_ϕT(bank) : nothing,
        w_2_02_ϕT_R1    = (have_RSD && have_PNG) ? w_2_02_ϕT_R1(bank) : nothing,
        w_2_20_ϕT       = (have_RSD && have_PNG) ? w_2_20_ϕT(bank) : nothing,
        w_2_20_ϕT_R1    = (have_RSD && have_PNG) ? w_2_20_ϕT_R1(bank) : nothing,
        w_2_00_ϕ        = have_PNG ? w_2_00_ϕ(bank) : nothing,
    )

    Plans = SpectralPlans(backend, plan_ϕTT, plan_ϕT, plan_ϕ, plan_limber, T_k, plan_ℓ, radial_plan)
    return W, Plans
end

function SetUp(C::CMB, G::GalaxyClustering; backend::AbstractNonLimberBackend=ChebyshevBackend(), radial_plan::RadialQuadraturePlan=LEGACY_RADIAL_PLAN)
    return SetUp(G, C; backend=backend,radial_plan)
end

function SetUp(L::WeakLensing, C::CMB; backend::AbstractNonLimberBackend=ChebyshevBackend(), radial_plan::RadialQuadraturePlan=LEGACY_RADIAL_PLAN)
    bank = tensor_bank(backend)
    _check_radial_bank(backend,radial_plan)
    plan_ϕTT = source_plan_ϕTT(backend)
    plan_limber, T_k, plan_ℓ = _setup_limber_plan(radial_plan)

    W = ProjectedMatterDensity(
        w_minus2_00_ϕTT = w_minus2_00_ϕTT(bank),  # γ (always); κ×γ cross uses this
    )

    Plans = SpectralPlans(backend, plan_ϕTT, nothing, nothing, plan_limber, T_k, plan_ℓ, radial_plan)
    return W, Plans
end

function SetUp(C::CMB, L::WeakLensing; backend::AbstractNonLimberBackend=ChebyshevBackend(), radial_plan::RadialQuadraturePlan=LEGACY_RADIAL_PLAN)
    return SetUp(L, C; backend=backend,radial_plan)
end

function SetUp(G::GalaxyClustering, L::WeakLensing, C::CMB; backend::AbstractNonLimberBackend=ChebyshevBackend(), radial_plan::RadialQuadraturePlan=LEGACY_RADIAL_PLAN)
    bank = tensor_bank(backend)
    _check_radial_bank(backend,radial_plan)
    plan_ϕTT = source_plan_ϕTT(backend)
    plan_limber, T_k, plan_ℓ = _setup_limber_plan(radial_plan)

    have_RSD = !isnothing(G.RSD)
    have_μ   = !isnothing(G.μ)
    have_PNG = !isnothing(G.PNG)
    have_γ   = !isnothing(L.γ)  # always true — L.γ required
    have_κ   = !isnothing(C.κ)  # always true — C.κ required

    plan_ϕT = have_PNG ? plan_ϕTT : nothing
    plan_ϕ  = have_PNG ? source_plan_ϕ(backend) : nothing

    W = ProjectedMatterDensity(
        w_2_00_ϕTT      = w_2_00_ϕTT(bank),             # δ
        w_minus2_00_ϕTT = w_minus2_00_ϕTT(bank),        # γ×γ, γ×κ
        w_0_00_ϕTT      = w_0_00_ϕTT(bank),             # δ×γ, δ×IA, δ×κ
        w_0_02_ϕTT      = (have_RSD && (have_μ || have_γ)) ? w_0_02_ϕTT(bank) : nothing,
        w_0_20_ϕTT      = (have_RSD && (have_μ || have_γ)) ? w_0_20_ϕTT(bank) : nothing,
        w_2_02_ϕTT      = have_RSD ? w_2_02_ϕTT(bank) : nothing,
        w_2_20_ϕTT      = have_RSD ? w_2_20_ϕTT(bank) : nothing,
        w_2_22_ϕTT      = have_RSD ? w_2_22_ϕTT(bank) : nothing,
        w_2_00_ϕT       = have_PNG ? w_2_00_ϕT(bank) : nothing,
        w_2_00_ϕT_R1    = have_PNG ? w_2_00_ϕT_R1(bank) : nothing,
        w_0_00_ϕT       = ((have_μ || have_γ) && have_PNG) ? w_0_00_ϕT(bank) : nothing,
        w_0_00_ϕT_R1    = ((have_μ || have_γ) && have_PNG) ? w_0_00_ϕT_R1(bank) : nothing,
        w_2_02_ϕT       = (have_RSD && have_PNG) ? w_2_02_ϕT(bank) : nothing,
        w_2_02_ϕT_R1    = (have_RSD && have_PNG) ? w_2_02_ϕT_R1(bank) : nothing,
        w_2_20_ϕT       = (have_RSD && have_PNG) ? w_2_20_ϕT(bank) : nothing,
        w_2_20_ϕT_R1    = (have_RSD && have_PNG) ? w_2_20_ϕT_R1(bank) : nothing,
        w_2_00_ϕ        = have_PNG ? w_2_00_ϕ(bank) : nothing,
    )

    Plans = SpectralPlans(backend, plan_ϕTT, plan_ϕT, plan_ϕ, plan_limber, T_k, plan_ℓ, radial_plan)
    return W, Plans
end

function SetUp(G::GalaxyClustering, C::CMB, L::WeakLensing; backend::AbstractNonLimberBackend=ChebyshevBackend(), radial_plan::RadialQuadraturePlan=LEGACY_RADIAL_PLAN)
    return SetUp(G, L, C; backend=backend,radial_plan)
end

function SetUp(L::WeakLensing, G::GalaxyClustering, C::CMB; backend::AbstractNonLimberBackend=ChebyshevBackend(), radial_plan::RadialQuadraturePlan=LEGACY_RADIAL_PLAN)
    return SetUp(G, L, C; backend=backend,radial_plan)
end

function SetUp(L::WeakLensing, C::CMB, G::GalaxyClustering; backend::AbstractNonLimberBackend=ChebyshevBackend(), radial_plan::RadialQuadraturePlan=LEGACY_RADIAL_PLAN)
    return SetUp(G, L, C; backend=backend,radial_plan)
end

function SetUp(C::CMB, G::GalaxyClustering, L::WeakLensing; backend::AbstractNonLimberBackend=ChebyshevBackend(), radial_plan::RadialQuadraturePlan=LEGACY_RADIAL_PLAN)
    return SetUp(G, L, C; backend=backend,radial_plan)
end

function SetUp(C::CMB, L::WeakLensing, G::GalaxyClustering; backend::AbstractNonLimberBackend=ChebyshevBackend(), radial_plan::RadialQuadraturePlan=LEGACY_RADIAL_PLAN)
    return SetUp(G, L, C; backend=backend,radial_plan)
end

"""
    PowerSpectrum

A mutable struct acting as a container for all processed spectral data required for the angular power spectrum measurements.

This struct caches both the Chebyshev coefficients used in the high-precision 
non-Limber integration and the pre-interpolated power spectra used for 
the Limber approximation.

# Fields
- `cϕTT`, `cϕT`, `cϕ`: Structs storing all the necessary Chebyshev coefficients.
- `ΔP_limber`: A 2D array of the non-linear correction to the matter power spectrum.  
- `Pδ_limber`: A 2D array of the total matter power spectrum.
- `signature`: [`BankSignature`](@ref) of the plans that produced the coefficients;
  checked against the workspace and plans at public boundaries.
"""
mutable struct PowerSpectrum{T<:Real}
    cϕTT::cϕTT{T}
    cϕT::Union{cϕT{T}, Nothing}
    cϕ::Union{cϕ{T}, Nothing}
    ΔP_limber::Matrix{T}
    Pδ_limber::Matrix{T}
    signature::BankSignature
    radial_plan::RadialQuadraturePlan
end

# Outer constructor: ONE common scalar type T promoted over every ACTIVE coefficient family and both
# Limber matrices. Fixed Float64 background with Dual Pmm (primordial coefficients stay Float64), or
# Limber-only Dual inputs with constant non-Limber coefficients, are valid public cases: constant
# families are promoted with `convert` (value preserved, zero partials), differentiated families are
# never stripped. Inactive (`nothing`) families stay `nothing`; a family/matrix already of type T is
# kept as the SAME object, so Float64 and reverse-mode paths do no extra work.
_coef_eltype(::Nothing) = Union{}
_coef_eltype(c::AbstractCoeff) = eltype(c.coefs)
_promote_coeff(::Nothing, ::Type) = nothing
# Three known containers, three typed methods (no type-name reflection); identity kept when already T.
_promote_coeff(c::cϕTT, ::Type{T}) where {T} = eltype(c.coefs) === T ? c : cϕTT(convert(Array{T, 3}, c.coefs))
_promote_coeff(c::cϕT,  ::Type{T}) where {T} = eltype(c.coefs) === T ? c : cϕT(convert(Array{T, 3}, c.coefs))
_promote_coeff(c::cϕ,   ::Type{T}) where {T} = eltype(c.coefs) === T ? c : cϕ(convert(Array{T, 1}, c.coefs))

function PowerSpectrum(c_TT::cϕTT,
                       c_T::Union{cϕT, Nothing},
                       c_ϕ::Union{cϕ, Nothing},
                       ΔP_limber::AbstractMatrix,
                       Pδ_limber::AbstractMatrix,
                       signature::BankSignature; radial_plan::RadialQuadraturePlan=LEGACY_RADIAL_PLAN)
    signature.geometry==radial_plan.geometry || throw(ArgumentError(
        "PowerSpectrum bank axes differ from its radial quadrature plan"))
    T = promote_type(_coef_eltype(c_TT), _coef_eltype(c_T), _coef_eltype(c_ϕ),
                     eltype(ΔP_limber), eltype(Pδ_limber))
    T <: Real || throw(ArgumentError("PowerSpectrum needs real scalar inputs; got $T"))
    return PowerSpectrum{T}(_promote_coeff(c_TT, T), _promote_coeff(c_T, T), _promote_coeff(c_ϕ, T),
                            convert(Matrix{T}, ΔP_limber), convert(Matrix{T}, Pδ_limber), signature, radial_plan)
end

# Legacy five-argument form: only valid for coefficients on the legacy
# Chebyshev sites; it never tags other shapes with the legacy signature.
function PowerSpectrum(cϕTT::cϕTT, cϕT, cϕ, ΔP_limber::AbstractMatrix, Pδ_limber::AbstractMatrix)
    size(cϕTT.coefs, 1) == length(k_cheb) || throw(ArgumentError(
        "legacy PowerSpectrum construction needs $(length(k_cheb)) Chebyshev coefficients; pass a BankSignature"))
    return PowerSpectrum(cϕTT, cϕT, cϕ, ΔP_limber, Pδ_limber, LEGACY_BANK_SIGNATURE)
end

"""
    get_PΦ(k::AbstractArray, cosmo::AbstractCosmology)

Computes the primordial power spectrum \$P_\\Phi(k)\$.

The implementation follows the standard inflationary prediction:
```math
P_\\Phi(k) = \\frac{9}{25} \\frac{2\\pi^2 A_s}{k^3} \\left( \\frac{k}{k_{pivot}}
```
where the pivot scale is fixed at `0.05 Mpc^{-1}``.

# Returns
- A 1D array of the primordial potential power spectrum at the requested scales. 
"""
function get_PΦ(k::AbstractArray{<:Any,1}, cosmo::AbstractCosmology)
    As = get_As(cosmo)
    ns = get_ns(cosmo)
    return @. 9/25 * 2 * π^2 * As / (k^3) * (k/0.05)^(ns - 1.)
end

# Background-dispatched overloads keep the primordial scalar type aligned with
# the prepared background. ACE 0.11 already parameterizes cosmology fields
# concretely; this dispatch remains useful when T is a ForwardDiff scalar.
function get_PΦ(k::AbstractArray{<:Any,1}, bg::Background{T}) where {T}
    As = get_As(bg)
    ns = get_ns(bg)
    return @. 9/25 * 2 * π^2 * As / (k^3) * (k/0.05)^(ns - 1.)
end

""" get_Tm(pk::AbstractArray{<:Any, 2}, k::AbstractArray, cosmo::AbstractCosmology)

Extracts the matter transfer function T_m(k,z) from a matter power spectrum.

The transfer function is defined as the square root of the ratio between the processed matter power spectrum and the primordial potential power spectrum:
```math
T_m(k, z) = \\sqrt{\\frac{P(k, z)}{P_\\Phi(k)}}
```
"""
function get_Tm(pk::AbstractArray{<:Any,2}, k::AbstractArray{<:Any, 1}, cosmo::AbstractCosmology)
    prim_pk = get_PΦ(k , cosmo)
    return sqrt.(pk ./ reshape(prim_pk, 1, :))
end

# Background-dispatched overload — see note on get_PΦ above.
function get_Tm(pk::AbstractArray{<:Any,2}, k::AbstractArray{<:Any, 1},
                bg::Background)
    prim_pk = get_PΦ(k, bg)
    return sqrt.(pk ./ reshape(prim_pk, 1, :))
end

"""
    transform_to_R_frame(matrix::AbstractArray{<:Any,2}, bg::Background)

Transforms a 2D array from standard coordinates `(z, k)` to a ratio-based 
coordinate system `(χ, Rχ, k)` using pre-built Background interpolators.
"""
function transform_to_R_frame(matrix::AbstractArray{<:Any,2}, bg::Background)
    new_χs = isnothing(bg.radial_plan) ? make_grid(bg.χ,Blast.R) : bg.radial_plan.chiR_flat
    # Pass the whole χ vector to the z_of_χ akima closure in ONE call rather
    # than scalar-broadcasting (`bg.z_of_χ.(new_χs)`). The dotted form
    # rebuilt the akima coefficient arrays (b, c, d of length ~1000) on every
    # single one of the 6144 scalar queries, costing ~100 ms and ~1 GiB per
    # Background evaluation. The vector form shares a single setup across
    # queries via `_akima_eval(..., tq::AbstractArray)` — ~1 ms, ~200 KiB.
    x = bg.z_of_χ(new_χs)
    interp = akima_interpolation(matrix, Blast.z_lin, x)
    nR = isnothing(bg.radial_plan) ? length(Blast.R) : length(bg.radial_plan.R)
    return reshape(interp, length(bg.χ), nR, size(interp, 2))
end

#TODO: this function will interface with mapse once ready.
"""
    prepare_pk_workspace(P::SpectralPlans,
                         pk::AbstractArray{<:Any, 2}, pk_limber_lin::AbstractArray{<:Any, 2}, 
                         pk_limber_nonlin::AbstractArray{<:Any, 2}, 
                         bg::Background)

Prepare total-matter `Pmm` for the complete Blast pipeline. `pk` must have
shape `(length(Blast.z_lin), length(nonlimber_sites(P.backend)))` (see [`nonlimber_grid`](@ref)); both Limber arrays must
have shape `(length(Blast.z_cheb), length(Blast.k_limber))`. All values must be
finite and strictly positive. Background and `P` must carry the same fixed
radial plan, and the backend's tensor bank must match that plan's axes exactly.
"""
function prepare_pk_workspace(P::SpectralPlans, pk::AbstractArray{<:Any, 2}, pk_limber_lin::AbstractArray{<:Any, 2}, pk_limber_nonlin::AbstractArray{<:Any, 2}, bg::Background;
                              pk_grid = nothing)
    _require_radial_background(bg,P.radial_plan)
    if !isnothing(pk_grid)    # exact coordinates supplied by the caller: must be the required sampling, never resampled
        required = nonlimber_grid(P)
        (pk_grid.z == required.z && pk_grid.k == required.k) ||
            throw(ArgumentError("pk_grid does not match the plans' non-Limber sampling grid (nonlimber_grid(plans)); Pmm is never resampled"))
    end
    k_sites = nonlimber_sites(P.backend)
    k_phys = nonlimber_k(P.backend)
    _validate_pmm("pk", pk, (length(Blast.z_lin), length(k_sites)))
    limber_shape = (length(Blast.z_cheb), length(Blast.k_limber))
    _validate_pmm("pk_limber_lin", pk_limber_lin, limber_shape)
    _validate_pmm("pk_limber_nonlin", pk_limber_nonlin, limber_shape)
    #Treating the non-limber power spectrum
    P_ϕ = get_PΦ(k_phys, bg)
    transfer_func = get_Tm(pk, k_phys, bg)
    
    transfer_func_χR = transform_to_R_frame(transfer_func, bg)
    transfer_func_χ1 = transfer_func_χR[:,end,:]

    P_ϕTT = _p_phi_TT_tullio(P_ϕ, transfer_func_χ1, transfer_func_χR)
    P_ϕT  = _p_phi_T_tullio(P_ϕ, transfer_func_χR)

    # Ratio backends divide the COMPLETE sources by the fixed q here; identity (no copy) otherwise.
    P_ϕTT, P_ϕT, P_ϕ = _normalize_sources(P.backend, P_ϕTT, P_ϕT, P_ϕ)

    c1 = build_coeff(cϕTT, P_ϕTT, P.plan_ϕTT)
    c2 = build_coeff(cϕT,  P_ϕT,  P.plan_ϕT)   
    c3 = build_coeff(cϕ,   P_ϕ,   P.plan_ϕ)    


    P_ϕ_limber = get_PΦ(10 .^ Blast.k_limber, bg)'
    
    T_m_limber_lin = get_Tm(pk_limber_lin, 10 .^ Blast.k_limber, bg)
    P_limber_lin = P_ϕ_limber .* T_m_limber_lin .^ 2.
    # permutedims materializes the transpose into a plain Matrix, avoiding
    # an Adjoint wrapper across the chebyshev_decomposition rrule boundary
    # (Mooncake's upstream rule in AbstractCosmologicalEmulators doesn't
    # support Adjoint tangents; plain Matrix works).
    c_lin = chebyshev_decomposition(P.plan_limber, permutedims(log10.(P_limber_lin)))

    T_m_limber_nonlin = get_Tm(pk_limber_nonlin, 10 .^ Blast.k_limber, bg)
    P_limber_nonlin = P_ϕ_limber .* T_m_limber_nonlin .^ 2.
    c_nonlin = chebyshev_decomposition(P.plan_limber, permutedims(log10.(P_limber_nonlin)))

    T_z_limber = get_limber_coords_polynomials(P.plan_limber, bg.z)

    P_lin_grid = 10.0 .^ limber_eval(c_lin, T_z_limber, P.T_k_limber)
    P_nonlin_grid = 10.0 .^ limber_eval(c_nonlin, T_z_limber, P.T_k_limber)

    ΔP_limber = P_nonlin_grid .- P_lin_grid
    Pδ_limber = P_nonlin_grid

    return PowerSpectrum(c1, c2, c3, ΔP_limber, Pδ_limber, bank_signature(P);radial_plan=P.radial_plan)
end

function _validate_pmm(name, pk::AbstractMatrix, expected_size)
    size(pk) == expected_size || throw(DimensionMismatch(
        "$name must have size $expected_size on Blast's production grids; got $(size(pk))"))
    all(isfinite, pk) || throw(ArgumentError("$name must contain only finite total-matter Pmm values"))
    all(x -> x > zero(x), pk) || throw(DomainError(
        minimum(pk), "$name must contain strictly positive total-matter Pmm values before square roots/logs"))
    return nothing
end

# ----------------------------------------------------------------------------
# prepare_pk_workspace tensor products.  Both are pure outer-products (no
# contraction), so these are really broadcasts dressed as @tullio.  We extract
# them as named helpers and register rrules + Mooncake primitives so the AD
# contract is explicit and stable under future Tullio/Mooncake changes.
#
#   P_ϕTT[k,i,j] = P_ϕ[k] · T_χ1[i,k] · T_χR[i,j,k]
#   P_ϕT[k,i,j]  = P_ϕ[k]              · T_χR[i,j,k]
# ----------------------------------------------------------------------------
function _p_phi_TT_tullio(P_ϕ::AbstractVector, T_χ1::AbstractMatrix,
                          T_χR::AbstractArray{<:Any, 3})
    @tullio out[k, i, j] := P_ϕ[k] * T_χ1[i, k] * T_χR[i, j, k]
    return out
end

function _p_phi_T_tullio(P_ϕ::AbstractVector, T_χR::AbstractArray{<:Any, 3})
    @tullio out[k, i, j] := P_ϕ[k] * T_χR[i, j, k]
    return out
end

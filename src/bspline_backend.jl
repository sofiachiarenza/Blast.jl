# Direct cubic B-spline non-Limber backend (phase 5). The bank is caller-owned; nothing is loaded here.
# Direct representation only: the same eight shared geometry tensors and eight-scan contractions as the
# Chebyshev path; source values are sampled at the bank's physical sites and turned into coefficients by
# a collocation SOLVE (public ACE `bspline_coefficients`), never read off as sample values.

abstract type AbstractSourcePolicy end
"Direct expansion: unweighted shared geometry tensors; sources are expanded as they are."
struct DirectExpansion <: AbstractSourcePolicy end
"Fiducial ratio: complete k-first sources are divided by the bank's FIXED q before the collocation solve."
struct FiducialRatio <: AbstractSourcePolicy end

_check_policy(::DirectSplineBank, ::DirectExpansion) = nothing
_check_policy(::RatioSplineBank, ::FiducialRatio) = nothing
_check_policy(bank, policy) = throw(ArgumentError(
    "policy $(typeof(policy)) is incompatible with a $(typeof(bank).name.name); use DirectExpansion for direct and FiducialRatio for ratio banks"))

_bank_signature(b::DirectSplineBank) = b.tensors.signature
_bank_signature(b::RatioSplineBank) = b.signature
_family_arrays(b::DirectSplineBank) = Tuple(getfield(b.tensors, g) for g in _BANK_GEOMETRIES)
_family_arrays(b::RatioSplineBank) = (values(b.TT)..., values(b.T)..., b.PHI)

"""
    CubicBSplineBackend(bank; policy = <matching policy>)

Explicit spline backend owning a validated caller-owned bank: a [`DirectSplineBank`](@ref)
(`DirectExpansion`) or a [`RatioSplineBank`](@ref) (`FiducialRatio`); an incompatible pair is rejected.
Pass it to `SetUp(...; backend = CubicBSplineBackend(bank))`. Sampling contract: the non-Limber `Pmm`
must be supplied on `nonlimber_grid(plans)`, whose `k` is the bank's stored physical sites (direct and
ratio banks of the same count use different sites). Source plans are created per `SetUp` call, so they are
private to that setup (one chain). For ratio banks the complete k-first sources are divided by the bank's
fixed `q` immediately before the collocation solve; `Pmm` inputs are never normalized and `q` is never
recomputed from the current cosmology.
"""
struct CubicBSplineBackend{B<:Union{DirectSplineBank,RatioSplineBank}, P<:AbstractSourcePolicy} <: AbstractNonLimberBackend
    bank::B
    policy::P
    function CubicBSplineBackend{B,P}(bank::B, policy::P) where {B<:Union{DirectSplineBank,RatioSplineBank}, P<:AbstractSourcePolicy}
        _check_policy(bank, policy)
        g, sig = bank.grid, _bank_signature(bank)
        sig.basis === :bspline_cubic_notaknot ||
            throw(ArgumentError("CubicBSplineBackend needs a bank with basis :bspline_cubic_notaknot; got $(sig.basis)"))
        n = length(g.logsites)
        (n ≥ 4 && length(g.sites) == n && length(g.knots) == n + 4 && sig.n_k == n) ||
            throw(ArgumentError("inconsistent spline bank: site/knot counts differ from the signature"))
        bank isa DirectSplineBank && bank.tensors.sites != g.logsites &&
            throw(ArgumentError("inconsistent spline bank: log₁₀ sites differ from the tensor bank"))
        all(a -> size(a, 4) == n, _family_arrays(bank)) ||
            throw(ArgumentError("inconsistent spline bank: a tensor's coefficient axis differs from the site count"))
        if bank isa RatioSplineBank
            all(v -> length(v) == n && all(x -> isfinite(x) && x > 0, v), (bank.q.qTT, bank.q.qT, bank.q.qPHI)) ||
                throw(ArgumentError("ratio bank normalizers must be finite and positive at every site"))
        end
        (all(isfinite, g.logsites) && issorted(g.logsites; lt = <=)) || throw(ArgumentError("spline sites must be finite and strictly increasing"))
        all(abs.(10 .^ g.logsites ./ g.sites .- 1) .≤ 4eps()) ||
            throw(ArgumentError("physical and log₁₀ sites disagree beyond roundoff"))
        knot_vector(bspline_basis(CubicBSplinePlan(g.logsites, g.logsites; extrapolation = :throw))) == g.knots ||
            throw(ArgumentError("bank knots differ from the not-a-knot basis of the runtime ACE"))
        return new{B,P}(bank, policy)
    end
end
CubicBSplineBackend(bank::DirectSplineBank; policy::AbstractSourcePolicy = DirectExpansion()) =
    CubicBSplineBackend{typeof(bank), typeof(policy)}(bank, policy)
CubicBSplineBackend(bank::RatioSplineBank; policy::AbstractSourcePolicy = FiducialRatio()) =
    CubicBSplineBackend{typeof(bank), typeof(policy)}(bank, policy)

# `tensor_bank` is what `SetUp` hands to the generated component constructors: a TensorBank (direct,
# all shared geometries) or a RatioSplineBank (source-aware family lookup, no dummy banks).
tensor_bank(b::CubicBSplineBackend{<:DirectSplineBank}) = b.bank.tensors
tensor_bank(b::CubicBSplineBackend{<:RatioSplineBank}) = b.bank
bank_signature(b::CubicBSplineBackend) = _bank_signature(b.bank)
nonlimber_k(b::CubicBSplineBackend) = b.bank.grid.sites          # physical sites, exact stored bits
nonlimber_sites(b::CubicBSplineBackend) = b.bank.grid.logsites   # fixed coordinate of the basis: not differentiated

# One plan per source-plan field, created per setup. Fixed log₁₀ sites double as the evaluation sites; the
# coefficient solve acts on the first axis, so the same plan type serves (k,χ,R) and k-only sources.
source_plan_ϕTT(b::CubicBSplineBackend) = (x = b.bank.grid.logsites; CubicBSplinePlan(x, x; extrapolation = :throw))
source_plan_ϕ(b::CubicBSplineBackend) = (x = b.bank.grid.logsites; CubicBSplinePlan(x, x; extrapolation = :throw))

# Complete k-first sources before the collocation solve (hook used by `prepare_pk_workspace`):
# identity (no copy) for Chebyshev and direct splines; ratio divides by the fixed q along the k axis.
_normalize_sources(::AbstractNonLimberBackend, P_ϕTT, P_ϕT, P_ϕ) = (P_ϕTT, P_ϕT, P_ϕ)
_normalize_sources(b::CubicBSplineBackend{<:RatioSplineBank}, P_ϕTT, P_ϕT, P_ϕ) =
    (P_ϕTT ./ reshape(b.bank.q.qTT, :, 1, 1), P_ϕT ./ reshape(b.bank.q.qT, :, 1, 1), P_ϕ ./ b.bank.q.qPHI)

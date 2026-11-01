# Non-Limber spectral backends: own a validated fixed bank and the sampling
# contract that goes with it.

"""
    AbstractNonLimberBackend

A concrete backend owns a validated [`TensorBank`](@ref) and defines the sites at
which the non-Limber power spectrum must be supplied. There is deliberately no
registry and no process-global "active backend": pass a backend to `SetUp`.
"""
abstract type AbstractNonLimberBackend end

# The source plans are rebuilt from the endpoints and count alone, so the
# declared sites must be exactly those nodes, in order. Only roundoff between
# native constructions is tolerated (a few ulp), never a physical perturbation.
function _is_canonical_chebyshev(sites::AbstractVector{Float64})
    n = length(sites)
    all(isfinite, sites) || return false
    minimum(sites) < maximum(sites) || return false
    canonical = chebpoints(n - 1, minimum(sites), maximum(sites))
    return maximum(abs.(sites .- canonical)) ≤ 8 * eps(maximum(abs, sites))
end

"""
    ChebyshevBackend(bank = legacy_tensor_bank())

Direct Chebyshev backend: Chebyshev-extrema sites in log₁₀ k and the matching
`T̃` bank. The default is the legacy 161-site production bank.
"""
struct ChebyshevBackend{B<:TensorBank} <: AbstractNonLimberBackend
    bank::B
    function ChebyshevBackend{B}(bank::B) where {B<:TensorBank}
        bank.signature.basis === :chebyshev ||
            throw(ArgumentError("ChebyshevBackend needs a bank with basis :chebyshev; got $(bank.signature.basis)"))
        n = length(bank.sites)
        n ≥ 2 || throw(ArgumentError("a Chebyshev bank needs at least two sites"))
        _is_canonical_chebyshev(bank.sites) ||
            throw(ArgumentError("bank sites are not the canonical Chebyshev extrema in log₁₀ k"))
        return new{B}(bank)
    end
end
ChebyshevBackend(bank::B) where {B<:TensorBank} = ChebyshevBackend{B}(bank)
ChebyshevBackend() = ChebyshevBackend(legacy_tensor_bank())

tensor_bank(b::ChebyshevBackend) = b.bank
bank_signature(b::AbstractNonLimberBackend) = tensor_bank(b).signature

"log₁₀ k sites at which the non-Limber power spectrum is sampled (bank order)."
nonlimber_sites(b::AbstractNonLimberBackend) = tensor_bank(b).sites

"Physical k [1/Mpc] sampling sites. Chebyshev: `10 .^ sites` (unchanged expression); splines: the stored physical sites."
nonlimber_k(b::ChebyshevBackend) = 10 .^ nonlimber_sites(b)

"""
    nonlimber_grid(plans)

Physical grid on which non-Limber `Pmm` must be supplied to
[`prepare_pk_workspace`](@ref): `(z, k)` with `k` in Mpc⁻¹ in the exact order of
the backend's sites; `Pmm` has shape `(length(z), length(k))`. A plain array carries no
coordinates, so only its shape can be checked unless the grid is passed back via
`prepare_pk_workspace(...; pk_grid = nonlimber_grid(plans))`, which is then compared exactly.
"""
nonlimber_grid(b::AbstractNonLimberBackend) = (z = Blast.z_lin, k = nonlimber_k(b))

# Source plans: the ϕTT/ϕT plans cover (k, χ, R), the ϕ plan covers k only.
function source_plan_ϕTT(b::ChebyshevBackend)
    s = nonlimber_sites(b)
    return prepare_chebyshev_plan(minimum(s), maximum(s), length(s) - 1;
                                  size_nd=(length(s), length(bank_axes(b.bank).chi), length(bank_axes(b.bank).R)), dim=1)
end
function source_plan_ϕ(b::ChebyshevBackend)
    s = nonlimber_sites(b)
    return prepare_chebyshev_plan(minimum(s), maximum(s), length(s) - 1;
                                  size_nd=(length(s),), dim=1)
end

# True iff `plan` is exactly this backend's source plan: transformed axis 1,
# degree and bounds from the sites, and FFT input of the expected array shape
# (the transformed-degree tuple `K` has one entry per *transformed* axis, not
# one per array dimension).
function _is_source_plan(b::ChebyshevBackend, plan, size_nd::Tuple)
    plan isa ChebyshevPlan || return false
    s = nonlimber_sites(b)
    return plan.K == (length(s) - 1,) && plan.dim == (1,) &&
           size(plan.fft_plan) == size_nd &&
           plan.nodes[1] == chebpoints(length(s) - 1, minimum(s), maximum(s))
end

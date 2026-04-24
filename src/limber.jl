"""
    _limber_contraction(P_term, K1, K2, weights, Δχ)

Internal helper that performs the χ-line-of-sight contraction for Limber terms.

Given a power-spectrum term sampled on `(ℓ, χ)` and two Limber kernels,
this computes binned angular spectra via Simpson-weighted integration.
"""
function _limber_contraction(P_term, K1, K2, weights, Δχ)
    @tullio Cℓ[idx_l, idx_i, idx_j] := P_term[idx_l, idx_m] * K1[idx_l, idx_m, idx_i] * K2[idx_l, idx_m, idx_j] * weights[idx_m] * Δχ
    return Cℓ
end

# Limber integration helpers: all based on fixed Blast.χ and ℓ grids.
# `LIMBER_INV_χ2_ROW` is pre-shaped as a `(1, nχ)` row so it can be broadcast
# against `pk.{ΔP,Pδ}_limber` rows without per-call `reshape` allocations.
const LIMBER_N_NONLIMBER = length(Blast.ℓ_nonlimber)
const LIMBER_INV_χ2_ROW = reshape(1.0 ./ (Blast.χ .^ 2), 1, :)
const LIMBER_Δχ = (last(Blast.χ) - first(Blast.χ)) / (length(Blast.χ) - 1)
const LIMBER_WEIGHTS = Blast.simpson_weights_array(length(Blast.χ))

"""
    get_limber_k_polynomials(plan::ChebyshevPlan, l::AbstractVector, chi::AbstractVector; is_log_k=true) -> Array

Precompute the Chebyshev basis polynomials for the Limber grid `k = (l .+ 0.5) ./ chi'`.
This part is constant if the multipole grid (ℓ) and distance grid (χ) are fixed.
"""
function get_limber_k_polynomials(plan::ChebyshevPlan, l::AbstractVector, chi::AbstractVector; is_log_k=true)
    K_k = plan.K[1]
    k_min = last(plan.nodes[1]); k_max = first(plan.nodes[1])
    n_l = length(l); n_chi = length(chi)

    k_mat = @. (l + 0.5) / chi'
    if is_log_k
        k_mat = log10.(k_mat)
    end
    
    T = eltype(plan.nodes[1])
    T_k = Array{T}(undef, n_l, K_k + 1, n_chi)
    for j in 1:n_chi
        T_k[:, :, j] = chebyshev_polynomials(vec(k_mat[:, j]), k_min, k_max, K_k)
    end
    return T_k
end

"""
    get_limber_coords_polynomials(plan::ChebyshevPlan, coords::AbstractVector) -> Matrix

Precompute the Chebyshev basis polynomials for the coordinate axis (usually redshift z).
"""
function get_limber_coords_polynomials(plan::ChebyshevPlan, coords::AbstractVector)
    K_chi = plan.K[2]
    chi_min = last(plan.nodes[2]); chi_max = first(plan.nodes[2])
    # ACE's `chebyshev_polynomials(::AbstractVector{T}, ::T, ::T, ::Int)` requires
    # all three first arguments to share a concrete type. Under ForwardDiff
    # through cosmology, `coords` (= bg.z) carries Duals while plan bounds are
    # Float64; promote the bounds to match.
    T = eltype(coords)
    return chebyshev_polynomials(coords, T(chi_min), T(chi_max), K_chi)
end

"""
    get_limber_polynomials(plan, l, chi; is_log_k=true) -> (Matrix, Array)

Legacy wrapper that returns both coordinate and k polynomials.
"""
function get_limber_polynomials(plan::ChebyshevPlan, l::AbstractVector, chi::AbstractVector; is_log_k=true)
    T_chi = get_limber_coords_polynomials(plan, chi)
    T_k = get_limber_k_polynomials(plan, l, chi; is_log_k=is_log_k)
    return T_chi, T_k
end

"""
    limber_eval(c, T_chi, T_k) -> Matrix

Evaluate a 2D Chebyshev expansion f(k, chi) on the Limber grid using precomputed polynomials.
"""
function limber_eval(c::AbstractMatrix, T_chi::AbstractMatrix, T_k::AbstractArray{<:Any,3})
    B = c * T_chi'
    @tullio result[i, j] := T_k[i, k, j] * B[k, j]
    return result
end

"""
    get_limber_kernel(Component::AbstractComponents)

Construct the Limber kernel associated with a single projected component.
Pure broadcast form: AD-native through both ForwardDiff and Mooncake without
needing a custom rrule, and fast because `prefactor .* kernel` fuses into a
single SIMD-friendly pass.
"""
function get_limber_kernel(Component::AbstractComponents)
    total_prefactor = Component.ell_prefactor .* Component.limber_factor
    kernel = Component.Kernel'
    kernel = reshape(kernel, 1, size(kernel, 1), size(kernel, 2))
    prefactor = reshape(total_prefactor, :, 1, 1)
    return prefactor .* kernel
end

function get_limber_kernel(Component::Nothing)
    return nothing
end

# Sum two same-shape Limber kernels, short-circuiting when a component is
# absent. Used by the per-probe aggregators below. Returning `nothing` when
# every component is missing lets downstream slicing / contraction helpers
# dispatch to their `::Nothing` fast paths instead of touching a scalar `0.`.
_combine_limber_kernels(K1::AbstractArray{<:Any,3}, K2::AbstractArray{<:Any,3}) = K1 .+ K2
_combine_limber_kernels(K1::AbstractArray{<:Any,3}, ::Nothing) = K1
_combine_limber_kernels(::Nothing, K2::AbstractArray{<:Any,3}) = K2
_combine_limber_kernels(::Nothing, ::Nothing) = nothing

"""
    get_limber_kernel(G::GalaxyClustering)

Sum Limber kernels over implemented `GalaxyClustering` Limber components.

Current implementation includes number counts `δ` and magnification `μ`.
"""
function get_limber_kernel(G::GalaxyClustering)
    return _combine_limber_kernels(get_limber_kernel(G.δ), get_limber_kernel(G.μ))
end

"""
    get_limber_kernel(L::WeakLensing)

Sum Limber kernels over all active `WeakLensing` components (γ, IA).
"""
function get_limber_kernel(L::WeakLensing)
    return _combine_limber_kernels(get_limber_kernel(L.γ), get_limber_kernel(L.IA))
end

"""
    get_limber_kernel(C::CMB)

Sum Limber kernels over all active `CMB` components (κ, ISW).
"""
function get_limber_kernel(C::CMB)
    return _combine_limber_kernels(get_limber_kernel(C.κ), get_limber_kernel(C.ISW))
end

function _low_ℓ_slice(K::AbstractArray{<:Any,3})
    return @view K[1:LIMBER_N_NONLIMBER, :, :]
end

function _high_ℓ_slice(K::AbstractArray{<:Any,3})
    return @view K[(LIMBER_N_NONLIMBER+1):end, :, :]
end

@doc raw"""
    get_limber_correction(Probe, pk)

Compute the low-ℓ Limber correction for auto-spectra of a single probe.

This function evaluates the correction term built from `pk.ΔP_limber` on
`Blast.ℓ_nonlimber` and contracts it with the probe Limber kernel:

```math
\Delta C_\ell^{ij} = \int d\chi\,\frac{\Delta P_\ell(\chi)}{\chi^2}
K_i(\ell,\chi)K_j(\ell,\chi).
```
"""
function get_limber_correction(Probe::Union{GalaxyClustering, WeakLensing, CMB}, pk::PowerSpectrum)
    return get_limber_correction(get_limber_kernel(Probe), pk)
end

get_limber_correction(::Nothing, pk::PowerSpectrum) = 0.

function get_limber_correction(K::AbstractArray{<:Any,3}, pk::PowerSpectrum)
    ΔP_over_χ2 = @views pk.ΔP_limber[1:LIMBER_N_NONLIMBER, :] .* LIMBER_INV_χ2_ROW
    K_low = _low_ℓ_slice(K)
    return _limber_contraction(ΔP_over_χ2, K_low, K_low, LIMBER_WEIGHTS, LIMBER_Δχ)
end

@doc raw"""
    get_limber_correction(ProbeA, ProbeB, pk)

Compute the low-ℓ Limber correction for cross-spectra between two probes.

This uses the same correction source term as the auto case (`pk.ΔP_limber`),
with one kernel from each probe:

```math
\Delta C_\ell^{ij} = \int d\chi\,\frac{\Delta P_\ell(\chi)}{\chi^2}
K_i^{A}(\ell,\chi)K_j^{B}(\ell,\chi).
```
"""
function get_limber_correction(ProbeA::Union{GalaxyClustering, WeakLensing, CMB}, ProbeB::Union{GalaxyClustering, WeakLensing, CMB}, pk::PowerSpectrum)
    return get_limber_correction(get_limber_kernel(ProbeA), get_limber_kernel(ProbeB), pk)
end

get_limber_correction(::Nothing, ::Nothing, pk::PowerSpectrum) = 0.
get_limber_correction(::Nothing, KB::AbstractArray{<:Any,3}, pk::PowerSpectrum) = 0.
get_limber_correction(KA::AbstractArray{<:Any,3}, ::Nothing, pk::PowerSpectrum) = 0.

function get_limber_correction(KA::AbstractArray{<:Any,3}, KB::AbstractArray{<:Any,3}, pk::PowerSpectrum)
    ΔP_over_χ2 = @views pk.ΔP_limber[1:LIMBER_N_NONLIMBER, :] .* LIMBER_INV_χ2_ROW
    return _limber_contraction(ΔP_over_χ2, _low_ℓ_slice(KA), _low_ℓ_slice(KB), LIMBER_WEIGHTS, LIMBER_Δχ)
end

@doc raw"""
    get_limber_Cℓ(Probe, pk)

Compute Limber auto-spectra for the high-ℓ regime.

This uses `pk.Pδ_limber` on multipoles above `Blast.ℓ_nonlimber` and performs
the same χ-contraction with probe kernels:

```math
C_\ell^{ij,\mathrm{Limber}} = \int d\chi\,\frac{P_\ell(\chi)}{\chi^2}
K_i(\ell,\chi)K_j(\ell,\chi).
```
"""
function get_limber_Cℓ(Probe::Union{GalaxyClustering, WeakLensing, CMB}, pk::PowerSpectrum)
    return get_limber_Cℓ(get_limber_kernel(Probe), pk)
end

get_limber_Cℓ(::Nothing, pk::PowerSpectrum) = 0.

function get_limber_Cℓ(K::AbstractArray{<:Any,3}, pk::PowerSpectrum)
    Pδ_over_χ2 = @views pk.Pδ_limber[(LIMBER_N_NONLIMBER+1):end, :] .* LIMBER_INV_χ2_ROW
    K_high = _high_ℓ_slice(K)
    return _limber_contraction(Pδ_over_χ2, K_high, K_high, LIMBER_WEIGHTS, LIMBER_Δχ)
end

@doc raw"""
    get_limber_Cℓ(ProbeA, ProbeB, pk)

Compute Limber cross-spectra for the high-ℓ regime between two probes:

```math
C_\ell^{ij,\mathrm{Limber}} = \int d\chi\,\frac{P_\ell(\chi)}{\chi^2}
K_i^{A}(\ell,\chi)K_j^{B}(\ell,\chi).
```
"""
function get_limber_Cℓ(ProbeA::Union{GalaxyClustering, WeakLensing, CMB}, ProbeB::Union{GalaxyClustering, WeakLensing, CMB}, pk::PowerSpectrum)
    return get_limber_Cℓ(get_limber_kernel(ProbeA), get_limber_kernel(ProbeB), pk)
end

get_limber_Cℓ(::Nothing, ::Nothing, pk::PowerSpectrum) = 0.
get_limber_Cℓ(::Nothing, KB::AbstractArray{<:Any,3}, pk::PowerSpectrum) = 0.
get_limber_Cℓ(KA::AbstractArray{<:Any,3}, ::Nothing, pk::PowerSpectrum) = 0.

function get_limber_Cℓ(KA::AbstractArray{<:Any,3}, KB::AbstractArray{<:Any,3}, pk::PowerSpectrum)
    Pδ_over_χ2 = @views pk.Pδ_limber[(LIMBER_N_NONLIMBER+1):end, :] .* LIMBER_INV_χ2_ROW
    return _limber_contraction(Pδ_over_χ2, _high_ℓ_slice(KA), _high_ℓ_slice(KB), LIMBER_WEIGHTS, LIMBER_Δχ)
end

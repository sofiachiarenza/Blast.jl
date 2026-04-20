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
"""
function get_limber_kernel(Component::AbstractComponents)
    total_prefactor = Component.ell_prefactor .* Component.limber_factor
    kernel = Component.Kernel'
    kernel = reshape(kernel, 1, size(kernel, 1), size(kernel, 2))  
    prefactor = reshape(total_prefactor, :, 1, 1)  
    return prefactor .* kernel  
end

function get_limber_kernel(Component::Nothing)
    return 0.
end

"""
    get_limber_kernel(G::GalaxyClustering)

Sum Limber kernels over implemented `GalaxyClustering` Limber components.

Current implementation includes number counts `δ` and magnification `μ`.
"""
function get_limber_kernel(G::GalaxyClustering)
    return get_limber_kernel(G.δ) .+ get_limber_kernel(G.μ) 
end

"""
    get_limber_kernel(L::WeakLensing)

Sum Limber kernels over all active `WeakLensing` components (γ, IA).
"""
function get_limber_kernel(L::WeakLensing)
    return get_limber_kernel(L.γ) .+ get_limber_kernel(L.IA)
end

"""
    get_limber_kernel(C::CMB)

Sum Limber kernels over all active `CMB` components (κ, ISW).
"""
function get_limber_kernel(C::CMB)
    return get_limber_kernel(C.κ) .+ get_limber_kernel(C.ISW)
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
    chi_grid = Blast.χ
    n = size(chi_grid, 1)
    ΔP_over_χ2 = pk.ΔP_limber[1:size(Blast.ℓ_nonlimber, 1), :] ./ reshape(chi_grid, 1, :) .^ 2

    Δχ = ((chi_grid[n]-chi_grid[1])/(n-1))
    weights = Blast.simpson_weights_array(n)

    K = get_limber_kernel(Probe)[1:size(Blast.ℓ_nonlimber, 1), :, :]

    return _limber_contraction(ΔP_over_χ2, K, K, weights, Δχ)
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
    chi_grid = Blast.χ
    n = size(chi_grid, 1)
    ΔP_over_χ2 = pk.ΔP_limber[1:size(Blast.ℓ_nonlimber, 1), :] ./ reshape(chi_grid, 1, :) .^ 2

    Δχ = ((chi_grid[n]-chi_grid[1])/(n-1))
    weights = Blast.simpson_weights_array(n)

    KA = get_limber_kernel(ProbeA)[1:size(Blast.ℓ_nonlimber, 1), :, :]
    KB = get_limber_kernel(ProbeB)[1:size(Blast.ℓ_nonlimber, 1), :, :]

    return _limber_contraction(ΔP_over_χ2, KA, KB, weights, Δχ)
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
    chi_grid = Blast.χ
    n = size(chi_grid, 1)
    Pδ_over_χ2 = pk.Pδ_limber[size(Blast.ℓ_nonlimber, 1)+1:end, :] ./ reshape(chi_grid, 1, :) .^ 2

    Δχ = ((chi_grid[n]-chi_grid[1])/(n-1))
    weights = Blast.simpson_weights_array(n)

    K = get_limber_kernel(Probe)[size(Blast.ℓ_nonlimber, 1)+1:end, :, :]

    return _limber_contraction(Pδ_over_χ2, K, K, weights, Δχ)
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
    chi_grid = Blast.χ
    n = size(chi_grid, 1)
    Pδ_over_χ2 = pk.Pδ_limber[size(Blast.ℓ_nonlimber, 1)+1:end, :] ./ reshape(chi_grid, 1, :) .^ 2

    Δχ = ((chi_grid[n]-chi_grid[1])/(n-1))
    weights = Blast.simpson_weights_array(n)

    KA = get_limber_kernel(ProbeA)[size(Blast.ℓ_nonlimber, 1)+1:end, :, :]
    KB = get_limber_kernel(ProbeB)[size(Blast.ℓ_nonlimber, 1)+1:end, :, :]

    return _limber_contraction(Pδ_over_χ2, KA, KB, weights, Δχ)
end

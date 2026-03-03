function _limber_contraction(P_term, K1, K2, weights, Δχ)
    @tullio Cℓ[idx_l, idx_i, idx_j] := P_term[idx_l, idx_m] * K1[idx_l, idx_m, idx_i] * K2[idx_l, idx_m, idx_j] * weights[idx_m] * Δχ
    return Cℓ
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

function get_limber_kernel(G::GalaxyClustering)
    return get_limber_kernel(G.δ) .+ get_limber_kernel(G.μ) 
end

function get_limber_kernel(L::WeakLensing)
    return get_limber_kernel(L.γ) .+ get_limber_kernel(L.IA)
end

function get_limber_kernel(C::CMB)
    return get_limber_kernel(C.κ) .+ get_limber_kernel(C.ISW)
end

"""
    get_limber_correction(Probe, pk)
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

"""
    get_limber_correction(ProbeA, ProbeB, pk)
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

"""
    get_limber_Cℓ(Probe, pk)
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

"""
    get_limber_Cℓ(ProbeA, ProbeB, pk)
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

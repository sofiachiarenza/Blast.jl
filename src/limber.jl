"""
    get_limber_kernel(Component::AbstractComponents)

Construct the Limber kernel associated with a single projected component.

The kernel is defined as the product of:
- an ℓ-dependent prefactor,
- a Limber projection factor,
- the line-of-sight kernel evaluated on the global comoving distance grid.

# Returns
A 3D array with dimensions: `(n\\ell, n\\chi, n_bins)`
where:
- `ℓ` runs over the multipole range,
- `χ` is the global comoving distance grid,
- `n_bins` labels tomographic bins.
"""
function get_limber_kernel(Component::AbstractComponents)
    total_prefactor = Component.ell_prefactor .* Component.limber_factor
    kernel = Component.Kernel'
    kernel = reshape(kernel, 1, size(kernel, 1), size(kernel, 2))  # (1, 200, nbins)
    prefactor = reshape(total_prefactor, :, 1, 1)  # (101, 1, 1)
    return prefactor .* kernel  # Result: (101, 200, nbins)
end

"""
    get_limber_kernel(::Nothing)

Return `nothing` for inactive components.
"""
function get_limber_kernel(Component::Nothing)
    return nothing
end

"""
    get_limber_kernel(G::GalaxyClustering)

Return the total Limber kernel for the galaxy clustering probe.

# Notes
Redshift-space distortions and PNG contributions are currently excluded
from the Limber correction as these effects are on very large scales.
"""
function get_limber_kernel(G::GalaxyClustering)
    #TODO: in the correction i'm currently excluding RSDs. The limber implementation still doesn't work, but also the contribution at the scales of interest is null basically.
    return get_limber_kernel(G.δ) .+ get_limber_kernel(G.μ) # .+ get_limber_kernel(G.PNG) .+ get_limber_kernel(G.RSD)
end

"""
    get_limber_kernel(L::WeakLensing)

Return the total Limber kernel for the weak lensing probe.
"""
function get_limber_kernel(L::WeakLensing)
    return get_limber_kernel(L.γ) .+ get_limber_kernel(L.IA)
end

"""
    get_limber_kernel(C::CMB)

Return the total Limber kernel for the CMB probe.
"""
function get_limber_kernel(C::CMB)
    return get_limber_kernel(C.κ) .+ get_limber_kernel(C.ISW)
end

"""
    get_limber_correction(Probe, pk)

Compute the Limber correction to the non-Limber angular power spectrum. 
For `\\ell < 215`, the non-Limber integral is performed with the linear power spectrum. The correction for the non-linear
power spectrum is computed according to Eq. (36) and (37) of https://arxiv.org/pdf/2410.03632.


# Arguments
- `Probe`: A cosmological probe (`GalaxyClustering`, `WeakLensing`, or `CMB`)
- `pk::PowerSpectrum`: Power spectrum object containing the linear and non-linear power spectra already on the correct grid.

# Returns
A 3D array `Cℓ[ℓ, i, j]` containing the correction to the angular power spectra.
"""
function get_limber_correction(Probe::Union{GalaxyClustering, WeakLensing, CMB}, pk::PowerSpectrum)
    ΔP_over_χ2 = pk.ΔP_limber[1:size(Blast.ℓ_nonlimber, 1), :] ./ reshape(Blast.χ, 1, :) .^ 2

    n = size(Blast.χ, 1)
    Δχ = ((χ[n]-χ[1])/(n-1))
    weights = Blast.simpson_weight_array(n)

    K = get_limber_kernel(Probe)[1:size(Blast.ℓ_nonlimber, 1), :, :]

    @tullio Cℓ[l,i,j] := ΔP_over_χ2[l,m]*K[l,m,i]*K[l,m,j]*weights[m]*Δχ
end

"""
    get_limber_correction(ProbeA, ProbeB, pk)

Compute the Limber correction to the non-Limber angular cross-power spectrum
between two probes.The correction for the non-linear power spectrum is computed according to Eq. (36) and (37) of https://arxiv.org/pdf/2410.03632.

# Arguments
- `Probe`: A cosmological probe (`GalaxyClustering`, `WeakLensing`, or `CMB`)
- `pk::PowerSpectrum`: Power spectrum object containing the linear and non-linear power spectra already on the correct grid.

# Returns
A 3D array `Cℓ[ℓ, i, j]` containing the correction to the angular power spectra.
"""
function get_limber_correction(ProbeA::Union{GalaxyClustering, WeakLensing, CMB}, ProbeB::Union{GalaxyClustering, WeakLensing, CMB}, pk::PowerSpectrum)
    ΔP_over_χ2 = pk.ΔP_limber[1:size(Blast.ℓ_nonlimber, 1), :] ./ reshape(Blast.χ, 1, :) .^ 2

    n = size(Blast.χ, 1)
    Δχ = ((χ[n]-χ[1])/(n-1))
    weights = Blast.simpson_weight_array(n)

    KA = get_limber_kernel(ProbeA)[1:size(Blast.ℓ_nonlimber, 1), :, :]
    KB = get_limber_kernel(ProbeB)[1:size(Blast.ℓ_nonlimber, 1), :, :]

    @tullio Cℓ[l,i,j] := ΔP_over_χ2[l,m]*KA[l,m,i]*KB[l,m,j]*weights[m]*Δχ
end

"""
    get_limber_Cℓ(Probe, pk)

Compute the angular power spectrum using the full Limber approximation
for multipoles `ℓ > 215`.

# Arguments
- `Probe`: A cosmological probe (`GalaxyClustering`, `WeakLensing`, or `CMB`)
- `pk::PowerSpectrum`: Power spectrum object containing the linear and non-linear power spectra already on the correct grid.

# Returns
A 3D array `Cℓ[ℓ, i, j]` containing the Limber angular power spectrum.
"""
function get_limber_Cℓ(Probe::Union{GalaxyClustering, WeakLensing, CMB}, pk::PowerSpectrum)
    Pδ_over_χ2 = pk.Pδ_limber[size(Blast.ℓ_nonlimber, 1)+1:end, :] ./ reshape(Blast.χ, 1, :) .^ 2

    n = size(Blast.χ, 1)
    Δχ = ((χ[n]-χ[1])/(n-1))
    weights = Blast.simpson_weight_array(n)

    K = get_limber_kernel(Probe)[size(Blast.ℓ_nonlimber, 1)+1:end, :, :]

    @tullio Cℓ[l,i,j] := Pδ_over_χ2[l,m]*K[l,m,i]*K[l,m,j]*weights[m]*Δχ
end

"""
    get_limber_Cℓ(Probe, pk)

Compute the angular cross-power spectrum using the full Limber approximation
for multipoles `ℓ > 215`.

# Arguments
- `Probe`: A cosmological probe (`GalaxyClustering`, `WeakLensing`, or `CMB`)
- `pk::PowerSpectrum`: Power spectrum object containing the linear and non-linear power spectra already on the correct grid.

# Returns
A 3D array `Cℓ[ℓ, i, j]` containing the Limber angular power spectrum.
"""
function get_limber_Cℓ(ProbeA::Union{GalaxyClustering, WeakLensing, CMB}, ProbeB::Union{GalaxyClustering, WeakLensing, CMB}, pk::PowerSpectrum)
    Pδ_over_χ2 = pk.Pδ_limber[size(Blast.ℓ_nonlimber, 1)+1:end, :] ./ reshape(Blast.χ, 1, :) .^ 2

    n = size(Blast.χ, 1)
    Δχ = ((χ[n]-χ[1])/(n-1))
    weights = Blast.simpson_weight_array(n)

    KA = get_limber_kernel(ProbeA)[size(Blast.ℓ_nonlimber, 1)+1:end, :, :]
    KB = get_limber_kernel(ProbeB)[size(Blast.ℓ_nonlimber, 1)+1:end, :, :]

    @tullio Cℓ[l,i,j] := Pδ_over_χ2[l,m]*KA[l,m,i]*KB[l,m,j]*weights[m]*Δχ
end



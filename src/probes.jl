"""
    AbstractCosmologicalProbes

Abstract supertype for cosmological probes in BLAST.
"""
abstract type AbstractCosmologicalProbes end

"""
    AbstractComponents

Abstract supertype for physical components entering cosmological probes.
"""
abstract type AbstractComponents end

"""
    NumberCounts <: AbstractComponents
"""
@kwdef mutable struct NumberCounts <: AbstractComponents
    nz::Array{<:Any, 2} = zeros(1, 1)
    z::Array{<:Any, 1} = zeros(1)
    bias::Array{<:Any, 2} = zeros(1, 1)
    Kernel::Array{<:Any, 2} = zeros(1, 1)
    ell_prefactor::AbstractVector = ones(size(Blast.full_ℓ_range, 1))
    limber_factor::AbstractVector = ones(size(Blast.full_ℓ_range, 1))
end

"""
    CosmicShear <: AbstractComponents
"""
@kwdef mutable struct CosmicShear <: AbstractComponents
    nz::Array{<:Any, 2} = zeros(1, 1)
    z::Array{<:Any, 1} = zeros(1)
    Kernel::Array{<:Any, 2} = zeros(1, 1)
    ell_prefactor::AbstractVector = @. sqrt(factorial_frac(Blast.full_ℓ_range))
    limber_factor::AbstractVector = (Blast.full_ℓ_range .+ 0.5) .^ (-2)
end

"""
    CMBLensing <: AbstractComponents
"""
@kwdef mutable struct CMBLensing <: AbstractComponents
    Kernel::Array{<:Any, 2} = zeros(1, 1)
    ell_prefactor::AbstractVector = @. Blast.full_ℓ_range * (Blast.full_ℓ_range + 1)
    limber_factor::AbstractVector = (Blast.full_ℓ_range .+ 0.5) .^ (-2)
end

"""
    RedshiftSpaceDistortions <: AbstractComponents
"""
@kwdef mutable struct RedshiftSpaceDistortions <: AbstractComponents
    nz::Array{<:Any, 2} = zeros(1,1)
    z::Array{<:Any, 1} = zeros(1)
    growth_rate::Array{<:Any, 1} = zeros(1)
    Kernel::Array{<:Any, 2} = zeros(1, 1)
    ell_prefactor::AbstractVector = ones(size(Blast.full_ℓ_range, 1))
    limber_factor::AbstractVector = ones(size(Blast.full_ℓ_range, 1))
end

"""
    MagnificationBias <: AbstractComponents
"""
@kwdef mutable struct MagnificationBias <: AbstractComponents
    nz::Array{<:Any, 2} = zeros(1, 1)
    z::Array{<:Any, 1} = zeros(1)
    s::Array{<:Any, 2} = zeros(1,1)
    Kernel::Array{<:Any, 2} = zeros(1, 1)
    ell_prefactor::AbstractVector = @. Blast.full_ℓ_range * (Blast.full_ℓ_range + 1)
    limber_factor::AbstractVector = (Blast.full_ℓ_range .+ 0.5) .^ (-2)
end

"""
    IntrinsicAlignment <: AbstractComponents
"""
@kwdef mutable struct IntrinsicAlignment <: AbstractComponents
    nz::Array{<:Any, 2} = zeros(1, 1)
    z::Array{<:Any, 1} = zeros(1)
    A_IA::Array{<:Any, 2} = zeros(1, 1)
    Kernel::Array{<:Any, 2} = zeros(1, 1)
    ell_prefactor::AbstractVector = @. sqrt(factorial_frac(Blast.full_ℓ_range))
    limber_factor::AbstractVector = (Blast.full_ℓ_range .+ 0.5) .^ (-2)
end

"""
    IntegratedSachsWolfe <: AbstractComponents
"""
@kwdef mutable struct IntegratedSachsWolfe <: AbstractComponents
    growth_rate::Array{<:Any, 1} = zeros(1)
    Kernel::Array{<:Any, 2} = zeros(1, 1)
    ell_prefactor::AbstractVector = ones(size(Blast.full_ℓ_range, 1))
    limber_factor::AbstractVector = ones(size(Blast.full_ℓ_range, 1))
end

"""
    PrimordialNonGaussianity <: AbstractComponents
"""
@kwdef mutable struct PrimordialNonGaussianity <: AbstractComponents
    nz::Array{<:Any, 2} = zeros(1,1)
    z::Array{<:Any, 1} = zeros(1)
    bias::Array{<:Any, 2} = zeros(1, 1)
    f_NL::Number = 0
    p::Number = 0
    Kernel::Array{<:Any, 2} = zeros(1, 1)
    ell_prefactor::AbstractVector = ones(size(Blast.full_ℓ_range, 1))
    limber_factor::AbstractVector = ones(size(Blast.full_ℓ_range, 1))
end

"""
    GalaxyClustering <: AbstractCosmologicalProbes
"""
@kwdef mutable struct GalaxyClustering <: AbstractCosmologicalProbes
    δ::NumberCounts
    RSD::Union{RedshiftSpaceDistortions, Nothing} = nothing
    μ::Union{MagnificationBias, Nothing} = nothing
    PNG::Union{PrimordialNonGaussianity, Nothing} = nothing
end

"""
    WeakLensing <: AbstractCosmologicalProbes
"""
@kwdef mutable struct WeakLensing <: AbstractCosmologicalProbes
    γ::CosmicShear
    IA::Union{IntrinsicAlignment, Nothing} = nothing
end

"""
    CMB <: AbstractCosmologicalProbes
"""
@kwdef mutable struct CMB <: AbstractCosmologicalProbes
    κ::CMBLensing
    ISW::Union{IntegratedSachsWolfe, Nothing} = nothing
end

# --- Kernel Computation Logic (Handling nz normalization manually inside) ---

function compute_kernel!(Component::NumberCounts, grid::CosmologicalGrid, bg::BackgroundQuantities, cosmo::AbstractCosmology) 
    n_bins = size(Component.nz, 1)
    kernel = zeros(n_bins, length(grid.z_range))
    for b in 1:n_bins
        nz_func = DataInterpolations.AkimaInterpolation(Component.nz[b,:], Component.z, extrapolation = ExtrapolationType.Extension)
        nz_norm, _ = quadgk(x->nz_func(x), first(grid.z_range), last(grid.z_range))
        kernel[b,:] = @. Component.bias[b,:] * (bg.Hz_array / C_LIGHT) * (nz_func.(grid.z_range) / nz_norm)
    end
    Component.Kernel = kernel
end

function compute_kernel_safe!(Component::CosmicShear, grid::CosmologicalGrid, bg::BackgroundQuantities, cosmo::AbstractCosmology) 
    n_bins = size(Component.nz, 1)
    kernel = zeros(n_bins, length(grid.z_range))
    prefac = 1.5 * cosmo.H0^2 * cosmo.Ωm / C_LIGHT^2
    χ_interp = DataInterpolations.AkimaInterpolation(bg.χz_array, grid.z_range, extrapolation=ExtrapolationType.Extension)
    for b in 1:n_bins
        nz_func = DataInterpolations.AkimaInterpolation(Component.nz[b,:], Component.z, extrapolation=ExtrapolationType.Extension)
        nz_norm, _ = quadgk(x->nz_func(x), first(grid.z_range), last(grid.z_range))
        for z_idx in 1:length(grid.z_range)
            integrand(x) = nz_func(x) * (1. - bg.χz_array[z_idx]/χ_interp(x))
            z_low = grid.z_range[z_idx]
            z_top = grid.z_range[end]
            int, _ = quadgk(x -> integrand(x), z_low, z_top) 
            kernel[b, z_idx] = prefac * bg.χz_array[z_idx] * (1. + grid.z_range[z_idx]) * int / nz_norm
        end
    end
    Component.Kernel = kernel
end

function compute_kernel!(Component::CosmicShear, grid::CosmologicalGrid, bg::BackgroundQuantities, cosmo::AbstractCosmology)
    n_bins = size(Component.nz, 1)
    kernel = zeros(n_bins, length(grid.z_range))
    n_z_array = zeros(n_bins, length(grid.z_range))
    χz_array = bg.χz_array
    z_range = grid.z_range
    prefac = 1.5 * cosmo.H0^2 * cosmo.Ωm / C_LIGHT^2
    simpson_matrix = simpson_weights_matrix(length(grid.z_range))
    Δχ = (χz_array[end] - χz_array[1]) / (length(χz_array) - 1)
    for b in 1:n_bins
        nz_func = DataInterpolations.AkimaInterpolation(Component.nz[b, :], Component.z, extrapolation=ExtrapolationType.Extension)
        nz_norm, _ = quadgk(x->nz_func(x), first(grid.z_range), last(grid.z_range))
        for (zidx, myz) in enumerate(grid.z_range)
            n_z_array[b, zidx] = nz_func(myz) / nz_norm
        end
    end
    @tullio kernel[idx_b, idx_zidx] := Δχ * (bg.Hz_array[idx_zp] / C_LIGHT) * simpson_matrix[idx_zidx, idx_zp] * prefac * χz_array[idx_zidx] * (1.0 + z_range[idx_zidx]) * n_z_array[idx_b, idx_zp] * (χz_array[idx_zp] - χz_array[idx_zidx]) / (χz_array[idx_zp] + 1e-18)
    Component.Kernel = kernel
end

function compute_kernel!(Component::CMBLensing, grid::CosmologicalGrid, bg::BackgroundQuantities, cosmo::AbstractCosmology) 
    prefac = 1.5 * cosmo.H0^2 * cosmo.Ωm / C_LIGHT^2
    χ_CMB = compute_χ(1090, cosmo) 
    kernel = @. prefac * bg.χz_array * (1. + grid.z_range) * (1 - bg.χz_array/χ_CMB)
    Component.Kernel = reshape(kernel, 1, size(kernel,1))
end

function compute_kernel!(Component::RedshiftSpaceDistortions, grid::CosmologicalGrid,  bg::BackgroundQuantities, cosmo::AbstractCosmology) 
    n_bins = size(Component.nz, 1)
    kernel = zeros(n_bins, length(grid.z_range))
    for b in 1:n_bins
        nz_func = DataInterpolations.AkimaInterpolation(Component.nz[b,:], Component.z, extrapolation=ExtrapolationType.Extension)
        nz_norm, _ = quadgk(x->nz_func(x), first(grid.z_range), last(grid.z_range))
        kernel[b,:] = @. Component.growth_rate * (bg.Hz_array / C_LIGHT) * (nz_func.(grid.z_range) / nz_norm) 
    end
    Component.Kernel = kernel
end

function compute_kernel_safe!(Component::MagnificationBias, grid::CosmologicalGrid, bg::BackgroundQuantities, cosmo::AbstractCosmology) 
    n_bins = size(Component.nz, 1)
    kernel = zeros(n_bins, length(grid.z_range))
    prefac = 1.5 * cosmo.H0^2 * cosmo.Ωm / C_LIGHT^2
    χ_interp = DataInterpolations.AkimaInterpolation(bg.χz_array, grid.z_range, extrapolation=ExtrapolationType.Extension)
    for b in 1:n_bins
        nz_func = DataInterpolations.AkimaInterpolation(Component.nz[b,:], Component.z, extrapolation=ExtrapolationType.Extension)
        nz_norm, _ = quadgk(x->nz_func(x), first(grid.z_range), last(grid.z_range))
        s_z = DataInterpolations.AkimaInterpolation(Component.s[b,:], grid.z_range, extrapolation=ExtrapolationType.Extension)
        for z_idx in 1:length(grid.z_range)
            integrand(x) = nz_func(x) * (1. - bg.χz_array[z_idx]/χ_interp(x)) * (5 .* s_z(x) .- 2)
            z_low = grid.z_range[z_idx]
            z_top =  grid.z_range[end]
            int, _ = quadgk(x -> integrand(x), z_low, z_top) 
            kernel[b, z_idx] = prefac * bg.χz_array[z_idx] * (1. + grid.z_range[z_idx]) * int / nz_norm
        end
    end
    Component.Kernel = kernel 
end

function compute_kernel!(Component::MagnificationBias, grid::CosmologicalGrid, bg::BackgroundQuantities, cosmo::AbstractCosmology)
    n_bins = size(Component.nz, 1)
    s_z_array = zeros(n_bins, length(grid.z_range))
    n_z_array = zeros(n_bins, length(grid.z_range))
    χz_array = bg.χz_array
    z_range = grid.z_range
    prefac = 1.5 * cosmo.H0^2 * cosmo.Ωm / C_LIGHT^2
    simpson_matrix = simpson_weights_matrix(length(grid.z_range))
    Δχ = (χz_array[end] - χz_array[1]) / (length(χz_array) - 1)
    for b in 1:n_bins
        nz_func = DataInterpolations.AkimaInterpolation(Component.nz[b, :], Component.z, extrapolation=ExtrapolationType.Extension)
        nz_norm, _ = quadgk(x->nz_func(x), first(grid.z_range), last(grid.z_range))
        s_z = DataInterpolations.AkimaInterpolation(Component.s[b, :], grid.z_range, extrapolation=ExtrapolationType.Extension)
        for (zidx, myz) in enumerate(grid.z_range)
            s_z_array[b, zidx] = s_z(myz)
            n_z_array[b, zidx] = nz_func(myz) / nz_norm
        end
    end
    @tullio kernel[idx_b, idx_zidx] := Δχ * (bg.Hz_array[idx_zp] / C_LIGHT) * simpson_matrix[idx_zidx, idx_zp] * prefac * χz_array[idx_zidx] * (1.0 + z_range[idx_zidx]) * n_z_array[idx_b, idx_zp] * (χz_array[idx_zp] - χz_array[idx_zidx]) / (χz_array[idx_zp] + 1e-18) * (5.0 * s_z_array[idx_b, idx_zp] - 2)
    Component.Kernel = kernel
end

function compute_kernel!(Component::IntrinsicAlignment, grid::CosmologicalGrid, bg::BackgroundQuantities, cosmo::AbstractCosmology) 
    n_bins = size(Component.nz, 1)
    kernel = zeros(n_bins, length(grid.z_range))
    for b in 1:n_bins
        nz_func = DataInterpolations.AkimaInterpolation(Component.nz[b,:], Component.z, extrapolation=ExtrapolationType.Extension)
        nz_norm, _ = quadgk(x->nz_func(x), first(grid.z_range), last(grid.z_range))
        kernel[b,:] = @. Component.A_IA[b,:] * (bg.Hz_array / C_LIGHT) * nz_func.(grid.z_range) / nz_norm 
    end
    Component.Kernel = kernel
end

function compute_kernel!(Component::IntegratedSachsWolfe, grid::CosmologicalGrid, bg::BackgroundQuantities, cosmo::AbstractCosmology) 
    T_CMB = 2.726
    prefac = 3T_CMB * cosmo.H0^2 * cosmo.Ωm / C_LIGHT^3
    kernel = @. prefac * bg.Hz_array * (1 - Component.growth_rate)
    Component.Kernel = reshape(kernel, 1, size(kernel, 1))
end

function compute_kernel!(Component::PrimordialNonGaussianity, grid::CosmologicalGrid, bg::BackgroundQuantities, cosmo::AbstractCosmology) 
    n_bins = size(Component.nz, 1)
    kernel = zeros(n_bins, length(grid.z_range))
    for b in 1:n_bins
        nz_func = DataInterpolations.AkimaInterpolation(Component.nz[b,:], Component.z, extrapolation = ExtrapolationType.Extension)
        nz_norm, _ = quadgk(x->nz_func(x), first(grid.z_range), last(grid.z_range))
        kernel[b,:] = (bg.Hz_array / C_LIGHT) .* Component.f_NL .* bΦ(Component.bias[b,:], Component.p) .* (nz_func.(grid.z_range) / nz_norm)
    end
    Component.Kernel = kernel
end

function compute_kernel!(Component::Nothing, grid::CosmologicalGrid, bg::BackgroundQuantities, cosmo::AbstractCosmology) 
    return nothing
end

function compute_kernel_safe!(Component::Nothing, grid::CosmologicalGrid, bg::BackgroundQuantities, cosmo::AbstractCosmology) 
    return nothing
end

function evaluate_components!(GC::GalaxyClustering, grid::CosmologicalGrid, bg::BackgroundQuantities, cosmo::AbstractCosmology) 
    compute_kernel!(GC.δ, grid, bg, cosmo)
    compute_kernel!(GC.RSD, grid, bg, cosmo)
    compute_kernel!(GC.μ, grid, bg, cosmo)
    compute_kernel!(GC.PNG, grid, bg, cosmo)
end

function evaluate_components!(WL::WeakLensing, grid::CosmologicalGrid, bg::BackgroundQuantities, cosmo::AbstractCosmology)
    compute_kernel!(WL.γ, grid, bg, cosmo)
    compute_kernel!(WL.IA, grid, bg, cosmo)
end

function evaluate_components!(cmb::CMB, grid::CosmologicalGrid, bg::BackgroundQuantities, cosmo::AbstractCosmology)
    compute_kernel!(cmb.κ, grid, bg, cosmo)
    compute_kernel!(cmb.ISW, grid, bg, cosmo)
end

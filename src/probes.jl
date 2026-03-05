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
    prepare_nz_matrix(nz::AbstractMatrix, z::AbstractVector, z_grid::AbstractVector)

Interpolate and normalize an n(z) matrix bin-by-bin onto the calculation grid.
"""
function prepare_nz_matrix(nz::AbstractMatrix, z::AbstractVector, z_grid::AbstractVector)
    n_bins = size(nz, 1)
    nz_normed = zeros(eltype(nz), n_bins, length(z_grid))

    for b in 1:n_bins
        nz_func = DataInterpolations.AkimaInterpolation(nz[b, :], z, extrapolation=ExtrapolationType.Extension)
        nz_norm_val, _ = quadgk(x -> nz_func(x), first(z_grid), last(z_grid))
        nz_normed[b, :] = nz_func.(z_grid) ./ nz_norm_val
    end
    return nz_normed
end

"""
    check_and_normalize!(Component, grid_z)

Internal helper: ensures nz_norm is populated for the current calculation grid.
"""
function check_and_normalize!(Component::AbstractComponents, z_grid::AbstractVector)
    if hasfield(typeof(Component), :nz) && hasfield(typeof(Component), :nz_norm)
        if size(Component.nz_norm) != (size(Component.nz, 1), length(z_grid))
            Component.nz_norm = prepare_nz_matrix(Component.nz, Component.z, z_grid)
        end
    end
end

function check_and_normalize!(Component::Nothing, z_grid::AbstractVector)
    return nothing
end

"""
    NumberCounts <: AbstractComponents
"""
@kwdef mutable struct NumberCounts <: AbstractComponents
    nz::Array{<:Any, 2} = zeros(1, 1)
    z::Array{<:Any, 1} = zeros(1)
    nz_norm::Array{<:Any, 2} = zeros(1, 1)
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
    nz_norm::Array{<:Any, 2} = zeros(1, 1)
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
    nz_norm::Array{<:Any, 2} = zeros(1, 1)
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
    nz_norm::Array{<:Any, 2} = zeros(1, 1)
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
    nz_norm::Array{<:Any, 2} = zeros(1, 1)
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
    nz_norm::Array{<:Any, 2} = zeros(1, 1)
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

# --- Simplified Kernel Computation Logic (Using Background object) ---

function compute_kernel!(Component::NumberCounts, bg::Background) 
    check_and_normalize!(Component, bg.z)
    Component.Kernel = @. Component.bias * (bg.H' / C_LIGHT) * Component.nz_norm
end

function compute_kernel_safe!(Component::CosmicShear, bg::Background) 
    check_and_normalize!(Component, bg.z)
    n_bins = size(Component.nz_norm, 1)
    kernel = zeros(n_bins, length(bg.z))
    
    H0 = get_H0(bg.cosmo)
    Ωm = get_Ωm(bg.cosmo)
    prefac = 1.5 * H0^2 * Ωm / C_LIGHT^2

    for b in 1:n_bins
        nz_interp = DataInterpolations.AkimaInterpolation(Component.nz_norm[b,:], bg.z, extrapolation=ExtrapolationType.Extension)
        for z_idx in 1:length(bg.z)
            integrand(x) = nz_interp(x) * (1. - bg.χ[z_idx]/bg.χ_of_z(x))
            z_low = bg.z[z_idx]
            z_top = bg.z[end]
            int, _ = quadgk(x -> integrand(x), z_low, z_top) 
            kernel[b, z_idx] = prefac * bg.χ[z_idx] * (1. + bg.z[z_idx]) * int
        end
    end
    Component.Kernel = kernel
end

function compute_kernel!(Component::CosmicShear, bg::Background)
    check_and_normalize!(Component, bg.z)
    n_bins = size(Component.nz_norm, 1)
    
    H0 = get_H0(bg.cosmo)
    Ωm = get_Ωm(bg.cosmo)
    prefac = 1.5 * H0^2 * Ωm / C_LIGHT^2
    
    simpson_matrix = simpson_weights_matrix(length(bg.z))
    Δχ = (bg.χ[end] - bg.χ[1]) / (length(bg.χ) - 1)

    @tullio kernel[idx_b, idx_zidx] := Δχ * (bg.H[idx_zp] / C_LIGHT) * simpson_matrix[idx_zidx, idx_zp] * prefac * bg.χ[idx_zidx] * (1.0 + bg.z[idx_zidx]) * Component.nz_norm[idx_b, idx_zp] * (bg.χ[idx_zp] - bg.χ[idx_zidx]) / (bg.χ[idx_zp] + 1e-18)
    Component.Kernel = kernel
end

function compute_kernel!(Component::CMBLensing, bg::Background) 
    H0 = get_H0(bg.cosmo)
    Ωm = get_Ωm(bg.cosmo)
    prefac = 1.5 * H0^2 * Ωm / C_LIGHT^2
    
    χ_CMB = compute_χ(1090, bg.cosmo) 
    kernel = @. prefac * bg.χ * (1. + bg.z) * (1 - bg.χ/χ_CMB)
    Component.Kernel = reshape(kernel, 1, size(kernel,1))
end

function compute_kernel!(Component::RedshiftSpaceDistortions, bg::Background) 
    check_and_normalize!(Component, bg.z)
    Component.Kernel = @. bg.f' * (bg.H' / C_LIGHT) * Component.nz_norm
end

function compute_kernel_safe!(Component::MagnificationBias, bg::Background) 
    check_and_normalize!(Component, bg.z)
    n_bins = size(Component.nz_norm, 1)
    kernel = zeros(n_bins, length(bg.z))
    
    H0 = get_H0(bg.cosmo)
    Ωm = get_Ωm(bg.cosmo)
    prefac = 1.5 * H0^2 * Ωm / C_LIGHT^2

    for b in 1:n_bins
        nz_interp = DataInterpolations.AkimaInterpolation(Component.nz_norm[b,:], bg.z, extrapolation=ExtrapolationType.Extension)
        s_vals = Component.s[b,:]
        for z_idx in 1:length(bg.z)
            integrand(x) = nz_interp(x) * (1. - bg.χ[z_idx]/bg.χ_of_z(x)) * (5 .* s_vals[z_idx] .- 2)
            z_low = bg.z[z_idx]
            z_top = bg.z[end]
            int, _ = quadgk(x -> integrand(x), z_low, z_top) 
            kernel[b, z_idx] = prefac * bg.χ[z_idx] * (1. + bg.z[z_idx]) * int
        end
    end
    Component.Kernel = kernel 
end

function compute_kernel!(Component::MagnificationBias, bg::Background)
    check_and_normalize!(Component, bg.z)
    n_bins = size(Component.nz_norm, 1)
    
    H0 = get_H0(bg.cosmo)
    Ωm = get_Ωm(bg.cosmo)
    prefac = 1.5 * H0^2 * Ωm / C_LIGHT^2
    
    simpson_matrix = simpson_weights_matrix(length(bg.z))
    Δχ = (bg.χ[end] - bg.χ[1]) / (length(bg.χ) - 1)

    @tullio kernel[idx_b, idx_zidx] := Δχ * (bg.H[idx_zp] / C_LIGHT) * simpson_matrix[idx_zidx, idx_zp] * prefac * bg.χ[idx_zidx] * (1.0 + bg.z[idx_zidx]) * Component.nz_norm[idx_b, idx_zp] * (bg.χ[idx_zp] - bg.χ[idx_zidx]) / (bg.χ[idx_zp] + 1e-18) * (5.0 * Component.s[idx_b, idx_zp] - 2)
    Component.Kernel = kernel
end

function compute_kernel!(Component::IntrinsicAlignment, bg::Background) 
    check_and_normalize!(Component, bg.z)
    Component.Kernel = @. Component.A_IA * (bg.H' / C_LIGHT) * Component.nz_norm
end

function compute_kernel!(Component::IntegratedSachsWolfe, bg::Background) 
    H0 = get_H0(bg.cosmo)
    Ωm = get_Ωm(bg.cosmo)
    T_CMB = 2.726
    prefac = 3T_CMB * H0^2 * Ωm / C_LIGHT^3
    kernel = @. prefac * bg.H * (1 - bg.f)
    Component.Kernel = reshape(kernel, 1, size(kernel, 1))
end

function compute_kernel!(Component::PrimordialNonGaussianity, bg::Background) 
    check_and_normalize!(Component, bg.z)
    b_phi_vals = bΦ(Component.bias, Component.p)
    Component.Kernel = @. (bg.H' / C_LIGHT) * Component.f_NL * b_phi_vals * Component.nz_norm
end

function compute_kernel!(Component::Nothing, bg::Background) 
    return nothing
end

function compute_kernel_safe!(Component::Nothing, bg::Background) 
    return nothing
end

function evaluate_components!(GC::GalaxyClustering, bg::Background) 
    compute_kernel!(GC.δ, bg)
    compute_kernel!(GC.RSD, bg)
    compute_kernel_safe!(GC.μ, bg)
    compute_kernel!(GC.PNG, bg)
end

function evaluate_components!(WL::WeakLensing, bg::Background)
    compute_kernel_safe!(WL.γ, bg)
    compute_kernel!(WL.IA, bg)
end

function evaluate_components!(cmb::CMB, bg::Background)
    compute_kernel!(cmb.κ, bg)
    compute_kernel!(cmb.ISW, bg)
end

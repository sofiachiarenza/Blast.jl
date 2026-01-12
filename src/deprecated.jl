#TODO I think this is not used anywhere anynore, might get rid of it.
"""
    resample_redshifts(bg::BackgroundQuantities, grid::AbstractCosmologicalGrid, new_χ::AbstractArray{T,1}) where T

Resamples redshift values corresponding to a new set of comoving distances using Akima interpolation.

# Arguments
- `bg::BackgroundQuantities`: An object containing background cosmological quantities, 
  including the mapping between redshift (`z`) and comoving distance (`χ`).
- `grid::AbstractCosmologicalGrid`: An object defining the range of redshifts (`z_range`) 
  and associated cosmological grid quantities.
- `new_χ::AbstractArray{T,1}`: A 1D array of comoving distances for which corresponding redshift values are desired.

# Returns
- `resampled_z::AbstractArray{T,1}`: A 1D array of resampled redshift values corresponding to the input comoving distances `new_χ`.
"""
function resample_redshifts(bg::BackgroundQuantities, grid::AbstractCosmologicalGrid, new_χ::AbstractArray{T,1}) where T
    z_of_χ = AkimaInterpolation(grid.z_range, bg.χz_array, extrapolation=ExtrapolationType.Extension)
    return z_of_χ.(new_χ)
end

function Cℓ_limber(pk, ℓ, bg::BackgroundQuantities, ProbeA::Blast.GalaxyKernel, ProbeB::Blast.GalaxyKernel )
    χ = bg.χz_array
    n = length(χ)
    
    F = 1

    Δχ = ((χ[n]-χ[1])/(n-1))
    pesi = Blast.simpson_weight_array(n)

    pk_over_chi = pk ./ (χ .^ 2)

    KA = ProbeA.Kernel
    KB = ProbeB.Kernel

    @tullio Cℓ[i,j] := Δχ*pk_over_chi[m]*KA[i,m]*KB[j,m]*pesi[m]*F
    return Cℓ
end

function Cℓ_limber(pk, ℓ, bg::BackgroundQuantities, ProbeA::Blast.GalaxyKernel, ProbeB::Blast.ShearKernel )
    χ = bg.χz_array
    n = length(χ)
    
    F = @. sqrt(Blast.factorial_frac(ℓ)) / (ℓ+0.5)^2

    Δχ = ((χ[n]-χ[1])/(n-1))
    pesi = Blast.simpson_weight_array(n)

    pk_over_chi = pk ./ (χ .^ 2)

    KA = ProbeA.Kernel
    KB = ProbeB.Kernel

    @tullio Cℓ[i,j] := Δχ*pk_over_chi[m]*KA[i,m]*KB[j,m]*pesi[m]*F
    return Cℓ
end

function Cℓ_limber(pk, ℓ, bg::BackgroundQuantities, ProbeA::Blast.ShearKernel, ProbeB::Blast.ShearKernel )
    χ = bg.χz_array
    n = length(χ)
    
    F = @. Blast.factorial_frac(ℓ) / (ℓ+0.5)^4

    Δχ = ((χ[n]-χ[1])/(n-1))
    pesi = Blast.simpson_weight_array(n)

    pk_over_chi = pk ./ (χ .^ 2)

    KA = ProbeA.Kernel
    KB = ProbeB.Kernel

    @tullio Cℓ[i,j] := Δχ*pk_over_chi[m]*KA[i,m]*KB[j,m]*pesi[m]*F
    return Cℓ
end

function Cℓ_limber(pk, ℓ, bg::BackgroundQuantities, ProbeA::Blast.CMBLensingKernel, ProbeB::Blast.ShearKernel )
    χ = bg.χz_array
    n = length(χ)
    
    F = @. sqrt(Blast.factorial_frac(ℓ)) * ℓ * (ℓ+1) / (ℓ+0.5)^4

    Δχ = ((χ[n]-χ[1])/(n-1))
    pesi = Blast.simpson_weight_array(n)

    pk_over_chi = pk ./ (χ .^ 2)

    KA = reshape(ProbeA.Kernel, 1, n)
    KB = ProbeB.Kernel

    @tullio Cℓ[i,j] := Δχ*pk_over_chi[m]*KA[i,m]*KB[j,m]*pesi[m]*F
    return Cℓ
end

function Cℓ_limber(pk, ℓ, bg::BackgroundQuantities, ProbeA::Blast.CMBLensingKernel, ProbeB::Blast.GalaxyKernel )
    χ = bg.χz_array
    n = length(χ)
    bins = size(ProbeB.Kernel)[1]
    F = @. ℓ * (ℓ + 1) / (ℓ+0.5)^2

    Δχ = ((χ[n]-χ[1])/(n-1))
    pesi = Blast.simpson_weight_array(n)

    pk_over_chi = pk ./ (χ .^ 2)

    KA = reshape(ProbeA.Kernel, 1, length(χ))
    KB = ProbeB.Kernel

    @tullio Cℓ[i,j] := Δχ*pk_over_chi[m]*KA[i,m]*KB[j,m]*pesi[m]*F
    return Cℓ
end

function Cℓ_limber(pk, ℓ, bg::BackgroundQuantities, ProbeA::Blast.CMBLensingKernel, ProbeB::Blast.RSDKernel )
    χ = bg.χz_array
    n = length(χ)
    bins = size(ProbeB.Kernel)[1]
    F = @. ℓ * (ℓ + 1) / (ℓ+0.5)^2

    Δχ = ((χ[n]-χ[1])/(n-1))
    pesi = Blast.simpson_weight_array(n)

    pk_over_chi = pk ./ (χ .^ 2)

    KA = reshape(ProbeA.Kernel, 1, length(χ))
    KB = reshape(limber_rsd_kernel(ℓ, bg, ProbeB), bins, n)

    @tullio Cℓ[i,j] := Δχ*pk_over_chi[m]*KA[i,m]*KB[j,m]*pesi[m]*F
    return Cℓ
end

function Cℓ_limber(pk, ℓ, bg::BackgroundQuantities, ProbeA::Blast.CMBLensingKernel, ProbeB::Blast.LensMagKernel )
    χ = bg.χz_array
    n = length(χ)
    F = @. (ℓ * (ℓ + 1))^2 / (ℓ+0.5)^4

    Δχ = ((χ[n]-χ[1])/(n-1))
    pesi = Blast.simpson_weight_array(n)

    pk_over_chi = pk ./ (χ .^ 2)

    KA = reshape(ProbeA.Kernel, 1, length(χ))
    KB = ProbeB.Kernel

    @tullio Cℓ[i,j] := Δχ*pk_over_chi[m]*KA[i,m]*KB[j,m]*pesi[m]*F
    return Cℓ
end

function Cℓ_limber(pk, ℓ, bg::BackgroundQuantities, ProbeA::Blast.GalaxyKernel, ProbeB::Blast.RSDKernel )
    χ = bg.χz_array
    n = length(χ)
    F = 1

    Δχ = ((χ[n]-χ[1])/(n-1))
    pesi = Blast.simpson_weight_array(n)

    pk_over_chi = pk ./ (χ .^ 2)

    KA = ProbeA.Kernel
    KB = reshape(limber_rsd_kernel(ℓ, bg, ProbeB), size(ProbeB.Kernel)[1], n)

    @tullio Cℓ[i,j] := Δχ*pk_over_chi[m]*KA[i,m]*KB[j,m]*pesi[m]*F
    return Cℓ
end

function Cℓ_limber(pk, ℓ, bg::BackgroundQuantities, ProbeA::Blast.RSDKernel, ProbeB::Blast.RSDKernel )
    χ = bg.χz_array
    n = length(χ)
    bins = size(ProbeA.Kernel)[1]
    F = 1

    Δχ = ((χ[n]-χ[1])/(n-1))
    pesi = Blast.simpson_weight_array(n)

    pk_over_chi = pk ./ (χ .^ 2)

    KA = reshape(limber_rsd_kernel(ℓ, bg, ProbeA), bins, n)
    KB = reshape(limber_rsd_kernel(ℓ, bg, ProbeB), bins, n)

    @tullio Cℓ[i,j] := Δχ*pk_over_chi[m]*KA[i,m]*KB[j,m]*pesi[m]*F
    return Cℓ
end

function Cℓ_limber(pk, ℓ, bg::BackgroundQuantities, ProbeA::Blast.GalaxyKernel, ProbeB::Blast.LensMagKernel )
    χ = bg.χz_array
    n = length(χ)
    
    F = @. ℓ * (ℓ + 1) / (ℓ+0.5)^2

    Δχ = ((χ[n]-χ[1])/(n-1))
    pesi = Blast.simpson_weight_array(n)

    pk_over_chi = pk ./ (χ .^ 2)

    KA = ProbeA.Kernel
    KB = ProbeB.Kernel

    @tullio Cℓ[i,j] := Δχ*pk_over_chi[m]*KA[i,m]*KB[j,m]*pesi[m]*F
    return Cℓ
end

function Cℓ_limber(pk, ℓ, bg::BackgroundQuantities, ProbeA::Blast.LensMagKernel, ProbeB::Blast.LensMagKernel )
    χ = bg.χz_array
    n = length(χ)
    
    F = @. (ℓ * (ℓ + 1))^2 / (ℓ+0.5)^4

    Δχ = ((χ[n]-χ[1])/(n-1))
    pesi = Blast.simpson_weight_array(n)

    pk_over_chi = pk ./ (χ .^ 2)

    KA = ProbeA.Kernel
    KB = ProbeB.Kernel

    @tullio Cℓ[i,j] := Δχ*pk_over_chi[m]*KA[i,m]*KB[j,m]*pesi[m]*F
    return Cℓ
end

function Cℓ_limber(pk, ℓ, bg::BackgroundQuantities, ProbeA::Blast.RSDKernel, ProbeB::Blast.LensMagKernel )
    χ = bg.χz_array
    n = length(χ)
    bins = size(ProbeA.Kernel)[1]
    F = @. ℓ * (ℓ + 1) / (ℓ+0.5)^2

    Δχ = ((χ[n]-χ[1])/(n-1))
    pesi = Blast.simpson_weight_array(n)

    pk_over_chi = pk ./ (χ .^ 2)

    KA = reshape(limber_rsd_kernel(ℓ, bg, ProbeA), bins, n)
    KB = ProbeB.Kernel

    @tullio Cℓ[i,j] := Δχ*pk_over_chi[m]*KA[i,m]*KB[j,m]*pesi[m]*F
    return Cℓ
end;

function Cℓ_limber(pk, ℓ, bg::BackgroundQuantities, ProbeA::Blast.ShearKernel, ProbeB::Blast.LensMagKernel )
    χ = bg.χz_array
    n = length(χ)
    bins = size(ProbeA.Kernel)[1]
    F = @. sqrt(factorial_frac(ℓ)) / (ℓ+0.5)^4 * ℓ * (ℓ + 1) 

    Δχ = ((χ[n]-χ[1])/(n-1))
    pesi = Blast.simpson_weight_array(n)

    pk_over_chi = pk ./ (χ .^ 2)

    KA = ProbeA.Kernel
    KB = ProbeB.Kernel

    @tullio Cℓ[i,j] := Δχ*pk_over_chi[m]*KA[i,m]*KB[j,m]*pesi[m]*F
    return Cℓ
end

function Cℓ_limber(pk, ℓ, bg::BackgroundQuantities, ProbeA::Blast.ShearKernel, ProbeB::Blast.RSDKernel )
    χ = bg.χz_array
    n = length(χ)
    F = @. sqrt(factorial_frac(ℓ)) / (ℓ+0.5)^2

    Δχ = ((χ[n]-χ[1])/(n-1))
    pesi = Blast.simpson_weight_array(n)

    pk_over_chi = pk ./ (χ .^ 2)

    KA = ProbeA.Kernel
    KB = reshape(limber_rsd_kernel(ℓ, bg, ProbeB), size(ProbeB.Kernel)[1], n)

    @tullio Cℓ[i,j] := Δχ*pk_over_chi[m]*KA[i,m]*KB[j,m]*pesi[m]*F
    return Cℓ
end

function load_Ts(folder, nχ, nR, nk)
    ell_vector = Blast.ℓ_nonlimber
    full_T = zeros(length(ell_vector), nχ, nR, nk)
    for i in 1:length(ell_vector)
        l_string = string(round(ell_vector[i]; digits=1))
        filename = folder * "/T_tilde_l_$l_string.npy"
        if isfile(filename)
            full_T[i,:,:,:] = npzread(filename)
        else
            println("Missing file!")
        end
    end
    return full_T
end

abstract type AbstractOLDCosmologicalProbes{T} end

@kwdef mutable struct GalaxyKernel{T} <: AbstractOLDCosmologicalProbes{T}
    Kernel::AbstractArray{T, 2} = zeros(1, 1)
end

GalaxyKernel{T}(n_bins::Int, nχ::Int) where T = GalaxyKernel{T}(Kernel = zeros(T, n_bins, nχ))

GalaxyKernel(n_bins::Int, nχ::Int) = GalaxyKernel{Float64}(n_bins, nχ)


"""
    ShearKernel{T} <: AbstractCosmologicalProbes{T}

Represents a shear kernel used in cosmological lensing calculations. The kernel is defined over multiple tomographic bins.

# Parameters
- `Kernel::AbstractArray{T, 2}`: A 2D array of type `T`, with dimensions `(n_bins, nχ)`. Stores the kernel values for each tomographic bin and a grid of χ values.

# Constructors
- `ShearKernel{T}(n_bins::Int, nχ::Int)`: Creates a `ShearKernel` with the specified `n_bins` and `nχ` values, initializing the kernel values to zeros of type `T`.
- `ShearKernel(n_bins::Int, nχ::Int)`: Creates a `ShearKernel` with the specified `n_bins` and `nχ` values, initializing the kernel values to zeros of type `Float64`.
"""
@kwdef mutable struct ShearKernel{T} <: AbstractOLDCosmologicalProbes{T}
    Kernel::AbstractArray{T, 2} = zeros(1, 1)
end

ShearKernel{T}(n_bins::Int, nχ::Int) where T = ShearKernel{T}(Kernel = zeros(T, n_bins, nχ))

ShearKernel(n_bins::Int, nχ::Int) = ShearKernel{Float64}(n_bins, nχ)

#TODO: missing docs
@kwdef mutable struct RSDKernel{T} <: AbstractOLDCosmologicalProbes{T}
    Kernel::AbstractArray{T, 2} = zeros(1, 1)
end

RSDKernel{T}(n_bins::Int, nχ::Int) where T = RSDKernel{T}(Kernel = zeros(T, n_bins, nχ))

RSDKernel(n_bins::Int, nχ::Int) = RSDKernel{Float64}(n_bins, nχ)

@kwdef mutable struct LensMagKernel{T} <: AbstractOLDCosmologicalProbes{T}
    Kernel::AbstractArray{T, 2} = zeros(1, 1)
end

LensMagKernel{T}(n_bins::Int, nχ::Int) where T = LensMagKernel{T}(Kernel = zeros(T, n_bins, nχ))

LensMagKernel(n_bins::Int, nχ::Int) = LensMagKernel{Float64}(n_bins, nχ)

"""
    CMBLensingKernel{T} <: AbstractCosmologicalProbes{T}

Represents a CMB lensing kernel.

# Parameters
- `Kernel::AbstractArray{T, 1}`: A 1D array of type `T`, with dimension `(nχ)`. Note that CMB Lensing by definition only has a single tomographic bin.

# Constructors
- `CMBLensingKernel{T}(nχ::Int)`: Creates a `CMBLensingKernel` with the specified `nχ` value, initializing the kernel values to zeros of type `T`.
- `CMBLensingKernel(nχ::Int)`: Creates a `CMBLensingKernel` with the specified `nχ` value, initializing the kernel values to zeros of type `Float64`.
"""
@kwdef mutable struct CMBLensingKernel{T} <: AbstractOLDCosmologicalProbes{T}
    Kernel::AbstractArray{T, 1} = zeros(1)
end

CMBLensingKernel{T}(nχ::Int) where T = CMBLensingKernel{T}(Kernel = zeros(T, nχ))

CMBLensingKernel(nχ::Int) = CMBLensingKernel{Float64}(nχ)

function compute_kernel!(nz::AbstractArray{T, 2}, z::AbstractArray{T, 1}, Probe::GalaxyKernel, 
                        grid::CosmologicalGrid, bg::BackgroundQuantities, 
                        cosmo::AbstractCosmology) where T

    #TODO: this test will suck for autodiff, will need fixing
    if all(iszero, bg.Hz_array) || all(iszero, bg.χz_array)
        evaluate_background_quantities!(grid, bg, cosmo)
    end

    n_bins = size(Probe.Kernel, 1)
    
    for b in 1:n_bins
        nz_func = DataInterpolations.AkimaInterpolation(nz[b,:], z, extrapolation=ExtrapolationType.Extension)
        nz_norm, _ = quadgk(x->nz_func(x), first(grid.z_range), last(grid.z_range))

        Probe.Kernel[b,:] = @. (bg.Hz_array / C_LIGHT) * (nz_func.(grid.z_range) / nz_norm)
    end
end

function compute_kernel!(nz::AbstractArray{T, 2}, z::AbstractArray{T, 1}, Probe::GalaxyKernel, 
    bias::AbstractArray{T,2}, grid::CosmologicalGrid, bg::BackgroundQuantities, 
    cosmo::AbstractCosmology) where T

    #TODO: this test will suck for autodiff, will need fixing
    if all(iszero, bg.Hz_array) || all(iszero, bg.χz_array)
    evaluate_background_quantities!(grid, bg, cosmo)
    end

    n_bins = size(Probe.Kernel, 1)

    for b in 1:n_bins
        nz_func = DataInterpolations.AkimaInterpolation(nz[b,:], z, extrapolation = ExtrapolationType.Extension)
        nz_norm, _ = quadgk(x->nz_func(x), first(grid.z_range), last(grid.z_range))

        Probe.Kernel[b,:] = @. bias[b,:] * (bg.Hz_array / C_LIGHT) * (nz_func.(grid.z_range) / nz_norm)
    end
end


function compute_kernel!(nz::AbstractArray{T, 2}, z::AbstractArray{T, 1}, Probe::ShearKernel, grid::CosmologicalGrid,
    bg::BackgroundQuantities, cosmo::AbstractCosmology) where T

    #TODO: this test will suck for autodiff, will need fixing
    if all(iszero, bg.Hz_array) || all(iszero, bg.χz_array)
        evaluate_background_quantities!(grid, bg, cosmo)
    end

    n_bins = size(Probe.Kernel, 1)

    for b in 1:n_bins
        nz_func = DataInterpolations.AkimaInterpolation(nz[b,:], z, extrapolation=ExtrapolationType.Extension)
        nz_norm, _ = quadgk(x->nz_func(x), first(grid.z_range), last(grid.z_range))

        prefac = 1.5 * cosmo.H0^2 * cosmo.Ωm / C_LIGHT^2

        for z_idx in 1:length(grid.z_range)
            integrand(x) = nz_func(x) * (1. - bg.χz_array[z_idx]/compute_χ(x, cosmo))
            z_low = grid.z_range[z_idx]
            z_top = grid.z_range[end]*1.1 #TODO: check max redshift, with n5k bins, lensing5 fallisce se uso valore diverso da 3.5
            int, err = quadgk(x -> integrand(x), z_low, z_top) 

            Probe.Kernel[b, z_idx] = prefac * bg.χz_array[z_idx] * (1. + grid.z_range[z_idx]) * int / nz_norm
        end
    end
end


function compute_kernel!(Probe::CMBLensingKernel, grid::CosmologicalGrid, 
    bg::BackgroundQuantities, cosmo::AbstractCosmology) 

    #TODO: this test will suck for autodiff, will need fixing
    if all(iszero, bg.Hz_array) || all(iszero, bg.χz_array)
        evaluate_background_quantities!(grid, bg, cosmo)
    end

    prefac = 1.5 * cosmo.H0^2 * cosmo.Ωm / C_LIGHT^2
    χ_CMB = compute_χ(1100., cosmo)

    Probe.Kernel = @. prefac * bg.χz_array * (1. + grid.z_range) * (1 - bg.χz_array/χ_CMB)
end

function compute_kernel!(nz::AbstractArray{T, 2}, z::AbstractArray{T, 1}, Probe::RSDKernel, 
    growth_factor::AbstractArray{T,1}, grid::CosmologicalGrid,  bg::BackgroundQuantities) where T

    n_bins = size(Probe.Kernel, 1)

    for b in 1:n_bins
        nz_func = DataInterpolations.AkimaInterpolation(nz[b,:], z, extrapolation=ExtrapolationType.Extension)
        nz_norm, _ = quadgk(x->nz_func(x), first(grid.z_range), last(grid.z_range))

        Probe.Kernel[b,:] = @. growth_factor * (bg.Hz_array / C_LIGHT) * (nz_func.(grid.z_range) / nz_norm) #TODO: might be missing C factors
    end
end


function compute_kernel!(nz::AbstractArray{T, 2}, z::AbstractArray{T, 1}, Probe::LensMagKernel, s::AbstractArray{T, 2}, grid::CosmologicalGrid,
    bg::BackgroundQuantities, cosmo::AbstractCosmology) where T

    #TODO: this test will suck for autodiff, will need fixing
    if all(iszero, bg.Hz_array) || all(iszero, bg.χz_array)
        evaluate_background_quantities!(grid, bg, cosmo)
    end

    n_bins = size(Probe.Kernel, 1)

    for b in 1:n_bins
        nz_func = DataInterpolations.AkimaInterpolation(nz[b,:], z, extrapolation=ExtrapolationType.Extension)
        nz_norm, _ = quadgk(x->nz_func(x), first(grid.z_range), last(grid.z_range))

        s_z = DataInterpolations.AkimaInterpolation(s[b,:], z, extrapolation=ExtrapolationType.Extension)

        prefac = 1.5 * cosmo.H0^2 * cosmo.Ωm / C_LIGHT^2

        for z_idx in 1:length(grid.z_range)
            integrand(x) = nz_func(x) * (1. - bg.χz_array[z_idx]/compute_χ(x, cosmo)) * (5 .* s_z(x) .- 2)
            z_low = grid.z_range[z_idx]
            z_top =  grid.z_range[end]*1.1 #TODO: check max redshift, with n5k bins, lensing5 fallisce se uso valore diverso da 3.5
            int, err = quadgk(x -> integrand(x), z_low, z_top) 

            Probe.Kernel[b, z_idx] = prefac * bg.χz_array[z_idx] * (1. + grid.z_range[z_idx]) * int / nz_norm
        end
    end
end

function P_phi(k::AbstractArray{<:Any,1}, cosmo::AbstractCosmology)
    return @. 9/25 * 2 * π^2 * cosmo.As / (k^3) * (k/0.05)^(cosmo.ns - 1.)
end

function extract_transfer_function(pk::AbstractArray{<:Any,2}, k::AbstractArray{<:Any, 1}, cosmo::AbstractCosmology)
    prim_pk = get_PΦ(k , cosmo)
    return sqrt.(pk ./ reshape(prim_pk, 1, :))
end

function old_to_χR_frame(matrix::AbstractArray{<:Any,2}, plan::FFTW.r2rFFTWPlan, bg::BackgroundQuantities, grid::AbstractCosmologicalGrid)
    coefs = fast_chebcoefs(matrix, plan)
    new_χs = make_grid(bg, R)
    x = resample_redshifts(bg, grid, new_χs)
    chebyshevs = chebyshev_polynomials(x, z_cheb)
    pk_interp = coefs' * chebyshevs  
    return reshape( pk_interp, size(k_cheb,1),  size(bg.χz_array, 1), size(R,1) )
end

"""
    chebyshev_polynomials(x::AbstractArray{T,1}, n_cheb::Int, z_min::T, z_max::T) where T

Computes the Chebyshev polynomials ``T_n(x)`` up to a specified order for a given range of `x`.

# Arguments
- `x::AbstractArray{T,1}`: An array of input values for which the Chebyshev polynomials will be evaluated.
- `n_cheb::Int`: The maximum order of Chebyshev polynomials to compute.
- `z_min::T`: The minimum value in the domain of `x`.
- `z_max::T`: The maximum value in the domain of `x`.

# Returns
A 2D array where each row corresponds to a Chebyshev polynomial ``T_n(x)``

# Notes
- Scales `x` to the Chebyshev domain ``[-1, 1]``.
- Recurrence relation:
  - ``T_0(x) = 1```
  - ``T_1(x) = x``
  - `` T_{n+1}(x) = 2xT_n(x) - T_{n-1}(x)`` for `` n \\geq 2``.
"""
function chebyshev_polynomials( x::AbstractArray{T,1}, n_cheb::Int, z_min::T, z_max::T) where T
    x_scaled = 2 .* (x .- z_min) ./ (z_max - z_min) .- 1
    all(abs.(x_scaled) .≤ 1)
    
    Tcheb = zeros(n_cheb, length(x_scaled))
    
    Tcheb[1, :] .= 1.0  # T0(x) = 1
    if n_cheb >= 2
        Tcheb[2, :] .= x_scaled  # T1(x) = x
    end
    
    for n in 2:n_cheb-1
        Tcheb[n+1, :] .= 2 .* x_scaled .* Tcheb[n, :] .- Tcheb[n-1, :]
    end
    
    return Tcheb
end

function chebyshev_polynomials( x::AbstractArray, cheb_nodes::AbstractArray )
    
    n_cheb = length(cheb_nodes)
    z_min = minimum(cheb_nodes)
    z_max = maximum(cheb_nodes)
    #AssertionError(maximum(x)<=z_max)

    app = LinRange(-1, 1, 1000)

    Tcheb = zeros(n_cheb, length(app))
    
    Tcheb[1, :] .= 1.0  # T0(x) = 1
    if n_cheb >= 2
        Tcheb[2, :] .= app  # T1(x) = x
    end
    
    for n in 2:n_cheb-1
        Tcheb[n+1, :] .= 2 .* app .* Tcheb[n, :] .- Tcheb[n-1, :]
    end

    T_cheb_return = zeros(n_cheb, length(x))
    x_scaled = 2 .* (x .-z_min) ./ (z_max - z_min) .- 1

    for i in 1:n_cheb
        interp = AkimaInterpolation(Tcheb[i,:], app, extrapolation=ExtrapolationType.Extension)
        T_cheb_return[i,:] = interp.(x_scaled)
    end
    
    return T_cheb_return
end

function chebyshev_frigo( x::AbstractArray, cheb_nodes::AbstractArray)

    n_cheb = length(cheb_nodes)
    Tcheb = zeros(n_cheb, length(x))

    c = FastChebInterp.ChebPoly(cheb_nodes, SA[minimum(cheb_nodes)], SA[maximum(cheb_nodes)])

    for i in 1:n_cheb
        copy_c = deepcopy(c) 
        copy_c.coefs .*= 0 
        copy_c.coefs[i] = 1.
        Tcheb[i,:] = copy_c.(x)
    end
    
    return Tcheb
end

"""
    interpolate_power_spectrum(pk::AbstractArray{T,2}, z_nodes::AbstractArray{T,1}, 
                               R::AbstractArray{T,1}, plan::FFTW.r2rFFTWPlan, 
                               bg::BackgroundQuantities, grid::AbstractCosmologicalGrid) where T

Interpolates the power spectrum  `P(z,k)` to put it on the ``\\chi-R`` grid optimal for the algorithm.
Returns the object ``P(k, \\chi, R)``. 

# Arguments
- `pk::AbstractArray{T,2}`: A 2D array of power spectrum values. The function expects the first axis to be `z`, and the second one to be `k`.
- `z_nodes::AbstractArray{T,1}`: Redshift values corresponding to the first axis of `pk`.
- `R::AbstractArray{T,1}`: Values of ``R \\equiv \\chi_2/\\chi_1``.
- `plan::FFTW.r2rFFTWPlan`: FFTW plan for computing Chebyshev coefficients.
- `bg::BackgroundQuantities`: Background cosmological quantities. Contains the comoving distance values.
- `grid::AbstractCosmologicalGrid`: Grid of cosmological quantities. Contains the redshift grid.

# Returns
A 3D array of interpolated power spectrum values with dimensions ``(k, \\chi, R)``
"""
function interpolate_power_spectrum(pk::AbstractArray{T,2}, z_nodes::AbstractArray{T,1}, R::AbstractArray{T,1},
    plan::FFTW.r2rFFTWPlan, bg::BackgroundQuantities, grid::AbstractCosmologicalGrid) where T

    coefs = fast_chebcoefs(pk, plan)
    new_χs = make_grid(bg, R)
    x = resample_redshifts(bg, grid, new_χs)
    chebyshevs = chebyshev_polynomials(x, z_nodes)
    pk_interp = coefs' * chebyshevs  #TODO: understand how to handle pk sizes
    return reshape(pk_interp, size(pk,2),  length(bg.χz_array), length(R))
end

"""
    unequal_time_power_spectrum(pk::AbstractArray{T,3}) where T

Takes in input the power spectrum on the ``(k, \\chi, R)`` grid and implements the equation: 
```math
P(k,\\chi, R\\chi)=\\sqrt{P(k,\\chi)P(k,R\\chi)},
``` 
which assumes that the quantities involved are perfectly correlated at different cosmic times.

# Arguments
- `pk::AbstractArray{T,3}`: A 3D array of power spectrum values on a grid ``(k, \\chi, R).``

# Returns
A 3D array with the same dimensions as `pk`.
"""
function unequal_time_power_spectrum(pk::AbstractArray{T,3}) where T
    pk_R1 = pk[:,:,end]
    @tullio final_pk[i,c,r] := sqrt(pk_R1[i,c] * pk[i,c,r])
    return final_pk
end

function Pmm_unequaltime(pk::AbstractArray{T,2}, k::AbstractArray{T,1}, z::AbstractArray{T,1}, R::AbstractArray{T,1}, bg::BackgroundQuantities, grid::AbstractCosmologicalGrid, cosmo::AbstractCosmology) where T
    primordial_pk = Blast.P_phi(k, cosmo)
    T_m = Blast.extract_transfer_function(pk, k, cosmo)
    plan = Blast.plan_fft(log10.(T_m), 1)
    T_m_interp = 10 .^ (Blast.interpolate_power_spectrum(log10.(T_m), z, R, plan, bg, grid))

    T_m_interp_R1 = T_m_interp[:,:,end]
    @tullio Pmm_unequaltime[k, i, j] := primordial_pk[k] * T_m_interp_R1[k,i] * T_m_interp[k, i, j]

    return Pmm_unequaltime
end

function Pgg_unequaltime(bias_kz::AbstractArray{T,2}, pk::AbstractArray{T,2}, k::AbstractArray{T,1}, z::AbstractArray{T,1}, R::AbstractArray{T,1}, bg::BackgroundQuantities, grid::AbstractCosmologicalGrid, cosmo::AbstractCosmology) where T
    primordial_pk = Blast.P_phi(k, cosmo)
    T_m = Blast.extract_transfer_function(pk, k, cosmo)
    plan = Blast.plan_fft(log10.(T_m), 1)
    T_m_interp = 10 .^ (Blast.interpolate_power_spectrum(log10.(T_m), z, R, plan, bg, grid))
    plan = Blast.plan_fft(bias_kz, 1)
    bias_interp = Blast.interpolate_power_spectrum(bias_kz, z, R, plan, bg, grid)

    T_m_interp_R1 = T_m_interp[:,:,end]
    bias_interp_R1 = bias_interp[:,:,end]

    return @tullio Pgg_unequaltime[k, i, j] := primordial_pk[k] * bias_interp_R1[k,i] * T_m_interp_R1[k,i] * bias_interp[k,i,j] * T_m_interp[k,i,j]
end


function Pgm_unequaltime(bias_kz::AbstractArray{T,2}, pk::AbstractArray{T,2}, k::AbstractArray{T,1}, z::AbstractArray{T,1}, R::AbstractArray{T,1}, bg::BackgroundQuantities, grid::AbstractCosmologicalGrid, cosmo::AbstractCosmology) where T
    primordial_pk = Blast.P_phi(k, cosmo)
    T_m = Blast.extract_transfer_function(pk, k, cosmo)
    plan = Blast.plan_fft(log10.(T_m), 1)
    T_m_interp = 10 .^ (Blast.interpolate_power_spectrum(log10.(T_m), z, R, plan, bg, grid))
    plan = Blast.plan_fft(bias_kz, 1)
    bias_interp = Blast.interpolate_power_spectrum(bias_kz, z, R, plan, bg, grid)

    T_m_interp_R1 = T_m_interp[:,:,end]
    bias_interp_R1 = bias_interp[:,:,end]

    return @tullio Pgm_unequaltime[k, i, j] := primordial_pk[k] * bias_interp_R1[k,i] * T_m_interp_R1[k,i] * T_m_interp[k,i,j]
end

function stoopid_2D_interpolator(x::AbstractArray{T,1}, y::AbstractArray{T,1}, f::AbstractArray{T,2}, logx::Bool, logy::Bool, logf::Bool ) where T
    if logx
        x_resc = LinRange(log10(first(x)),log10(last(x)), length(x))
    else
        x_resc = LinRange(first(x),last(x), length(x))
    end 

    if logy
        y_resc = LinRange(log10(first(y)),log10(last(y)), length(y))
    else
        y_resc = LinRange(first(y),last(y), length(y))
    end

    if logf 
        Interpolator = Interpolations.interpolate(log10.(f),BSpline(Cubic(Line(OnGrid()))))
        Interpolator = scale(Interpolator, (x_resc, y_resc))
        Interpolator = Interpolations.extrapolate(Interpolator, Line());
    else
        Interpolator = Interpolations.interpolate(f,BSpline(Cubic(Line(OnGrid()))))
        Interpolator = scale(Interpolator, (x_resc, y_resc))
        Interpolator = Interpolations.extrapolate(Interpolator, Line());
    end

    return Interpolator
end

function stoopid_unequal_time_pk(interpolator, k::AbstractArray{T,1}, z1::Number, z2::Number) where T
    #TODO: this only works if the interpolator is created with false true true
    return @. sqrt(10^interpolator(z1,log10(k)) * 10^interpolator(z2,log10(k)))
end

function limber_rsd_kernel(ℓ::Number, bg::BackgroundQuantities, RSDK::Blast.RSDKernel)
    χ = bg.χz_array
    nbins = size(RSDK.Kernel)[1]
    rds_kernels = zeros( nbins, length(χ) )

    for b in 1:nbins
        kernel_interp = AkimaInterpolation(RSDK.Kernel[b, :], χ, extrapolation=ExtrapolationType.Extension)
        piece1 = @. (2*ℓ^2 + 2*ℓ - 1) / ((2*ℓ - 1)*(2*ℓ + 3)) * RSDK.Kernel[b, :]
        piece2 = @. (ℓ - 1)*ℓ / ((2*ℓ - 1) * sqrt(2*ℓ - 3)*(2*ℓ + 1)) * kernel_interp.((2*ℓ - 3)/(2*ℓ + 1) * χ)
        piece3 = @. (ℓ + 1)*(ℓ + 2) / ((2*ℓ + 3) * sqrt((2*ℓ + 1)*(2*ℓ + 5))) * kernel_interp.((2*ℓ + 5)/(2*ℓ + 1) * χ)
        rds_kernels[b, :] .= piece1 .- piece2 .- piece3
    end

    return rds_kernels
end

function bΦ(z, bias_model, p::Number)
    return 2 * 1.686 * (bias_model(z) .- p)
end

function scale_dependent_bias(f_NL::Number, pk::AbstractArray{T,2}, k::AbstractArray{T,1}, z::AbstractArray{T,1}, bias_model, p::Number, cosmo::AbstractCosmology) where T
    transfer_func = Blast.extract_transfer_function(pk, k, cosmo)
    return (bΦ(z, bias_model, p) .* f_NL) ./ transfer_func
end


## old integrals.jl dropped here for now
"""
    make_grid(BackgroundQuantities::BackgroundQuantities, R::Vector{T}) where T

Constructs a grid by multiplying the `χz_array` from `BackgroundQuantities` with the vector `R`.

# Arguments
- `BackgroundQuantities::BackgroundQuantities`: An instance of the `BackgroundQuantities` type that contains the `χz_array`.
- `R::Vector{T}`: A vector of values to be used in the grid construction, where `T` can be any type.
"""
function make_grid(bg::BackgroundQuantities, R::Vector{T}) where T
    return vec(bg.χz_array * R')
end

"""
    grid_interpolator(Probe::Union{GalaxyKernel, ShearKernel}, 
        BackgroundQuantities::BackgroundQuantities, grid::Vector{T}) where T

Interpolates the kernel values for a given grid based on the specified cosmological probes. 
Returns a 2D array of interpolated kernel values, where rows correspond to the number of bins and columns correspond to the grid points.

# Arguments
- `Probe::Union{GalaxyKernel, ShearKernel}`: The kernel data to be interpolated.
- `BackgroundQuantities::BackgroundQuantities`: An instance of the `BackgroundQuantities` type that contains the `χz_array`.
- `grid::Vector{T}`: A vector of values where the interpolated kernel values will be evaluated.
"""
function grid_interpolator(Probe::Union{GalaxyKernel, ShearKernel, RSDKernel, LensMagKernel}, 
    bg::BackgroundQuantities, grid::Vector{T}) where T

    n_bins = size(Probe.Kernel, 1)
    kernel_interpolated = zeros(n_bins, length(grid))

    for b in 1:n_bins
        interp = AkimaInterpolation(Probe.Kernel[b,:], bg.χz_array, extrapolation=ExtrapolationType.Extension)
        kernel_interpolated[b, :] = interp.(grid)
    end

    return kernel_interpolated
end

"""
    grid_interpolator(Probe::CMBLensingKernel, 
        bg::BackgroundQuantities, grid::Vector{T}) where T

Interpolates the kernel values for a given grid based on the specified cosmological probes. 
Returns a 2D array of interpolated kernel values, where rows correspond to the number of bins and columns correspond to the grid points.

# Arguments
- `Probe::CMBLensingKernel`: The kernel data to be interpolated.
- `bg::BackgroundQuantities`: An instance of the `BackgroundQuantities` type that contains the `χz_array`.
- `grid::Vector{T}`: A vector of values where the interpolated kernel values will be evaluated.
"""
function grid_interpolator(Probe::CMBLensingKernel, 
    bg::BackgroundQuantities, grid::Vector{T}) where T

    kernel_interpolated = zeros(1, length(grid))

    interp = AkimaInterpolation(Probe.Kernel, bg.χz_array, extrapolation=ExtrapolationType.Extension)
    kernel_interpolated[1, :] = interp.(grid)

    return kernel_interpolated
end

"""
    get_kernel_array(AbstractCosmologicalProbes::GalaxyKernel, 
        BackgroundQuantities::BackgroundQuantities, R::Vector{T}) where T

Obtains the kernel array for the `GalaxyKernel` probe, with dimensions (bins, nχ, nR). 


# Arguments
- `AbstractCosmologicalProbes::GalaxyKernel`: An instance of the `GalaxyKernel` type.
- `BackgroundQuantities::BackgroundQuantities`: An instance of the `BackgroundQuantities` type.
- `R::Vector{T}`: A vector of values for which the kernel array is to be computed.
"""
function get_kernel_array(Probe::GalaxyKernel, 
    bg::BackgroundQuantities, R::Vector{T}) where T

    n_bins = size(Probe.Kernel, 1)
    nχ = size(Probe.Kernel, 2)
    nR = length(R)
    
    W_array = reshape(grid_interpolator(Probe, bg, make_grid(bg, R)), n_bins, nχ, nR)

    return W_array
end

function get_kernel_array(Probe::RSDKernel, 
    bg::BackgroundQuantities, R::Vector{T}) where T

    n_bins = size(Probe.Kernel, 1)
    nχ = size(Probe.Kernel, 2)
    nR = length(R)
    
    W_array = reshape(grid_interpolator(Probe, bg, make_grid(bg, R)), n_bins, nχ, nR)

    return W_array
end

"""
    get_kernel_array(Probe::ShearKernel, 
        bg::BackgroundQuantities, R::Vector{T}) where T

Obtains the kernel array for the weak lensing probe, with dimensions (bins, nχ, nR)
The difference with respect to the galaxy case is that these kernels are divided by ``\\chi^2``.

# Arguments
- `Probe::ShearKernel`: An instance of `ShearKernel`.
- `bg::BackgroundQuantities`: An instance of the `BackgroundQuantities` type.
- `R::Vector{T}`: A vector of values for which the kernel array is to be computed.
"""
function get_kernel_array(Probe::Union{ShearKernel, LensMagKernel}, 
    bg::BackgroundQuantities, R::Vector{T}) where T

    n_bins = size(Probe.Kernel, 1)
    nχ = size(Probe.Kernel, 2)
    nR = length(R)
    
    W_L = grid_interpolator(Probe, bg, make_grid(bg, R))

    χ2_app = zeros(n_bins, nχ*nR)
    for i in 1:n_bins
        χ2_app[i,:] = make_grid(bg, R) .^ 2
    end
    
    W_array = reshape( W_L./χ2_app , n_bins, nχ, nR)

    return W_array
end

"""
    get_kernel_array(Probe::CMBLensingKernel, 
        bg::BackgroundQuantities, R::Vector{T}) where T

Obtains the kernel array for the CMB lensing probe, with dimensions (1, nχ, nR)
Similarly to the weak lensing case, the kernel is divided by ``\\chi^2``.

# Arguments
- `Probe::CMBLensingKernel`: An instance of `CMBLensingKernel`.
- `bg::BackgroundQuantities`: An instance of the `BackgroundQuantities` type.
- `R::Vector{T}`: A vector of values for which the kernel array is to be computed.
"""
function get_kernel_array(Probe::CMBLensingKernel, 
    bg::BackgroundQuantities, R::Vector{T}) where T

    nχ = size(Probe.Kernel, 1)
    nR = length(R)
    
    W_L = grid_interpolator(Probe, bg, make_grid(bg, R))

    W_L[1,:] = W_L[1,:] ./ (make_grid(bg, R) .^ 2)
    
    W_array = reshape(W_L, 1, nχ, nR)

    return W_array
end

"""
    combine_kernels(ProbeA::Union{GalaxyKernel, ShearKernel, CMBLensingKernel}, 
        ProbeB::Union{GalaxyKernel, ShearKernel, CMBLensingKernel}, 
        BackgroundQuantities::BackgroundQuantities, R::Vector{T}) where T

Combines the kernels from two different cosmological probes into a single array. 
This is needed to perform the integration in the χ-R coordinates.

# Arguments
- `ProbeA::Union{GalaxyKernel, ShearKernel, CMBLensingKernel}`: The first cosmological probe to combine.
- `ProbeB::Union{GalaxyKernel, ShearKernel, CMBLensingKernel}`: The second cosmological probe to combine.
- `BackgroundQuantities::BackgroundQuantities`: An instance of the `BackgroundQuantities` type.
- `R::Vector{T}`: A vector of values used in the combination.
"""
function combine_kernels(ProbeA::Union{GalaxyKernel, ShearKernel, CMBLensingKernel, RSDKernel, LensMagKernel}, 
    ProbeB::Union{GalaxyKernel, ShearKernel, CMBLensingKernel, RSDKernel, LensMagKernel},
    bg::BackgroundQuantities, R::Vector{T}) where T

    W_A = get_kernel_array(ProbeA, bg, R)
    W_A_r1 = W_A[:,:,end]

    W_B = get_kernel_array(ProbeB, bg, R)
    W_B_r1 = W_B[:,:,end]

    @tullio K[i,j,c,r] := W_A_r1[i,c] * W_B[j,c,r] + W_A[i,c,r]*W_B_r1[j,c]

    return K
end


"""
    get_ell_prefactor(ProbeA::GalaxyKernel, ProbeB::GalaxyKernel, ℓ_list::Vector)

Calculates the prefactor for the angular power spectrum when both probes are `GalaxyKernel`.

# Arguments
- `ℓ_list::Vector`: A vector of angular multipole values.
"""
function get_ell_prefactor(ProbeA::GalaxyKernel, ProbeB::GalaxyKernel, ℓ_list::Vector)
    return 2 / π * ones(length(ℓ_list))
end

function get_ell_prefactor(ProbeA::RSDKernel, ProbeB::RSDKernel, ℓ_list::Vector)
    return 2 / π * ones(length(ℓ_list))
end

"""
    get_ell_prefactor(ProbeA::GalaxyKernel, ProbeB::ShearKernel, ℓ_list::Vector)

Calculates the prefactor for the angular power spectrum when probes are `GalaxyKernel` and `ShearKernel`.

# Arguments
- `ℓ_list::Vector`: A vector of angular multipole values.
"""
function get_ell_prefactor(ProbeA::GalaxyKernel, ProbeB::ShearKernel, ℓ_list::Vector)
    return 2 / π * sqrt.(factorial_frac(ℓ_list))
end

# Define the same prefactor for ShearKernel and GalaxyKernel order
get_ell_prefactor(ProbeA::ShearKernel, ProbeB::GalaxyKernel, ℓ_list::Vector) = 
    get_ell_prefactor(ProbeB, ProbeA, ℓ_list)

"""
    get_ell_prefactor(ProbeA::ShearKernel, ProbeB::ShearKernel, ℓ_list::Vector)

Calculates the prefactor for the angular power spectrum when both probes are `ShearKernel`.

# Arguments
- `ℓ_list::Vector`: A vector of angular multipole values.
"""
function get_ell_prefactor(ProbeA::ShearKernel, ProbeB::ShearKernel, ℓ_list::Vector)
    return 2 / π * factorial_frac(ℓ_list)
end

"""
    get_ell_prefactor(ProbeA::CMBLensingKernel, ProbeB::ShearKernel, ℓ_list::Vector)

Calculates the prefactor for the angular power spectrum when probes are `CMBLensingKernel` and `ShearKernel`.

# Arguments
- `ℓ_list::Vector`: A vector of angular multipole values.
"""
function get_ell_prefactor(ProbeA::CMBLensingKernel, ProbeB::ShearKernel, ℓ_list::Vector)
    return @. 2 * ℓ_list * (ℓ_list +1) * sqrt(factorial_frac(ℓ_list)) / π #TODO: double check  prefactor
end

# Define the same prefactor for ShearKernel and GalaxyKernel order
get_ell_prefactor(ProbeA::ShearKernel, ProbeB::CMBLensingKernel, ℓ_list::Vector) = 
    get_ell_prefactor(ProbeB, ProbeA, ℓ_list)

"""
    get_ell_prefactor(ProbeA::CMBLensingKernel, ProbeB::CMBLensingKernel, ℓ_list::Vector)

Calculates the prefactor for the angular power spectrum when both probes are `CMBLensingKernel`.

# Arguments
- `ℓ_list::Vector`: A vector of angular multipole values.
"""
function get_ell_prefactor(ProbeA::CMBLensingKernel, ProbeB::CMBLensingKernel, ℓ_list::Vector)
    return @. π * ℓ_list^2 * (ℓ_list + 1)^2  #TODO: double check  prefactor
end

"""
    get_ell_prefactor(ProbeA::CMBLensingKernel, ProbeB::GalaxyKernel, ℓ_list::Vector)

Calculates the prefactor for the angular power spectrum when probes are `CMBLensingKernel` and `GalaxyKernel`.

# Arguments
- `ℓ_list::Vector`: A vector of angular multipole values.
"""
function get_ell_prefactor(ProbeA::CMBLensingKernel, ProbeB::GalaxyKernel, ℓ_list::Vector)
    return @. 2/π * ℓ_list * (ℓ_list + 1)  #TODO: double check prefactor
end

# Define the same prefactor for ShearKernel and GalaxyKernel order
get_ell_prefactor(ProbeA::GalaxyKernel, ProbeB::CMBLensingKernel, ℓ_list::Vector) = 
    get_ell_prefactor(ProbeB, ProbeA, ℓ_list)

function get_ell_prefactor(ProbeA::CMBLensingKernel, ProbeB::RSDKernel, ℓ_list::Vector)
    return @. 2 * ℓ_list * (ℓ_list + 1) / π 
end

get_ell_prefactor(ProbeA::RSDKernel, ProbeB::CMBLensingKernel, ℓ_list::Vector) = 
    get_ell_prefactor(ProbeB, ProbeA, ℓ_list)

function get_ell_prefactor(ProbeA::GalaxyKernel, ProbeB::RSDKernel, ℓ_list::Vector)
    return  2 / π * ones(length(ℓ_list)) 
end
    
get_ell_prefactor(ProbeA::RSDKernel, ProbeB::GalaxyKernel, ℓ_list::Vector) = 
    get_ell_prefactor(ProbeB, ProbeA, ℓ_list)

function get_ell_prefactor(ProbeA::ShearKernel, ProbeB::RSDKernel, ℓ_list::Vector)
    return  2 / π * sqrt.(factorial_frac(ℓ_list)) 
end
        
get_ell_prefactor(ProbeA::RSDKernel, ProbeB::ShearKernel, ℓ_list::Vector) = 
    get_ell_prefactor(ProbeB, ProbeA, ℓ_list)

function get_ell_prefactor(ProbeA::CMBLensingKernel, ProbeB::LensMagKernel, ℓ_list::Vector)
    return @.  2 * (ℓ_list * (ℓ_list + 1))^2 / π 
end

get_ell_prefactor(ProbeA::LensMagKernel, ProbeB::CMBLensingKernel, ℓ_list::Vector) = 
    get_ell_prefactor(ProbeB, ProbeA, ℓ_list)

function get_ell_prefactor(ProbeA::RSDKernel, ProbeB::LensMagKernel, ℓ_list::Vector)
    return @. 2 / π * ℓ_list * (ℓ_list + 1)
end
    
get_ell_prefactor(ProbeA::LensMagKernel, ProbeB::RSDKernel, ℓ_list::Vector) = 
    get_ell_prefactor(ProbeB, ProbeA, ℓ_list)

function get_ell_prefactor(ProbeA::GalaxyKernel, ProbeB::LensMagKernel, ℓ_list::Vector)
    return @. 2 / π * ℓ_list * (ℓ_list + 1) 
end
        
get_ell_prefactor(ProbeA::LensMagKernel, ProbeB::GalaxyKernel, ℓ_list::Vector) = 
    get_ell_prefactor(ProbeB, ProbeA, ℓ_list)


function get_ell_prefactor(ProbeA::LensMagKernel, ProbeB::LensMagKernel, ℓ_list::Vector)
    return @. 2 / π * (ℓ_list * (ℓ_list + 1))^2  
end

function get_ell_prefactor(ProbeA::ShearKernel, ProbeB::LensMagKernel, ℓ_list::Vector)
    return @. 2 / π * ℓ_list * (ℓ_list + 1) * sqrt(factorial_frac(ℓ_list))
end

get_ell_prefactor(ProbeA::LensMagKernel, ProbeB::ShearKernel, ℓ_list::Vector) = 
    get_ell_prefactor(ProbeB, ProbeA, ℓ_list)


"""
    compute_Cℓ(w::AbstractArray{T, 3}, 
               ProbeA::Union{GalaxyKernel, ShearKernel, CMBLensingKernel}, 
               ProbeB::Union{GalaxyKernel, ShearKernel, CMBLensingKernel}, 
               BackgroundQuantities::BackgroundQuantities, 
               R::AbstractVector, 
               ℓ_list::AbstractArray{T,1} = Blast.ℓ) where T

Computes the angular power spectrum `Cℓ` by performing two outer integrals over `χ` and `R`. 
The Simpson quadrature rule is used for integration over `χ`, while Clenshaw-Curtis quadrature is used for `R`.

# Arguments
- `w::AbstractArray{T, 3}`: A 3D array representing the projected matter densities, containing the inner integrals over `k`. The array dimensions are (ℓ, χ, R).
- `ProbeA::Union{GalaxyKernel, ShearKernel, CMBLensingKernel}`: The first cosmological probe kernel type.
- `ProbeB::Union{GalaxyKernel, ShearKernel, CMBLensingKernel}`: The second cosmological probe kernel type.
- `BackgroundQuantities::BackgroundQuantities`: Contains background information, including `χz_array` for distances in the cosmology.
- `R::AbstractVector`: A 1D array representing the radial grid values for integration.
- `ℓ_list::AbstractArray{T,1}`: An optional 1D array of angular multipole values, defaulting to the global variable Blast.ℓ, the list of ℓ values used for the precomputed part.

# Returns
- A 3D array `Cℓ` with dimensions (ℓ, i, j), where `i` and `j` represent the tomographic bins. The array contains the computed angular power spectrum coefficients for each combination of `ℓ` values and tomographic bins.
"""
function compute_Cℓ(w::AbstractArray{T, 3}, ProbeA::Union{GalaxyKernel, ShearKernel, CMBLensingKernel, RSDKernel, LensMagKernel}, 
    ProbeB::Union{GalaxyKernel, ShearKernel, CMBLensingKernel, RSDKernel, LensMagKernel}, bkgq::BackgroundQuantities, R::AbstractVector, ℓ_list::AbstractArray{T,1} = Blast.ℓ_nonlimber) where T

    nχ = length(bkgq.χz_array)
    nR = length(R)

    K = combine_kernels(ProbeA, ProbeB, bkgq, R)

    #Integration in χ is performed using the Simpson quadrature rule
    Δχ = ((last(bkgq.χz_array)-first(bkgq.χz_array))/(nχ-1))
    w_χ = simpson_weight_array(nχ)

    #Integration in R is performed using the Clenshaw-Curtis quadrature rule
    w_R = get_clencurt_weights_R_integration(2*nR+1)

    ell_prefactor = get_ell_prefactor(ProbeA, ProbeB, ℓ_list)

    @tullio Cℓ[l,i,j] := ell_prefactor[l]*bkgq.χz_array[n]*K[i,j,n,m]*w[l,n,m]*w_χ[n]*w_R[m]*Δχ

    return Cℓ 
end

"""
    compute_Cℓ(w::AbstractArray{T, 3}, 
               K::AbstractArray{T, 4}, 
               bkgq::BackgroundQuantities, 
               weights_χ::AbstractArray{T, 1}, 
               weights_R::AbstractArray{T, 1}, 
               ell_prefactor::AbstractArray{T,1}) where T

Computes the angular power spectrum `Cℓ` by performing integrals over `χ` and `R` for a given set of inputs.

# Arguments
- `w::AbstractArray{T, 3}`: A 3D array representing the projected matter densities, containing the inner integrals over `k`. The array dimensions are (ℓ, χ, R).
- `K::AbstractArray{T, 4}`: A 4D array representing the kernel product `K[i, j, χ, R]` for the given cosmological probes. It combines the effects of the probes and tomographic bin combinations.
- `bkgq::BackgroundQuantities`: Contains background information, including `χz_array` for distances in the cosmology.
- `weights_χ::AbstractArray{T, 1}`: A 1D array of weights for Simpson quadrature integration over `χ`.
- `weights_R::AbstractArray{T, 1}`: A 1D array of weights for Clenshaw-Curtis quadrature integration over `R`.
- `ell_prefactor::AbstractArray{T,1}`: A 1D array of prefactors dependent on angular multipole values `ℓ`, combining with other components to form the power spectrum.

# Returns
- A 3D array `Cℓ` with dimensions (ℓ, i, j), where `i` and `j` represent the tomographic bins. The array contains the computed angular power spectrum coefficients for each combination of `ℓ` values and tomographic bins.
"""
function compute_Cℓ(w::AbstractArray{T, 3}, K::AbstractArray{T, 4}, bkgq::BackgroundQuantities, weights_χ::AbstractArray{T, 1},
                    weights_R::AbstractArray{T, 1}, ell_prefactor::AbstractArray{T,1}) where T
    
    nχ = length(bkgq.χz_array)
    Δχ = ((last(bkgq.χz_array)-first(bkgq.χz_array))/(nχ-1))

    @tullio Cℓ[l,i,j] := ell_prefactor[l]*bkgq.χz_array[n]*K[i,j,n,m]*w[l,n,m]*weights_χ[n]*weights_R[m]*Δχ

    return Cℓ 
end

#TODO: this is not nice, need to think of a better way to use multiple dispatch or someting like that here
function compute_Cℓ_rsd(w_02::AbstractArray{T, 3}, w_20::AbstractArray{T, 3}, ProbeA::Union{GalaxyKernel, ShearKernel, CMBLensingKernel, RSDKernel, LensMagKernel}, 
    ProbeB::RSDKernel, bkgq::BackgroundQuantities, R::AbstractVector, ℓ_list::AbstractArray{T,1} = Blast.ℓ_nonlimber) where T

    nχ = length(bkgq.χz_array)
    nR = length(R)

    W_A = get_kernel_array(ProbeA, bkgq, R)
    W_A_r1 = W_A[:,:,end]

    W_B = get_kernel_array(ProbeB, bkgq, R)
    W_B_r1 = W_B[:,:,end]

    @tullio K[l,i,j,c,r] := W_A_r1[i,c] * W_B[j,c,r] * w_02[l,c,r] + W_A[i,c,r] * W_B_r1[j,c] * w_20[l,c,r]

    #Integration in χ is performed using the Simpson quadrature rule
    Δχ = ((last(bkgq.χz_array)-first(bkgq.χz_array))/(nχ-1))
    w_χ = simpson_weight_array(nχ)

    #Integration in R is performed using the Clenshaw-Curtis quadrature rule
    w_R = get_clencurt_weights_R_integration(2*nR+1)

    ell_prefactor = get_ell_prefactor(ProbeA, ProbeB, ℓ_list)

    @tullio Cℓ[l,i,j] := ell_prefactor[l]*bkgq.χz_array[n]*K[l,i,j,n,m]*w_χ[n]*w_R[m]*Δχ
end
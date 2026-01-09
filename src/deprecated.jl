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
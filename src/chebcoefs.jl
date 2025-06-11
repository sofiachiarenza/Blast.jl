# FUNCTIONS FOR THE DECOMPOSITION OF THE POWER SPECTRUM ON THE BASIS OF THE CHEBYSHEV POLYNOMIALS, CENTRAL PART OF BLAST
"""
    plan_fft(vals::AbstractArray{<:Number, N}, axis::Int)

Create an FFTW real-to-real (R2R) transformation plan for a specified axis of a given multidimensional array `vals`. 
For example, if the `vals` array is the power spectrum P(k,z), one can set `axis=1` and perform the FFT in `k`, or `axis=2` if the FFT is to be performed along `z`.

# Arguments
- `vals::AbstractArray{<:Number, N}`: The input array of any numerical type with `N` dimensions.
- `axis::Int`: The axis along which the FFT transformation will be applied (e.g., `1` for the first axis, `2` for the second axis, etc.).

# Returns
- `p::FFTW.rFFTWPlan`: An FFTW plan object for transforming `vals` with the appropriate real-to-real transformations. This plan can be applied using the `*` operator (e.g., `transformed_vals = p * vals`).

"""
function plan_fft(vals::AbstractArray{<:Number, N}, axis::Int) where {N}
    kind = map(n -> n > 1 ? FFTW.REDFT00 : FFTW.DHT, size(vals)[axis])
    p = FFTW.plan_r2r(deepcopy(vals), kind, [axis]; flags=FFTW.PATIENT, timelimit=Inf)   
                                                                                    
    return p 
end



"""
    fast_chebcoefs(vals::AbstractArray{<:Number,N}, plan::FFTW.r2rFFTWPlan)

Efficiently compute the Chebyshev coefficients of a multidimensional array `vals` using an O(n log n) method. This method leverages FFT-based type-I Discrete Cosine Transform (DCT-I).

Arguments:
- `vals::AbstractArray{<:Number,N}`: A multidimensional array of values for which to compute the Chebyshev coefficients.

- `plan::FFTW.r2rFFTWPlan`: A FFTW plan object for transforming `vals` with the appropriate real to real transformations. This plan is applied using the `*` operator (e.g., `transformed_vals = p * vals`) and performs the DCT of the `vals` array along the first axis.

Returns:
- `coefs`: An array of the same size as `vals`, containing the computed Chebyshev coefficients.
"""
function fast_chebcoefs(vals::AbstractArray, plan::FFTW.r2rFFTWPlan)
    coefs = plan * vals

    s = size(coefs)
    coefs ./= 2*(s[1]-1)
    
    N = length(s)
    coefs[CartesianIndices(ntuple(i -> i == 1 ? (2:s[1]-1) : (1:s[i]), Val{N}()))] *= 2

    return coefs
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

abstract type AbstractCoeff end
abstract type AbstractCoeffComponents end

@kwdef mutable struct NullCoeff <: AbstractCoeffComponents
    coefs::AbstractArray{<:Any, 3} = zeros(1, 1, 1)
end

@kwdef mutable struct cϕTT <: AbstractCoeffComponents
    coefs::AbstractArray{<:Any, 3} = zeros(1, 1, 1)
end

@kwdef mutable struct cϕT <: AbstractCoeffComponents
    coefs::AbstractArray{<:Any, 3} = zeros(1, 1, 1)
end

@kwdef mutable struct cϕ <: AbstractCoeffComponents
    coefs::AbstractArray{<:Any, 1} = zeros(1)
end

@kwdef mutable struct ChebCoefs <: AbstractCoeff
    cϕTT::cϕTT = cϕTT()
    cϕT::Union{cϕT, NullCoeff} = NullCoeff()
    cϕ::Union{cϕ, NullCoeff} = NullCoeff()
end

#TODO: power spectrum should be treated differently, i'm repeating the operations here!! 

function evaluate_coefs!(c::cϕTT, pk::AbstractArray{<:Any, 2}, bg::BackgroundQuantities, grid::AbstractCosmologicalGrid, cosmo::AbstractCosmology)
    primordial_pk = Blast.P_phi( 10 .^ (k_cheb), cosmo)
    T_m = Blast.extract_transfer_function(pk, 10 .^ (k_cheb), cosmo)
    plan = Blast.plan_fft(log10.(T_m), 1)
    T_m_interp = 10 .^ (Blast.interpolate_power_spectrum(log10.(T_m), Blast.z_cheb, Blast.R, plan, bg, grid))

    T_m_interp_R1 = T_m_interp[:,:,end]
    @tullio P_ϕTT[k, i, j] := primordial_pk[k] * T_m_interp_R1[k,i] * T_m_interp[k, i, j]

    plan = Blast.plan_fft(P_ϕTT,1)
    cheb_coeff = Blast.fast_chebcoefs(P_ϕTT, plan)
    c.coefs = permutedims(cheb_coeff, (2,3,1))
end

function evaluate_coefs!(c::cϕT, pk::AbstractArray{<:Any, 2}, bg::BackgroundQuantities, grid::AbstractCosmologicalGrid, cosmo::AbstractCosmology)
    primordial_pk = Blast.P_phi( 10 .^ (k_cheb), cosmo)
    T_m = Blast.extract_transfer_function(pk, 10 .^ (k_cheb), cosmo)
    plan = Blast.plan_fft(log10.(T_m), 1)
    T_m_interp = 10 .^ (Blast.interpolate_power_spectrum(log10.(T_m), Blast.z_cheb, Blast.R, plan, bg, grid))

    @tullio P_ϕT[k, i, j] := primordial_pk[k]* T_m_interp[k, i, j]

    plan = Blast.plan_fft(P_ϕT,1)
    cheb_coeff = Blast.fast_chebcoefs(P_ϕT, plan)
    c.coefs = permutedims(cheb_coeff, (2,3,1))
end

function evaluate_coefs!(c::cϕ, pk::AbstractArray{<:Any, 2}, bg::BackgroundQuantities, grid::AbstractCosmologicalGrid, cosmo::AbstractCosmology)
    primordial_pk = Blast.P_phi( 10 .^ (k_cheb), cosmo)
    
    plan = Blast.plan_fft(primordial_pk,1)
    cheb_coeff = Blast.fast_chebcoefs(primordial_pk, plan)
    c.coefs = cheb_coeff
end

function evaluate_coefs!(c::NullCoeff, pk::AbstractArray{<:Any, 2}, bg::BackgroundQuantities, grid::AbstractCosmologicalGrid, cosmo::AbstractCosmology)
    c.coefs = zeros(size(k_cheb), size(χ), size(R))
end

function evaluate_coefs!(c::ChebCoefs, pk::AbstractArray{<:Any, 2}, bg::BackgroundQuantities, grid::AbstractCosmologicalGrid, cosmo::AbstractCosmology)
    evaluate_coefs!(c.cϕTT, pk, bg, grid, cosmo)
    evaluate_coefs!(c.cϕT, pk, bg, grid, cosmo)
    evaluate_coefs!(c.cϕ, pk, bg, grid, cosmo)
end
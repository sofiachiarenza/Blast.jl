# FUNCTIONS FOR THE DECOMPOSITION OF THE POWER SPECTRUM ON THE BASIS OF THE CHEBYSHEV POLYNOMIALS
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

function fast_chebcoefs(vals::AbstractArray, plan::Nothing)
    return nothing
end

abstract type AbstractCoeff end

@kwdef mutable struct cϕTT <: AbstractCoeff
    coefs::AbstractArray{<:Any, 3} 
end

@kwdef mutable struct cϕT <: AbstractCoeff
    coefs::AbstractArray{<:Any, 3}
end

@kwdef mutable struct cϕ <: AbstractCoeff
    coefs::AbstractArray{<:Any, 1}
end

#Constructors
make_coeff(::Type{T}, c) where {T<:AbstractCoeff} = T(c)
make_coeff(::Type{<:AbstractCoeff}, ::Nothing) = nothing

@inline function build_coeff(::Type{T}, vals::AbstractArray, plan::Union{FFTW.r2rFFTWPlan, Nothing}) where {T<:AbstractCoeff}
    c = Blast.fast_chebcoefs(vals, plan)
    return make_coeff(T, c)
end



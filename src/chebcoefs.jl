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

"""
    AbstractCoeff

Abstract supertype for Chebyshev-coefficient containers in BLAST.

Concrete subtypes of `AbstractCoeff` store Chebyshev expansion coefficients.

All coefficient types are constructed uniformly via and `build_coeff`.
"""
abstract type AbstractCoeff end

"""
    cϕTT <: AbstractCoeff

Container for Chebyshev coefficients of the unequal time power spectrum.

This type stores a three-dimensional array of the Chebyshev coefficients, where the decomposition is performed on the wakenumbers k.

Fields:
- `coefs::AbstractArray{<:Any,3}`: Chebyshev coefficients array of size (nk, nχ, nR).
"""
@kwdef mutable struct cϕTT <: AbstractCoeff
    coefs::AbstractArray{<:Any, 3} 
end

"""
    cϕT <: AbstractCoeff

Container for Chebyshev coefficients of the power spectrum built with only one transfer function (needed when PNGs are active).

This type stores a three-dimensional array of the Chebyshev coefficients, where the decomposition is performed on the wakenumbers k.

Fields:
- `coefs::AbstractArray{<:Any,3}`: Chebyshev coefficients array of size (nk, nχ, nR).
"""
@kwdef mutable struct cϕT <: AbstractCoeff
    coefs::AbstractArray{<:Any, 3}
end

"""
    cϕ <: AbstractCoeff

Container for Chebyshev coefficients of the primordial power spectrum (needed when PNGs are active).

This type stores a one-dimensional array of the Chebyshev coefficients, where the decomposition is performed on the wakenumbers k.

Fields:
- `coefs::AbstractArray{<:Any,3}`: Chebyshev coefficients array of size (nk).
"""
@kwdef mutable struct cϕ <: AbstractCoeff
    coefs::AbstractArray{<:Any, 1}
end

#Constructors
"""
    make_coeff(::Type{T}, c) where {T<:AbstractCoeff}

Construct a Chebyshev-coefficient container of type `T` from a coefficient array.

This is a lightweight constructor wrapper that standardizes the creation of
`AbstractCoeff` subtypes.

Arguments:
- `T<:AbstractCoeff`: The concrete coefficient type to construct (either cϕTT, cϕT, cϕ).
- `c`: Array of Chebyshev coefficients.

Returns:
- An instance of type `T` containing the provided coefficients.
"""
make_coeff(::Type{T}, c) where {T<:AbstractCoeff} = T(c)
make_coeff(::Type{<:AbstractCoeff}, ::Nothing) = nothing

"""
    build_coeff(::Type{T}, vals::AbstractArray, plan::Union{FFTW.r2rFFTWPlan,Nothing}) where {T<:AbstractCoeff}

Compute Chebyshev coefficients from input values and wrap them in a coefficient container.

This function applies an FFT-based Chebyshev transform to `vals` using the provided
FFTW plan, then constructs a coefficient object of type `T`. If `plan` is `nothing`,
no computation is performed and `nothing` is returned.

Arguments:
- `T<:AbstractCoeff`: The concrete coefficient type to construct.
- `vals::AbstractArray`: Input array to be decomposed on the Chebyshev polynomials (i.e. the power spectrum).
- `plan::Union{FFTW.r2rFFTWPlan,Nothing}`: FFTW plan used to compute Chebyshev coefficients.

Returns:
- An instance of type `T` containing the Chebyshev coefficients, or `nothing` if
  `plan` is `nothing` (i.e. if the coefficient is not active.).
"""
@inline function build_coeff(::Type{T}, vals::AbstractArray, plan::Union{FFTW.r2rFFTWPlan, Nothing}) where {T<:AbstractCoeff}
    c = Blast.fast_chebcoefs(vals, plan)
    return make_coeff(T, c)
end



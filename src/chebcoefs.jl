"""
   plan_fft(vals::AbstractArray{<:Number,N})

Create an FFTW real-to-real (R2R) transformation plan for the first axis of a given multidimensional array `vals`. 
In practise, the `vals` array is the power spectrum P(k,χ). The first axis should then contains the wavenumbers `k`, while the second axis contains the `χ` information.

# Arguments
- `vals::AbstractArray{<:Number, N}`: The input array of any numerical type with `N` dimensions.

# Returns
- `p::FFTW.rFFTWPlan`: A FFTW plan object for transforming `vals` with the appropriate real to real transformations. This plan can be applied using the `*` operator (e.g., `transformed_vals = p * vals`).

"""
function plan_fft(vals::AbstractArray{<:Number,N}) where {N}
    kind = map(n -> n > 1 ? FFTW.REDFT00 : FFTW.DHT, size(vals)[1])
    p = FFTW.plan_r2r(deepcopy(vals), kind, [1]; flags=FFTW.PATIENT, timelimit=Inf)   
                                                                                    
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

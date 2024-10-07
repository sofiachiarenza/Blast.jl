"""
   plan_fft(vals::AbstractArray{<:Number,N})

Create an FFTW real-to-real (R2R) transformation plan for a given multidimensional array `vals`.

# Arguments
- `vals::AbstractArray{<:Number, N}`: The input array of any numerical type with `N` dimensions.

# Returns
- `p::FFTW.rFFTWPlan`: An FFTW plan object for transforming `vals` with the appropriate R2R transformations. This plan can be applied using the `*` operator (e.g., `transformed_vals = p * vals`).

"""
function plan_fft(vals::AbstractArray{<:Number,N}) where {N}
    kind = map(n -> n > 1 ? FFTW.REDFT00 : FFTW.DHT, size(vals))
    p = FFTW.plan_r2r(copy(vals), kind; flags=FFTW.PATIENT, timelimit=Inf) #TODO: FFTW:ESTIMATE is the default. FFTW:PATIENT spend several seconds (or more) benchmarking different possible FFT algorithms and picking the fastest one. 
                                                                    # in this case, i think it's worth it to use patient but let's discuss. 
    return p 
end


"""
    fast_chebcoefs(vals::AbstractArray{<:Number,N})

Efficiently compute the Chebyshev coefficients of a multidimensional array `vals` using an O(n log n) method. This method leverages FFT-based type-I Discrete Cosine Transform (DCT-I).

Arguments:
- `vals::AbstractArray{<:Number,N}`: A multidimensional array of values for which to compute the Chebyshev coefficients.

Returns:
- `coefs`: An array of the same size as `vals`, containing the computed Chebyshev coefficients.
"""
function fast_chebcoefs(vals::AbstractArray{<:Number,N}, plan::FFTW.r2rFFTWPlan) where {N}
   coefs = plan*vals

   s = size(coefs)
   coefs ./= prod(map(n -> n > 1 ? 2(n-1) : 1, s))
   for dim = 1:N
       if size(coefs, dim) > 1
           coefs[CartesianIndices(ntuple(i -> i == dim ? (2:s[i]-1) : (1:s[i]), Val{N}()))] .*= 2
       end
   end

   return coefs
end
# FUNCTIONS FOR THE DECOMPOSITION OF THE POWER SPECTRUM ON THE BASIS OF THE CHEBYSHEV POLYNOMIALS, CENTRAL PART OF BLAST
"""
    plan_fft(vals::AbstractArray{<:Number, N}, axis::Int)

Create an FFTW real-to-real (R2R) transformation plan for a specified axis of a given multidimensional array `vals`. 
In practice, the `vals` array is often the power spectrum P(k,χ): the specified axis should contain the wavenumbers `k`. So if the power spectrum is given in a matrix of shape (nk, nχ), axis should be `1`. Instead, axis = `2` should be used for a matrix of shape (nχ, nk).

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

# FUNCTIONS FOR THE INTERPOLATION OF THE POWER SPECTRUM GIVEN BY THE EMULATOR 
function chebyshev_polynomials( x::AbstractArray{T,1}, n_cheb::Int, z_min::T, z_max::T) where T
    x_scaled = 2 .* (x .- z_min) ./ (z_max - z_min) .- 1.0
    
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

#TODO: this works for pk in shape (nz, nk), so axis in plan should be 1. Fix this!
function interpolate_power_spectrum(pk::AbstractArray{T,2}, z_nodes::AbstractArray{T,1}, R::AbstractArray{T,1},
    plan::FFTW.r2rFFTWPlan, bg::BackgroundQuantities, grid::AbstractCosmologicalGrid) where T

    coefs = fast_chebcoefs(pk, plan)
    new_χs = make_grid(bg, R)
    x = resample_redshifts(bg, grid, new_χs)
    chebyshevs = chebyshev_polynomials(x, length(z_nodes), minimum(z_nodes), maximum(z_nodes))
    pk_interp = coefs' * chebyshevs  #TODO: understand how to handle pk sizes
    return reshape(pk_interp, size(pk,2),  length(bg.χz_array), length(R))
end

function power_spectrum(pk::AbstractArray{T,3}) where T
    pk_R1 = pk[:,:,end]
    @tullio final_pk[i,c,r] := sqrt(pk_R1[i,c] * pk[i,c,r])
    return final_pk
end

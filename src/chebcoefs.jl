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
function cheb_poly(z_cheb::AbstractArray{T,1}, bg::BackgroundQuantities, grid::AbstractCosmologicalGrid, new_χ::AbstractArray{T,1}) where T
    nz = length(z_cheb)
    cheb_interp = FastChebInterp.ChebPoly(z_cheb, SA[minimum(z_cheb)], SA[maximum(z_cheb)])

    cheb_poly = zeros(nz, length(new_χ))

    for i in 1:nz
        copy_cheb = deepcopy(cheb_interp)
        copy_cheb.coefs .= 0
        copy_cheb.coefs[i] = 1.0
        cheb_poly[i,:] =  copy_cheb.(resample_redshifts(bg, grid, new_χ))
    end
    return cheb_poly
end

function make_grid_chebinterp(bg::BackgroundQuantities, R::Vector{T}) where T
    return vec(reverse(bg.χz_array) * R')
end

#TODO: this works for pk in shape (nz, nχ), so axis in plan should be 1. Fix this!
function interpolate_power_spectrum(pk::AbstractArray{T,2}, z_grid::AbstractArray{T,1}, plan::FFTW.r2rFFTWPlan, 
    bg::BackgroundQuantities, grid::AbstractCosmologicalGrid, R::AbstractArray{T,1}) where T

    coefs = fast_chebcoefs(pk, plan)
    χR_grid = make_grid_chebinterp(bg, R)
    chebyshevs = cheb_poly(z_grid, bg, grid, χR_grid)
    pk_interp = zeros(size(pk,2),length(χR_grid)) #TODO: understand how to handle pk sizes
    @tullio pk_interp[i,j] = coefs[k,i] * chebyshevs[k,j]

    return reshape(pk_interp, size(pk,2),  length(bg.χz_array), length(R))
end
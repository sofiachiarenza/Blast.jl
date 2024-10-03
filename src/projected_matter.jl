"""
    w_ell_tullio(c::AbstractArray,T::AbstractArray)
Compute the tensor contraction of the chebyshev coefficients of the power spectrum 'c' and the precomputed
integrals 'T' to obtain the projected matter densities.
"""
function w_ell_tullio(c::AbstractArray,T::AbstractArray)
    return @tullio w[i,j,k] := c[j,k,l] * T[i,j,k,l]
end

"""
    get_clencurt_grid(kmin::Number, kmax::Number, N::Number)
Return the integration points in k. They are a set of 'N' Chebyshev points rescaled between 'kmin' and 'kmax'.
"""
function get_clencurt_grid(kmin::Number, kmax::Number, N::Number)
    CC_obj = FastTransforms.chebyshevmoments1(Float64, N)
    x = FastTransforms.clenshawcurtisnodes(Float64, N)
    x = (kmax - kmin) / 2 * x .+ (kmin + kmax) / 2 

    x[1] *= (1-1e-8)
    x[end] *= (1+1e-8) #TODO: this is just a quick patch, need to figure this out properly.

    return x
end

"""
    get_clencurt_weights(kmin::Number, kmax::Number, N::Number)
Return the set of 'N' weights needed to perform the integration with the Clenshaw-Curtis quadrature rule.
The weights are rescaled between 'kmin' and 'kmax'.  
"""
function get_clencurt_weights(kmin::Number, kmax::Number, N::Number)
    CC_obj = FastTransforms.chebyshevmoments1(Float64, N)
    w = FastTransforms.clenshawcurtisweights(CC_obj)
    w = (kmax - kmin) / 2 * w

    return w
end

"""
    bessel_cheb_eval(ℓ::Number, kmin::Number, kmax::Number, χ::AbstractArray, n_cheb::Int, N::Number)
Return the Chebyshev polynomials up to order 'n_cheb+1' and the Bessel function of order 'ℓ' evaluated on the grid of 'N' Chebyshev points in the interval ['kmin', 'kmax'] and on the specified 'χ' points. 
"""
function bessel_cheb_eval(ℓ::Number, kmin::Number, kmax::Number, χ::AbstractArray, n_cheb::Int, N::Number)

    nχ = length(χ)
    x = get_clencurt_grid(kmin, kmax, N)

    k_cheb = chebpoints(n_cheb, log10(kmin), log10(kmax)) 
    c = FastChebInterp.ChebPoly(k_cheb, SA[log10(kmin)], SA[log10(kmax)])

    T = zeros(n_cheb+1,N) 
    Threads.@threads for i in 1:n_cheb+1
        copy_c = deepcopy(c) 
        copy_c.coefs .*= 0 
        copy_c.coefs[i] = 1.
        T[i,:] = copy_c.(log10.(x))
    end

    Bessel = zeros(nχ, N)
    Threads.@threads for i in 1:nχ
            Bessel[i,:] = @views SpecialFunctions.sphericalbesselj.(ℓ, χ[i] * x)
    end

    return T, Bessel

end

"""
    compute_T̃(ℓ::Number, χ::AbstractArray, R::AbstractArray, kmin::Number, kmax::Number, β::Number; n_cheb = 119, N=2^(15)+1)
Compute integrals of the Bessels function and the Chebyshev polynomials. This is the precomputation part of the code.
The parameters are:

    - ℓ: Multipole order

    - χ: Array containing values of the comoving distance. 

    - R: Array containing values for the R=χ₁/χ₂ variable.

    - kmin-kmax: Integration range in k.

    - β: Exponent of the k dependence of the integral. This parameter depends on the combination of tracers: β=2,-2,0 for clustering, cosmic shear and the cross-correlation respectively.

    - n_cheb: Number of chebyshev polynomials used in the approximation of the power spectra.

    - N: Number of integration points in k.
"""
function compute_T̃(ℓ::Number, χ::AbstractArray, R::AbstractArray, kmin::Number, kmax::Number, β::Number; n_cheb = 119, N=2^(15)+1)
    @assert kmin < kmax "The integration range is unphysical. Make sure kmin < kmax." #TODO: added assert because the test doesn't even work if they are the same.
    nχ = length(χ)
    nR = length(R)

    x = get_clencurt_grid(kmin, kmax, N)
    w = get_clencurt_weights(kmin, kmax, N)
    T, Bessel1 = Bessel_Cheb_eval(ℓ, kmin, kmax, χ, n_cheb, N)

    T_tilde = zeros(1, nχ, nR, n_cheb+1)
    
    for (ridx, r) in enumerate(R)
        Bessel2 = zeros(nχ, N)
        
        Threads.@threads for i in 1:nχ
            Bessel2[i,:] = @views SpecialFunctions.sphericalbesselj.(ℓ, r*χ[i] * x)
        end

        α = w .* (x .^ β) #β = 2 for CC, -2 for LL and 0 for CL.
         
        @tturbo for l in 1:n_cheb+1, i in 1:nχ
            Cij = zero(eltype(w))
            for k in 1:N
                Cij +=  T[l,k] * Bessel1[i,k] * Bessel2[i,k] * α[k]
            end
            T_tilde[1,i,ridx,l] = Cij
        end
    end

    return T_tilde

end
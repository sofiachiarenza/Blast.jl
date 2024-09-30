"""
    w_ell_tullio(c,T)
Compute the tensor contraction of the chebyshev coefficients of the power spectrum and the precomputed
integrals to obtain the projected matter densities (the inner integrals in k).
"""
function w_ell_tullio(c,T)
    return @tullio w[i,j,k] := c[j,k,l] * T[i,j,k,l]
end

function get_clencurt_grid(kmin::Number, kmax::Number, N::Number)
    CC_obj = FastTransforms.chebyshevmoments1(Float64, N)
    x = FastTransforms.clenshawcurtisnodes(Float64, N)
    x = (kmax - kmin) / 2 * x .+ (kmin + kmax) / 2 

    x[1] *= (1-1e-8)
    x[end] *= (1+1e-8) #TODO: this is just a quick patch, need to figure this out properly.

    return x
end

function get_clencurt_weights(kmin::Number, kmax::Number, N::Number)
    CC_obj = FastTransforms.chebyshevmoments1(Float64, N)
    w = FastTransforms.clenshawcurtisweights(CC_obj)
    w = (kmax - kmin) / 2 * w

    return w
end

function Bessel_Cheb_eval( ℓ::Number, kmin::Number, kmax::Number, χ::AbstractArray, n_cheb::Int, N::Number)

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
    compute_T̃(ℓ, χ, R, kmin, kmax, tracers; n_cheb = 119, N=2^(15)+1)
Compute integrals of the Bessels function and the Chebyshev polynomials. 
This is the precomputation part of the code.
"""

function compute_T̃(ℓ::Number, χ::AbstractArray, R::AbstractArray, kmin::Number, kmax::Number, β::Number; n_cheb = 119, N=2^(15)+1)
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
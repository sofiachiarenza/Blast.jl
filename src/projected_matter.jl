"""
    w_ell_tullio(c,T)
Compute the tensor contraction of the chebyshev coefficients of the power spectrum and the precomputed
integrals to obtain the projected matter densities (the inner integrals in k).
"""
function w_ell_tullio(c,T)
    return @tullio w[i,j,k] := c[j,k,l] * T[i,j,k,l]
end

"""
    compute_T̃(f, ℓ, χ, R, kmin, kmax, tracers; n_cheb = 119, N=2^(15)+1)
Compute integrals of the Bessels function and the Chebyshev polynomials. 
This is the precomputation part of the code.
"""


#TODO: the power spectrum is passed, but then I only use it to get the chebyshev polynomials. 
#I think this should be done internally, without the need to pass the power spectrum.
function compute_T̃(f, ℓ, χ, R, kmin, kmax, tracers; n_cheb = 119, N=2^(15)+1)
    nχ = length(χ)
    nR = length(R)

    CC_obj = FastTransforms.chebyshevmoments1(Float64, N)
    w = FastTransforms.clenshawcurtisweights(CC_obj)
    x = FastTransforms.clenshawcurtisnodes(Float64, N)
    
    #rescaling to my interval
    x = (kmax - kmin) / 2 * x .+ (kmin + kmax) / 2 
    w = (kmax - kmin) / 2 * w

    k_cheb = chebpoints(n_cheb, log10(minimum(minimum.(x))), log10(maximum(maximum.(x))))
    c = chebinterp(f(10. .^ k_cheb,χ[1],χ[1]), log10(minimum(minimum.(x))), log10(maximum(maximum.(x))))

    T = zeros(n_cheb+1,N)
    Threads.@threads for i in 1:n_cheb+1
        copy_c = deepcopy(c) #copio l'interpolante 
        copy_c.coefs .*= 0 #azzero i coeff del polinomio
        copy_c.coefs[i] = 1.
        T[i,:] = copy_c.(log10.(x))
    end

    f_bessel_χ1 = zeros(nχ, N)
    Threads.@threads for i in 1:nχ
            f_bessel_χ1[i,:] = @views SpecialFunctions.sphericalbesselj.(ℓ, χ[i] * x)
    end

    T_tilde = zeros(1, nχ, nR, n_cheb+1)
    
    for (ridx, r) in enumerate(R)
        f_bessel_χ2 = zeros(nχ, N)
        
        Threads.@threads for i in 1:nχ
            f_bessel_χ2[i,:] = @views SpecialFunctions.sphericalbesselj.(ℓ, r*χ[i] * x)
        end
    
        α = zeros(N)
    
        # TODO: update tracers names, this is not generalizable.
        if tracers == "CC"
            α = w .* (x .^ 2) 
        elseif tracers == "LL"
            α = w ./ (x .^ 2)
        elseif tracers == "CL"
            α = w
        end
         
        @tturbo for l in 1:n_cheb+1, i in 1:nχ
            Cij = zero(eltype(w))
            for k in 1:N
                Cij +=  T[l,k] * f_bessel_χ1[i,k] * f_bessel_χ2[i,k] * α[k]
            end
            T_tilde[1,i,ridx,l] = Cij
        end
    end

    return T_tilde

end
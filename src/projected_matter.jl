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

#TODO: missing documentation
function bessel_second_derivative(ℓ::Number, x::AbstractArray)
    return @views @. 2/x * SpecialFunctions.sphericalbesselj.(ℓ+1, x) + (ℓ^2-ℓ-x^2)/x^2 * SpecialFunctions.sphericalbesselj.(ℓ, x)
end

"""
    bessel_cheb_eval(ℓ::Number, kmin::Number, kmax::Number, χ::AbstractArray, n_cheb::Int, N::Number)
Return the Chebyshev polynomials up to order 'n_cheb+1' and the Bessel function of order 'ℓ' evaluated on the grid of 'N' Chebyshev points in the interval ['kmin', 'kmax'] and on the specified 'χ' points. 
"""
function bessel_cheb_eval(ℓ::Number, kmin::Number, kmax::Number, χ::AbstractArray, n_cheb::Int, N::Number, deriv_order::Int)

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
        if deriv_order == 0
            Bessel[i, :] = @views SpecialFunctions.sphericalbesselj.(ℓ, χ[i] * x)
        elseif deriv_order == 2
            Bessel[i, :] = @views bessel_second_derivative(ℓ, χ[i] * x)
        else
            error("Invalid derivative order: $deriv_order. Expected 0 or 2.")
        end
    end

    return T, Bessel

end



"""
    compute_T̃(ℓ::Number, χ::AbstractArray, R::AbstractArray, kmin::Number, kmax::Number, β::Number; n_cheb::Int = 119, N::Int = 2^(15)+1)
Compute integrals of the Bessels function and the Chebyshev polynomials. This is the precomputation part of the code.

# Arguments
- `ℓ::Number`: Multipole order

- `χ::AbstractArray`: Array containing values of the comoving distance. 

- `R::AbstractArray`: Array containing values for the R=χ₁/χ₂ variable.

- `kmin::Number` and `kmax::Number`: Integration range in k.

- `β::Number`: Exponent of the k dependence of the integral. This parameter depends on the combination of tracers: β=2,-2,0 for clustering, cosmic shear and the cross-correlation respectively.

- `n_cheb::Int`: Number of chebyshev polynomials used in the approximation of the power spectra.

- `N::Int`: Number of integration points in k.
"""
function compute_T̃(ℓ::Number, χ::AbstractArray, R::AbstractArray, kmin::Number, kmax::Number, β::Number, der_A::Int, der_B::Int; n_cheb::Int = 119, N::Int = 2^(15)+1)
    if kmin >= kmax 
        throw(DomainError("The integration range is unphysical. Make sure kmin < kmax.")) 
    end
    
    nχ = length(χ)
    nR = length(R)

    x = get_clencurt_grid(kmin, kmax, N)
    w = get_clencurt_weights(kmin, kmax, N)

    T, Bessel1 = bessel_cheb_eval(ℓ, kmin, kmax, χ, n_cheb, N, der_A)


    T_tilde = zeros(1, nχ, nR, n_cheb+1)
    
    for (ridx, r) in enumerate(R)
        Bessel2 = zeros(nχ, N)
        
        if der_B == 0
            Threads.@threads for i in 1:nχ
                Bessel2[i,:] = @views SpecialFunctions.sphericalbesselj.(ℓ, r*χ[i] * x)
            end
        elseif der_B == 2
            Threads.@threads for i in 1:nχ
                Bessel2[i,:] = bessel_second_derivative(ℓ, r*χ[i] * x)
            end
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

## Handling of the various combinations of T̃

function w_ell_tullio(c::AbstractCoeffComponents, T::AbstractArray{<:Any, 4})
    cheb_coeffs = c.coefs
    return @tullio w[i,j,k] := cheb_coeffs[j,k,l] * T[i,j,k,l]
end

abstract type AbstractProjectedMatterDensity end

@kwdef mutable struct w_2_00_ϕTT <: AbstractProjectedMatterDensity
    active::Bool = false
    T̃ = T_tildes.T_2_00
    w::AbstractArray{<:Any, 3}
end

function compute_w!(w::w_2_00_ϕTT, c::ChebCoefs)
    w.w = w_ell_tullio(c.cϕTT, w.T̃)
end

@kwdef mutable struct w_minus2_00_ϕTT <: AbstractProjectedMatterDensity
    active::Bool = false
    T̃ = T_tildes.T_minus2_00
    w::AbstractArray{<:Any, 3}
end

function compute_w!(w::w_minus2_00_ϕTT, c::ChebCoefs)
    w.w = w_ell_tullio(c.cϕTT, w.T̃)
end

@kwdef mutable struct w_0_00_ϕTT <: AbstractProjectedMatterDensity
    active::Bool = false
    T̃ = T_tildes.T_0_00
    w::AbstractArray{<:Any, 3}
end

function compute_w!(w::w_0_00_ϕTT, c::ChebCoefs)
    w.w = w_ell_tullio(c.cϕTT, w.T̃)
end

@kwdef mutable struct w_0_02_ϕTT <: AbstractProjectedMatterDensity
    active::Bool = false
    T̃ = T_tildes.T_0_02
    w::AbstractArray{<:Any, 3}
end

function compute_w!(w::w_0_02_ϕTT, c::ChebCoefs)
    w.w = w_ell_tullio(c.cϕTT, w.T̃)
end

@kwdef mutable struct w_0_20_ϕTT <: AbstractProjectedMatterDensity
    active::Bool = false
    T̃ = T_tildes.T_0_20
    w::AbstractArray{<:Any, 3}
end

function compute_w!(w::w_0_20_ϕTT, c::ChebCoefs)
    w.w = w_ell_tullio(c.cϕTT, w.T̃)
end

@kwdef mutable struct w_2_02_ϕTT <: AbstractProjectedMatterDensity
    active::Bool = false
    T̃ = T_tildes.T_2_02
    w::AbstractArray{<:Any, 3}
end

function compute_w!(w::w_2_02_ϕTT, c::ChebCoefs)
    w.w = w_ell_tullio(c.cϕTT, w.T̃)
end

@kwdef mutable struct w_2_20_ϕTT <: AbstractProjectedMatterDensity
    active::Bool = false
    T̃ = T_tildes.T_2_20
    w::AbstractArray{<:Any, 3}
end

function compute_w!(w::w_2_20_ϕTT, c::ChebCoefs)
    w.w = w_ell_tullio(c.cϕTT, w.T̃)
end

@kwdef mutable struct w_2_22_ϕTT <: AbstractProjectedMatterDensity
    active::Bool = false
    T̃ = T_tildes.T_2_22
    w::AbstractArray{<:Any, 3}
end

function compute_w!(w::w_2_22_ϕTT, c::ChebCoefs)
    w.w = w_ell_tullio(c.cϕTT, w.T̃)
end

@kwdef mutable struct w_2_00_ϕT <: AbstractProjectedMatterDensity
    active::Bool = false
    T̃ = T_tildes.T_2_00
    w::AbstractArray{<:Any, 3}
end

function compute_w!(w::w_2_00_ϕT, c::ChebCoefs)
    w.w = w_ell_tullio(c.cϕT, w.T̃)
end

@kwdef mutable struct w_0_00_ϕT <: AbstractProjectedMatterDensity
    active::Bool = false
    T̃ = T_tildes.T_0_00
    w::AbstractArray{<:Any, 3}
end

function compute_w!(w::w_0_00_ϕT, c::ChebCoefs)
    w.w = w_ell_tullio(c.cϕT, w.T̃)
end

@kwdef mutable struct w_2_02_ϕT <: AbstractProjectedMatterDensity
    active::Bool = false
    T̃ = T_tildes.T_2_02
    w::AbstractArray{<:Any, 3}
end

function compute_w!(w::w_2_02_ϕT, c::ChebCoefs)
    w.w = w_ell_tullio(c.cϕT, w.T̃)
end

@kwdef mutable struct w_2_20_ϕT <: AbstractProjectedMatterDensity
    active::Bool = false
    T̃ = T_tildes.T_2_20
    w::AbstractArray{<:Any, 3}
end

function compute_w!(w::w_2_20_ϕT, c::ChebCoefs)
    w.w = w_ell_tullio(c.cϕT, w.T̃)
end

@kwdef mutable struct w_2_00_ϕ <: AbstractProjectedMatterDensity
    active::Bool = false
    T̃ = T_tildes.T_2_00
    w::AbstractArray{<:Any, 3}
end

function compute_w!(w::w_2_00_ϕ, c::ChebCoefs)
    w.w = w_ell_tullio(c.cϕ, w.T̃)
end

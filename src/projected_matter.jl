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
    bessel_second_derivative(ℓ::Number, x::AbstractArray)

Compute the second derivative with respect to the argument of the spherical
Bessel function of order `ℓ`.

Specifically, this returns 
```math
\\frac{d^2 j_\\ell(x)}{dx^2}
```
using the analytic recurrence relation in terms of `j_ℓ(x)` and `j_{ℓ+1}(x)`.

# Arguments
- `ℓ::Number`: Multipole order.
- `x::AbstractArray`: Array of arguments at which the derivative is evaluated.

# Returns
An array of the same shape as `x`, containing the second derivative of
`j_ℓ(x)`.
"""
function bessel_second_derivative(ℓ::Number, x::AbstractArray)
    return @views @. 2/x * SpecialFunctions.sphericalbesselj.(ℓ+1, x) + (ℓ^2-ℓ-x^2)/x^2 * SpecialFunctions.sphericalbesselj.(ℓ, x)
end

"""
    bessel_cheb_eval(ℓ::Number, kmin::Number, kmax::Number, χ::AbstractArray, n_cheb::Int, N::Number)

Evaluate spherical Bessel functions and Chebyshev polynomials on a common
Clenshaw–Curtis integration grid.

This function precomputes:

1. The values of Chebyshev polynomials up to order `n_cheb`,
   evaluated at `N` Clenshaw–Curtis nodes between `kmin` and `kmax`.
2. The spherical Bessel function `j_ℓ(kχ)` (or its second derivative)
   evaluated at the same `k` nodes for all `χ`.

These quantities are used to build the precomputed kernels `T̃` appearing
in the projected matter integrals.

# Arguments
- `ℓ::Number`: Multipole order.
- `kmin::Number`, `kmax::Number`: Integration range in wavenumber `k`.
- `χ::AbstractArray`: Comoving distance grid.
- `n_cheb::Int`: Maximum Chebyshev polynomial order.
- `N::Number`: Number of Clenshaw–Curtis integration nodes.
- `deriv_order::Int`: Derivative order of the Bessel function (`0` or `2`).

# Returns
- `T::Array`: Chebyshev polynomials evaluated on the `k` grid, with shape `(n_cheb+1, N)`.
- `Bessel::Array`: Bessel function values with shape `(length(χ), N)`.
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

Compute the precomputed kernel `T̃` entering the projected matter density
integrals.

This function evaluates integrals of the form:
```math
\\tilde{T}(\\ell, \\chi, R, n) = \\int dk k^\\beta j_\\ell^{A}(k\\chi) j_\\ell^{B}(k R\\chi) T_n(k)
```
where `T_n` are Chebyshev polynomials and `j_ℓ^{A,B}` denote spherical
Bessel functions or their second derivatives.

These kernels are **cosmology-independent** once the global `χ` grid is fixed,
and are therefore precomputed and stored as artifacts.

# Arguments
- `ℓ::Number`: Multipole order.
- `χ::AbstractArray`: Global comoving distance grid.
- `R::AbstractArray`: Ratio grid `R = χ₁ / χ₂`.
- `kmin::Number`, `kmax::Number`: Integration limits in wavenumber.
- `β::Number`: Power-law exponent of `k` in the integrand.
- `der_A::Int`, `der_B::Int`: Derivative orders (`0` or `2`) for the two Bessel factors.
- `n_cheb::Int`: Number of Chebyshev polynomials used to approximate the power spectrum.
- `N::Int`: Number of Clenshaw–Curtis integration nodes.

# Returns
A 4D array with shape `(1, length(χ), length(R), n_cheb+1)` containing the
precomputed kernel `T̃`.
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

"""
    w_ell_tullio(c, T)

Contract Chebyshev coefficients of the power spectrum with the precomputed `T̃`
to form the projected matter density `w_\ell(χ_1, χ_2)`.

This performs the sum over Chebyshev indices using `Tullio` for efficient
tensor contraction. Multiple methods are provided depending on the dimensionality
of the coefficient array `c`.

# Arguments
- `c`: Chebyshev coefficients of the power spectrum.
- `T`: Precomputed kernel `T̃`.

# Returns
An array containing the projected matter density `w`.
"""
function w_ell_tullio(c::AbstractArray{<:Any, 3}, T::AbstractArray{<:Any, 4})
    return @tullio w[i,j,k] := c[l,j,k] * T[i,j,k,l]
end

function w_ell_tullio(c::AbstractArray{<:Any, 2}, T::AbstractArray{<:Any, 4})
    return @tullio w[i,j,k] := c[l,j] * T[i,j,k,l]
end

function w_ell_tullio(c::AbstractArray{<:Any, 1}, T::AbstractArray{<:Any, 4})
    return @tullio w[i,j,k] := c[l] * T[i,j,k,l]
end


abstract type AbstractProjectedMatterDensity end

"""
Abstract supertype for projected matter density components.

Each concrete subtype represents a specific combination, which depends on:
- The power of k in the precomputed `\\tilde{T}`.
- The order of the derivatives of the Bessel functions in the precomputed `\\tilde{T}`.
- The power spectrum (i.e. the full unequal time, the transfer function of the primordial power spectrum.)

Each component stores:
- A reference to the relevant `T̃`,
- The corresponding projected weight array `w`.
"""
abstract type ProjectedMatterDensityComponent end

function compute_w!(w::Nothing, c::PowerSpectrum)
    return nothing
end

@kwdef mutable struct w_2_00_ϕTT <: ProjectedMatterDensityComponent
    T̃ = T_tildes.T_2_00
    w::AbstractArray{<:Any, 3} = zeros(1, 1, 1)
end

function compute_w!(w::w_2_00_ϕTT, c::PowerSpectrum)
    w.w = w_ell_tullio(c.cϕTT.coefs, w.T̃)
end

@kwdef mutable struct w_minus2_00_ϕTT <: ProjectedMatterDensityComponent
    T̃ = T_tildes.T_minus2_00
    w::AbstractArray{<:Any, 3} = zeros(1, 1, 1)
end

function compute_w!(w::w_minus2_00_ϕTT, c::PowerSpectrum)
    w.w = w_ell_tullio(c.cϕTT.coefs, w.T̃)
end

@kwdef mutable struct w_0_00_ϕTT <: ProjectedMatterDensityComponent
    T̃ = T_tildes.T_0_00
    w::AbstractArray{<:Any, 3} = zeros(1, 1, 1)
end

function compute_w!(w::w_0_00_ϕTT, c::PowerSpectrum)
    w.w = w_ell_tullio(c.cϕTT.coefs, w.T̃)
end

@kwdef mutable struct w_0_02_ϕTT <: ProjectedMatterDensityComponent
    T̃ = T_tildes.T_0_02
    w::AbstractArray{<:Any, 3} = zeros(1, 1, 1)
end

function compute_w!(w::w_0_02_ϕTT, c::PowerSpectrum)
    w.w = w_ell_tullio(c.cϕTT.coefs, w.T̃)
end

@kwdef mutable struct w_0_20_ϕTT <: ProjectedMatterDensityComponent
    T̃ = T_tildes.T_0_20
    w::AbstractArray{<:Any, 3} = zeros(1, 1, 1)
end

function compute_w!(w::w_0_20_ϕTT, c::PowerSpectrum)
    w.w = w_ell_tullio(c.cϕTT.coefs, w.T̃)
end

@kwdef mutable struct w_2_02_ϕTT <: ProjectedMatterDensityComponent
    T̃ = T_tildes.T_2_02
    w::AbstractArray{<:Any, 3} = zeros(1, 1, 1)
end

function compute_w!(w::w_2_02_ϕTT, c::PowerSpectrum)
    w.w = w_ell_tullio(c.cϕTT.coefs, w.T̃)
end

@kwdef mutable struct w_2_20_ϕTT <: ProjectedMatterDensityComponent
    T̃ = T_tildes.T_2_20
    w::AbstractArray{<:Any, 3} = zeros(1, 1, 1)
end

function compute_w!(w::w_2_20_ϕTT, c::PowerSpectrum)
    w.w = w_ell_tullio(c.cϕTT.coefs, w.T̃)
end

@kwdef mutable struct w_2_22_ϕTT <: ProjectedMatterDensityComponent
    T̃ = T_tildes.T_2_22
    w::AbstractArray{<:Any, 3} = zeros(1, 1, 1)
end

function compute_w!(w::w_2_22_ϕTT, c::PowerSpectrum)
    w.w = w_ell_tullio(c.cϕTT.coefs, w.T̃)
end

@kwdef mutable struct w_2_00_ϕT <: ProjectedMatterDensityComponent
    T̃ = T_tildes.T_2_00
    w::AbstractArray{<:Any, 3} = zeros(1, 1, 1)
end

function compute_w!(w::w_2_00_ϕT, c::PowerSpectrum)
    w.w = w_ell_tullio(c.cϕT.coefs, w.T̃)
end

@kwdef mutable struct w_2_00_ϕT_R1 <: ProjectedMatterDensityComponent
    T̃ = T_tildes.T_2_00
    w::AbstractArray{<:Any, 3} = zeros(1, 1, 1)
end

function compute_w!(w::w_2_00_ϕT_R1, c::PowerSpectrum)
    coefs_R1 = c.cϕT.coefs[:,:, end]
    w.w = w_ell_tullio(coefs_R1, w.T̃)
end

@kwdef mutable struct w_0_00_ϕT <: ProjectedMatterDensityComponent
    T̃ = T_tildes.T_0_00
    w::AbstractArray{<:Any, 3} = zeros(1, 1, 1)
end

function compute_w!(w::w_0_00_ϕT, c::PowerSpectrum)
    w.w = w_ell_tullio(c.cϕT.coefs, w.T̃)
end

@kwdef mutable struct w_0_00_ϕT_R1 <: ProjectedMatterDensityComponent
    T̃ = T_tildes.T_0_00
    w::AbstractArray{<:Any, 3} = zeros(1, 1, 1)
end

function compute_w!(w::w_0_00_ϕT_R1, c::PowerSpectrum)
    coefs_R1 = c.cϕT.coefs[:,:, end]
    w.w = w_ell_tullio(coefs_R1, w.T̃)
end

@kwdef mutable struct w_2_02_ϕT <: ProjectedMatterDensityComponent
    T̃ = T_tildes.T_2_02
    w::AbstractArray{<:Any, 3} = zeros(1, 1, 1)
end

function compute_w!(w::w_2_02_ϕT, c::PowerSpectrum)
    w.w = w_ell_tullio(c.cϕT.coefs, w.T̃)
end

@kwdef mutable struct w_2_02_ϕT_R1 <: ProjectedMatterDensityComponent
    T̃ = T_tildes.T_2_02
    w::AbstractArray{<:Any, 3} = zeros(1, 1, 1)
end

function compute_w!(w::w_2_02_ϕT_R1, c::PowerSpectrum)
    coefs_R1 = c.cϕT.coefs[:,:, end]
    w.w = w_ell_tullio(coefs_R1, w.T̃)
end

@kwdef mutable struct w_2_20_ϕT <: ProjectedMatterDensityComponent
    T̃ = T_tildes.T_2_20
    w::AbstractArray{<:Any, 3} = zeros(1, 1, 1)
end

function compute_w!(w::w_2_20_ϕT, c::PowerSpectrum)
    w.w = w_ell_tullio(c.cϕT.coefs, w.T̃)
end

@kwdef mutable struct w_2_20_ϕT_R1 <: ProjectedMatterDensityComponent
    T̃ = T_tildes.T_2_20
    w::AbstractArray{<:Any, 3} = zeros(1, 1, 1)
end

function compute_w!(w::w_2_20_ϕT_R1, c::PowerSpectrum)
    coefs_R1 = c.cϕT.coefs[:,:,end]
    w.w = w_ell_tullio(coefs_R1, w.T̃)
end

@kwdef mutable struct w_2_00_ϕ <: ProjectedMatterDensityComponent
    T̃ = T_tildes.T_2_00
    w::AbstractArray{<:Any, 3} = zeros(1, 1, 1)
end

function compute_w!(w::w_2_00_ϕ, c::PowerSpectrum)
    w.w = w_ell_tullio(c.cϕ.coefs, w.T̃)
end

"""
    ProjectedMatterDensity

Container holding all projected matter density components required for all the active observables and components.

Each field corresponds to a specific kernel contribution (e.g. `w_2_00_ϕTT`,
`w_0_02_ϕT_R1`). Fields may be set to `nothing` if the corresponding contribution
is not required.

The container is populated during setup and filled by calling `compute_w!`
with a `PowerSpectrum` object.
"""
@kwdef mutable struct ProjectedMatterDensity <: AbstractProjectedMatterDensity
    w_2_00_ϕTT::Union{w_2_00_ϕTT, Nothing} = nothing
    w_minus2_00_ϕTT::Union{w_minus2_00_ϕTT, Nothing} = nothing
    w_0_00_ϕTT::Union{w_0_00_ϕTT, Nothing} = nothing
    w_0_02_ϕTT::Union{w_0_02_ϕTT, Nothing} = nothing
    w_0_20_ϕTT::Union{w_0_20_ϕTT, Nothing} = nothing
    w_2_02_ϕTT::Union{w_2_02_ϕTT, Nothing} = nothing
    w_2_20_ϕTT::Union{w_2_20_ϕTT, Nothing} = nothing
    w_2_22_ϕTT::Union{w_2_22_ϕTT, Nothing} = nothing
    w_2_00_ϕT::Union{w_2_00_ϕT, Nothing} = nothing
    w_2_00_ϕT_R1::Union{w_2_00_ϕT_R1, Nothing} = nothing
    w_0_00_ϕT::Union{w_0_00_ϕT, Nothing} = nothing
    w_0_00_ϕT_R1::Union{w_0_00_ϕT_R1, Nothing} = nothing
    w_2_02_ϕT::Union{w_2_02_ϕT, Nothing} = nothing
    w_2_02_ϕT_R1::Union{w_2_02_ϕT_R1, Nothing} = nothing
    w_2_20_ϕT::Union{w_2_20_ϕT, Nothing} = nothing
    w_2_20_ϕT_R1::Union{w_2_20_ϕT_R1, Nothing} = nothing
    w_2_00_ϕ::Union{w_2_00_ϕ, Nothing} = nothing
end 

"""
    compute_w!(W::ProjectedMatterDensity, c::PowerSpectrum)

Compute all active projected matter density components.
"""
function compute_w!(W::ProjectedMatterDensity, c::PowerSpectrum)
    compute_w!(W.w_2_00_ϕTT, c)
    compute_w!(W.w_minus2_00_ϕTT, c)
    compute_w!(W.w_0_00_ϕTT, c)
    compute_w!(W.w_0_02_ϕTT, c)
    compute_w!(W.w_0_20_ϕTT, c)
    compute_w!(W.w_2_02_ϕTT, c)
    compute_w!(W.w_2_20_ϕTT, c)
    compute_w!(W.w_2_22_ϕTT, c)
    compute_w!(W.w_2_00_ϕT, c)
    compute_w!(W.w_2_00_ϕT_R1, c)
    compute_w!(W.w_0_00_ϕT, c)
    compute_w!(W.w_0_00_ϕT_R1, c)
    compute_w!(W.w_2_02_ϕT, c)
    compute_w!(W.w_2_02_ϕT_R1, c)
    compute_w!(W.w_2_20_ϕT, c)
    compute_w!(W.w_2_20_ϕT_R1, c)
    compute_w!(W.w_2_00_ϕ, c)
end
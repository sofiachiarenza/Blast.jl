using AbstractCosmologicalEmulators
using ChainRulesCore
using ForwardDiff: Dual, value, partials, Partials, tagtype

import AbstractCosmologicalEmulators: ChebyshevPlan, chebpoints, prepare_chebyshev_plan, chebyshev_decomposition, chebyshev_polynomials

"""
    chebinterp_native(c::AbstractVector, x_grid, x_min, x_max)

Evaluates a 1D Chebyshev expansion on a grid.
"""
function chebinterp_native(c::AbstractVector, x_grid, x_min, x_max)
    xx = float.(x_grid)
    T_type = eltype(xx)
    T = chebyshev_polynomials(xx, T_type(x_min), T_type(x_max), length(c) - 1)
    return T * c
end

"""
    AbstractCoeff

Abstract supertype for Chebyshev-coefficient containers in BLAST.
"""
abstract type AbstractCoeff end

"""
    cϕTT <: AbstractCoeff
Container for Chebyshev coefficients of the unequal time power spectrum.
"""
@kwdef mutable struct cϕTT <: AbstractCoeff
    coefs::AbstractArray{<:Any, 3} 
end

"""
    cϕT <: AbstractCoeff
Container for Chebyshev coefficients of the power spectrum built with only one transfer function.
"""
@kwdef mutable struct cϕT <: AbstractCoeff
    coefs::AbstractArray{<:Any, 3}
end

"""
    cϕ <: AbstractCoeff
Container for Chebyshev coefficients of the primordial power spectrum.
"""
@kwdef mutable struct cϕ <: AbstractCoeff
    coefs::AbstractArray{<:Any, 1}
end

# Constructors
make_coeff(::Type{T}, c) where {T<:AbstractCoeff} = T(c)
make_coeff(::Type{<:AbstractCoeff}, ::Nothing) = nothing

"""
    build_coeff(::Type{T}, vals::AbstractArray, plan::Union{ChebyshevPlan,Nothing}) where {T<:AbstractCoeff}

Compute Chebyshev coefficients from input values and wrap them in a coefficient container.
"""
@inline function build_coeff(::Type{T}, vals::AbstractArray, plan::Union{ChebyshevPlan, Nothing}) where {T<:AbstractCoeff}
    c = chebyshev_decomposition(plan, vals)
    return make_coeff(T, c)
end

@inline function build_coeff(::Type{T}, vals::AbstractArray, plan::Nothing) where {T<:AbstractCoeff}
    return nothing
end

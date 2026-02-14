import ChainRulesCore: rrule, NoTangent, unthunk, ProjectTo
import Tullio: @tullio
import FFTW
import Mooncake
import Mooncake: @from_chainrules, MinimalCtx
import Mooncake: tangent_type, fdata_type, rdata_type, zero_tangent_internal, fdata, rdata, NoFData, NoRData, increment_rdata!!

# w_ell_tullio rrules 
# NOTE: In Blast, T is treated as constant 

function rrule(::typeof(w_ell_tullio),
               c::AbstractArray{<:Any, 3},
               T::AbstractArray{<:Any, 4})
    y  = w_ell_tullio(c, T)
    pc = ProjectTo(c)

    function w_ell_tullio_pullback(ȳ)
        ȳ = unthunk(ȳ)
        @tullio ∂c[l,j,k] := ȳ[i,j,k] * T[i,j,k,l]
        return (NoTangent(), pc(∂c), NoTangent())
    end

    return (y, w_ell_tullio_pullback)
end

@from_chainrules MinimalCtx Tuple{typeof(w_ell_tullio), AbstractArray{<:Any, 3}, AbstractArray{<:Any, 4}}

function rrule(::typeof(w_ell_tullio),
               c::AbstractArray{<:Any, 2},
               T::AbstractArray{<:Any, 4})
    y  = w_ell_tullio(c, T)
    pc = ProjectTo(c)

    function w_ell_tullio_pullback(ȳ)
        ȳ = unthunk(ȳ)
        @tullio ∂c[l,j] := ȳ[i,j,k] * T[i,j,k,l]
        return (NoTangent(), pc(∂c), NoTangent())
    end

    return (y, w_ell_tullio_pullback)
end

@from_chainrules MinimalCtx Tuple{typeof(w_ell_tullio), AbstractArray{<:Any, 2}, AbstractArray{<:Any, 4}}

function rrule(::typeof(w_ell_tullio),
               c::AbstractArray{<:Any, 1},
               T::AbstractArray{<:Any, 4})
    y  = w_ell_tullio(c, T)
    pc = ProjectTo(c)

    function w_ell_tullio_pullback(ȳ)
        ȳ = unthunk(ȳ)
        @tullio ∂c[l] := ȳ[i,j,k] * T[i,j,k,l]
        return (NoTangent(), pc(∂c), NoTangent())
    end

    return (y, w_ell_tullio_pullback)
end

@from_chainrules MinimalCtx Tuple{typeof(w_ell_tullio), AbstractArray{<:Any, 1}, AbstractArray{<:Any, 4}}

# fast_chebcoefs rrule 
# - `plan` is treated as constant

function rrule(::typeof(fast_chebcoefs), vals::AbstractArray, plan::FFTW.r2rFFTWPlan)
    coefs = fast_chebcoefs(vals, plan)

    s = size(coefs)
    N = length(s)

    pv = ProjectTo(vals)

    function fast_chebcoefs_pullback(ȳ)
        dy = unthunk(ȳ)

        K = 2 * (s[1] - 1)
        z̄ = dy ./ K

        vals̄ = plan * z̄

        indices = CartesianIndices(ntuple(i -> i == 1 ? (2:s[1]-1) : (1:s[i]), Val{N}()))
        vals̄[indices] .*= 2

        return (NoTangent(), pv(vals̄), NoTangent())
    end

    return (coefs, fast_chebcoefs_pullback)
end

# Mooncake patches for FFTW plans, makes Mooncake see FFTW plans as constants
Mooncake.@from_rrule Mooncake.DefaultCtx Tuple{typeof(fast_chebcoefs), AbstractArray, FFTW.r2rFFTWPlan}

Mooncake.tangent_type(::Type{P}) where {P<:FFTW.FFTWPlan} = P
Mooncake.fdata_type(::Type{P})   where {P<:FFTW.FFTWPlan} = NoFData
Mooncake.rdata_type(::Type{P})   where {P<:FFTW.FFTWPlan} = NoRData

function Mooncake.zero_tangent_internal(p::P, ::IdDict{Any, Any}) where {P<:FFTW.FFTWPlan}
    return p
end

Mooncake.fdata(p::FFTW.FFTWPlan) = NoFData()
Mooncake.rdata(p::FFTW.FFTWPlan) = NoRData()
Mooncake.increment_rdata!!(x::FFTW.FFTWPlan, ::NoRData) = x

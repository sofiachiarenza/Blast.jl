module MooncakeExt

import Mooncake
import Mooncake:
    @from_chainrules, @is_primitive, CoDual, MinimalCtx,
    NoFData, NoRData, primal, rrule!!, tangent

import Blast:
    get_clencurt_weights, simpson_weights_array,
    get_clencurt_weights_R_integration, simpson_weights_matrix,
    _combine_kernels_tullio, _compute_Cℓ_tullio, _compute_Cℓ_rsd_tullio,
    _compute_Cℓ_fused_tullio, _compute_Cℓ_rsd_fused_tullio,
    _limber_contraction, w_ell_tullio, limber_eval,
    _w_single_inplace!, _w_fused3_inplace!, _w_fused4_inplace!,
    _p_phi_TT_tullio, _p_phi_T_tullio,
    _cosmic_shear_kernel_tullio, _magnification_bias_kernel_tullio

# Note: Akima interpolation and chebyshev_decomposition Mooncake rules live in
# AbstractCosmologicalEmulators.jl's own MooncakeExt — we rely on those.

# =============================================================================
# Quadrature weights (treat as constants — zero-gradient pass-through)
# =============================================================================

@from_chainrules MinimalCtx Tuple{typeof(get_clencurt_weights), Any, Any, Any}
@from_chainrules MinimalCtx Tuple{typeof(simpson_weights_array), Int}
@from_chainrules MinimalCtx Tuple{typeof(get_clencurt_weights_R_integration), Int}
@from_chainrules MinimalCtx Tuple{typeof(simpson_weights_matrix), Int}

# =============================================================================
# Kernel combination and Cℓ integration
# =============================================================================

@from_chainrules MinimalCtx Tuple{typeof(_combine_kernels_tullio), AbstractArray{<:Any, 3}, AbstractArray{<:Any, 3}}

@from_chainrules MinimalCtx Tuple{typeof(_compute_Cℓ_tullio), AbstractArray, AbstractArray, AbstractVector, AbstractVector, AbstractVector, Number, AbstractVector}
@from_chainrules MinimalCtx Tuple{typeof(_compute_Cℓ_rsd_tullio), AbstractArray, AbstractArray, AbstractArray, AbstractArray, AbstractArray, AbstractArray, AbstractVector, AbstractVector, AbstractVector, Number, AbstractVector}

@from_chainrules MinimalCtx Tuple{typeof(_compute_Cℓ_fused_tullio), AbstractArray{<:Any, 3}, AbstractArray{<:Any, 3}, AbstractArray, AbstractVector, AbstractVector, AbstractVector, Number, AbstractVector}
@from_chainrules MinimalCtx Tuple{typeof(_compute_Cℓ_rsd_fused_tullio), AbstractArray{<:Any, 3}, AbstractArray{<:Any, 3}, AbstractArray, AbstractArray, AbstractVector, AbstractVector, AbstractVector, Number, AbstractVector}

@from_chainrules MinimalCtx Tuple{typeof(_limber_contraction), AbstractArray, AbstractArray, AbstractArray, AbstractVector, Number}
@from_chainrules MinimalCtx Tuple{typeof(limber_eval), AbstractMatrix, AbstractMatrix, AbstractArray{<:Any, 3}}

# =============================================================================
# Projected matter density contraction
# =============================================================================

@from_chainrules MinimalCtx Tuple{typeof(w_ell_tullio), AbstractArray{<:Any, 3}, AbstractArray{<:Any, 4}}
@from_chainrules MinimalCtx Tuple{typeof(w_ell_tullio), AbstractArray{<:Any, 2}, AbstractArray{<:Any, 4}}
@from_chainrules MinimalCtx Tuple{typeof(w_ell_tullio), AbstractArray{<:Any, 1}, AbstractArray{<:Any, 4}}

# The mutating kernels cannot use `@from_chainrules`: they overwrite their
# output arguments, and their result arrays alias those arguments. Native
# Mooncake rules are required to preserve mutation-restoration semantics.

function _increment_coefficient_adjoint!(dc::AbstractArray{<:Any,3}, dw, T)
    Threads.@threads for k in axes(T, 3)
        @inbounds for j in axes(T, 2), l in axes(T, 4)
            value = zero(eltype(dc))
            @simd for i in axes(T, 1)
                value += dw[i, j, k] * T[i, j, k, l]
            end
            dc[l, j, k] += value
        end
    end
    return nothing
end

function _increment_coefficient_adjoint!(dc::AbstractArray{<:Any,2}, dw, T)
    Threads.@threads for j in axes(T, 2)
        @inbounds for l in axes(T, 4)
            value = zero(eltype(dc))
            for k in axes(T, 3)
                @simd for i in axes(T, 1)
                    value += dw[i, j, k] * T[i, j, k, l]
                end
            end
            dc[l, j] += value
        end
    end
    return nothing
end

function _increment_coefficient_adjoint!(dc::AbstractArray{<:Any,1}, dw, T)
    Threads.@threads for l in axes(T, 4)
        value = zero(eltype(dc))
        @inbounds for k in axes(T, 3), j in axes(T, 2)
            @simd for i in axes(T, 1)
                value += dw[i, j, k] * T[i, j, k, l]
            end
        end
        dc[l] += value
    end
    return nothing
end

_array_tangent(x::CoDual) = last(Mooncake.arrayify(x))

function _increment_fused3_adjoint!(dc_tt, dc_t, dc_t_r1,
                                    dw_tt, dw_t, dw_t_r1, T)
    Threads.@threads for l in axes(T, 4)
        @inbounds for j in axes(T, 2)
            t_r1_increment = zero(eltype(dc_t_r1))
            for k in axes(T, 3)
                tt_increment = zero(eltype(dc_tt))
                t_increment = zero(eltype(dc_t))
                @simd for i in axes(T, 1)
                    value = T[i, j, k, l]
                    tt_increment += dw_tt[i, j, k] * value
                    t_increment += dw_t[i, j, k] * value
                    t_r1_increment += dw_t_r1[i, j, k] * value
                end
                dc_tt[l, j, k] += tt_increment
                dc_t[l, j, k] += t_increment
            end
            # `dc_t_r1` can alias the final slice of `dc_t`; increment only
            # after all rank-3 writes so neither contribution is overwritten.
            dc_t_r1[l, j] += t_r1_increment
        end
    end
    return nothing
end

function _increment_fused4_adjoint!(dc_tt, dc_t, dc_t_r1, dc_phi,
                                    dw_tt, dw_t, dw_t_r1, dw_phi, T)
    Threads.@threads for l in axes(T, 4)
        phi_increment = zero(eltype(dc_phi))
        @inbounds for j in axes(T, 2)
            t_r1_increment = zero(eltype(dc_t_r1))
            for k in axes(T, 3)
                tt_increment = zero(eltype(dc_tt))
                t_increment = zero(eltype(dc_t))
                @simd for i in axes(T, 1)
                    value = T[i, j, k, l]
                    tt_increment += dw_tt[i, j, k] * value
                    t_increment += dw_t[i, j, k] * value
                    t_r1_increment += dw_t_r1[i, j, k] * value
                    phi_increment += dw_phi[i, j, k] * value
                end
                dc_tt[l, j, k] += tt_increment
                dc_t[l, j, k] += t_increment
            end
            dc_t_r1[l, j] += t_r1_increment
        end
        dc_phi[l] += phi_increment
    end
    return nothing
end

@is_primitive MinimalCtx Tuple{
    typeof(_w_single_inplace!),
    Array{<:Base.IEEEFloat,3},
    AbstractArray{<:Base.IEEEFloat},
    Array{Float64,4},
}

function rrule!!(
    ::CoDual{typeof(_w_single_inplace!)},
    w::CoDual{<:Array{<:Base.IEEEFloat,3}},
    c::CoDual{<:AbstractArray{<:Base.IEEEFloat}},
    T::CoDual{<:Array{Float64,4}},
)
    pw = primal(w)
    pc = primal(c)
    pT = primal(T)
    old_w = copy(pw)
    _w_single_inplace!(pw, pc, pT)

    function _w_single_pullback!!(::NoRData)
        dw = tangent(w)
        _increment_coefficient_adjoint!(_array_tangent(c), dw, pT)
        fill!(dw, zero(eltype(dw)))
        copyto!(pw, old_w)
        return NoRData(), NoRData(), NoRData(), NoRData()
    end
    return CoDual(nothing, NoFData()), _w_single_pullback!!
end

@is_primitive MinimalCtx Tuple{
    typeof(_w_fused3_inplace!),
    Array{<:Base.IEEEFloat,3}, Array{<:Base.IEEEFloat,3},
    Array{<:Base.IEEEFloat,3},
    Array{<:Base.IEEEFloat,3}, Array{<:Base.IEEEFloat,3},
    AbstractArray{<:Base.IEEEFloat,2},
    Array{Float64,4},
}

function rrule!!(
    ::CoDual{typeof(_w_fused3_inplace!)},
    w_tt::CoDual{<:Array{<:Base.IEEEFloat,3}},
    w_t::CoDual{<:Array{<:Base.IEEEFloat,3}},
    w_t_r1::CoDual{<:Array{<:Base.IEEEFloat,3}},
    c_tt::CoDual{<:Array{<:Base.IEEEFloat,3}},
    c_t::CoDual{<:Array{<:Base.IEEEFloat,3}},
    c_t_r1::CoDual{<:AbstractArray{<:Base.IEEEFloat,2}},
    T::CoDual{<:Array{Float64,4}},
)
    p_w_tt, p_w_t, p_w_t_r1 = primal(w_tt), primal(w_t), primal(w_t_r1)
    p_c_tt, p_c_t, p_c_t_r1 = primal(c_tt), primal(c_t), primal(c_t_r1)
    pT = primal(T)
    old_w_tt, old_w_t, old_w_t_r1 = copy(p_w_tt), copy(p_w_t), copy(p_w_t_r1)
    _w_fused3_inplace!(p_w_tt, p_w_t, p_w_t_r1, p_c_tt, p_c_t, p_c_t_r1, pT)

    function _w_fused3_pullback!!(::NoRData)
        d_w_tt, d_w_t, d_w_t_r1 = tangent(w_tt), tangent(w_t), tangent(w_t_r1)
        _increment_fused3_adjoint!(
            _array_tangent(c_tt), _array_tangent(c_t), _array_tangent(c_t_r1),
            d_w_tt, d_w_t, d_w_t_r1, pT,
        )
        fill!(d_w_tt, zero(eltype(d_w_tt)))
        fill!(d_w_t, zero(eltype(d_w_t)))
        fill!(d_w_t_r1, zero(eltype(d_w_t_r1)))
        copyto!(p_w_tt, old_w_tt)
        copyto!(p_w_t, old_w_t)
        copyto!(p_w_t_r1, old_w_t_r1)
        return ntuple(_ -> NoRData(), 8)
    end
    return CoDual(nothing, NoFData()), _w_fused3_pullback!!
end

@is_primitive MinimalCtx Tuple{
    typeof(_w_fused4_inplace!),
    Array{<:Base.IEEEFloat,3}, Array{<:Base.IEEEFloat,3},
    Array{<:Base.IEEEFloat,3}, Array{<:Base.IEEEFloat,3},
    Array{<:Base.IEEEFloat,3}, Array{<:Base.IEEEFloat,3},
    AbstractArray{<:Base.IEEEFloat,2}, AbstractVector{<:Base.IEEEFloat},
    Array{Float64,4},
}

function rrule!!(
    ::CoDual{typeof(_w_fused4_inplace!)},
    w_tt::CoDual{<:Array{<:Base.IEEEFloat,3}},
    w_t::CoDual{<:Array{<:Base.IEEEFloat,3}},
    w_t_r1::CoDual{<:Array{<:Base.IEEEFloat,3}},
    w_phi::CoDual{<:Array{<:Base.IEEEFloat,3}},
    c_tt::CoDual{<:Array{<:Base.IEEEFloat,3}},
    c_t::CoDual{<:Array{<:Base.IEEEFloat,3}},
    c_t_r1::CoDual{<:AbstractArray{<:Base.IEEEFloat,2}},
    c_phi::CoDual{<:AbstractVector{<:Base.IEEEFloat}},
    T::CoDual{<:Array{Float64,4}},
)
    outputs = (w_tt, w_t, w_t_r1, w_phi)
    p_outputs = map(primal, outputs)
    old_outputs = map(copy, p_outputs)
    pT = primal(T)
    _w_fused4_inplace!(
        p_outputs...,
        primal(c_tt), primal(c_t), primal(c_t_r1), primal(c_phi), pT,
    )

    function _w_fused4_pullback!!(::NoRData)
        d_outputs = map(tangent, outputs)
        _increment_fused4_adjoint!(
            _array_tangent(c_tt), _array_tangent(c_t),
            _array_tangent(c_t_r1), _array_tangent(c_phi),
            d_outputs..., pT,
        )
        for (p_output, old_output, d_output) in zip(p_outputs, old_outputs, d_outputs)
            fill!(d_output, zero(eltype(d_output)))
            copyto!(p_output, old_output)
        end
        return ntuple(_ -> NoRData(), 10)
    end
    return CoDual(nothing, NoFData()), _w_fused4_pullback!!
end

# =============================================================================
# prepare_pk_workspace tensor products and probe-kernel contractions
# =============================================================================

@from_chainrules MinimalCtx Tuple{typeof(_p_phi_TT_tullio), AbstractVector, AbstractMatrix, AbstractArray{<:Any, 3}}
@from_chainrules MinimalCtx Tuple{typeof(_p_phi_T_tullio),  AbstractVector, AbstractArray{<:Any, 3}}

@from_chainrules MinimalCtx Tuple{typeof(_cosmic_shear_kernel_tullio),
    AbstractVector, AbstractVector, AbstractVector,
    AbstractMatrix, AbstractMatrix, Number, Number}
@from_chainrules MinimalCtx Tuple{typeof(_magnification_bias_kernel_tullio),
    AbstractVector, AbstractVector, AbstractVector,
    AbstractMatrix, AbstractMatrix, AbstractMatrix, Number, Number}

end # module MooncakeExt

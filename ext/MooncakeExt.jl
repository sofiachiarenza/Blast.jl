module MooncakeExt

import Mooncake: @from_chainrules, MinimalCtx

import Blast:
    get_clencurt_weights, simpson_weights_array,
    get_clencurt_weights_R_integration, simpson_weights_matrix,
    _combine_kernels_tullio, _compute_Cℓ_tullio, _compute_Cℓ_rsd_tullio,
    _limber_contraction, w_ell_tullio

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

@from_chainrules MinimalCtx Tuple{typeof(_limber_contraction), AbstractArray, AbstractArray, AbstractArray, AbstractVector, Number}

# =============================================================================
# Projected matter density contraction
# =============================================================================

@from_chainrules MinimalCtx Tuple{typeof(w_ell_tullio), AbstractArray{<:Any, 3}, AbstractArray{<:Any, 4}}
@from_chainrules MinimalCtx Tuple{typeof(w_ell_tullio), AbstractArray{<:Any, 2}, AbstractArray{<:Any, 4}}
@from_chainrules MinimalCtx Tuple{typeof(w_ell_tullio), AbstractArray{<:Any, 1}, AbstractArray{<:Any, 4}}

end # module MooncakeExt

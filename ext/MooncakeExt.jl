module MooncakeExt

import Mooncake: @from_chainrules, MinimalCtx

import Blast:
    get_clencurt_weights, simpson_weights_array,
    get_clencurt_weights_R_integration, simpson_weights_matrix,
    _akima_slopes, _akima_coefficients, _akima_eval, _akima_interpolation,
    _combine_kernels_tullio, _compute_Cℓ_tullio, _compute_Cℓ_rsd_tullio,
    _limber_contraction, w_ell_tullio

import AbstractCosmologicalEmulators: chebyshev_decomposition

# =============================================================================
# Chebyshev decomposition (from AbstractCosmologicalEmulators)
# =============================================================================

@from_chainrules MinimalCtx Tuple{typeof(chebyshev_decomposition), Any, AbstractArray}

# =============================================================================
# Quadrature weights (treat as constants — zero-gradient pass-through)
# =============================================================================

@from_chainrules MinimalCtx Tuple{typeof(get_clencurt_weights), Any, Any, Any}
@from_chainrules MinimalCtx Tuple{typeof(simpson_weights_array), Int}
@from_chainrules MinimalCtx Tuple{typeof(get_clencurt_weights_R_integration), Int}
@from_chainrules MinimalCtx Tuple{typeof(simpson_weights_matrix), Int}

# =============================================================================
# Akima interpolation
# =============================================================================

@from_chainrules MinimalCtx Tuple{typeof(_akima_slopes), AbstractVector, AbstractVector}
@from_chainrules MinimalCtx Tuple{typeof(_akima_coefficients), Any, Any}
@from_chainrules MinimalCtx Tuple{typeof(_akima_eval), Any, Any, Any, Any, Any, AbstractArray}
@from_chainrules MinimalCtx Tuple{typeof(_akima_slopes), AbstractMatrix, Any}
@from_chainrules MinimalCtx Tuple{typeof(_akima_coefficients), Any, AbstractMatrix}
@from_chainrules MinimalCtx Tuple{typeof(_akima_eval), AbstractMatrix, Any, AbstractMatrix, AbstractMatrix, AbstractMatrix, Any}
@from_chainrules MinimalCtx Tuple{typeof(_akima_eval), AbstractMatrix, Any, AbstractMatrix, AbstractMatrix, AbstractMatrix, AbstractArray}
@from_chainrules MinimalCtx Tuple{typeof(_akima_interpolation), AbstractMatrix, Any, AbstractArray}

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

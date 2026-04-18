import ChainRulesCore: rrule, NoTangent, unthunk, ProjectTo
import Tullio: @tullio
import FFTW

# =============================================================================
# 1. CONSTANT PATCHES
# =============================================================================

function rrule(::typeof(get_clencurt_weights), kmin, kmax, N)
    return get_clencurt_weights(kmin, kmax, N), _ -> (NoTangent(), NoTangent(), NoTangent(), NoTangent())
end
function rrule(::typeof(simpson_weights_array), n)
    return simpson_weights_array(n), _ -> (NoTangent(), NoTangent())
end
function rrule(::typeof(get_clencurt_weights_R_integration), n)
    return get_clencurt_weights_R_integration(n), _ -> (NoTangent(), NoTangent())
end
function rrule(::typeof(simpson_weights_matrix), n)
    return simpson_weights_matrix(n), _ -> (NoTangent(), NoTangent())
end

# =============================================================================
# 3. KERNEL COMBINATION
# =============================================================================

function rrule(::typeof(_combine_kernels_tullio), W_A::AbstractArray{<:Any, 3}, W_B::AbstractArray{<:Any, 3})
    K = _combine_kernels_tullio(W_A, W_B); W_A_r1 = W_A[:, :, end]; W_B_r1 = W_B[:, :, end]
    function _combine_kernels_tullio_pullback(ȳ)
        ȳ = unthunk(ȳ)
        @tullio ∂W_A[idx_i, idx_c, idx_r] := ȳ[idx_i, idx_j, idx_c, idx_r] * W_B_r1[idx_j, idx_c]
        @tullio ∂W_A_r1[idx_i, idx_c] := ȳ[idx_i, idx_j, idx_c, idx_r] * W_B[idx_j, idx_c, idx_r]
        ∂W_A[:, :, end] .+= ∂W_A_r1
        @tullio ∂W_B[idx_j, idx_c, idx_r] := ȳ[idx_i, idx_j, idx_c, idx_r] * W_A_r1[idx_i, idx_c]
        @tullio ∂W_B_r1[idx_j, idx_c] := ȳ[idx_i, idx_j, idx_c, idx_r] * W_A[idx_i, idx_c, idx_r]
        ∂W_B[:, :, end] .+= ∂W_B_r1
        return (NoTangent(), ∂W_A, ∂W_B)
    end
    return K, _combine_kernels_tullio_pullback
end


# =============================================================================
# 4. FINAL INTEGRATION
# =============================================================================

function rrule(::typeof(_compute_Cℓ_tullio), K::AbstractArray, pmj::AbstractArray, w_χ::AbstractVector, w_R::AbstractVector, prefactor::AbstractVector, Δχ::Number, χ_grid::AbstractVector)
    Cℓ = _compute_Cℓ_tullio(K, pmj, w_χ, w_R, prefactor, Δχ, χ_grid)
    project_K = ProjectTo(K); project_pmj = ProjectTo(pmj)
    function _compute_Cℓ_tullio_pullback(ȳ)
        ȳ = unthunk(ȳ)
        @tullio ∂K[idx_i, idx_j, idx_n, idx_m] := ȳ[idx_l, idx_i, idx_j] * prefactor[idx_l] * χ_grid[idx_n] * pmj[idx_l, idx_n, idx_m] * w_χ[idx_n] * w_R[idx_m] * Δχ
        @tullio ∂pmj[idx_l, idx_n, idx_m] := ȳ[idx_l, idx_i, idx_j] * prefactor[idx_l] * χ_grid[idx_n] * K[idx_i, idx_j, idx_n, idx_m] * w_χ[idx_n] * w_R[idx_m] * Δχ
        return (NoTangent(), project_K(∂K), project_pmj(∂pmj), NoTangent(), NoTangent(), NoTangent(), NoTangent(), NoTangent())
    end
    return Cℓ, _compute_Cℓ_tullio_pullback
end

function rrule(::typeof(_compute_Cℓ_rsd_tullio), W_A_r1::AbstractArray, W_B::AbstractArray, pmj02::AbstractArray, W_A::AbstractArray, W_B_r1::AbstractArray, pmj20::AbstractArray, w_χ::AbstractVector, w_R::AbstractVector, prefactor::AbstractVector, Δχ::Number, χ_grid::AbstractVector)
    Cℓ = _compute_Cℓ_rsd_tullio(W_A_r1, W_B, pmj02, W_A, W_B_r1, pmj20, w_χ, w_R, prefactor, Δχ, χ_grid)
    proj_W_A_r1 = ProjectTo(W_A_r1); proj_W_B = ProjectTo(W_B); proj_pmj02 = ProjectTo(pmj02)
    proj_W_A = ProjectTo(W_A); proj_W_B_r1 = ProjectTo(W_B_r1); proj_pmj20 = ProjectTo(pmj20)
    function _compute_Cℓ_rsd_tullio_pullback(ȳ)
        ȳ = unthunk(ȳ)
        @tullio ∂K[idx_l, idx_i, idx_j, idx_n, idx_m] := ȳ[idx_l, idx_i, idx_j] * prefactor[idx_l] * χ_grid[idx_n] * w_χ[idx_n] * w_R[idx_m] * Δχ
        @tullio ∂W_A_r1[idx_i, idx_c] := ∂K[idx_l, idx_i, idx_j, idx_c, idx_r] * W_B[idx_j, idx_c, idx_r] * pmj02[idx_l, idx_c, idx_r]
        @tullio ∂W_B[idx_j, idx_c, idx_r] := ∂K[idx_l, idx_i, idx_j, idx_c, idx_r] * W_A_r1[idx_i, idx_c] * pmj02[idx_l, idx_c, idx_r]
        @tullio ∂pmj02[idx_l, idx_c, idx_r] := ∂K[idx_l, idx_i, idx_j, idx_c, idx_r] * W_A_r1[idx_i, idx_c] * W_B[idx_j, idx_c, idx_r]
        @tullio ∂W_A[idx_i, idx_c, idx_r] := ∂K[idx_l, idx_i, idx_j, idx_c, idx_r] * W_B_r1[idx_j, idx_c] * pmj20[idx_l, idx_c, idx_r]
        @tullio ∂W_B_r1[idx_j, idx_c] := ∂K[idx_l, idx_i, idx_j, idx_c, idx_r] * W_A[idx_i, idx_c, idx_r] * pmj20[idx_l, idx_c, idx_r]
        @tullio ∂pmj20[idx_l, idx_c, idx_r] := ∂K[idx_l, idx_i, idx_j, idx_c, idx_r] * W_A[idx_i, idx_c, idx_r] * W_B_r1[idx_j, idx_c]
        return (NoTangent(), proj_W_A_r1(∂W_A_r1), proj_W_B(∂W_B), proj_pmj02(∂pmj02), proj_W_A(∂W_A), proj_W_B_r1(∂W_B_r1), proj_pmj20(∂pmj20), NoTangent(), NoTangent(), NoTangent(), NoTangent(), NoTangent())
    end
    return Cℓ, _compute_Cℓ_rsd_tullio_pullback
end


# =============================================================================
# 5. LIMBER integration
# =============================================================================

function rrule(::typeof(_limber_contraction), P_term::AbstractArray, K1::AbstractArray, K2::AbstractArray, weights::AbstractVector, Δχ::Number)
    Cℓ = _limber_contraction(P_term, K1, K2, weights, Δχ)
    proj_P = ProjectTo(P_term); proj_K1 = ProjectTo(K1); proj_K2 = ProjectTo(K2)
    function _limber_contraction_pullback(ȳ)
        ȳ = unthunk(ȳ)
        @tullio ∂P_term[idx_l, idx_m] := ȳ[idx_l, idx_i, idx_j] * K1[idx_l, idx_m, idx_i] * K2[idx_l, idx_m, idx_j] * weights[idx_m] * Δχ
        @tullio ∂K1[idx_l, idx_m, idx_i] := ȳ[idx_l, idx_i, idx_j] * P_term[idx_l, idx_m] * K2[idx_l, idx_m, idx_j] * weights[idx_m] * Δχ
        @tullio ∂K2[idx_l, idx_m, idx_j] := ȳ[idx_l, idx_i, idx_j] * P_term[idx_l, idx_m] * K1[idx_l, idx_m, idx_i] * weights[idx_m] * Δχ
        return (NoTangent(), proj_P(∂P_term), proj_K1(∂K1), proj_K2(∂K2), NoTangent(), NoTangent())
    end
    return Cℓ, _limber_contraction_pullback
end


# w_ell_tullio
function rrule(::typeof(w_ell_tullio), c::AbstractArray{<:Any, 3}, T::AbstractArray{<:Any, 4})
    y  = w_ell_tullio(c, T); pc = ProjectTo(c)
    function w_ell_tullio_pullback_3(ȳ)
        ȳ = unthunk(ȳ); @tullio ∂c[idx_l,idx_j,idx_k] := ȳ[idx_i,idx_j,idx_k] * T[idx_i,idx_j,idx_k,idx_l]
        return (NoTangent(), pc(∂c), NoTangent())
    end
    return (y, w_ell_tullio_pullback_3)
end
function rrule(::typeof(w_ell_tullio), c::AbstractArray{<:Any, 2}, T::AbstractArray{<:Any, 4})
    y  = w_ell_tullio(c, T); pc = ProjectTo(c)
    function w_ell_tullio_pullback_2(ȳ)
        ȳ = unthunk(ȳ); @tullio ∂c[idx_l,idx_j] := ȳ[idx_i,idx_j,idx_k] * T[idx_i,idx_j,idx_k,idx_l]
        return (NoTangent(), pc(∂c), NoTangent())
    end
    return (y, w_ell_tullio_pullback_2)
end
function rrule(::typeof(w_ell_tullio), c::AbstractArray{<:Any, 1}, T::AbstractArray{<:Any, 4})
    y  = w_ell_tullio(c, T); pc = ProjectTo(c)
    function w_ell_tullio_pullback_1(ȳ)
        ȳ = unthunk(ȳ); @tullio ∂c[idx_l] := ȳ[idx_i,idx_j,idx_k] * T[idx_i,idx_j,idx_k,idx_l]
        return (NoTangent(), pc(∂c), NoTangent())
    end
    return (y, w_ell_tullio_pullback_1)
end


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

# limber_eval(c, T_chi, T_k):
#   B[k,j]      = Σ_n c[k,n] * T_chi[j,n]                      (= c * T_chi')
#   result[i,j] = Σ_k T_k[i,k,j] * B[k,j]                      (@tullio contraction)
#
# Mooncake's auto-generated rule on the Tullio tiled-execution path throws a
# MooncakeRuleCompilationError. Registering this ChainRules rrule as a
# primitive via @from_chainrules lets Mooncake skip the body trace entirely.
function rrule(::typeof(limber_eval), c::AbstractMatrix, T_chi::AbstractMatrix, T_k::AbstractArray{<:Any,3})
    B = c * T_chi'
    @tullio result[i, j] := T_k[i, k, j] * B[k, j]
    proj_c = ProjectTo(c); proj_T_chi = ProjectTo(T_chi); proj_T_k = ProjectTo(T_k)
    function limber_eval_pullback(ȳ)
        ȳ = unthunk(ȳ)
        # ∂L/∂B[k,j] = Σ_i T_k[i,k,j] * ȳ[i,j]
        @tullio ∂B[k, j] := T_k[i, k, j] * ȳ[i, j]
        # ∂L/∂T_k[i,k,j] = ȳ[i,j] * B[k,j]
        @tullio ∂T_k[i, k, j] := ȳ[i, j] * B[k, j]
        # ∂L/∂c = ∂B * T_chi   (from B = c * T_chi')
        ∂c = ∂B * T_chi
        # ∂L/∂T_chi = ∂B' * c
        ∂T_chi = ∂B' * c
        return (NoTangent(), proj_c(∂c), proj_T_chi(∂T_chi), proj_T_k(∂T_k))
    end
    return result, limber_eval_pullback
end

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

# =============================================================================
# prepare_pk_workspace tensor products (pure outer products, no contraction).
#
#   P_ϕTT[k,i,j] = P_ϕ[k] * T_χ1[i,k] * T_χR[i,j,k]
#   P_ϕT[k,i,j]  = P_ϕ[k] * T_χR[i,j,k]
# =============================================================================

function rrule(::typeof(_p_phi_TT_tullio),
               P_ϕ::AbstractVector, T_χ1::AbstractMatrix,
               T_χR::AbstractArray{<:Any, 3})
    y = _p_phi_TT_tullio(P_ϕ, T_χ1, T_χR)
    proj_P = ProjectTo(P_ϕ); proj_T1 = ProjectTo(T_χ1); proj_TR = ProjectTo(T_χR)
    function _p_phi_TT_pullback(ȳ)
        ȳ = unthunk(ȳ)
        # ∂L/∂P_ϕ[k]    = Σ_{i,j} ȳ[k,i,j] * T_χ1[i,k] * T_χR[i,j,k]
        @tullio ∂P_ϕ[k]       := ȳ[k, i, j] * T_χ1[i, k] * T_χR[i, j, k]
        # ∂L/∂T_χ1[i,k]  = Σ_j   ȳ[k,i,j] * P_ϕ[k]    * T_χR[i,j,k]
        @tullio ∂T_χ1[i, k]   := ȳ[k, i, j] * P_ϕ[k]    * T_χR[i, j, k]
        # ∂L/∂T_χR[i,j,k] =      ȳ[k,i,j] * P_ϕ[k]    * T_χ1[i,k]
        @tullio ∂T_χR[i, j, k] := ȳ[k, i, j] * P_ϕ[k]    * T_χ1[i, k]
        return (NoTangent(), proj_P(∂P_ϕ), proj_T1(∂T_χ1), proj_TR(∂T_χR))
    end
    return y, _p_phi_TT_pullback
end

function rrule(::typeof(_p_phi_T_tullio),
               P_ϕ::AbstractVector, T_χR::AbstractArray{<:Any, 3})
    y = _p_phi_T_tullio(P_ϕ, T_χR)
    proj_P = ProjectTo(P_ϕ); proj_TR = ProjectTo(T_χR)
    function _p_phi_T_pullback(ȳ)
        ȳ = unthunk(ȳ)
        # ∂L/∂P_ϕ[k]      = Σ_{i,j} ȳ[k,i,j] * T_χR[i,j,k]
        @tullio ∂P_ϕ[k]        := ȳ[k, i, j] * T_χR[i, j, k]
        # ∂L/∂T_χR[i,j,k] =        ȳ[k,i,j] * P_ϕ[k]
        @tullio ∂T_χR[i, j, k] := ȳ[k, i, j] * P_ϕ[k]
        return (NoTangent(), proj_P(∂P_ϕ), proj_TR(∂T_χR))
    end
    return y, _p_phi_T_pullback
end

# =============================================================================
# CosmicShear / MagnificationBias kernel contractions.
#
# Let α = Δχ · prefac / C_LIGHT.  Then
#
#   K[b,i] = α · Σ_p H[p] · S[i,p] · χ[i] · (1+z[i]) · nz[b,p] · (χ[p]-χ[i])/χ[p]
#   K̃[b,i] = α · Σ_p H[p] · S[i,p] · χ[i] · (1+z[i]) · nz[b,p] · (χ[p]-χ[i])/χ[p]
#                                                                 · (5·s[b,p] - 2)
#
# χ appears in two roles (χ[i] as "observer", χ[p] as "source"), so ∂L/∂χ has
# two additive contributions.  `simpson_matrix` is treated as NoTangent.
# =============================================================================

function rrule(::typeof(_cosmic_shear_kernel_tullio),
               H::AbstractVector, χ::AbstractVector, z::AbstractVector,
               nz_norm::AbstractMatrix, simpson_matrix::AbstractMatrix,
               Δχ::Number, prefac::Number)
    K = _cosmic_shear_kernel_tullio(H, χ, z, nz_norm, simpson_matrix, Δχ, prefac)
    proj_H = ProjectTo(H); proj_χ = ProjectTo(χ); proj_z = ProjectTo(z)
    proj_nz = ProjectTo(nz_norm)
    function _cs_kernel_pullback(ȳ)
        ȳ = unthunk(ȳ)
        α = Δχ * prefac / C_LIGHT

        # ∂L/∂nz[b,p] = α · Σ_i ȳ[b,i] · H[p] · S[i,p] · χ[i] · (1+z[i]) · (χ[p]-χ[i])/χ[p]
        @tullio ∂nz[b, p] := α * ȳ[b, i] * H[p] * simpson_matrix[i, p] *
                             χ[i] * (1.0 + z[i]) * (χ[p] - χ[i]) / χ[p]

        # ∂L/∂H[p] = α · Σ_{b,i} ȳ[b,i] · S[i,p] · χ[i] · (1+z[i]) · nz[b,p] · (χ[p]-χ[i])/χ[p]
        @tullio ∂H[p]     := α * ȳ[b, i] * simpson_matrix[i, p] *
                             χ[i] * (1.0 + z[i]) * nz_norm[b, p] *
                             (χ[p] - χ[i]) / χ[p]

        # ∂L/∂z[i] = α · χ[i] · Σ_{b,p} ȳ[b,i] · H[p] · S[i,p] · nz[b,p] · (χ[p]-χ[i])/χ[p]
        @tullio ∂z[i]     := α * ȳ[b, i] * H[p] * simpson_matrix[i, p] *
                             χ[i] * nz_norm[b, p] * (χ[p] - χ[i]) / χ[p]

        # χ dual role. Observer (q=i): d/dχ[q] {χ[q] · (χ[p] - χ[q])/χ[p]}
        #                              = (χ[p] - 2·χ[q]) / χ[p]
        @tullio ∂χ_obs[q] := α * ȳ[b, q] * H[p] * simpson_matrix[q, p] *
                             (1.0 + z[q]) * nz_norm[b, p] *
                             (χ[p] - 2.0 * χ[q]) / χ[p]
        # Source (q=p):    d/dχ[q] {(χ[q] - χ[i])/χ[q]} = χ[i] / χ[q]^2
        @tullio ∂χ_src[q] := α * ȳ[b, i] * H[q] * simpson_matrix[i, q] *
                             χ[i] * (1.0 + z[i]) * nz_norm[b, q] *
                             χ[i] / (χ[q] * χ[q])
        ∂χ = ∂χ_obs .+ ∂χ_src

        # K is linear in Δχ and in prefac → use K itself.
        s_yK = sum(ȳ .* K)
        ∂Δχ     = s_yK / Δχ
        ∂prefac = s_yK / prefac

        return (NoTangent(), proj_H(∂H), proj_χ(∂χ), proj_z(∂z),
                proj_nz(∂nz), NoTangent(), ∂Δχ, ∂prefac)
    end
    return K, _cs_kernel_pullback
end

function rrule(::typeof(_magnification_bias_kernel_tullio),
               H::AbstractVector, χ::AbstractVector, z::AbstractVector,
               nz_norm::AbstractMatrix, s::AbstractMatrix,
               simpson_matrix::AbstractMatrix, Δχ::Number, prefac::Number)
    K = _magnification_bias_kernel_tullio(H, χ, z, nz_norm, s, simpson_matrix, Δχ, prefac)
    proj_H = ProjectTo(H); proj_χ = ProjectTo(χ); proj_z = ProjectTo(z)
    proj_nz = ProjectTo(nz_norm); proj_s = ProjectTo(s)
    function _mb_kernel_pullback(ȳ)
        ȳ = unthunk(ȳ)
        α = Δχ * prefac / C_LIGHT

        # Same pattern as CosmicShear but every partial (except ∂s itself)
        # picks up the extra factor (5·s[b,p] - 2).
        @tullio ∂nz[b, p] := α * ȳ[b, i] * H[p] * simpson_matrix[i, p] *
                             χ[i] * (1.0 + z[i]) * (χ[p] - χ[i]) / χ[p] *
                             (5.0 * s[b, p] - 2)

        @tullio ∂H[p]     := α * ȳ[b, i] * simpson_matrix[i, p] *
                             χ[i] * (1.0 + z[i]) * nz_norm[b, p] *
                             (χ[p] - χ[i]) / χ[p] * (5.0 * s[b, p] - 2)

        @tullio ∂z[i]     := α * ȳ[b, i] * H[p] * simpson_matrix[i, p] *
                             χ[i] * nz_norm[b, p] * (χ[p] - χ[i]) / χ[p] *
                             (5.0 * s[b, p] - 2)

        # ∂L/∂s[b,p]: factor of 5 from d(5·s - 2)/ds.
        @tullio ∂s[b, p]  := α * ȳ[b, i] * H[p] * simpson_matrix[i, p] *
                             χ[i] * (1.0 + z[i]) * nz_norm[b, p] *
                             (χ[p] - χ[i]) / χ[p] * 5.0

        @tullio ∂χ_obs[q] := α * ȳ[b, q] * H[p] * simpson_matrix[q, p] *
                             (1.0 + z[q]) * nz_norm[b, p] *
                             (χ[p] - 2.0 * χ[q]) / χ[p] *
                             (5.0 * s[b, p] - 2)
        @tullio ∂χ_src[q] := α * ȳ[b, i] * H[q] * simpson_matrix[i, q] *
                             χ[i] * (1.0 + z[i]) * nz_norm[b, q] *
                             χ[i] / (χ[q] * χ[q]) *
                             (5.0 * s[b, q] - 2)
        ∂χ = ∂χ_obs .+ ∂χ_src

        s_yK = sum(ȳ .* K)
        ∂Δχ     = s_yK / Δχ
        ∂prefac = s_yK / prefac

        return (NoTangent(), proj_H(∂H), proj_χ(∂χ), proj_z(∂z),
                proj_nz(∂nz), proj_s(∂s), NoTangent(), ∂Δχ, ∂prefac)
    end
    return K, _mb_kernel_pullback
end


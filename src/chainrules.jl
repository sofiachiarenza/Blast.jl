import ChainRulesCore: rrule, NoTangent, unthunk, ProjectTo
import Tullio: @tullio
import FFTW
import Mooncake
import Mooncake: @from_chainrules, MinimalCtx
import Mooncake: tangent_type, fdata_type, rdata_type, zero_tangent_internal, fdata, rdata, NoFData, NoRData, increment_rdata!!

# =============================================================================
# 1. CONSTANT PATCHES
# =============================================================================

# These patches (for FFTWPlan and ChebyshevPlan) are provided by AbstractCosmologicalEmulators.
# Removal fixes method overwriting error.
Mooncake.@from_chainrules MinimalCtx Tuple{typeof(chebyshev_decomposition), Any, AbstractArray}

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

@from_chainrules MinimalCtx Tuple{typeof(get_clencurt_weights), Any, Any, Any}
@from_chainrules MinimalCtx Tuple{typeof(simpson_weights_array), Int}
@from_chainrules MinimalCtx Tuple{typeof(get_clencurt_weights_R_integration), Int}
@from_chainrules MinimalCtx Tuple{typeof(simpson_weights_matrix), Int}

# =============================================================================
# 2. AKIMA INTERPOLATION
# =============================================================================

function rrule(::typeof(_akima_slopes), u::AbstractVector, t::AbstractVector)
    n = length(u); dt = diff(t); m = zeros(eltype(u), n + 3)
    m[3:(n+1)] .= diff(u) ./ dt; m[2] = 2m[3] - m[4]; m[1] = 2m[2] - m[3]
    m[n+2] = 2m[n+1] - m[n]; m[n+3] = 2m[n+2] - m[n+1]
    function _akima_slopes_pullback(Δm)
        gm = collect(unthunk(Δm))
        gm[n+2] += 2gm[n+3]; gm[n+1] -= gm[n+3]; gm[n+1] += 2gm[n+2]; gm[n] -= gm[n+2]
        gm[2] += 2gm[1]; gm[3] -= gm[1]; gm[3] += 2gm[2]; gm[4] -= gm[2]
        sm_bar = gm[3:(n+1)]; δu = zero(u); δt = zero(t)
        @inbounds for i in 1:n-1
            g = sm_bar[i]; invdt = 1 / dt[i]
            δu[i] -= g * invdt; δu[i+1] += g * invdt
            diffu = u[i+1] - u[i]; invdt2 = invdt^2
            δt[i] += g * diffu * invdt2; δt[i+1] -= g * diffu * invdt2
        end
        return (NoTangent(), δu, δt)
    end
    return m, _akima_slopes_pullback
end

function rrule(::typeof(_akima_coefficients), t, m)
    n = length(t); dt = diff(t); dm = abs.(diff(m)); f1 = dm[3:(n+2)]; f2 = dm[1:n]; f12 = f1 + f2
    b = (m[4:end] .+ m[1:(end-3)]) ./ 2; eps_akima = eps(eltype(f12)) * 100; use_weighted = f12 .> eps_akima
    for i in eachindex(f12)
        if use_weighted[i]
            b[i] = (f1[i] * m[i+1] + f2[i] * m[i+2]) / f12[i]
        end
    end
    c = (3 .* m[3:(end-2)] .- 2 .* b[1:(end-1)] .- b[2:end]) ./ dt
    d = (b[1:(end-1)] .+ b[2:end] .- 2 .* m[3:(end-2)]) ./ dt .^ 2
    function _akima_coefficients_pullback(Δ)
        Δb, Δc, Δd = unthunk(Δ); ∂t = zeros(eltype(t), length(t)); ∂m = zeros(eltype(m), length(m)); ∂b_acc = zeros(eltype(b), length(b))
        dt_inv = 1.0 ./ dt; dt_inv_sq = dt_inv .^ 2
        if Δd !== nothing
            @. ∂b_acc[1:(end-1)] += Δd * dt_inv_sq; @. ∂b_acc[2:end] += Δd * dt_inv_sq; @. ∂m[3:(end-2)] -= 2.0 * Δd * dt_inv_sq
            ∂dt_d = @. -2.0 * Δd * (b[1:(end-1)] + b[2:end] - 2.0 * m[3:(end-2)]) * dt_inv_sq * dt_inv
            @. ∂t[1:(end-1)] -= ∂dt_d; @. ∂t[2:end] += ∂dt_d
        end
        if Δc !== nothing
            @. ∂m[3:(end-2)] += 3.0 * Δc * dt_inv; @. ∂b_acc[1:(end-1)] -= 2.0 * Δc * dt_inv; @. ∂b_acc[2:end] -= Δc * dt_inv
            ∂dt_c = @. -Δc * (3.0 * m[3:(end-2)] - 2.0 * b[1:(end-1)] - b[2:end]) * dt_inv^2
            @. ∂t[1:(end-1)] -= ∂dt_c; @. ∂t[2:end] += ∂dt_c
        end
        if Δb !== nothing; @. ∂b_acc += Δb; end
        if any(!iszero, ∂b_acc)
            ∂f1 = zeros(length(f1)); ∂f2 = zeros(length(f2)); ∂f12 = zeros(length(f12))
            for i in eachindex(use_weighted)
                if use_weighted[i]
                    f12_inv = 1.0 / f12[i]
                    ∂f1[i] += ∂b_acc[i] * m[i+1] * f12_inv; ∂f2[i] += ∂b_acc[i] * m[i+2] * f12_inv
                    ∂m[i+1] += ∂b_acc[i] * f1[i] * f12_inv; ∂m[i+2] += ∂b_acc[i] * f2[i] * f12_inv
                    ∂f12[i] += -∂b_acc[i] * (f1[i] * m[i+1] + f2[i] * m[i+2]) * f12_inv^2
                else
                    ∂m[i+3] += ∂b_acc[i] / 2; ∂m[i] += ∂b_acc[i] / 2
                end
            end
            @. ∂f1 += ∂f12; @. ∂f2 += ∂f12; ∂dm = zeros(length(dm)); @. ∂dm[3:(n+2)] += ∂f1; @. ∂dm[1:n] += ∂f2
            ∂diff_m = ∂dm .* sign.(diff(m)); @. ∂m[1:(end-1)] -= ∂diff_m; @. ∂m[2:end] += ∂diff_m
        end
        return (NoTangent(), ∂t, ∂m)
    end
    return (b, c, d), _akima_coefficients_pullback
end

function rrule(::typeof(_akima_eval), u, t, b, c, d, tq::AbstractArray)
    n_q = length(tq); T = promote_type(eltype(u), eltype(t), eltype(b), eltype(c), eltype(d), eltype(tq))
    res = zeros(T, n_q)
    @inbounds for i in eachindex(tq)
        idx = _akima_find_interval(t, tq[i]); w_val = tq[i] - t[idx]
        res[i] = ((d[idx] * w_val + c[idx]) * w_val + b[idx]) * w_val + u[idx]
    end
    function _akima_eval_pullback(ȳ)
        ȳ_un = unthunk(ȳ); ū_tot = zero(u); t̄_tot = zero(t); b̄_tot = zero(b); c̄_tot = zero(c); d̄_tot = zero(d)
        tq̄ = similar(tq, promote_type(eltype(ȳ_un), eltype(tq)))
        @inbounds for i in eachindex(tq)
            ȳ_i = ȳ_un[i]
            if !iszero(ȳ_i)
                idx = _akima_find_interval(t, tq[i]); w_val = tq[i] - t[idx]; w_val_sq = w_val * w_val
                dw_val = 3 * d[idx] * w_val_sq + 2 * c[idx] * w_val + b[idx]
                ū_tot[idx] += ȳ_i; t̄_tot[idx] -= ȳ_i * dw_val; tq̄[i] = ȳ_i * dw_val
                b̄_tot[idx] += ȳ_i * w_val; c̄_tot[idx] += ȳ_i * w_val_sq; d̄_tot[idx] += ȳ_i * w_val * w_val_sq
            else
                tq̄[i] = zero(eltype(tq̄))
            end
        end
        return NoTangent(), ū_tot, t̄_tot, b̄_tot, c̄_tot, d̄_tot, tq̄
    end
    return res, _akima_eval_pullback
end

# Fused Matrix Rule
function rrule(::typeof(_akima_interpolation), u::AbstractMatrix, t, t_new::AbstractArray)
    n, n_cols = size(u); dt = diff(t); m = zeros(promote_type(eltype(u), eltype(t)), n + 3, n_cols)
    for col in 1:n_cols
        m[3:(end-2), col] .= diff(view(u, :, col)) ./ dt
        m[2, col] = 2m[3, col] - m[4, col]; m[1, col] = 2m[2, col] - m[3, col]
        m[n+2, col] = 2m[n+1, col] - m[n, col]; m[n+3, col] = 2m[n+2, col] - m[n+1, col]
    end
    n_q = length(t_new); n_c = n_cols
    b = zeros(eltype(m), n, n_c); c = zeros(eltype(m), n-1, n_c); d = zeros(eltype(m), n-1, n_c)
    use_w = falses(n, n_c); eps_akima = eps(eltype(m)) * 100
    for col in 1:n_c
        b[:, col] = (view(m, 4:(n+3), col) .+ view(m, 1:n, col)) ./ 2
        dm = abs.(diff(view(m, :, col))); f1 = view(dm, 3:(n+2)); f2 = view(dm, 1:n); f12 = f1 .+ f2
        for i in 1:n
            if f12[i] > eps_akima
                b[i, col] = (f1[i] * m[i+1, col] + f2[i] * m[i+2, col]) / f12[i]
                use_w[i, col] = true
            end
        end
        c[:, col] = (3 .* view(m, 3:(n+1), col) .- 2 .* view(b, 1:(n-1), col) .- view(b, 2:n, col)) ./ dt
        d[:, col] = (view(b, 1:(n-1), col) .+ view(b, 2:n, col) .- 2 .* view(m, 3:(n+1), col)) ./ dt .^ 2
    end
    res = zeros(promote_type(eltype(u), eltype(t_new)), n_q, n_c)
    @inbounds for i in 1:n_q
        idx = _akima_find_interval(t, t_new[i]); w_val = t_new[i] - t[idx]
        for col in 1:n_c
            res[i, col] = ((d[idx, col] * w_val + c[idx, col]) * w_val + b[idx, col]) * w_val + u[idx, col]
        end
    end
    function _akima_interpolation_matrix_fused_pullback(ȳ)
        ȳ_un = unthunk(ȳ); ∂u = zero(u); ∂t = zero(t); ∂t_new = zeros(n_q)
        ∂b = zero(b); ∂c = zero(c); ∂d = zero(d); ∂m = zeros(n+3, n_c)
        @inbounds for i in 1:n_q
            idx = _akima_find_interval(t, t_new[i]); w_val = t_new[i] - t[idx]; w_val_sq = w_val * w_val; w_val_cb = w_val * w_val_sq
            for col in 1:n_c
                ȳ_ic = ȳ_un[i, col]
                if !iszero(ȳ_ic)
                    dw_val = 3 * d[idx, col] * w_val_sq + 2 * c[idx, col] * w_val + b[idx, col]
                    ∂u[idx, col] += ȳ_ic; ∂t[idx] -= ȳ_ic * dw_val; ∂t_new[i] += ȳ_ic * dw_val
                    ∂b[idx, col] += ȳ_ic * w_val; ∂c[idx, col] += ȳ_ic * w_val_sq; ∂d[idx, col] += ȳ_ic * w_val_cb
                end
            end
        end
        dt_inv = 1.0 ./ dt; dt_inv_sq = dt_inv .^ 2
        for col in 1:n_c
            @. ∂b[1:(end-1), col] += ∂d[:, col] * dt_inv_sq; @. ∂b[2:end, col] += ∂d[:, col] * dt_inv_sq
            @. ∂m[3:(end-2), col] -= 2.0 * ∂d[:, col] * dt_inv_sq
            ∂dt_d = @. -2.0 * ∂d[:, col] * (b[1:(end-1), col] + b[2:end, col] - 2.0 * m[3:(end-2), col]) * dt_inv_sq * dt_inv
            @. ∂t[1:(end-1)] -= ∂dt_d; @. ∂t[2:end] += ∂dt_d
            @. ∂m[3:(end-2), col] += 3.0 * ∂c[:, col] * dt_inv
            @. ∂b[1:(end-1), col] -= 2.0 * ∂c[:, col] * dt_inv; @. ∂b[2:end, col] -= ∂c[:, col] * dt_inv
            ∂dt_c = @. -∂c[:, col] * (3.0 * m[3:(end-2), col] - 2.0 * b[1:(end-1), col] - b[2:end, col]) * dt_inv^2
            @. ∂t[1:(end-1)] -= ∂dt_c; @. ∂t[2:end] += ∂dt_c
            diff_m_col = diff(view(m, :, col)); dm_col = abs.(diff_m_col)
            f1_c = view(dm_col, 3:(n+2)); f2_c = view(dm_col, 1:n); f12_c = f1_c .+ f2_c
            ∂f1 = zeros(n); ∂f2 = zeros(n); ∂f12 = zeros(n)
            for i in 1:n
                if use_w[i, col]
                    f12_inv = 1.0 / f12_c[i]
                    ∂f1[i] += ∂b[i, col] * m[i+1, col] * f12_inv; ∂f2[i] += ∂b[i, col] * m[i+2, col] * f12_inv
                    ∂m[i+1, col] += ∂b[i, col] * f1_c[i] * f12_inv; ∂m[i+2, col] += ∂b[i, col] * f2_c[i] * f12_inv
                    ∂f12[i] += -∂b[i, col] * (f1_c[i] * m[i+1, col] + f2_c[i] * m[i+2, col]) * f12_inv^2
                else
                    ∂m[i+3, col] += ∂b[i, col] / 2; ∂m[i, col] += ∂b[i, col] / 2
                end
            end
            @. ∂f1 += ∂f12; @. ∂f2 += ∂f12; ∂dm = zeros(n+2)
            @. ∂dm[3:(n+2)] += ∂f1; @. ∂dm[1:n] += ∂f2
            ∂diff_m = ∂dm .* sign.(diff_m_col)
            @. ∂m[1:(end-1), col] -= ∂diff_m; @. ∂m[2:end, col] += ∂diff_m
            ∂m[n+2, col] += 2∂m[n+3, col]; ∂m[n+1, col] -= ∂m[n+3, col]; ∂m[n+1, col] += 2∂m[n+2, col]; ∂m[n, col] -= ∂m[n+2, col]
            ∂m[2, col] += 2∂m[1, col]; ∂m[3, col] -= ∂m[1, col]; ∂m[3, col] += 2∂m[2, col]; ∂m[4, col] -= ∂m[2, col]
            sm_b = view(∂m, 3:(n+1), col)
            @inbounds for i in 1:n-1
                g = sm_b[i]; invdt = 1 / dt[i]; ∂u[i, col] -= g * invdt; ∂u[i+1, col] += g * invdt
                diffu = u[i+1, col] - u[i, col]; invdt2 = invdt^2
                ∂t[i] += g * diffu * invdt2; ∂t[i+1] -= g * diffu * invdt2
            end
        end
        return (NoTangent(), ∂u, ∂t, ∂t_new)
    end
    return res, _akima_interpolation_matrix_fused_pullback
end

@from_chainrules MinimalCtx Tuple{typeof(_akima_slopes), AbstractVector, AbstractVector}
@from_chainrules MinimalCtx Tuple{typeof(_akima_coefficients), Any, Any}
@from_chainrules MinimalCtx Tuple{typeof(_akima_eval), Any, Any, Any, Any, Any, AbstractArray}
@from_chainrules MinimalCtx Tuple{typeof(_akima_slopes), AbstractMatrix, Any}
@from_chainrules MinimalCtx Tuple{typeof(_akima_coefficients), Any, AbstractMatrix}
@from_chainrules MinimalCtx Tuple{typeof(_akima_eval), AbstractMatrix, Any, AbstractMatrix, AbstractMatrix, AbstractMatrix, Any}
@from_chainrules MinimalCtx Tuple{typeof(_akima_eval), AbstractMatrix, Any, AbstractMatrix, AbstractMatrix, AbstractMatrix, AbstractArray}
@from_chainrules MinimalCtx Tuple{typeof(_akima_interpolation), AbstractMatrix, Any, AbstractArray}

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

@from_chainrules MinimalCtx Tuple{typeof(_combine_kernels_tullio), AbstractArray{<:Any, 3}, AbstractArray{<:Any, 3}}

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

@from_chainrules MinimalCtx Tuple{typeof(_compute_Cℓ_tullio), AbstractArray, AbstractArray, AbstractVector, AbstractVector, AbstractVector, Number, AbstractVector}
@from_chainrules MinimalCtx Tuple{typeof(_compute_Cℓ_rsd_tullio), AbstractArray, AbstractArray, AbstractArray, AbstractArray, AbstractArray, AbstractArray, AbstractVector, AbstractVector, AbstractVector, Number, AbstractVector}

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

@from_chainrules MinimalCtx Tuple{typeof(_limber_contraction), AbstractArray, AbstractArray, AbstractArray, AbstractVector, Number}

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

@from_chainrules MinimalCtx Tuple{typeof(w_ell_tullio), AbstractArray{<:Any, 3}, AbstractArray{<:Any, 4}}
@from_chainrules MinimalCtx Tuple{typeof(w_ell_tullio), AbstractArray{<:Any, 2}, AbstractArray{<:Any, 4}}
@from_chainrules MinimalCtx Tuple{typeof(w_ell_tullio), AbstractArray{<:Any, 1}, AbstractArray{<:Any, 4}}

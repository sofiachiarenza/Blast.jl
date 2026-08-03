using Blast
using Enzyme
using ForwardDiff
using Reactant
using Random

function pullback_contract(P, K1, K2, upstream, weights, Δχ)
    loss(P_, K1_, K2_, U_, w_) = sum(
        Blast.reactant_limber_contraction(P_, K1_, K2_, w_, Δχ) .* U_,
    )
    grad_fun(P_, K1_, K2_, U_, w_) = Enzyme.gradient(
        Reverse, loss, P_, Const(K1_), Const(K2_), Const(U_), Const(w_),
    )[1]
    args = map(Reactant.to_rarray, (P, K1, K2, upstream, weights))
    compiled = Reactant.@compile sync=true grad_fun(args...)
    result = compiled(args...)
    Reactant.synchronize(result)
    return Array(result)
end

function pullback_grid(coefficients, T_z, T_k, upstream)
    loss(c, Tz, Tk, U) = sum(
        Blast.reactant_limber_grid_from_coefficients(c, Tz, Tk) .* U,
    )
    grad_fun(c, Tz, Tk, U) = Enzyme.gradient(
        Reverse, loss, c, Const(Tz), Const(Tk), Const(U),
    )[1]
    args = map(Reactant.to_rarray, (coefficients, T_z, T_k, upstream))
    compiled = Reactant.@compile sync=true grad_fun(args...)
    result = compiled(args...)
    Reactant.synchronize(result)
    return Array(result)
end

function pullback_coeff(logP, Mk, Mz, upstream)
    loss(v, K, Z, U) = sum(Blast.reactant_limber_chebyshev_coefficients(v, K, Z) .* U)
    grad_fun(v, K, Z, U) = Enzyme.gradient(
        Reverse, loss, v, Const(K), Const(Z), Const(U),
    )[1]
    args = map(Reactant.to_rarray, (logP, Mk, Mz, upstream))
    compiled = Reactant.@compile sync=true grad_fun(args...)
    result = compiled(args...)
    Reactant.synchronize(result)
    return Array(result)
end

function main()
    Reactant.set_default_backend("cpu")
    Random.seed!(4567)
    nk, nzc, nz, nl, nA, nB = 3, 2, 2, 3, 2, 2
    log_lin = randn(nk, nzc)
    log_non = randn(nk, nzc)
    Mk, Mz = randn(2, nk), randn(2, nzc)
    Tz, Tk = randn(nz, 2), randn(nl, 2, nz)
    weights, invχ2 = rand(nz), rand(nz)
    K_low = randn(1, nz, nA)
    K_high = randn(2, nz, nA)
    U_low = randn(1, nA, nB)
    U_high = randn(2, nA, nB)
    Δχ = 0.4

    function primal(a, b)
        c_lin = Blast.reactant_limber_chebyshev_coefficients(a, Mk, Mz)
        c_non = Blast.reactant_limber_chebyshev_coefficients(b, Mk, Mz)
        P_lin = Blast.reactant_limber_grid_from_coefficients(c_lin, Tz, Tk)
        P_non = Blast.reactant_limber_grid_from_coefficients(c_non, Tz, Tk)
        ΔP = P_non .- P_lin
        C_low = Blast.reactant_limber_contraction(ΔP[1:1, :] .* reshape(invχ2, 1, :), K_low, K_low, weights, Δχ)
        C_high = Blast.reactant_limber_contraction(P_non[2:3, :] .* reshape(invχ2, 1, :), K_high, K_high, weights, Δχ)
        return sum(C_low .* U_low) + sum(C_high .* U_high)
    end
    ref_a = ForwardDiff.gradient(a -> primal(a, log_non), log_lin)
    ref_b = ForwardDiff.gradient(b -> primal(log_lin, b), log_non)

    c_lin = Blast.reactant_limber_chebyshev_coefficients(log_lin, Mk, Mz)
    c_non = Blast.reactant_limber_chebyshev_coefficients(log_non, Mk, Mz)
    P_lin = Blast.reactant_limber_grid_from_coefficients(c_lin, Tz, Tk)
    P_non = Blast.reactant_limber_grid_from_coefficients(c_non, Tz, Tk)
    dC_low = pullback_contract(P_non[1:1, :] .* reshape(invχ2, 1, :), K_low, K_low, U_low, weights, Δχ)
    dC_high = pullback_contract(P_non[2:3, :] .* reshape(invχ2, 1, :), K_high, K_high, U_high, weights, Δχ)
    dP_non = vcat(dC_low .* reshape(invχ2, 1, :), dC_high .* reshape(invχ2, 1, :))
    dP_lin = vcat(-dC_low .* reshape(invχ2, 1, :), zeros(2, nz))
    dc_non = pullback_grid(c_non, Tz, Tk, dP_non)
    dc_lin = pullback_grid(c_lin, Tz, Tk, dP_lin)
    dlog_non = pullback_coeff(log_non, Mk, Mz, dc_non)
    dlog_lin = pullback_coeff(log_lin, Mk, Mz, dc_lin)

    println("Reactant=", pkgversion(Reactant))
    println("staged_log_lin_gradient_error=", maximum(abs, dlog_lin .- ref_a))
    println("staged_log_non_gradient_error=", maximum(abs, dlog_non .- ref_b))
end

main()

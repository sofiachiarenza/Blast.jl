"""Staged scalar-loss pullbacks for the Reactant Blast path.

This file deliberately keeps the reverse driver outside Blast's core package:
Enzyme is an optional AD dependency, while the forward Reactant extension must
remain usable without it. Each pullback is compiled independently and the
returned cotangents are ordinary arrays passed to the next stage.
"""

using Blast
using Enzyme
using Reactant

struct LimberPullbackState{A,B,C,D,E,F,G,H,I,J}
    logP_lin::A
    logP_nonlin::B
    coeff_lin::C
    coeff_nonlin::D
    P_lin::E
    P_nonlin::F
    dlogP_lin::G
    dlogP_nonlin::H
    dP_lin::I
    dP_nonlin::J
end

const _PULLBACK_CACHE = Dict{Any,Any}()

function _compile_pullback(f, args; key=nothing)
    if key !== nothing && haskey(_PULLBACK_CACHE, key)
        return _PULLBACK_CACHE[key]
    end
    compiled = Reactant.@compile sync=true f(args...)
    key === nothing || (_PULLBACK_CACHE[key] = compiled)
    return compiled
end

function _pullback_contract(P, K1, K2, upstream, weights, Δχ)
    loss(P_, K1_, K2_, U_, weights_) = sum(
        Blast.reactant_limber_contraction(P_, K1_, K2_, weights_, Δχ) .* U_,
    )
    grad(P_, K1_, K2_, U_, weights_) = Enzyme.gradient(
        Reverse, loss, P_, Const(K1_), Const(K2_), Const(U_), Const(weights_),
    )[1]
    args = map(Reactant.to_rarray, (P, K1, K2, upstream, weights))
    compiled = _compile_pullback(grad, args; key=(:contract, size(P), size(K1), size(K2)))
    return Array(compiled(args...))
end

function _pullback_grid(coefficients, T_z, T_k, upstream)
    loss(c, Tz, Tk, U) = sum(
        Blast.reactant_limber_grid_from_coefficients(c, Tz, Tk) .* U,
    )
    grad(c, Tz, Tk, U) = Enzyme.gradient(
        Reverse, loss, c, Const(Tz), Const(Tk), Const(U),
    )[1]
    args = map(Reactant.to_rarray, (coefficients, T_z, T_k, upstream))
    compiled = _compile_pullback(grad, args; key=(:grid, size(coefficients), size(T_z), size(T_k)))
    return Array(compiled(args...))
end

function _pullback_coeff(logP, transform_k, transform_z, upstream)
    loss(values, Mk, Mz, U) = sum(
        Blast.reactant_limber_chebyshev_coefficients(values, Mk, Mz) .* U,
    )
    grad(values, Mk, Mz, U) = Enzyme.gradient(
        Reverse, loss, values, Const(Mk), Const(Mz), Const(U),
    )[1]
    args = map(Reactant.to_rarray, (logP, transform_k, transform_z, upstream))
    compiled = _compile_pullback(grad, args; key=(:coeff, size(logP), size(transform_k), size(transform_z)))
    return Array(compiled(args...))
end

"""Reverse the split Limber branch for a scalar loss.

`dC_correction` and `dC_limber` are cotangents supplied by the final
`C_l`/likelihood pullback. The returned pair is the cotangent with respect to
the logarithmic linear and nonlinear Limber power grids. The caller then uses
the power-product pullback to obtain cotangents with respect to the two input
`P(k,z)` arrays.
"""
function staged_limber_pullback(
    state::LimberPullbackState,
    dC_correction,
    dC_limber,
    K1_low,
    K2_low,
    K1_high,
    K2_high,
    weights,
    inv_χ2,
    Δχ,
    transform_k,
    transform_z,
    T_z,
    T_k,
)
    nlow = size(K1_low, 1)
    dP_correction = _pullback_contract(
        (state.P_nonlin .- state.P_lin)[1:nlow, :] .* reshape(inv_χ2, 1, :),
        K1_low, K2_low, dC_correction, weights, Δχ,
    )
    dP_limber = _pullback_contract(
        state.P_nonlin[(nlow + 1):end, :] .* reshape(inv_χ2, 1, :),
        K1_high, K2_high, dC_limber, weights, Δχ,
    )

    dP_nonlin = vcat(
        dP_correction .* reshape(inv_χ2, 1, :),
        dP_limber .* reshape(inv_χ2, 1, :),
    )
    dP_lin = vcat(
        -dP_correction .* reshape(inv_χ2, 1, :),
        zeros(size(dP_limber)),
    )
    dcoeff_lin = _pullback_grid(state.coeff_lin, T_z, T_k, dP_lin)
    dcoeff_nonlin = _pullback_grid(state.coeff_nonlin, T_z, T_k, dP_nonlin)
    dlogP_lin = _pullback_coeff(state.logP_lin, transform_k, transform_z, dcoeff_lin)
    dlogP_nonlin = _pullback_coeff(state.logP_nonlin, transform_k, transform_z, dcoeff_nonlin)
    return dlogP_lin, dlogP_nonlin
end

"""Chain the Limber endpoint cotangent to both Limber input spectra."""
function staged_limber_gradient(
    pk_lin,
    pk_nonlin,
    background,
    state,
    dC_correction,
    dC_limber,
    K1_low,
    K2_low,
    K1_high,
    K2_high,
    weights,
    inv_χ2,
    Δχ,
    transform_k,
    transform_z,
    T_z,
    T_k,
)
    dlog_lin, dlog_nonlin = staged_limber_pullback(
        state, dC_correction, dC_limber, K1_low, K2_low, K1_high, K2_high,
        weights, inv_χ2, Δχ, transform_k, transform_z, T_z, T_k,
    )
    return _pullback_power_products(
        pk_lin, pk_nonlin, background, dlog_lin, dlog_nonlin,
    )
end

function _pullback_compute_w(cϕTT, cϕT, cϕ, T_tuple, dW_tuple)
    loss(cTT, cT, c0, T1, T2, T3, T4, T5, T6, T7, T8,
         d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12, d13, d14, d15, d16, d17) = begin
        W = Blast.reactant_compute_w(cTT, cT, c0, T1, T2, T3, T4, T5, T6, T7, T8)
        sum(sum, ntuple(i -> W[i] .* (i == 1 ? d1 : i == 2 ? d2 : i == 3 ? d3 :
            i == 4 ? d4 : i == 5 ? d5 : i == 6 ? d6 : i == 7 ? d7 : i == 8 ? d8 :
            i == 9 ? d9 : i == 10 ? d10 : i == 11 ? d11 : i == 12 ? d12 :
            i == 13 ? d13 : i == 14 ? d14 : i == 15 ? d15 : i == 16 ? d16 : d17), 17))
    end
    grad(cTT, cT, c0, T1, T2, T3, T4, T5, T6, T7, T8,
         d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12, d13, d14, d15, d16, d17) = Enzyme.gradient(
        Reverse, loss, cTT, cT, c0,
        Const(T1), Const(T2), Const(T3), Const(T4), Const(T5), Const(T6), Const(T7), Const(T8),
        Const(d1), Const(d2), Const(d3), Const(d4), Const(d5), Const(d6), Const(d7), Const(d8),
        Const(d9), Const(d10), Const(d11), Const(d12), Const(d13), Const(d14), Const(d15), Const(d16), Const(d17),
    )
    args = (
        Reactant.to_rarray(cϕTT), Reactant.to_rarray(cϕT), Reactant.to_rarray(cϕ),
        map(Reactant.to_rarray, T_tuple)..., map(Reactant.to_rarray, dW_tuple)...,
    )
    compiled = _compile_pullback(grad, args; key=(:compute_w, size(cϕTT), size(cϕT), size(cϕ)))
    result = compiled(args...)
    return Array(result[1]), Array(result[2]), Array(result[3])
end

function _pullback_nonlimber_prepare(pk, background, transform, dcoeffs)
    loss(pk_, M, d1, d2, d3) = begin
        c = Blast.reactant_prepare_nonlimber_spectrum(pk_, background, M)
        sum(c[1] .* d1) + sum(c[2] .* d2) + sum(c[3] .* d3)
    end
    grad(pk_, M, d1, d2, d3) = Enzyme.gradient(
        Reverse, loss, pk_, Const(M), Const(d1), Const(d2), Const(d3),
    )[1]
    args = (
        Reactant.to_rarray(pk), Reactant.to_rarray(transform),
        map(Reactant.to_rarray, dcoeffs)...,
    )
    compiled = _compile_pullback(grad, args; key=(:prepare, size(pk), size(transform)))
    return Array(compiled(args...))
end

function _pullback_power_products(pk_lin, pk_nonlin, background, dlog_lin, dlog_nonlin)
    loss(pk1, pk2, U1, U2) = begin
        log1, log2 = Blast.reactant_limber_power_products(pk1, pk2, background)
        sum(log1 .* U1) + sum(log2 .* U2)
    end
    grad(pk1, pk2, U1, U2) = Enzyme.gradient(
        Reverse, loss, pk1, pk2, Const(U1), Const(U2),
    )
    args = (
        Reactant.to_rarray(pk_lin), Reactant.to_rarray(pk_nonlin),
        Reactant.to_rarray(dlog_lin), Reactant.to_rarray(dlog_nonlin),
    )
    compiled = _compile_pullback(grad, args; key=(:power, size(pk_lin), size(pk_nonlin)))
    result = compiled(args...)
    return Array(result[1]), Array(result[2])
end

function _endpoint_value(w_tuple, config)
    return Blast.reactant_full_3x2pt(
        w_tuple...,
        config.kernels..., config.integration...,
        config.C_terms..., config.finalization..., config.Δχ, config.pref,
    )
end

function _nonlimber_value(w_tuple, config)
    return Blast.reactant_nonlimber_3x2pt(
        w_tuple...,
        config.kernels..., config.integration..., config.Δχ, config.pref,
    )
end

"""Return all 17 W cotangents for a scalar C_l/likelihood cotangent.

Each W component is active in its own compiled Enzyme pullback. This is
deliberately staged: one reverse graph with all 17 active W tensors is much
larger and provides no useful performance benefit for the scalar likelihood
case.
"""
function _staged_endpoint_pullback(w_tuple, config, upstream, evaluator)
    w_r = map(Reactant.to_rarray, w_tuple)
    cfg_r = merge(config, (
        kernels=map(Reactant.to_rarray, config.kernels),
        integration=map(Reactant.to_rarray, config.integration),
        C_terms=map(Reactant.to_rarray, config.C_terms),
        finalization=map(Reactant.to_rarray, config.finalization),
    ))
    upstream_r = map(Reactant.to_rarray, upstream)
    result = ntuple(17) do active
        other = ntuple(i -> w_r[i], 17)
        loss(w_active, other_, cfg_, upstream_) = begin
            w = ntuple(i -> i == active ? w_active : other_[i], 17)
            y = evaluator(w, cfg_)
            sum(y[1] .* upstream_[1]) + sum(y[2] .* upstream_[2]) + sum(y[3] .* upstream_[3])
        end
        grad(w_active, other_, cfg_, upstream_) = Enzyme.gradient(
            Reverse, loss, w_active, Const(other_), Const(cfg_), Const(upstream_),
        )[1]
        compiled = _compile_pullback(
            grad, (w_r[active], other, cfg_r, upstream_r);
            key=(:endpoint, evaluator, active, size(w_r[active])),
        )
        Array(compiled(w_r[active], other, cfg_r, upstream_r))
    end
    return result
end

function staged_endpoint_pullback(w_tuple, config, upstream)
    return _staged_endpoint_pullback(w_tuple, config, upstream, _endpoint_value)
end

function staged_nonlimber_endpoint_pullback(w_tuple, config, upstream)
    return _staged_endpoint_pullback(w_tuple, config, upstream, _nonlimber_value)
end

function staged_finalization_pullback(
    C_nonlimber,
    C_correction,
    C_limber,
    finalization,
    upstream,
)
    ell2, transform, T_eval, inv_ell2 = finalization
    # The endpoint has three probe outputs; invoke the finalizer three times
    # so the cotangents remain separate and retain their probe shapes.
    function loss_all(Cn1, Cc1, Cl1, Cn2, Cc2, Cl2, Cn3, Cc3, Cl3,
                      ell2_, M_, Te_, inv2_, U1, U2, U3)
        y1 = Blast.reactant_finalize_c_ell(Cn1, Cc1, Cl1, ell2_, M_, Te_, inv2_)
        y2 = Blast.reactant_finalize_c_ell(Cn2, Cc2, Cl2, ell2_, M_, Te_, inv2_)
        y3 = Blast.reactant_finalize_c_ell(Cn3, Cc3, Cl3, ell2_, M_, Te_, inv2_)
        return sum(y1 .* U1) + sum(y2 .* U2) + sum(y3 .* U3)
    end
    function grad(Cn1, Cc1, Cl1, Cn2, Cc2, Cl2, Cn3, Cc3, Cl3,
                  ell2_, M_, Te_, inv2_, U1, U2, U3)
        return Enzyme.gradient(
            Reverse, loss_all,
            Cn1, Cc1, Cl1, Cn2, Cc2, Cl2, Cn3, Cc3, Cl3,
            Const(ell2_), Const(M_), Const(Te_), Const(inv2_),
            Const(U1), Const(U2), Const(U3),
        )
    end
    args = (
        Reactant.to_rarray(C_nonlimber[1]), Reactant.to_rarray(C_correction[1]), Reactant.to_rarray(C_limber[1]),
        Reactant.to_rarray(C_nonlimber[2]), Reactant.to_rarray(C_correction[2]), Reactant.to_rarray(C_limber[2]),
        Reactant.to_rarray(C_nonlimber[3]), Reactant.to_rarray(C_correction[3]), Reactant.to_rarray(C_limber[3]),
        map(Reactant.to_rarray, finalization)...,
        map(Reactant.to_rarray, upstream)...,
    )
    compiled = _compile_pullback(grad, args; key=(:finalization, size(C_nonlimber[1]), size(C_correction[1]), size(C_limber[1])))
    result = compiled(args...)
    return (
        (Array(result[1]), Array(result[4]), Array(result[7])),
        (Array(result[2]), Array(result[5]), Array(result[8])),
        (Array(result[3]), Array(result[6]), Array(result[9])),
    )
end

"""Chain the non-Limber endpoint cotangent back to the input `pk`."""
function staged_nonlimber_pullback(
    pk,
    background,
    transform,
    coefficients,
    w_tuple,
    T_tuple,
    endpoint_config,
    upstream,
)
    dW = staged_nonlimber_endpoint_pullback(w_tuple, endpoint_config, upstream)
    dcTT, dcT, dcϕ = _pullback_compute_w(
        coefficients[1], coefficients[2], coefficients[3], T_tuple, dW,
    )
    # The W stage uses a rank-2 broadcast representation of primordial
    # coefficients; preparation returns one coefficient per k node.
    dcϕ_vector = vec(sum(dcϕ, dims=2))
    return _pullback_nonlimber_prepare(pk, background, transform, (dcTT, dcT, dcϕ_vector))
end

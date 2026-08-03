module ReactantExt

using LinearAlgebra
using Reactant
import Blast
import AbstractCosmologicalEmulators: ChebyshevPlan, chebyshev_decomposition

function Blast.reactant_p_phi_TT(P_ϕ, T_χ1, T_χR)
    return Blast._p_phi_TT_broadcast(P_ϕ, T_χ1, T_χR)
end

function Blast.reactant_p_phi_T(P_ϕ, T_χR)
    return Blast._p_phi_T_broadcast(P_ϕ, T_χR)
end

# Dispatch the existing production helpers to the validated Reactant forms
# once tracing has produced an RArray. Plain Julia keeps the original Tullio
# methods above in src/setup.jl and src/projected_matter.jl.
function Blast._p_phi_TT_tullio(P_ϕ::AbstractVector,
                                T_χ1::Reactant.TracedRArray,
                                T_χR::AbstractArray{<:Any, 3})
    return Blast._p_phi_TT_broadcast(P_ϕ, T_χ1, T_χR)
end

function Blast._p_phi_T_tullio(
    P_ϕ::AbstractVector,
    T_χR::Union{
        Reactant.TracedRArray,
        Base.ReshapedArray{<:Any, 3, <:Reactant.TracedRArray},
    },
)
    return Blast._p_phi_T_broadcast(P_ϕ, T_χR)
end

"""Build the exact one-dimensional Chebyshev transform matrix for `plan`.

The matrix is host-side setup data. Applying it is a regular matrix
multiplication and is the Reactant-compatible replacement for the FFTW plan in
the hot traced path.
"""
function Blast.reactant_chebyshev_matrix(plan::ChebyshevPlan{1})
    n = first(plan.K) + 1
    basis = Matrix{Float64}(I, n, n)
    return reduce(hcat, (
        chebyshev_decomposition(plan, basis[:, i]) for i in axes(basis, 2)
    ))
end

function Blast.reactant_chebyshev_matmul(vals, transform)
    return reshape(
        transform * reshape(vals, size(vals, 1), :),
        size(vals),
    )
end

function Blast.reactant_chebyshev_2d_matmul(vals, transform_k, transform_z)
    return transform_k * vals * transform_z'
end

_stablehlo_eltype(::Type{Float64}) = "f64"
_stablehlo_eltype(::Type{Float32}) = "f32"
_stablehlo_eltype(::Type{Float16}) = "f16"
_stablehlo_eltype(::Type{<:Reactant.TracedRNumber{T}}) where {T} = _stablehlo_eltype(T)

function _w_hlo_code(::Type{T}, nl, ni, nj, nk, c_rank) where {T}
    ty = _stablehlo_eltype(T)
    einsum, c_type = if c_rank == 3
        ("lijk,ljk->ijk", "tensor<$(nl)x$(nj)x$(nk)x$(ty)>")
    elseif c_rank == 2
        ("lijk,lj->ijk", "tensor<$(nl)x$(nj)x$(ty)>")
    elseif c_rank == 1
        ("lijk,l->ijk", "tensor<$(nl)x$(ty)>")
    else
        throw(ArgumentError("unsupported coefficient rank: $c_rank"))
    end
    return """
module {
  func.func @main(%arg0: $(c_type), %arg1: tensor<$(nl)x$(ni)x$(nj)x$(nk)x$(ty)>) -> tensor<$(ni)x$(nj)x$(nk)x$(ty)> {
    %0 = "stablehlo.einsum"(%arg1, %arg0) {einsum_config = "$(einsum)"} : (tensor<$(nl)x$(ni)x$(nj)x$(nk)x$(ty)>, $(c_type)) -> tensor<$(ni)x$(nj)x$(nk)x$(ty)>
    return %0 : tensor<$(ni)x$(nj)x$(nk)x$(ty)>
  }
}
"""
end

function Blast.reactant_w_ell_hlo(c, T_lijk)
    if Reactant.within_compile()
        nl, ni, nj, nk = size(T_lijk)
        code = _w_hlo_code(eltype(c), nl, ni, nj, nk, ndims(c))
        out_jki = only(Reactant.Ops.hlo_call(code, c, T_lijk))
        return out_jki
    end

    T_ijkl = permutedims(T_lijk, (2, 3, 4, 1))
    return Blast.w_ell_tullio(c, T_ijkl)
end

function Blast.w_ell_tullio(c::Reactant.TracedRArray{<:Any, 3}, T::Array{Float64, 4})
    return Blast.reactant_w_ell_hlo(c, permutedims(T, (4, 1, 2, 3)))
end

function Blast.w_ell_tullio(c::Reactant.TracedRArray{<:Any, 2}, T::Array{Float64, 4})
    return Blast.reactant_w_ell_hlo(c, permutedims(T, (4, 1, 2, 3)))
end

function Blast.w_ell_tullio(c::Reactant.TracedRArray{<:Any, 1}, T::Array{Float64, 4})
    return Blast.reactant_w_ell_hlo(c, permutedims(T, (4, 1, 2, 3)))
end

function Blast.reactant_chebyshev_w_ell_hlo(vals, transform, T_lijk)
    c = Blast.reactant_chebyshev_matmul(vals, transform)
    return Blast.reactant_w_ell_hlo(c, T_lijk)
end

"""Prepare the non-Limber coefficient arrays without host-side wrappers.

`transform` is host setup data represented by a dynamic array argument during
compilation so tests can prove that the compiled program does not capture it.
The background and Blast grids are fixed configuration for this first port.
"""
function Blast.reactant_prepare_nonlimber_spectrum(pk, bg, transform)
    k = 10 .^ Blast.k_cheb
    P_ϕ = Blast.get_PΦ(k, bg)
    transfer_func = Blast.get_Tm(pk, k, bg)
    transfer_func_χR = Blast.transform_to_R_frame(transfer_func, bg)
    transfer_func_χ1 = transfer_func_χR[:, end, :]

    P_ϕTT = Blast._p_phi_TT_broadcast(P_ϕ, transfer_func_χ1, transfer_func_χR)
    P_ϕT = Blast._p_phi_T_broadcast(P_ϕ, transfer_func_χR)

    return (
        Blast.reactant_chebyshev_matmul(P_ϕTT, transform),
        Blast.reactant_chebyshev_matmul(P_ϕT, transform),
        Blast.reactant_chebyshev_matmul(P_ϕ, transform),
    )
end

"""Compute all active projected-matter arrays as a flat tuple.

Every `T_lijk` argument is the corresponding precomputed kernel with the
Chebyshev index first. Returning arrays rather than `ProjectedMatterDensity`
avoids converting traced arrays into the host-only concrete fields of the
ordinary Blast structs.
"""
function Blast.reactant_compute_w(
    cϕTT, cϕT, cϕ_2d,
    T_2_00, T_minus2_00, T_0_00, T_0_02, T_0_20,
    T_2_02, T_2_20, T_2_22,
)
    cϕT_R1 = cϕT[:, :, end]
    cϕ_for_w = if ndims(cϕ_2d) == 1
        reshape(cϕ_2d, size(cϕ_2d, 1), 1) .* ones(Float64, 1, size(T_2_00, 3))
    else
        cϕ_2d
    end
    return (
        Blast.reactant_w_ell_hlo(cϕTT, T_2_00),
        Blast.reactant_w_ell_hlo(cϕTT, T_minus2_00),
        Blast.reactant_w_ell_hlo(cϕTT, T_0_00),
        Blast.reactant_w_ell_hlo(cϕTT, T_0_02),
        Blast.reactant_w_ell_hlo(cϕTT, T_0_20),
        Blast.reactant_w_ell_hlo(cϕTT, T_2_02),
        Blast.reactant_w_ell_hlo(cϕTT, T_2_20),
        Blast.reactant_w_ell_hlo(cϕTT, T_2_22),
        Blast.reactant_w_ell_hlo(cϕT, T_2_00),
        Blast.reactant_w_ell_hlo(cϕT_R1, T_2_00),
        Blast.reactant_w_ell_hlo(cϕT, T_0_00),
        Blast.reactant_w_ell_hlo(cϕT_R1, T_0_00),
        Blast.reactant_w_ell_hlo(cϕT, T_2_02),
        Blast.reactant_w_ell_hlo(cϕT_R1, T_2_02),
        Blast.reactant_w_ell_hlo(cϕT, T_2_20),
        Blast.reactant_w_ell_hlo(cϕT_R1, T_2_20),
        Blast.reactant_w_ell_hlo(cϕ_for_w, T_2_00),
    )
end

function _einsum_code(::Type{T}, shape_a, shape_b, shape_out, equation) where {T}
    ty = _stablehlo_eltype(T)
    dims(shape) = join(("$(d)" for d in shape), "x")
    return """
module {
  func.func @main(%arg0: tensor<$(dims(shape_a))x$(ty)>, %arg1: tensor<$(dims(shape_b))x$(ty)>) -> tensor<$(dims(shape_out))x$(ty)> {
    %0 = "stablehlo.einsum"(%arg0, %arg1) {einsum_config = "$(equation)"} : (tensor<$(dims(shape_a))x$(ty)>, tensor<$(dims(shape_b))x$(ty)>) -> tensor<$(dims(shape_out))x$(ty)>
    return %0 : tensor<$(dims(shape_out))x$(ty)>
  }
}
"""
end

function _reactant_einsum(a, b, shape_a, shape_b, shape_out, equation)
    a = a isa Reactant.TracedRArray ? a : Reactant.to_rarray(a)
    b = b isa Reactant.TracedRArray ? b : Reactant.to_rarray(b)
    code = _einsum_code(eltype(a), shape_a, shape_b, shape_out, equation)
    return only(Reactant.Ops.hlo_call(code, a, b))
end

function Blast.reactant_nonlimber_c_ell(
    W_A, W_B, pmj, prefactor, w_χ, w_R, χ_grid, Δχ,
)
    if !Reactant.within_compile()
        return Blast._compute_Cℓ_fused_tullio(W_A, W_B, pmj, w_χ, w_R, prefactor, Δχ, χ_grid)
    end

    nA, nχ, nR = size(W_A)
    nB = size(W_B, 1)
    nℓ = size(pmj, 1)
    W_A_r1 = W_A[:, :, end]
    W_B_r1 = W_B[:, :, end]
    q = pmj .* reshape(prefactor, nℓ, 1, 1) .*
        reshape(χ_grid .* w_χ .* Δχ, 1, nχ, 1) .*
        reshape(w_R, 1, 1, nR)

    # Split the fused host Tullio expression into two contractions per term.
    # This avoids materializing K[i,j,n,m] while keeping every operation a
    # dense tensor operation representable by StableHLO einsum.
    q_b = _reactant_einsum(q, W_B, (nℓ, nχ, nR), (nB, nχ, nR), (nℓ, nB, nχ), "lnm,jnm->ljn")
    term_a = _reactant_einsum(q_b, W_A_r1, (nℓ, nB, nχ), (nA, nχ), (nℓ, nA, nB), "ljn,in->lij")
    q_a = _reactant_einsum(q, W_A, (nℓ, nχ, nR), (nA, nχ, nR), (nℓ, nA, nχ), "lnm,inm->lin")
    term_b = _reactant_einsum(q_a, W_B_r1, (nℓ, nA, nχ), (nB, nχ), (nℓ, nA, nB), "lin,jn->lij")
    return term_a + term_b
end

function Blast.reactant_nonlimber_rsd_c_ell(
    W_A, W_B, pmj02, pmj20, prefactor, w_χ, w_R, χ_grid, Δχ,
)
    if !Reactant.within_compile()
        return Blast._compute_Cℓ_rsd_fused_tullio(W_A, W_B, pmj02, pmj20, w_χ, w_R, prefactor, Δχ, χ_grid)
    end

    nA, nχ, nR = size(W_A)
    nB = size(W_B, 1)
    nℓ = size(pmj02, 1)
    W_A_r1 = W_A[:, :, end]
    W_B_r1 = W_B[:, :, end]
    scale = reshape(prefactor, nℓ, 1, 1) .*
        reshape(χ_grid .* w_χ .* Δχ, 1, nχ, 1) .*
        reshape(w_R, 1, 1, nR)
    q02 = pmj02 .* scale
    q20 = pmj20 .* scale

    q02_b = _reactant_einsum(q02, W_B, (nℓ, nχ, nR), (nB, nχ, nR), (nℓ, nB, nχ), "lnm,jnm->ljn")
    term02 = _reactant_einsum(q02_b, W_A_r1, (nℓ, nB, nχ), (nA, nχ), (nℓ, nA, nB), "ljn,in->lij")
    q20_a = _reactant_einsum(q20, W_A, (nℓ, nχ, nR), (nA, nχ, nR), (nℓ, nA, nχ), "lnm,inm->lin")
    term20 = _reactant_einsum(q20_a, W_B_r1, (nℓ, nA, nχ), (nB, nχ), (nℓ, nA, nB), "lin,jn->lij")
    return term02 + term20
end

function Blast.reactant_limber_contraction(P_term, K1, K2, weights, Δχ)
    if !Reactant.within_compile()
        return Blast._limber_contraction(P_term, K1, K2, weights, Δχ)
    end

    nℓ, nχ = size(P_term)
    nA = size(K1, 3)
    nB = size(K2, 3)
    weighted_P = P_term .* reshape(weights .* Δχ, 1, nχ)
    first_contraction = _reactant_einsum(
        weighted_P, K1,
        (nℓ, nχ), (nℓ, nχ, nA), (nℓ, nA, nχ),
        "ln,lni->lin",
    )
    return _reactant_einsum(
        first_contraction, K2,
        (nℓ, nA, nχ), (nℓ, nχ, nB), (nℓ, nA, nB),
        "lin,lnj->lij",
    )
end

function Blast.reactant_limber_eval(c, T_z, T_k)
    if !Reactant.within_compile()
        return Blast.limber_eval(c, T_z, T_k)
    end
    B = c * T_z'
    nl, nk, nz = size(T_k)
    return _reactant_einsum(
        T_k, B,
        (nl, nk, nz), (nk, nz), (nl, nz), "lkn,kn->ln",
    )
end

function Blast.reactant_prepare_limber(pk_lin, pk_nonlin, bg, transform_k, transform_z, T_z, T_k)
    k = 10 .^ Blast.k_limber
    P_ϕ = reshape(Blast.get_PΦ(k, bg), 1, :)
    T_lin = Blast.get_Tm(pk_lin, k, bg)
    T_nonlin = Blast.get_Tm(pk_nonlin, k, bg)
    P_lin = permutedims(P_ϕ .* T_lin .^ 2)
    P_nonlin = permutedims(P_ϕ .* T_nonlin .^ 2)
    c_lin = Blast.reactant_chebyshev_2d_matmul(log10.(P_lin), transform_k, transform_z)
    c_nonlin = Blast.reactant_chebyshev_2d_matmul(log10.(P_nonlin), transform_k, transform_z)
    P_lin_grid = 10.0 .^ Blast.reactant_limber_eval(c_lin, T_z, T_k)
    P_nonlin_grid = 10.0 .^ Blast.reactant_limber_eval(c_nonlin, T_z, T_k)
    return P_nonlin_grid .- P_lin_grid, P_nonlin_grid
end

function Blast.reactant_limber_power_products(pk_lin, pk_nonlin, bg)
    k = 10 .^ Blast.k_limber
    P_ϕ = reshape(Blast.get_PΦ(k, bg), 1, :)
    T_lin = Blast.get_Tm(pk_lin, k, bg)
    T_nonlin = Blast.get_Tm(pk_nonlin, k, bg)
    return log10.(permutedims(P_ϕ .* T_lin .^ 2)), log10.(permutedims(P_ϕ .* T_nonlin .^ 2))
end

function Blast.reactant_limber_chebyshev_coefficients(logP, transform_k, transform_z)
    return Blast.reactant_chebyshev_2d_matmul(logP, transform_k, transform_z)
end

function Blast.reactant_limber_grid_from_coefficients(coefficients, T_z, T_k)
    return 10.0 .^ Blast.reactant_limber_eval(coefficients, T_z, T_k)
end

function Blast.reactant_limber_grid_difference(P_lin, P_nonlin)
    return P_nonlin .- P_lin, P_nonlin
end

function Blast.reactant_compute_w_from_spectrum(cϕTT, cϕT, cϕ, T_tuple...)
    nχ = size(cϕTT, 2)
    cϕTT_ = cϕTT .* ones(eltype(cϕTT), size(cϕTT))
    cϕT_ = cϕT .* ones(eltype(cϕT), size(cϕT))
    cϕ_2d = reshape(cϕ, :, 1) .* ones(eltype(cϕ), 1, nχ)
    return Blast.reactant_compute_w(cϕTT_, cϕT_, cϕ_2d, T_tuple...)
end

function Blast.reactant_limber_terms_from_prepared(
    ΔP, P_nonlin, KGG_low, KGG_high, KG_low, KL_low,
    KG_high, KL_high, KLL_low, KLL_high, inv_χ2, weights, Δχ,
)
    nlow = size(KGG_low, 1)
    low = ΔP[1:nlow, :] .* reshape(inv_χ2, 1, :)
    high = P_nonlin[(nlow + 1):end, :] .* reshape(inv_χ2, 1, :)
    return (
        Blast.reactant_limber_contraction(low, KGG_low, KGG_low, weights, Δχ),
        Blast.reactant_limber_contraction(high, KGG_high, KGG_high, weights, Δχ),
        Blast.reactant_limber_contraction(low, KG_low, KL_low, weights, Δχ),
        Blast.reactant_limber_contraction(high, KG_high, KL_high, weights, Δχ),
        Blast.reactant_limber_contraction(low, KLL_low, KLL_low, weights, Δχ),
        Blast.reactant_limber_contraction(high, KLL_high, KLL_high, weights, Δχ),
    )
end

function Blast.reactant_finalize_c_ell(
    C_nonlimber, C_correction, C_limber,
    ell2_reversed, transform, T_eval, inv_ell2,
)
    full = vcat(C_nonlimber .+ C_correction, C_limber)
    weighted = reverse(full; dims=1) .* reshape(ell2_reversed, :, 1, 1)
    nfull = size(weighted, 1)
    coefficients = Blast.reactant_chebyshev_matmul(
        weighted,
        transform,
    )
    evaluated = T_eval * reshape(coefficients, nfull, :)
    return reshape(evaluated, length(inv_ell2), size(full, 2), size(full, 3)) .* reshape(inv_ell2, :, 1, 1)
end

function Blast.reactant_nonlimber_3x2pt(
    w_2_00_ϕTT, w_minus2_00_ϕTT, w_0_00_ϕTT,
    w_0_02_ϕTT, w_0_20_ϕTT, w_2_02_ϕTT, w_2_20_ϕTT, w_2_22_ϕTT,
    w_2_00_ϕT, w_2_00_ϕT_R1, w_0_00_ϕT, w_0_00_ϕT_R1,
    w_2_02_ϕT, w_2_02_ϕT_R1, w_2_20_ϕT, w_2_20_ϕT_R1, w_2_00_ϕ,
    Kδ, KRSD, Kμ, KfNL, Kγ, KIA,
    wχ, wR, χ_grid, Δχ, pref,
)
    reactant_nonlimber_c_ell = Blast.reactant_nonlimber_c_ell
    reactant_nonlimber_rsd_c_ell = Blast.reactant_nonlimber_rsd_c_ell

    # Galaxy auto.
    gg = reactant_nonlimber_c_ell(Kδ, Kδ, w_2_00_ϕTT, pref.δδ, wχ, wR, χ_grid, Δχ)
    gg -= reactant_nonlimber_rsd_c_ell(Kδ, KRSD, w_2_02_ϕTT, w_2_20_ϕTT, pref.δRSD, wχ, wR, χ_grid, Δχ)
    gg -= reactant_nonlimber_rsd_c_ell(KRSD, Kδ, w_2_20_ϕTT, w_2_02_ϕTT, pref.δRSD, wχ, wR, χ_grid, Δχ)
    gg += reactant_nonlimber_c_ell(KRSD, KRSD, w_2_22_ϕTT, pref.RSDRSD, wχ, wR, χ_grid, Δχ)
    gg += reactant_nonlimber_c_ell(Kδ, Kμ, w_0_00_ϕTT, pref.δμ, wχ, wR, χ_grid, Δχ)
    gg += reactant_nonlimber_c_ell(Kμ, Kδ, w_0_00_ϕTT, pref.δμ, wχ, wR, χ_grid, Δχ)
    gg += reactant_nonlimber_c_ell(Kμ, Kμ, w_minus2_00_ϕTT, pref.μμ, wχ, wR, χ_grid, Δχ)
    gg -= reactant_nonlimber_rsd_c_ell(Kμ, KRSD, w_0_02_ϕTT, w_0_20_ϕTT, pref.μRSD, wχ, wR, χ_grid, Δχ)
    gg -= reactant_nonlimber_rsd_c_ell(KRSD, Kμ, w_0_20_ϕTT, w_0_02_ϕTT, pref.μRSD, wχ, wR, χ_grid, Δχ)
    gg += reactant_nonlimber_c_ell(Kδ, KfNL, w_2_00_ϕT_R1, pref.δfNL, wχ, wR, χ_grid, Δχ)
    gg += reactant_nonlimber_c_ell(KfNL, Kδ, w_2_00_ϕT, pref.fNLδ, wχ, wR, χ_grid, Δχ)
    gg -= reactant_nonlimber_rsd_c_ell(KfNL, KRSD, w_2_02_ϕT, w_2_20_ϕT, pref.fNLRSD, wχ, wR, χ_grid, Δχ)
    gg -= reactant_nonlimber_rsd_c_ell(KRSD, KfNL, w_2_20_ϕT_R1, w_2_02_ϕT_R1, pref.RSDfNL, wχ, wR, χ_grid, Δχ)
    gg += reactant_nonlimber_c_ell(Kμ, KfNL, w_0_00_ϕT_R1, pref.μfNL, wχ, wR, χ_grid, Δχ)
    gg += reactant_nonlimber_c_ell(KfNL, Kμ, w_0_00_ϕT, pref.fNLμ, wχ, wR, χ_grid, Δχ)
    gg += reactant_nonlimber_c_ell(KfNL, KfNL, w_2_00_ϕ, pref.fNLfNL, wχ, wR, χ_grid, Δχ)

    # Weak-lensing auto.
    ss = reactant_nonlimber_c_ell(Kγ, Kγ, w_minus2_00_ϕTT, pref.γγ, wχ, wR, χ_grid, Δχ)
    ss += reactant_nonlimber_c_ell(Kγ, KIA, w_minus2_00_ϕTT, pref.γIA, wχ, wR, χ_grid, Δχ)
    ss += reactant_nonlimber_c_ell(KIA, Kγ, w_minus2_00_ϕTT, pref.IAγ, wχ, wR, χ_grid, Δχ)
    ss += reactant_nonlimber_c_ell(KIA, KIA, w_minus2_00_ϕTT, pref.IAIA, wχ, wR, χ_grid, Δχ)

    # Galaxy-lensing cross.
    gs = reactant_nonlimber_c_ell(Kδ, Kγ, w_0_00_ϕTT, pref.δγ, wχ, wR, χ_grid, Δχ)
    gs += reactant_nonlimber_c_ell(Kδ, KIA, w_0_00_ϕTT, pref.δIA, wχ, wR, χ_grid, Δχ)
    gs -= reactant_nonlimber_rsd_c_ell(KRSD, Kγ, w_0_20_ϕTT, w_0_02_ϕTT, pref.RSDγ, wχ, wR, χ_grid, Δχ)
    gs -= reactant_nonlimber_rsd_c_ell(KRSD, KIA, w_0_20_ϕTT, w_0_02_ϕTT, pref.RSDIA, wχ, wR, χ_grid, Δχ)
    gs += reactant_nonlimber_c_ell(Kμ, Kγ, w_minus2_00_ϕTT, pref.μγ, wχ, wR, χ_grid, Δχ)
    gs += reactant_nonlimber_c_ell(Kμ, KIA, w_minus2_00_ϕTT, pref.μIA, wχ, wR, χ_grid, Δχ)
    gs += reactant_nonlimber_c_ell(KfNL, Kγ, w_0_00_ϕT, pref.fNLγ, wχ, wR, χ_grid, Δχ)
    gs += reactant_nonlimber_c_ell(KfNL, KIA, w_0_00_ϕT, pref.fNLIA, wχ, wR, χ_grid, Δχ)
    return (; gg, gs, ss)
end

function Blast.reactant_full_3x2pt(
    w_2_00_ϕTT, w_minus2_00_ϕTT, w_0_00_ϕTT,
    w_0_02_ϕTT, w_0_20_ϕTT, w_2_02_ϕTT, w_2_20_ϕTT, w_2_22_ϕTT,
    w_2_00_ϕT, w_2_00_ϕT_R1, w_0_00_ϕT, w_0_00_ϕT_R1,
    w_2_02_ϕT, w_2_02_ϕT_R1, w_2_20_ϕT, w_2_20_ϕT_R1, w_2_00_ϕ,
    Kδ, KRSD, Kμ, KfNL, Kγ, KIA, wχ, wR, χ_grid,
    Cgg_correction, Cgg_limber, Cgs_correction, Cgs_limber,
    Css_correction, Css_limber,
    ell2_reversed, transform, T_eval, inv_ell2, Δχ, pref,
)
    nonlimber = Blast.reactant_nonlimber_3x2pt(
        w_2_00_ϕTT, w_minus2_00_ϕTT, w_0_00_ϕTT,
        w_0_02_ϕTT, w_0_20_ϕTT, w_2_02_ϕTT, w_2_20_ϕTT, w_2_22_ϕTT,
        w_2_00_ϕT, w_2_00_ϕT_R1, w_0_00_ϕT, w_0_00_ϕT_R1,
        w_2_02_ϕT, w_2_02_ϕT_R1, w_2_20_ϕT, w_2_20_ϕT_R1, w_2_00_ϕ,
        Kδ, KRSD, Kμ, KfNL, Kγ, KIA, wχ, wR, χ_grid,
        Δχ, pref,
    )
    return (
        Blast.reactant_finalize_c_ell(nonlimber.gg, Cgg_correction, Cgg_limber, ell2_reversed, transform, T_eval, inv_ell2),
        Blast.reactant_finalize_c_ell(nonlimber.gs, Cgs_correction, Cgs_limber, ell2_reversed, transform, T_eval, inv_ell2),
        Blast.reactant_finalize_c_ell(nonlimber.ss, Css_correction, Css_limber, ell2_reversed, transform, T_eval, inv_ell2),
    )
end

end

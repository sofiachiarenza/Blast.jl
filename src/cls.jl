"""
    make_grid(χ::Array{<:Any, 1}, R::Array{<:Any, 1}) 

Construct a flattened χ–R grid corresponding to all combinations of comoving distance `χ` and ratio `R`.
"""
function make_grid(χ::Array{<:Any, 1}, R::Array{<:Any, 1}) 
    return vec(χ * R') 
end

# χ–R grid helpers are global invariants (Blast.χ and Blast.R are constants).
# Cache once to avoid rebuilding in every get_Cℓ path.
const χR_GRID_FLAT = make_grid(Blast.χ, Blast.R)
const INV_χ2_GRID_FLAT = 1.0 ./ (χR_GRID_FLAT .^ 2)
const INV_χ2_GRID_ROW = reshape(INV_χ2_GRID_FLAT, 1, :)
const FULL_ℓ2 = Blast.full_ℓ_range .^ 2
const FULL_ℓ2_REVERSED = reverse(FULL_ℓ2)

# Interpolated probe kernels can pick up numerically useless tiny tails from the
# tomographic-overlap geometry. Those values do not make the contraction sparse,
# but they do drive products into denormal territory and slow the hot non-Limber
# contractions down badly. Clip them once, right after interpolation.
const KERNEL_ARRAY_ABS_CUTOFF = 1e-100

function _apply_kernel_array_cutoff!(A::AbstractArray)
    z = zero(eltype(A))
    @inbounds for idx in eachindex(A)
        if abs(A[idx]) < KERNEL_ARRAY_ABS_CUTOFF
            A[idx] = z
        end
    end
    return A
end

"""
    grid_interpolator(Probe, bg)

Interpolate a probe kernel from the χ grid onto a χ–R grid using Background interpolators.
"""
function grid_interpolator(Probe::AbstractComponents, bg::Background) 
    kernel = Probe.Kernel
    chi_grid = bg.χ
    return collect(akima_interpolation(collect(kernel'), chi_grid, χR_GRID_FLAT)')
end

"""
    get_kernel_array(Probe, bg)

Return the probe kernel evaluated on the χ–R grid.
"""
function get_kernel_array(Probe::Union{NumberCounts, RedshiftSpaceDistortions, PrimordialNonGaussianity, IntegratedSachsWolfe}, bg::Background)
    n_bins = size(Probe.Kernel, 1)
    nχ = length(bg.χ)
    nR = length(Blast.R)
    W_array = reshape(grid_interpolator(Probe, bg), n_bins, nχ, nR)
    _apply_kernel_array_cutoff!(W_array)
    return W_array
end

function get_kernel_array(Probe::Union{CosmicShear, IntrinsicAlignment, MagnificationBias, CMBLensing}, bg::Background)
    n_bins = size(Probe.Kernel, 1)
    nχ = length(bg.χ)
    nR = length(Blast.R)
    W_L = grid_interpolator(Probe, bg)
    W_L .*= INV_χ2_GRID_ROW
    _apply_kernel_array_cutoff!(W_L)
    return reshape(W_L, n_bins, nχ, nR)
end

function _combine_kernels_tullio(W_A, W_B)
    W_A_r1 = W_A[:, :, end]
    W_B_r1 = W_B[:, :, end]
    @tullio K[idx_i, idx_j, idx_c, idx_r] := W_A_r1[idx_i, idx_c] * W_B[idx_j, idx_c, idx_r] + W_A[idx_i, idx_c, idx_r] * W_B_r1[idx_j, idx_c]
    return K
end

"""
    combine_kernels(ProbeA, ProbeB, bg)
"""
function combine_kernels(ProbeA::AbstractComponents, ProbeB::AbstractComponents, bg::Background)
    W_A = get_kernel_array(ProbeA, bg)
    W_B = get_kernel_array(ProbeB, bg)
    return _combine_kernels_tullio(W_A, W_B)
end

function _compute_Cℓ_tullio(K, pmj, w_χ, w_R, prefactor, Δχ, χ_grid)
    @tullio Cℓ[idx_l, idx_i, idx_j] := prefactor[idx_l] * χ_grid[idx_n] * K[idx_i, idx_j, idx_n, idx_m] * pmj[idx_l, idx_n, idx_m] * w_χ[idx_n] * w_R[idx_m] * Δχ
    return Cℓ
end

function _compute_Cℓ_rsd_tullio(W_A_r1, W_B, pmj02, W_A, W_B_r1, pmj20, w_χ, w_R, prefactor, Δχ, χ_grid)
    @tullio K[idx_l, idx_i, idx_j, idx_c, idx_r] := W_A_r1[idx_i, idx_c] * W_B[idx_j, idx_c, idx_r] * pmj02[idx_l, idx_c, idx_r] + W_A[idx_i, idx_c, idx_r] * W_B_r1[idx_j, idx_c] * pmj20[idx_l, idx_c, idx_r]
    @tullio Cℓ[idx_l, idx_i, idx_j] := prefactor[idx_l] * χ_grid[idx_n] * K[idx_l, idx_i, idx_j, idx_n, idx_m] * w_χ[idx_n] * w_R[idx_m] * Δχ
    return Cℓ
end

# -----------------------------------------------------------------------------
# Fused variants — used on the hot path.
#
# The pair (_combine_kernels_tullio, _compute_Cℓ_tullio) was responsible for a
# 4-D (nbins_A, nbins_B, nχ, nR) intermediate `K[i,j,c,r]` that was allocated
# then immediately consumed.  Substituting the K-expression into the Cℓ sum
# lets Tullio contract straight to the 3-D result, eliminating the biggest
# non-AD allocation in the get_Cℓ pipeline (≈67% of GG non-AD allocation).
#
# Mathematically identical to the two-step form (see tests and rrules).
# -----------------------------------------------------------------------------

function _compute_Cℓ_fused_tullio(W_A, W_B, pmj, w_χ, w_R, prefactor, Δχ, χ_grid)
    W_A_r1 = W_A[:, :, end]
    W_B_r1 = W_B[:, :, end]
    @tullio Cℓ[idx_l, idx_i, idx_j] := prefactor[idx_l] * χ_grid[idx_n] * w_χ[idx_n] * w_R[idx_m] * Δχ * pmj[idx_l, idx_n, idx_m] *
                                       (W_A_r1[idx_i, idx_n] * W_B[idx_j, idx_n, idx_m] +
                                        W_A[idx_i, idx_n, idx_m] * W_B_r1[idx_j, idx_n])
    return Cℓ
end

function _compute_Cℓ_rsd_fused_tullio(W_A, W_B, pmj02, pmj20, w_χ, w_R, prefactor, Δχ, χ_grid)
    W_A_r1 = W_A[:, :, end]
    W_B_r1 = W_B[:, :, end]
    @tullio Cℓ[idx_l, idx_i, idx_j] := prefactor[idx_l] * χ_grid[idx_n] * w_χ[idx_n] * w_R[idx_m] * Δχ *
                                       (W_A_r1[idx_i, idx_n] * W_B[idx_j, idx_n, idx_m] * pmj02[idx_l, idx_n, idx_m] +
                                        W_A[idx_i, idx_n, idx_m] * W_B_r1[idx_j, idx_n] * pmj20[idx_l, idx_n, idx_m])
    return Cℓ
end

function _prepare_nonlimber_integration(bg::Background)
    _require_production_background(bg)
    nχ = size(bg.χ, 1)
    nR = size(Blast.R, 1)
    Δχ = ((last(bg.χ)-first(bg.χ))/(nχ-1))
    w_χ = simpson_weights_array(nχ)
    w_R = get_clencurt_weights_R_integration(2*nR+1)
    χ_grid = bg.χ
    return (; Δχ, w_χ, w_R, χ_grid)
end

function _require_production_background(bg::Background)
    length(bg.χ) == length(Blast.χ) || throw(DimensionMismatch(
        "the full Blast pipeline requires the production χ grid of length $(length(Blast.χ)); " *
        "got $(length(bg.χ))"))
    all(eachindex(Blast.χ)) do i
        bg.χ[i] == Blast.χ[i]
    end || throw(ArgumentError(
        "the full Blast pipeline currently requires the exact production Blast.χ grid; " *
        "custom Background grids are supported only for isolated background/kernel calculations"))
    return nothing
end

function _validate_multipoles(ℓ)
    isempty(ℓ) && throw(ArgumentError("requested multipoles must not be empty"))
    all(isfinite, ℓ) || throw(ArgumentError("requested multipoles must be finite"))
    all(x -> 2 <= x <= 2000, ℓ) || throw(DomainError(
        extrema(ℓ), "requested multipoles must lie within Blast's interpolation domain [2, 2000]"))
    return nothing
end

function _pair_prefactor(Component1::AbstractComponents, Component2::AbstractComponents)
    return 2/π .* Component1.ell_prefactor[1:length(ℓ_nonlimber)] .* Component2.ell_prefactor[1:length(ℓ_nonlimber)]
end

function _compute_Cℓ_from_kernels(W_A, W_B, pmj, prefactor, integ)
    return _compute_Cℓ_fused_tullio(W_A, W_B, pmj, integ.w_χ, integ.w_R, prefactor, integ.Δχ, integ.χ_grid)
end

function _compute_Cℓ_rsd_from_kernels(W_A, W_B, pmj02, pmj20, prefactor, integ)
    return _compute_Cℓ_rsd_fused_tullio(W_A, W_B, pmj02, pmj20, integ.w_χ, integ.w_R, prefactor, integ.Δχ, integ.χ_grid)
end

_get_kernel_or_nothing(Component::Nothing, bg::Background) = nothing
_get_kernel_or_nothing(Component::AbstractComponents, bg::Background) = get_kernel_array(Component, bg)

function _compute_Cℓ_cached(Component1::Union{AbstractComponents, Nothing}, Component2::Union{AbstractComponents, Nothing}, W_A, W_B, w::Union{ProjectedMatterDensityComponent, Nothing}, integ)
    if isnothing(Component1) || isnothing(Component2) || isnothing(w)
        return 0.
    end
    return _compute_Cℓ_from_kernels(W_A, W_B, w.w, _pair_prefactor(Component1, Component2), integ)
end

function _compute_Cℓ_rsd_cached(Component1::Union{AbstractComponents, Nothing}, Component2::Union{AbstractComponents, Nothing}, W_A, W_B, w02::Union{ProjectedMatterDensityComponent, Nothing}, w20::Union{ProjectedMatterDensityComponent, Nothing}, integ)
    if isnothing(Component1) || isnothing(Component2) || isnothing(w02) || isnothing(w20)
        return 0.
    end
    return _compute_Cℓ_rsd_from_kernels(W_A, W_B, w02.w, w20.w, _pair_prefactor(Component1, Component2), integ)
end

# Shared finalize step: given `full_Cℓ` stacked on `full_ℓ_range`, interpolate
# onto the requested ℓ grid via a Chebyshev decomposition of ℓ² · C_ℓ (the
# ℓ² weighting smooths the spectrum before decomposition; we divide it back
# out after chebinterp_native). Called by every `get_Cℓ(ℓ, ...)` mixing the
# non-Limber and Limber paths.
#
# `eltype(full_Cℓ)` is used for the zeros(...) initializer so ForwardDiff.Dual
# threads through the loop assignment; plain `zeros(...)` would silently strip
# Duals at the `Cℓ_final[:, i, j] = ...` store.
function _finalize_Cℓ(full_Cℓ, ℓ, nbins_A::Integer, nbins_B::Integer, P::Union{FFTPlans, Nothing})
    _validate_multipoles(ℓ)
    nℓ = size(ℓ, 1)
    nℓ_full = size(full_Cℓ, 1)
    Cℓ_final = zeros(eltype(full_Cℓ), nℓ, nbins_A, nbins_B)
    plan_interp = isnothing(P) ? prepare_chebyshev_plan(2.0, 2000.0, 100) : P.plan_ℓ
    inv_ℓ2 = inv.(ℓ .^ 2)
    ℓ_eval = float.(ℓ)
    T_eval = chebyshev_polynomials(ℓ_eval, eltype(ℓ_eval)(2.0), eltype(ℓ_eval)(2000.0), nℓ_full - 1)

    tmp_weighted = Vector{eltype(full_Cℓ)}(undef, nℓ_full)
    tmp_interp = Vector{promote_type(eltype(T_eval), eltype(full_Cℓ))}(undef, nℓ)

    for i in 1:nbins_A, j in 1:nbins_B
        @inbounds for k in 1:nℓ_full
            tmp_weighted[k] = full_Cℓ[nℓ_full - k + 1, i, j] * FULL_ℓ2_REVERSED[k]
        end
        c_coeffs = chebyshev_decomposition(plan_interp, tmp_weighted)
        mul!(tmp_interp, T_eval, c_coeffs)
        @views Cℓ_final[:, i, j] .= tmp_interp .* inv_ℓ2
    end
    return Cℓ_final
end

function _finalize_Cℓ_parts(Cℓ_nonlimber::AbstractArray{<:Any, 3}, Cℓ_correction::AbstractArray{<:Any, 3}, Cℓ_limber::AbstractArray{<:Any, 3}, ℓ, nbins_A::Integer, nbins_B::Integer, P::Union{FFTPlans, Nothing})
    _validate_multipoles(ℓ)
    nℓ = size(ℓ, 1)
    n_nonlimber = size(Cℓ_nonlimber, 1)
    n_limber = size(Cℓ_limber, 1)
    nℓ_full = n_nonlimber + n_limber

    Cℓ_final = zeros(promote_type(eltype(Cℓ_nonlimber), eltype(Cℓ_correction), eltype(Cℓ_limber)), nℓ, nbins_A, nbins_B)
    plan_interp = isnothing(P) ? prepare_chebyshev_plan(2.0, 2000.0, 100) : P.plan_ℓ
    inv_ℓ2 = inv.(ℓ .^ 2)
    ℓ_eval = float.(ℓ)
    T_eval = chebyshev_polynomials(ℓ_eval, eltype(ℓ_eval)(2.0), eltype(ℓ_eval)(2000.0), nℓ_full - 1)

    tmp_weighted = Vector{eltype(Cℓ_final)}(undef, nℓ_full)
    tmp_interp = Vector{eltype(Cℓ_final)}(undef, nℓ)
    ℓ2_rev = @view FULL_ℓ2_REVERSED[1:nℓ_full]

    for i in 1:nbins_A, j in 1:nbins_B
        @inbounds for k in 1:nℓ_full
            idx = nℓ_full - k + 1
            if idx <= n_nonlimber
                tmp_weighted[k] = (Cℓ_nonlimber[idx, i, j] + Cℓ_correction[idx, i, j]) * ℓ2_rev[k]
            else
                tmp_weighted[k] = Cℓ_limber[idx - n_nonlimber, i, j] * ℓ2_rev[k]
            end
        end
        c_coeffs = chebyshev_decomposition(plan_interp, tmp_weighted)
        mul!(tmp_interp, T_eval, c_coeffs)
        @views Cℓ_final[:, i, j] .= tmp_interp .* inv_ℓ2
    end
    return Cℓ_final
end

function _finalize_Cℓ_parts(Cℓ_nonlimber::AbstractArray{<:Any, 3}, Cℓ_correction::Number, Cℓ_limber::AbstractArray{<:Any, 3}, ℓ, nbins_A::Integer, nbins_B::Integer, P::Union{FFTPlans, Nothing})
    _validate_multipoles(ℓ)
    nℓ = size(ℓ, 1)
    n_nonlimber = size(Cℓ_nonlimber, 1)
    n_limber = size(Cℓ_limber, 1)
    nℓ_full = n_nonlimber + n_limber

    Cℓ_final = zeros(promote_type(eltype(Cℓ_nonlimber), typeof(Cℓ_correction), eltype(Cℓ_limber)), nℓ, nbins_A, nbins_B)
    plan_interp = isnothing(P) ? prepare_chebyshev_plan(2.0, 2000.0, 100) : P.plan_ℓ
    inv_ℓ2 = inv.(ℓ .^ 2)
    ℓ_eval = float.(ℓ)
    T_eval = chebyshev_polynomials(ℓ_eval, eltype(ℓ_eval)(2.0), eltype(ℓ_eval)(2000.0), nℓ_full - 1)

    tmp_weighted = Vector{eltype(Cℓ_final)}(undef, nℓ_full)
    tmp_interp = Vector{eltype(Cℓ_final)}(undef, nℓ)
    ℓ2_rev = @view FULL_ℓ2_REVERSED[1:nℓ_full]

    for i in 1:nbins_A, j in 1:nbins_B
        @inbounds for k in 1:nℓ_full
            idx = nℓ_full - k + 1
            if idx <= n_nonlimber
                tmp_weighted[k] = (Cℓ_nonlimber[idx, i, j] + Cℓ_correction) * ℓ2_rev[k]
            else
                tmp_weighted[k] = Cℓ_limber[idx - n_nonlimber, i, j] * ℓ2_rev[k]
            end
        end
        c_coeffs = chebyshev_decomposition(plan_interp, tmp_weighted)
        mul!(tmp_interp, T_eval, c_coeffs)
        @views Cℓ_final[:, i, j] .= tmp_interp .* inv_ℓ2
    end
    return Cℓ_final
end

"""
    compute_Cℓ(Component1, Component2, w, bg)
"""
function compute_Cℓ(Component1::AbstractComponents, Component2::AbstractComponents, w::ProjectedMatterDensityComponent, bg::Background) 
    integ = _prepare_nonlimber_integration(bg)
    W_A = get_kernel_array(Component1, bg)
    W_B = get_kernel_array(Component2, bg)
    prefactor = _pair_prefactor(Component1, Component2)
    return _compute_Cℓ_from_kernels(W_A, W_B, w.w, prefactor, integ)
end

"""
    compute_Cℓ(Component1, Component2, w02, w20, bg)
"""
function compute_Cℓ(Component1::AbstractComponents, Component2::AbstractComponents, 
    w02::ProjectedMatterDensityComponent, w20::ProjectedMatterDensityComponent, bg::Background) 
    integ = _prepare_nonlimber_integration(bg)
    W_A = get_kernel_array(Component1, bg)
    W_B = get_kernel_array(Component2, bg)
    prefactor = _pair_prefactor(Component1, Component2)
    return _compute_Cℓ_rsd_from_kernels(W_A, W_B, w02.w, w20.w, prefactor, integ)
end

"""
    get_Cℓ(...)

Compute angular spectra at finite multipoles in Blast's supported interpolation
domain `2 ≤ ℓ ≤ 2000`. The full pipeline requires the exact production
`Blast.χ` background grid.
"""
function get_Cℓ(Component1::Union{AbstractComponents, Nothing}, Component2::Union{AbstractComponents, Nothing}, W::ProjectedMatterDensity, bg::Background)
    return 0.
end

## Galaxy clustering auto
function get_Cℓ(Component1::NumberCounts, Component2::NumberCounts, W::ProjectedMatterDensity, bg::Background)
    return compute_Cℓ(Component1, Component2, W.w_2_00_ϕTT, bg)
end

function get_Cℓ(Component1::NumberCounts, Component2::RedshiftSpaceDistortions, W::ProjectedMatterDensity, bg::Background)
    return compute_Cℓ(Component1, Component2, W.w_2_02_ϕTT, W.w_2_20_ϕTT, bg)
end

function get_Cℓ(Component1::RedshiftSpaceDistortions, Component2::NumberCounts, W::ProjectedMatterDensity, bg::Background)
    return compute_Cℓ(Component1, Component2, W.w_2_20_ϕTT, W.w_2_02_ϕTT, bg)
end

function get_Cℓ(Component1::RedshiftSpaceDistortions, Component2::RedshiftSpaceDistortions, W::ProjectedMatterDensity, bg::Background)
    return compute_Cℓ(Component1, Component2, W.w_2_22_ϕTT, bg)
end

function get_Cℓ(Component1::NumberCounts, Component2::MagnificationBias, W::ProjectedMatterDensity, bg::Background)
    return compute_Cℓ(Component1, Component2, W.w_0_00_ϕTT, bg)
end

function get_Cℓ(Component1::MagnificationBias, Component2::NumberCounts, W::ProjectedMatterDensity, bg::Background)
    return compute_Cℓ(Component1, Component2, W.w_0_00_ϕTT, bg)
end

function get_Cℓ(Component1::MagnificationBias, Component2::MagnificationBias, W::ProjectedMatterDensity, bg::Background)
    return compute_Cℓ(Component1, Component2, W.w_minus2_00_ϕTT, bg)
end

function get_Cℓ(Component1::MagnificationBias, Component2::RedshiftSpaceDistortions, W::ProjectedMatterDensity, bg::Background)
    return compute_Cℓ(Component1, Component2, W.w_0_02_ϕTT, W.w_0_20_ϕTT, bg)
end

function get_Cℓ(Component1::RedshiftSpaceDistortions, Component2::MagnificationBias, W::ProjectedMatterDensity, bg::Background)
    return compute_Cℓ(Component1, Component2, W.w_0_20_ϕTT, W.w_0_02_ϕTT, bg)
end

function get_Cℓ(Component1::NumberCounts, Component2::PrimordialNonGaussianity, W::ProjectedMatterDensity, bg::Background)
    return compute_Cℓ(Component1, Component2, W.w_2_00_ϕT_R1, bg)
end

function get_Cℓ(Component1::PrimordialNonGaussianity, Component2::NumberCounts, W::ProjectedMatterDensity, bg::Background)
    return compute_Cℓ(Component1, Component2, W.w_2_00_ϕT, bg)
end

function get_Cℓ(Component1::PrimordialNonGaussianity, Component2::RedshiftSpaceDistortions, W::ProjectedMatterDensity, bg::Background)
    return compute_Cℓ(Component1, Component2, W.w_2_02_ϕT, W.w_2_20_ϕT, bg)
end

function get_Cℓ(Component1::RedshiftSpaceDistortions, Component2::PrimordialNonGaussianity, W::ProjectedMatterDensity, bg::Background)
    return compute_Cℓ(Component1, Component2, W.w_2_20_ϕT_R1, W.w_2_02_ϕT_R1, bg)
end

function get_Cℓ(Component1::MagnificationBias, Component2::PrimordialNonGaussianity, W::ProjectedMatterDensity, bg::Background)
    return compute_Cℓ(Component1, Component2, W.w_0_00_ϕT_R1, bg)
end

function get_Cℓ(Component1::PrimordialNonGaussianity, Component2::MagnificationBias, W::ProjectedMatterDensity, bg::Background)
    return compute_Cℓ(Component1, Component2, W.w_0_00_ϕT, bg)
end

function get_Cℓ(Component1::PrimordialNonGaussianity, Component2::PrimordialNonGaussianity, W::ProjectedMatterDensity, bg::Background)
    return compute_Cℓ(Component1, Component2, W.w_2_00_ϕ, bg)
end

# ----------------------------------------------------------------------------
# Pair-spectrum transpose helper.
#
# When two clustering components A,B share the same projected-matter weight
# tensor, the fused-Tullio expression for `Cℓ_AB[l,i,j]` is invariant under
# the simultaneous swaps  i↔j  and  A↔B,  so
#
#       Cℓ_BA[l,i,j]  ≡  Cℓ_AB[l,j,i]   (verified empirically to FP roundoff,
#                                        rel ≈ 1e-16; see
#                                        benchmark/perf_GG_perpair_symmetry.jl)
#
# We exploit this in `get_Cℓ(::GalaxyClustering, ...)` to compute one pair
# of each (δ↔μ, δ↔RSD, μ↔RSD) and recover its transpose for free.
#
# `_compute_Cℓ_cached` / `_compute_Cℓ_rsd_cached` short-circuit to the scalar
# `0.` when any component or weight is nothing; the helper passes scalars
# through unchanged so the auto-only / partial-probe paths still work.
# ----------------------------------------------------------------------------
_transpose_cℓ_pair(c::AbstractArray{<:Any, 3}) = permutedims(c, (1, 3, 2))
_transpose_cℓ_pair(c::Number) = c
_cℓ_term_eltype(c::AbstractArray) = eltype(c)
_cℓ_term_eltype(c::Number) = typeof(c)

function _combine_gc_nonlimber(
    Cℓ_δδ, Cℓ_δRSD, Cℓ_RSDδ, Cℓ_RSDRSD,
    Cℓ_δμ, Cℓ_μδ, Cℓ_μμ, Cℓ_μRSD, Cℓ_RSDμ,
)
    T = promote_type(
        _cℓ_term_eltype(Cℓ_δδ), _cℓ_term_eltype(Cℓ_δRSD),
        _cℓ_term_eltype(Cℓ_RSDδ), _cℓ_term_eltype(Cℓ_RSDRSD),
        _cℓ_term_eltype(Cℓ_δμ), _cℓ_term_eltype(Cℓ_μδ),
        _cℓ_term_eltype(Cℓ_μμ), _cℓ_term_eltype(Cℓ_μRSD),
        _cℓ_term_eltype(Cℓ_RSDμ),
    )
    result = similar(Cℓ_δδ, T)
    @. result = Cℓ_δδ - Cℓ_δRSD - Cℓ_RSDδ + Cℓ_RSDRSD +
                Cℓ_δμ + Cℓ_μδ + Cℓ_μμ - Cℓ_μRSD - Cℓ_RSDμ
    return result
end

function _combine_gc_nonlimber(
    Cℓ_δδ, Cℓ_δRSD, Cℓ_RSDδ, Cℓ_RSDRSD,
    Cℓ_δμ, Cℓ_μδ, Cℓ_μμ, Cℓ_μRSD, Cℓ_RSDμ,
    Cℓ_δfNL, Cℓ_fNLδ, Cℓ_fNLRSD, Cℓ_RSDfNL,
    Cℓ_μfNL, Cℓ_fNLμ, Cℓ_fNLfNL,
)
    T = promote_type(
        _cℓ_term_eltype(Cℓ_δδ), _cℓ_term_eltype(Cℓ_δRSD),
        _cℓ_term_eltype(Cℓ_RSDδ), _cℓ_term_eltype(Cℓ_RSDRSD),
        _cℓ_term_eltype(Cℓ_δμ), _cℓ_term_eltype(Cℓ_μδ),
        _cℓ_term_eltype(Cℓ_μμ), _cℓ_term_eltype(Cℓ_μRSD),
        _cℓ_term_eltype(Cℓ_RSDμ), _cℓ_term_eltype(Cℓ_δfNL),
        _cℓ_term_eltype(Cℓ_fNLδ), _cℓ_term_eltype(Cℓ_fNLRSD),
        _cℓ_term_eltype(Cℓ_RSDfNL), _cℓ_term_eltype(Cℓ_μfNL),
        _cℓ_term_eltype(Cℓ_fNLμ), _cℓ_term_eltype(Cℓ_fNLfNL),
    )
    result = similar(Cℓ_δδ, T)
    @. result = Cℓ_δδ - Cℓ_δRSD - Cℓ_RSDδ + Cℓ_RSDRSD +
                Cℓ_δμ + Cℓ_μδ + Cℓ_μμ - Cℓ_μRSD - Cℓ_RSDμ +
                Cℓ_δfNL + Cℓ_fNLδ - Cℓ_fNLRSD - Cℓ_RSDfNL +
                Cℓ_μfNL + Cℓ_fNLμ + Cℓ_fNLfNL
    return result
end

function _gc_nonlimber_base_pairs(G, W, integ, W_δ, W_RSD, W_μ)
    Cℓ_δδ = _compute_Cℓ_cached(G.δ, G.δ, W_δ, W_δ, W.w_2_00_ϕTT, integ)
    Cℓ_RSDRSD = _compute_Cℓ_cached(G.RSD, G.RSD, W_RSD, W_RSD, W.w_2_22_ϕTT, integ)
    Cℓ_μμ = _compute_Cℓ_cached(G.μ, G.μ, W_μ, W_μ, W.w_minus2_00_ϕTT, integ)
    Cℓ_δμ = _compute_Cℓ_cached(G.δ, G.μ, W_δ, W_μ, W.w_0_00_ϕTT, integ)
    Cℓ_δRSD = _compute_Cℓ_rsd_cached(G.δ, G.RSD, W_δ, W_RSD, W.w_2_02_ϕTT, W.w_2_20_ϕTT, integ)
    Cℓ_μRSD = _compute_Cℓ_rsd_cached(G.μ, G.RSD, W_μ, W_RSD, W.w_0_02_ϕTT, W.w_0_20_ϕTT, integ)
    return (
        Cℓ_δδ, Cℓ_δRSD, _transpose_cℓ_pair(Cℓ_δRSD), Cℓ_RSDRSD,
        Cℓ_δμ, _transpose_cℓ_pair(Cℓ_δμ), Cℓ_μμ,
        Cℓ_μRSD, _transpose_cℓ_pair(Cℓ_μRSD),
    )
end

function _gc_nonlimber_png_pairs(G, W, integ, W_δ, W_RSD, W_μ, W_fNL)
    return (
        _compute_Cℓ_cached(G.δ, G.PNG, W_δ, W_fNL, W.w_2_00_ϕT_R1, integ),
        _compute_Cℓ_cached(G.PNG, G.δ, W_fNL, W_δ, W.w_2_00_ϕT, integ),
        _compute_Cℓ_rsd_cached(G.PNG, G.RSD, W_fNL, W_RSD, W.w_2_02_ϕT, W.w_2_20_ϕT, integ),
        _compute_Cℓ_rsd_cached(G.RSD, G.PNG, W_RSD, W_fNL, W.w_2_20_ϕT_R1, W.w_2_02_ϕT_R1, integ),
        _compute_Cℓ_cached(G.μ, G.PNG, W_μ, W_fNL, W.w_0_00_ϕT_R1, integ),
        _compute_Cℓ_cached(G.PNG, G.μ, W_fNL, W_μ, W.w_0_00_ϕT, integ),
        _compute_Cℓ_cached(G.PNG, G.PNG, W_fNL, W_fNL, W.w_2_00_ϕ, integ),
    )
end

function get_Cℓ(ℓ::AbstractArray{<:Any,1}, G::GalaxyClustering, Pk::PowerSpectrum, W::ProjectedMatterDensity, bg::Background, P::Union{FFTPlans, Nothing}=nothing)
    integ = _prepare_nonlimber_integration(bg)
    W_δ = _get_kernel_or_nothing(G.δ, bg)
    W_RSD = _get_kernel_or_nothing(G.RSD, bg)
    W_μ = _get_kernel_or_nothing(G.μ, bg)
    W_fNL = _get_kernel_or_nothing(G.PNG, bg)

    base_pairs = _gc_nonlimber_base_pairs(G, W, integ, W_δ, W_RSD, W_μ)

    Cℓ_nonlimber = if isnothing(G.PNG)
        # Keep the common no-PNG path free of scalar 0. sentinel terms.  Besides
        # avoiding unnecessary calls, this keeps the fused broadcast small enough
        # for Julia inference to materialize a concrete Array return type.
        _combine_gc_nonlimber(base_pairs...)
    else
        # PNG cross-pairs use _R1-sliced vs full weight tensors that are NOT
        # equal under index transpose, so they cannot be folded; computed
        # independently.
        png_pairs = _gc_nonlimber_png_pairs(G, W, integ, W_δ, W_RSD, W_μ, W_fNL)
        _combine_gc_nonlimber(base_pairs..., png_pairs...)
    end
    K_limber_G = get_limber_kernel(G)
    Cℓ_correction = get_limber_correction(K_limber_G, Pk)
    Cℓ_limber = get_limber_Cℓ(K_limber_G, Pk)

    nbins = size(G.δ.Kernel, 1)
    return _finalize_Cℓ_parts(Cℓ_nonlimber, Cℓ_correction, Cℓ_limber, ℓ, nbins, nbins, P)
end

## Weak lensing auto
function get_Cℓ(Component1::CosmicShear, Component2::CosmicShear, W::ProjectedMatterDensity, bg::Background)
    return compute_Cℓ(Component1, Component2, W.w_minus2_00_ϕTT, bg)
end

function get_Cℓ(Component1::CosmicShear, Component2::IntrinsicAlignment, W::ProjectedMatterDensity, bg::Background)
    return compute_Cℓ(Component1, Component2, W.w_minus2_00_ϕTT, bg)
end

function get_Cℓ(Component1::IntrinsicAlignment, Component2::CosmicShear, W::ProjectedMatterDensity, bg::Background)
    return compute_Cℓ(Component1, Component2, W.w_minus2_00_ϕTT, bg)
end

function get_Cℓ(Component1::IntrinsicAlignment, Component2::IntrinsicAlignment, W::ProjectedMatterDensity, bg::Background)
    return compute_Cℓ(Component1, Component2, W.w_minus2_00_ϕTT, bg)
end

function get_Cℓ(ℓ::AbstractArray{<:Any, 1}, L::WeakLensing, Pk::PowerSpectrum, W::ProjectedMatterDensity, bg::Background, P::Union{FFTPlans, Nothing}=nothing)
    integ = _prepare_nonlimber_integration(bg)
    W_γ = _get_kernel_or_nothing(L.γ, bg)
    W_IA = _get_kernel_or_nothing(L.IA, bg)

    Cℓ_γγ = _compute_Cℓ_cached(L.γ, L.γ, W_γ, W_γ, W.w_minus2_00_ϕTT, integ)
    Cℓ_γI = _compute_Cℓ_cached(L.γ, L.IA, W_γ, W_IA, W.w_minus2_00_ϕTT, integ)
    Cℓ_Iγ = _compute_Cℓ_cached(L.IA, L.γ, W_IA, W_γ, W.w_minus2_00_ϕTT, integ)
    Cℓ_II = _compute_Cℓ_cached(L.IA, L.IA, W_IA, W_IA, W.w_minus2_00_ϕTT, integ)

    Cℓ_nonlimber = @. Cℓ_γγ + Cℓ_γI + Cℓ_Iγ + Cℓ_II
    K_limber_L = get_limber_kernel(L)
    Cℓ_correction = get_limber_correction(K_limber_L, Pk)
    Cℓ_limber = get_limber_Cℓ(K_limber_L, Pk)

    nbins = size(L.γ.Kernel, 1)
    return _finalize_Cℓ_parts(Cℓ_nonlimber, Cℓ_correction, Cℓ_limber, ℓ, nbins, nbins, P)
end

## Cross clustering-lensing
function get_Cℓ(Component1::NumberCounts, Component2::CosmicShear, W::ProjectedMatterDensity, bg::Background)
    return compute_Cℓ(Component1, Component2, W.w_0_00_ϕTT, bg)
end

function get_Cℓ(Component1::NumberCounts, Component2::IntrinsicAlignment, W::ProjectedMatterDensity, bg::Background)
    return compute_Cℓ(Component1, Component2, W.w_0_00_ϕTT, bg)
end

function get_Cℓ(Component1::RedshiftSpaceDistortions, Component2::CosmicShear, W::ProjectedMatterDensity, bg::Background)
    return compute_Cℓ(Component1, Component2, W.w_0_20_ϕTT, W.w_0_02_ϕTT, bg)
end

function get_Cℓ(Component1::RedshiftSpaceDistortions, Component2::IntrinsicAlignment, W::ProjectedMatterDensity, bg::Background)
    return compute_Cℓ(Component1, Component2, W.w_0_20_ϕTT, W.w_0_02_ϕTT, bg)
end

function get_Cℓ(Component1::MagnificationBias, Component2::CosmicShear, W::ProjectedMatterDensity, bg::Background)
    return compute_Cℓ(Component1, Component2, W.w_minus2_00_ϕTT, bg)
end

function get_Cℓ(Component1::MagnificationBias, Component2::IntrinsicAlignment, W::ProjectedMatterDensity, bg::Background)
    return compute_Cℓ(Component1, Component2, W.w_minus2_00_ϕTT, bg)
end

function get_Cℓ(Component1::PrimordialNonGaussianity, Component2::CosmicShear, W::ProjectedMatterDensity, bg::Background)
    return compute_Cℓ(Component1, Component2, W.w_0_00_ϕT, bg)
end

function get_Cℓ(Component1::PrimordialNonGaussianity, Component2::IntrinsicAlignment, W::ProjectedMatterDensity, bg::Background)
    return compute_Cℓ(Component1, Component2, W.w_0_00_ϕT, bg)
end

function get_Cℓ(ℓ::AbstractArray{<:Any, 1}, G::GalaxyClustering, L::WeakLensing, Pk::PowerSpectrum, W::ProjectedMatterDensity, bg::Background, P::Union{FFTPlans, Nothing}=nothing)
    integ = _prepare_nonlimber_integration(bg)
    W_δ = _get_kernel_or_nothing(G.δ, bg)
    W_RSD = _get_kernel_or_nothing(G.RSD, bg)
    W_μ = _get_kernel_or_nothing(G.μ, bg)
    W_fNL = _get_kernel_or_nothing(G.PNG, bg)
    W_γ = _get_kernel_or_nothing(L.γ, bg)
    W_IA = _get_kernel_or_nothing(L.IA, bg)

    Cℓ_δγ = _compute_Cℓ_cached(G.δ, L.γ, W_δ, W_γ, W.w_0_00_ϕTT, integ)
    Cℓ_δI = _compute_Cℓ_cached(G.δ, L.IA, W_δ, W_IA, W.w_0_00_ϕTT, integ)
    Cℓ_RSDγ = _compute_Cℓ_rsd_cached(G.RSD, L.γ, W_RSD, W_γ, W.w_0_20_ϕTT, W.w_0_02_ϕTT, integ)
    Cℓ_RSDI = _compute_Cℓ_rsd_cached(G.RSD, L.IA, W_RSD, W_IA, W.w_0_20_ϕTT, W.w_0_02_ϕTT, integ)
    Cℓ_μγ = _compute_Cℓ_cached(G.μ, L.γ, W_μ, W_γ, W.w_minus2_00_ϕTT, integ)
    Cℓ_μI = _compute_Cℓ_cached(G.μ, L.IA, W_μ, W_IA, W.w_minus2_00_ϕTT, integ)
    Cℓ_fNLγ = _compute_Cℓ_cached(G.PNG, L.γ, W_fNL, W_γ, W.w_0_00_ϕT, integ)
    Cℓ_fNLI = _compute_Cℓ_cached(G.PNG, L.IA, W_fNL, W_IA, W.w_0_00_ϕT, integ)

    Cℓ_nonlimber = @. Cℓ_δγ + Cℓ_δI - Cℓ_RSDγ - Cℓ_RSDI + Cℓ_μγ + Cℓ_μI + Cℓ_fNLγ + Cℓ_fNLI
    K_limber_G = get_limber_kernel(G)
    K_limber_L = get_limber_kernel(L)
    Cℓ_correction = get_limber_correction(K_limber_G, K_limber_L, Pk)
    Cℓ_limber = get_limber_Cℓ(K_limber_G, K_limber_L, Pk)

    nbins_A = size(G.δ.Kernel, 1)
    nbins_B = size(L.γ.Kernel, 1)
    return _finalize_Cℓ_parts(Cℓ_nonlimber, Cℓ_correction, Cℓ_limber, ℓ, nbins_A, nbins_B, P)
end

function get_Cℓ(ℓ::AbstractArray{<:Any, 1}, L::WeakLensing, G::GalaxyClustering, Pk::PowerSpectrum, W::ProjectedMatterDensity, bg::Background, P::Union{FFTPlans, Nothing}=nothing)
    return get_Cℓ(ℓ, G, L, Pk, W, bg, P)
end

## Cross clustering - CMB 
function get_Cℓ(Component1::CMBLensing, Component2::NumberCounts, W::ProjectedMatterDensity, bg::Background)
    return compute_Cℓ(Component1, Component2, W.w_0_00_ϕTT, bg)
end

function get_Cℓ(Component1::CMBLensing, Component2::RedshiftSpaceDistortions, W::ProjectedMatterDensity, bg::Background)
    return compute_Cℓ(Component1, Component2, W.w_0_02_ϕTT, W.w_0_20_ϕTT, bg)
end

function get_Cℓ(Component1::CMBLensing, Component2::MagnificationBias, W::ProjectedMatterDensity, bg::Background)
    return compute_Cℓ(Component1, Component2, W.w_minus2_00_ϕTT, bg)
end

function get_Cℓ(Component1::CMBLensing, Component2::PrimordialNonGaussianity, W::ProjectedMatterDensity, bg::Background)
    return compute_Cℓ(Component1, Component2, W.w_0_00_ϕT_R1, bg)
end

function get_Cℓ(Component1::IntegratedSachsWolfe, Component2::NumberCounts, W::ProjectedMatterDensity, bg::Background)
    return compute_Cℓ(Component1, Component2, W.w_0_00_ϕTT, bg)
end

function get_Cℓ(Component1::IntegratedSachsWolfe, Component2::RedshiftSpaceDistortions, W::ProjectedMatterDensity, bg::Background)
    return compute_Cℓ(Component1, Component2, W.w_0_02_ϕTT, W.w_0_20_ϕTT, bg)
end

function get_Cℓ(Component1::IntegratedSachsWolfe, Component2::MagnificationBias, W::ProjectedMatterDensity, bg::Background)
    return compute_Cℓ(Component1, Component2, W.w_minus2_00_ϕTT, bg)
end

function get_Cℓ(Component1::IntegratedSachsWolfe, Component2::PrimordialNonGaussianity, W::ProjectedMatterDensity, bg::Background)
    return compute_Cℓ(Component1, Component2, W.w_0_00_ϕT_R1, bg)
end

function get_Cℓ(ℓ::AbstractArray{<:Any, 1}, K::CMB, G::GalaxyClustering, Pk::PowerSpectrum, W::ProjectedMatterDensity, bg::Background, P::Union{FFTPlans, Nothing}=nothing)
    Cℓ_κδ = get_Cℓ(K.κ, G.δ, W, bg)
    Cℓ_κRSD = get_Cℓ(K.κ, G.RSD, W, bg)
    Cℓ_κμ = get_Cℓ(K.κ, G.μ, W, bg)
    Cℓ_κfNL = get_Cℓ(K.κ, G.PNG, W, bg)
    Cℓ_Tδ = get_Cℓ(K.ISW, G.δ, W, bg)
    Cℓ_TRSD = get_Cℓ(K.ISW, G.RSD, W, bg)
    Cℓ_Tμ = get_Cℓ(K.ISW, G.μ, W, bg)
    Cℓ_TfNL = get_Cℓ(K.ISW, G.PNG, W, bg)

    Cℓ_nonlimber = @. Cℓ_κδ - Cℓ_κRSD + Cℓ_κμ + Cℓ_κfNL + Cℓ_Tδ - Cℓ_TRSD + Cℓ_Tμ + Cℓ_TfNL
    Cℓ_correction = get_limber_correction(K, G, Pk)
    Cℓ_limber = get_limber_Cℓ(K, G, Pk)
    full_Cℓ = cat(Cℓ_nonlimber .+ Cℓ_correction, Cℓ_limber; dims=1) 

    nbins_A = size(K.κ.Kernel, 1)
    nbins_B = size(G.δ.Kernel, 1)
    return _finalize_Cℓ(full_Cℓ, ℓ, nbins_A, nbins_B, P)
end

function get_Cℓ(ℓ::AbstractArray{<:Any, 1}, G::GalaxyClustering, K::CMB, Pk::PowerSpectrum, W::ProjectedMatterDensity, bg::Background, P::Union{FFTPlans, Nothing}=nothing)
    return get_Cℓ(ℓ, K, G, Pk, W, bg, P)
end

## Cross lensing - CMB
function get_Cℓ(Component1::CMBLensing, Component2::CosmicShear, W::ProjectedMatterDensity, bg::Background)
    return compute_Cℓ(Component1, Component2, W.w_minus2_00_ϕTT, bg)
end

function get_Cℓ(Component1::CMBLensing, Component2::IntrinsicAlignment, W::ProjectedMatterDensity, bg::Background)
    return compute_Cℓ(Component1, Component2, W.w_minus2_00_ϕTT, bg)
end

function get_Cℓ(Component1::IntegratedSachsWolfe, Component2::CosmicShear, W::ProjectedMatterDensity, bg::Background)
    return compute_Cℓ(Component1, Component2, W.w_minus2_00_ϕTT, bg)
end

function get_Cℓ(Component1::IntegratedSachsWolfe, Component2::IntrinsicAlignment, W::ProjectedMatterDensity, bg::Background)
    return compute_Cℓ(Component1, Component2, W.w_minus2_00_ϕTT, bg)
end

function get_Cℓ(ℓ::AbstractArray{<:Any, 1}, K::CMB, L::WeakLensing, Pk::PowerSpectrum, W::ProjectedMatterDensity, bg::Background, P::Union{FFTPlans, Nothing}=nothing)
    Cℓ_κγ = get_Cℓ(K.κ, L.γ, W, bg)
    Cℓ_κI = get_Cℓ(K.κ, L.IA, W, bg)
    Cℓ_Tγ = get_Cℓ(K.ISW, L.γ, W, bg)
    Cℓ_TI = get_Cℓ(K.ISW, L.IA, W, bg)

    Cℓ_nonlimber = @. Cℓ_κγ + Cℓ_κI + Cℓ_Tγ + Cℓ_TI
    Cℓ_correction = get_limber_correction(K, L, Pk)
    Cℓ_limber = get_limber_Cℓ(K, L, Pk)
    full_Cℓ = cat(Cℓ_nonlimber .+ Cℓ_correction, Cℓ_limber; dims=1) 

    nbins_A = size(K.κ.Kernel, 1)
    nbins_B = size(L.γ.Kernel, 1)
    return _finalize_Cℓ(full_Cℓ, ℓ, nbins_A, nbins_B, P)
end

function get_Cℓ(ℓ::AbstractArray{<:Any, 1}, S::WeakLensing, K::CMB, Pk::PowerSpectrum, W::ProjectedMatterDensity, bg::Background, P::Union{FFTPlans, Nothing}=nothing)
    return get_Cℓ(ℓ, K, S, Pk, W, bg, P)
end

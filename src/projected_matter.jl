"""
    w_ell_tullio(c, T)

Contract Chebyshev coefficients of the power spectrum with the precomputed `T̃`
to form the projected matter density `w_ℓ(χ_1, χ_2)`.

This performs the sum over Chebyshev indices using `Tullio` for efficient
tensor contraction. Multiple methods are provided depending on the dimensionality
of the coefficient array `c`.

# Arguments
- `c`: Chebyshev coefficients of the power spectrum.
- `T`: Precomputed kernel `T̃`.

# Returns
An array containing the projected matter density `w`.
"""
function w_ell_tullio(c::AbstractArray{<:Any, 3}, T::AbstractArray{<:Any, 4})
    return @tullio w[i,j,k] := c[l,j,k] * T[i,j,k,l]
end

function w_ell_tullio(c::AbstractArray{<:Any, 2}, T::AbstractArray{<:Any, 4})
    return @tullio w[i,j,k] := c[l,j] * T[i,j,k,l]
end

function w_ell_tullio(c::AbstractArray{<:Any, 1}, T::AbstractArray{<:Any, 4})
    return @tullio w[i,j,k] := c[l] * T[i,j,k,l]
end


# Mutating kernels used by `compute_w!`. These loops deliberately partition the
# output along `k`, so every thread owns disjoint output slices and no atomic
# accumulation is required. Keeping `i` innermost follows the contiguous memory
# layout of both `w` and `T̃`.
function _w_single_inplace!(w::AbstractArray{<:Any,3},
                            c::AbstractArray{<:Any,3},
                            T::AbstractArray{<:Any,4})
    fill!(w, zero(eltype(w)))
    Threads.@threads for k in axes(T, 3)
        @inbounds for l in axes(T, 4), j in axes(T, 2)
            coefficient = c[l, j, k]
            @simd for i in axes(T, 1)
                w[i, j, k] += coefficient * T[i, j, k, l]
            end
        end
    end
    return nothing
end

function _w_single_inplace!(w::AbstractArray{<:Any,3},
                            c::AbstractArray{<:Any,2},
                            T::AbstractArray{<:Any,4})
    fill!(w, zero(eltype(w)))
    Threads.@threads for k in axes(T, 3)
        @inbounds for l in axes(T, 4), j in axes(T, 2)
            coefficient = c[l, j]
            @simd for i in axes(T, 1)
                w[i, j, k] += coefficient * T[i, j, k, l]
            end
        end
    end
    return nothing
end

function _w_single_inplace!(w::AbstractArray{<:Any,3},
                            c::AbstractArray{<:Any,1},
                            T::AbstractArray{<:Any,4})
    fill!(w, zero(eltype(w)))
    Threads.@threads for k in axes(T, 3)
        @inbounds for l in axes(T, 4), j in axes(T, 2)
            coefficient = c[l]
            @simd for i in axes(T, 1)
                w[i, j, k] += coefficient * T[i, j, k, l]
            end
        end
    end
    return nothing
end

function _w_fused3_inplace!(w_tt, w_t, w_t_r1, c_tt, c_t, c_t_r1, T)
    fill!(w_tt, zero(eltype(w_tt)))
    fill!(w_t, zero(eltype(w_t)))
    fill!(w_t_r1, zero(eltype(w_t_r1)))
    Threads.@threads for k in axes(T, 3)
        @inbounds for l in axes(T, 4), j in axes(T, 2)
            tt_ljk = c_tt[l, j, k]
            t_ljk = c_t[l, j, k]
            t_r1_lj = c_t_r1[l, j]
            @simd for i in axes(T, 1)
                value = T[i, j, k, l]
                w_tt[i, j, k] += tt_ljk * value
                w_t[i, j, k] += t_ljk * value
                w_t_r1[i, j, k] += t_r1_lj * value
            end
        end
    end
    return nothing
end

function _w_fused4_inplace!(w_tt, w_t, w_t_r1, w_phi,
                            c_tt, c_t, c_t_r1, c_phi, T)
    fill!(w_tt, zero(eltype(w_tt)))
    fill!(w_t, zero(eltype(w_t)))
    fill!(w_t_r1, zero(eltype(w_t_r1)))
    fill!(w_phi, zero(eltype(w_phi)))
    Threads.@threads for k in axes(T, 3)
        @inbounds for l in axes(T, 4), j in axes(T, 2)
            tt_ljk = c_tt[l, j, k]
            t_ljk = c_t[l, j, k]
            t_r1_lj = c_t_r1[l, j]
            phi_l = c_phi[l]
            @simd for i in axes(T, 1)
                value = T[i, j, k, l]
                w_tt[i, j, k] += tt_ljk * value
                w_t[i, j, k] += t_ljk * value
                w_t_r1[i, j, k] += t_r1_lj * value
                w_phi[i, j, k] += phi_l * value
            end
        end
    end
    return nothing
end


abstract type AbstractProjectedMatterDensity end

"""
Abstract supertype for projected matter density components.

Each concrete subtype represents a specific combination, which depends on:
- The power of k in the precomputed `\\tilde{T}`.
- The order of the derivatives of the Bessel functions in the precomputed `\\tilde{T}`.
- The power spectrum (i.e. the full unequal time, the transfer function of the primordial power spectrum.)

Each component stores:
- A reference to the relevant `T̃`,
- The corresponding projected weight array `w`.
"""
abstract type ProjectedMatterDensityComponent end

# ────────────────────────────────────────────────────────────────────────────
# Parametric w_* structs (Phase B).
#
# Each struct carries a type parameter T matching the element type of `w`.
# `T̃` is always `Array{Float64, 4}` — a reference to one of the precomputed
# `T_tildes.T_*` constants, never differentiated wrt.
#
# The functional `compute_w(w, c)` returns a fresh struct whose array eltype
# follows `c`; this is the ForwardDiff path because Dual-valued coefficients
# require Dual-valued outputs. The separate `compute_w!` path reuses an
# explicitly allocated workspace for Float64 primal and Mooncake execution.
#
# All 17 struct definitions share the exact same layout and the exact same
# compute_w body, varying only in:
#   - which T_tildes.T_* constant to wrap,
#   - which c.c{ϕTT,ϕT,ϕ}.coefs to contract,
#   - whether to pre-slice those coefs to the last-R column (the `_R1`
#     variants reduce a 3D coefficient tensor to 2D, dispatching into the
#     rank-2 w_ell_tullio method).
# The table below drives a `@eval` loop that generates all 17 (struct,
# default-constructor, compute_w) triples. Adding a new component means
# appending one line to the table.
# ────────────────────────────────────────────────────────────────────────────

# (struct_name, T_tildes field, coef source ∈ {:ϕTT,:ϕT,:ϕ}, slice to R1 column?)
const _PROJECTED_MATTER_COMPONENTS = (
    (:w_2_00_ϕTT,      :T_2_00,      :ϕTT, false),
    (:w_minus2_00_ϕTT, :T_minus2_00, :ϕTT, false),
    (:w_0_00_ϕTT,      :T_0_00,      :ϕTT, false),
    (:w_0_02_ϕTT,      :T_0_02,      :ϕTT, false),
    (:w_0_20_ϕTT,      :T_0_20,      :ϕTT, false),
    (:w_2_02_ϕTT,      :T_2_02,      :ϕTT, false),
    (:w_2_20_ϕTT,      :T_2_20,      :ϕTT, false),
    (:w_2_22_ϕTT,      :T_2_22,      :ϕTT, false),
    (:w_2_00_ϕT,       :T_2_00,      :ϕT,  false),
    (:w_2_00_ϕT_R1,    :T_2_00,      :ϕT,  true),
    (:w_0_00_ϕT,       :T_0_00,      :ϕT,  false),
    (:w_0_00_ϕT_R1,    :T_0_00,      :ϕT,  true),
    (:w_2_02_ϕT,       :T_2_02,      :ϕT,  false),
    (:w_2_02_ϕT_R1,    :T_2_02,      :ϕT,  true),
    (:w_2_20_ϕT,       :T_2_20,      :ϕT,  false),
    (:w_2_20_ϕT_R1,    :T_2_20,      :ϕT,  true),
    (:w_2_00_ϕ,        :T_2_00,      :ϕ,   false),
)

for (_name, _T_tilde_field, _coef_source, _r1_slice) in _PROJECTED_MATTER_COMPONENTS
    _coef_field = Symbol(:c, _coef_source)  # :cϕTT, :cϕT, or :cϕ
    _coefs_expr = _r1_slice ? :(c.$_coef_field.coefs[:, :, end]) : :(c.$_coef_field.coefs)
    @eval begin
        mutable struct $_name{T<:Real} <: ProjectedMatterDensityComponent
            T̃::Array{Float64, 4}
            w::Array{T, 3}
        end
        $_name() = $_name{Float64}(T_tildes.$_T_tilde_field, zeros(Float64, 1, 1, 1))
        function compute_w(w::$_name, c::PowerSpectrum)
            result = w_ell_tullio($_coefs_expr, w.T̃)
            return $_name{eltype(result)}(w.T̃, result)
        end
        function allocate_compute_w(w::$_name, c::PowerSpectrum)
            output_type = promote_type(eltype(c.$_coef_field.coefs), eltype(w.T̃))
            output = zeros(output_type, size(w.T̃)[1:3])
            return $_name{output_type}(w.T̃, output)
        end
    end
end

# Nothing-carrying components: absent contribution, no-op.
compute_w(::Nothing, ::PowerSpectrum) = nothing
allocate_compute_w(::Nothing, ::PowerSpectrum) = nothing

"""
    ProjectedMatterDensity

Container holding all projected matter density components required for all the active observables and components.

Each field corresponds to a specific kernel contribution (e.g. `w_2_00_ϕTT`,
`w_0_02_ϕT_R1`). Fields may be set to `nothing` if the corresponding contribution
is not required.

The container is populated during setup and filled by calling `compute_w!`
with a `PowerSpectrum` object.
"""
@kwdef mutable struct ProjectedMatterDensity{
    T1, T2, T3, T4, T5, T6, T7, T8, T9,
    T10, T11, T12, T13, T14, T15, T16, T17
} <: AbstractProjectedMatterDensity
    w_2_00_ϕTT::T1         = nothing
    w_minus2_00_ϕTT::T2    = nothing
    w_0_00_ϕTT::T3         = nothing
    w_0_02_ϕTT::T4         = nothing
    w_0_20_ϕTT::T5         = nothing
    w_2_02_ϕTT::T6         = nothing
    w_2_20_ϕTT::T7         = nothing
    w_2_22_ϕTT::T8         = nothing
    w_2_00_ϕT::T9          = nothing
    w_2_00_ϕT_R1::T10      = nothing
    w_0_00_ϕT::T11         = nothing
    w_0_00_ϕT_R1::T12      = nothing
    w_2_02_ϕT::T13         = nothing
    w_2_02_ϕT_R1::T14      = nothing
    w_2_20_ϕT::T15         = nothing
    w_2_20_ϕT_R1::T16      = nothing
    w_2_00_ϕ::T17          = nothing
end

"""
    allocate_compute_w(W, c)

Allocate a `ProjectedMatterDensity` workspace with the same active components
as `W`, sized for the precomputed `T̃` tensors and with an element type inferred
from `c`. Reuse the returned workspace with [`compute_w!`](@ref).

The returned workspace is mutable and belongs to one inference execution. Do
not share it between concurrently evaluated chains or tasks; allocate one
workspace per concurrent execution.
"""
function allocate_compute_w(W::ProjectedMatterDensity, c::PowerSpectrum)
    return ProjectedMatterDensity(
        w_2_00_ϕTT      = allocate_compute_w(W.w_2_00_ϕTT,      c),
        w_minus2_00_ϕTT = allocate_compute_w(W.w_minus2_00_ϕTT, c),
        w_0_00_ϕTT      = allocate_compute_w(W.w_0_00_ϕTT,      c),
        w_0_02_ϕTT      = allocate_compute_w(W.w_0_02_ϕTT,      c),
        w_0_20_ϕTT      = allocate_compute_w(W.w_0_20_ϕTT,      c),
        w_2_02_ϕTT      = allocate_compute_w(W.w_2_02_ϕTT,      c),
        w_2_20_ϕTT      = allocate_compute_w(W.w_2_20_ϕTT,      c),
        w_2_22_ϕTT      = allocate_compute_w(W.w_2_22_ϕTT,      c),
        w_2_00_ϕT       = allocate_compute_w(W.w_2_00_ϕT,       c),
        w_2_00_ϕT_R1    = allocate_compute_w(W.w_2_00_ϕT_R1,    c),
        w_0_00_ϕT       = allocate_compute_w(W.w_0_00_ϕT,       c),
        w_0_00_ϕT_R1    = allocate_compute_w(W.w_0_00_ϕT_R1,    c),
        w_2_02_ϕT       = allocate_compute_w(W.w_2_02_ϕT,       c),
        w_2_02_ϕT_R1    = allocate_compute_w(W.w_2_02_ϕT_R1,    c),
        w_2_20_ϕT       = allocate_compute_w(W.w_2_20_ϕT,       c),
        w_2_20_ϕT_R1    = allocate_compute_w(W.w_2_20_ϕT_R1,    c),
        w_2_00_ϕ        = allocate_compute_w(W.w_2_00_ϕ,        c),
    )
end

"""
    compute_w(W::ProjectedMatterDensity, c::PowerSpectrum)

Compute all active projected matter density components and return a fresh
`ProjectedMatterDensity` whose field eltypes track the eltype of `c` (Float64
for primal, Dual under ForwardDiff). This is functional (no mutation): the
input `W` is unchanged. Required because `ProjectedMatterDensity` is
parameterized on each field's concrete type, so an AD eltype change cannot
be absorbed by in-place reassignment.

Nothing-valued fields stay `nothing` via the `compute_w(::Nothing, ...)`
method above.
"""
function compute_w(W::ProjectedMatterDensity, c::PowerSpectrum)
    # Keep the container-level path sequential. Each component contraction is
    # already handled by the lower-level Tullio kernel; wrapping all fields in
    # Threads.@spawn adds scheduler overhead, gives no meaningful overlap on the
    # realistic 3×2pt workload, and introduces task/foreigncall nodes that AD
    # backends such as Mooncake cannot differentiate through.
    return ProjectedMatterDensity(
        w_2_00_ϕTT      = compute_w(W.w_2_00_ϕTT,      c),
        w_minus2_00_ϕTT = compute_w(W.w_minus2_00_ϕTT, c),
        w_0_00_ϕTT      = compute_w(W.w_0_00_ϕTT,      c),
        w_0_02_ϕTT      = compute_w(W.w_0_02_ϕTT,      c),
        w_0_20_ϕTT      = compute_w(W.w_0_20_ϕTT,      c),
        w_2_02_ϕTT      = compute_w(W.w_2_02_ϕTT,      c),
        w_2_20_ϕTT      = compute_w(W.w_2_20_ϕTT,      c),
        w_2_22_ϕTT      = compute_w(W.w_2_22_ϕTT,      c),
        w_2_00_ϕT       = compute_w(W.w_2_00_ϕT,       c),
        w_2_00_ϕT_R1    = compute_w(W.w_2_00_ϕT_R1,    c),
        w_0_00_ϕT       = compute_w(W.w_0_00_ϕT,       c),
        w_0_00_ϕT_R1    = compute_w(W.w_0_00_ϕT_R1,    c),
        w_2_02_ϕT       = compute_w(W.w_2_02_ϕT,       c),
        w_2_02_ϕT_R1    = compute_w(W.w_2_02_ϕT_R1,    c),
        w_2_20_ϕT       = compute_w(W.w_2_20_ϕT,       c),
        w_2_20_ϕT_R1    = compute_w(W.w_2_20_ϕT_R1,    c),
        w_2_00_ϕ        = compute_w(W.w_2_00_ϕ,        c),
    )
end

_compute_w_single_inplace!(::Nothing, c) = nothing
function _compute_w_single_inplace!(component, c)
    _w_single_inplace!(component.w, c, component.T̃)
    return nothing
end

"""
    compute_w!(workspace, c)

Overwrite a preallocated projected-matter workspace. Components sharing a
`T̃` tensor are contracted together, reducing the realistic 17 tensor scans to
eight. Construct `workspace` once with [`allocate_compute_w`](@ref).

The workspace element type must match the coefficient element type. Use the
functional [`compute_w`](@ref) path when differentiating with ForwardDiff,
because its output arrays must carry `Dual` values.

This function mutates `workspace` and is not safe for concurrent calls on the
same object. Independent tasks/chains must use independent workspaces.
"""
function compute_w!(W::ProjectedMatterDensity, c::PowerSpectrum)
    c_tt = c.cϕTT.coefs
    c_t = isnothing(c.cϕT) ? nothing : c.cϕT.coefs
    c_t_r1 = isnothing(c_t) ? nothing : @view(c_t[:, :, end])
    c_phi = isnothing(c.cϕ) ? nothing : c.cϕ.coefs

    if !isnothing(W.w_2_00_ϕTT) && !isnothing(W.w_2_00_ϕT) &&
       !isnothing(W.w_2_00_ϕT_R1) && !isnothing(W.w_2_00_ϕ)
        _w_fused4_inplace!(
            W.w_2_00_ϕTT.w, W.w_2_00_ϕT.w,
            W.w_2_00_ϕT_R1.w, W.w_2_00_ϕ.w,
            c_tt, c_t, c_t_r1, c_phi, W.w_2_00_ϕTT.T̃,
        )
    else
        _compute_w_single_inplace!(W.w_2_00_ϕTT, c_tt)
        _compute_w_single_inplace!(W.w_2_00_ϕT, c_t)
        _compute_w_single_inplace!(W.w_2_00_ϕT_R1, c_t_r1)
        _compute_w_single_inplace!(W.w_2_00_ϕ, c_phi)
    end

    if !isnothing(W.w_0_00_ϕTT) && !isnothing(W.w_0_00_ϕT) &&
       !isnothing(W.w_0_00_ϕT_R1)
        _w_fused3_inplace!(
            W.w_0_00_ϕTT.w, W.w_0_00_ϕT.w, W.w_0_00_ϕT_R1.w,
            c_tt, c_t, c_t_r1, W.w_0_00_ϕTT.T̃,
        )
    else
        _compute_w_single_inplace!(W.w_0_00_ϕTT, c_tt)
        _compute_w_single_inplace!(W.w_0_00_ϕT, c_t)
        _compute_w_single_inplace!(W.w_0_00_ϕT_R1, c_t_r1)
    end

    if !isnothing(W.w_2_02_ϕTT) && !isnothing(W.w_2_02_ϕT) &&
       !isnothing(W.w_2_02_ϕT_R1)
        _w_fused3_inplace!(
            W.w_2_02_ϕTT.w, W.w_2_02_ϕT.w, W.w_2_02_ϕT_R1.w,
            c_tt, c_t, c_t_r1, W.w_2_02_ϕTT.T̃,
        )
    else
        _compute_w_single_inplace!(W.w_2_02_ϕTT, c_tt)
        _compute_w_single_inplace!(W.w_2_02_ϕT, c_t)
        _compute_w_single_inplace!(W.w_2_02_ϕT_R1, c_t_r1)
    end

    if !isnothing(W.w_2_20_ϕTT) && !isnothing(W.w_2_20_ϕT) &&
       !isnothing(W.w_2_20_ϕT_R1)
        _w_fused3_inplace!(
            W.w_2_20_ϕTT.w, W.w_2_20_ϕT.w, W.w_2_20_ϕT_R1.w,
            c_tt, c_t, c_t_r1, W.w_2_20_ϕTT.T̃,
        )
    else
        _compute_w_single_inplace!(W.w_2_20_ϕTT, c_tt)
        _compute_w_single_inplace!(W.w_2_20_ϕT, c_t)
        _compute_w_single_inplace!(W.w_2_20_ϕT_R1, c_t_r1)
    end

    _compute_w_single_inplace!(W.w_minus2_00_ϕTT, c_tt)
    _compute_w_single_inplace!(W.w_0_02_ϕTT, c_tt)
    _compute_w_single_inplace!(W.w_0_20_ϕTT, c_tt)
    _compute_w_single_inplace!(W.w_2_22_ϕTT, c_tt)
    return W
end

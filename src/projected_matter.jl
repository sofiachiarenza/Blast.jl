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
# `compute_w!(w::w_*, c::PowerSpectrum)` is replaced by a functional
# `compute_w(w, c)` that returns a fresh struct of the appropriate eltype. The
# container-level `compute_w!(W::ProjectedMatterDensity, c)` reassigns fields
# with the result — no in-place mutation of `w.w`. This eliminates the
# Dual→Float64 strip when `c.cϕTT.coefs` carries Duals (ForwardDiff through
# cosmology parameters).
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
    end
end

# Nothing-carrying components: absent contribution, no-op.
compute_w(::Nothing, ::PowerSpectrum) = nothing

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

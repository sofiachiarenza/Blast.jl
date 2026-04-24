"""
    AbstractCosmologicalProbes

Abstract supertype for cosmological probes in BLAST.
"""
abstract type AbstractCosmologicalProbes end

"""
    AbstractComponents

Abstract supertype for physical components entering cosmological probes.
"""
abstract type AbstractComponents end

"""
    _trapezoidal_norms(nz::AbstractMatrix, z::AbstractVector) -> Vector

Per-bin trapezoidal integral ``∫ n(z)\\,dz`` over the user's ``z`` grid.

Used as a cosmology-independent normalization constant: since realistic
photometric ``n(z)`` distributions taper to zero at both ends of the user's
``z`` support, the integral computed on the user's grid matches the
integral over any wider grid (e.g. `bg.z`) to far better than per-mille.
This lets us precompute norms **once at component construction** and skip
QuadGK entirely on every `evaluate_components!` call.

Trapezoidal is deliberate: for densely-sampled ``n(z)`` (typical photo-z
bins have ≥ 100 samples), the trapezoidal error on a smooth ``n(z)`` is
``O((dz)^2)`` — negligible against the per-cent-level systematic budget
of realistic surveys. It is pure broadcast arithmetic so it propagates
ForwardDiff `Dual`s and Mooncake tangents natively, with no sensealg
needed.

Compare: the old path built an adaptive QuadGK + akima_interpolation
spline for every quadrature sample, rebuilt the spline's slope &
coefficient arrays from scratch on each call, and leaned on an
`Integrals.jl + ReCallVJP{MooncakeVJP}` wrapper — ~1 s / 11 GiB per
`f_full` evaluation on a 10-bin × 500-knot problem. Trapezoidal on the
same data is ~10 µs / ~1 KiB.
"""
function _trapezoidal_norms(nz::AbstractMatrix, z::AbstractVector)
    n_bins, n_zpts = size(nz)
    T = promote_type(eltype(nz), eltype(z))
    if n_zpts < 2
        # Degenerate default-constructor case: return zeros so `zeros(T, n_bins)`
        # doesn't crash the promotion logic; real `nz` inputs always have
        # ≥2 knots and take the loop below.
        return zeros(T, n_bins)
    end
    norms = Vector{T}(undef, n_bins)
    @inbounds for b in 1:n_bins
        s = zero(T)
        for i in 1:n_zpts-1
            s += (nz[b, i] + nz[b, i+1]) * (z[i+1] - z[i]) * 0.5
        end
        norms[b] = s
    end
    return norms
end

"""
    prepare_nz_matrix(nz, z, z_grid, nz_norms)
    prepare_nz_matrix(nz, z, z_grid)

Interpolate an n(z) matrix bin-by-bin onto the calculation grid and
normalize by the supplied `nz_norms` vector (one entry per bin).

The 4-arg form takes precomputed norms — this is the hot path used by
`check_and_normalize!`. Every component with `has_nz == true` carries
`nz_norms` precomputed at construction time (see `_trapezoidal_norms`),
so this form only ever does a single batched akima interpolation plus
one broadcast division.

The 3-arg form (kept for tests and direct user calls) computes
trapezoidal norms internally and forwards.
"""
function prepare_nz_matrix(nz::AbstractMatrix, z::AbstractVector,
                           z_grid::AbstractVector, nz_norms::AbstractVector)
    # Vectorized Akima interpolation: ACE's matrix dispatch expects
    # (length(z), n_bins) layout and returns (length(z_grid), n_bins).
    interpolated = akima_interpolation(permutedims(nz), z, z_grid)
    # Normalize each bin, then return to (n_bins, length(z_grid)) layout.
    return permutedims(interpolated ./ reshape(nz_norms, 1, :))
end

function prepare_nz_matrix(nz::AbstractMatrix, z::AbstractVector,
                           z_grid::AbstractVector)
    return prepare_nz_matrix(nz, z, z_grid, _trapezoidal_norms(nz, z))
end

# NOTE: The previous `_compute_nz_norms` — which wrapped QuadGK adaptive
# quadrature through `Integrals.jl` + `ReCallVJP{MooncakeVJP}` sensealg and
# rebuilt an akima spline from scratch on every quadrature sample — has been
# replaced by `_trapezoidal_norms`. See that docstring for the reasoning.
# The old path cost ~1 s / 11 GiB per `evaluate_components!(GC, bg)` call
# (87 % of the whole `f_full` pipeline, ≥ 2 M allocations) on a 10-bin ×
# 500-knot realistic fixture; trapezoidal on the same inputs is microseconds
# and handful-of-KiB, at an accuracy far below any realistic systematic
# floor.

"""
    smooth_nz(nz::AbstractMatrix, z::AbstractVector; z_out::AbstractVector=z, span::Real=0.02)

Smooth n(z) bin-by-bin with `Loess.jl` and evaluate on `z_out`.
"""
function smooth_nz(
    nz::AbstractMatrix,
    z::AbstractVector,
    ;
    z_out::AbstractVector=z,
    span::Real=0.02
)
    size(nz, 2) == length(z) || throw(DimensionMismatch("size(nz, 2) must match length(z)."))
    0 < span <= 1 || throw(ArgumentError("span must satisfy 0 < span <= 1."))

    n_bins = size(nz, 1)
    nz_smooth = Matrix{Float64}(undef, n_bins, length(z_out))

    for b in 1:n_bins
        model = Loess.loess(z, nz[b, :], span=span)
        vs = collect(Loess.predict(model, z_out))
        vs[vs .< 0] .= 0.0
        nz_smooth[b, :] = vs
    end

    return nz_smooth
end

"""
    has_nz(::AbstractComponents) → Bool

Trait indicating whether a component carries an n(z) redshift distribution
that must be normalised before kernel computation. Defaults to `false`.
"""
has_nz(::AbstractComponents) = false
has_nz(::Nothing) = false

"""
    check_and_normalize!(Component, grid_z)

Internal helper: ensures nz_norm is populated for the current calculation grid.
Uses the `has_nz` trait to avoid `hasfield` runtime overhead.
"""
function check_and_normalize!(Component::AbstractComponents, z_grid::AbstractVector)
    if has_nz(Component)
        if size(Component.nz_norm) != (size(Component.nz, 1), length(z_grid))
            Component.nz_norm = prepare_nz_matrix(Component.nz, Component.z,
                                                  z_grid, Component.nz_norms)
        end
    end
end

function check_and_normalize!(Component::Nothing, z_grid::AbstractVector)
    return nothing
end

"""
    NLA_model(bg::Background; A=1.72, C1=0.0134)

Computes the Non-Linear Alignment (NLA) model for intrinsic alignments.
Returns an array evaluated on the Background grid.
"""
function NLA_model(bg::Background; A=1.72, C1=0.0134)
    Ωm = get_Ωm(bg)
    return @. - A * C1 * Ωm / bg.D
end

@doc raw"""
    NumberCounts{T<:Real} <: AbstractComponents

Galaxy number-counts density component.

This component stores the redshift distribution and linear bias entering

```math
K_i^{\delta}(z) = \frac{H(z)}{c}\, b_i(z)\, \hat n_i(z).
```

The type parameter `T` is the element type of the Dual-carrying fields
(`nz`, `nz_norm`, `bias`, `Kernel`). It defaults to `Float64` and is
inferred automatically from the input arrays by the keyword constructor,
so `NumberCounts(nz=..., z=..., bias=...)` works as before. When any input
is a `ForwardDiff.Dual`, `T` promotes accordingly and ForwardDiff can
propagate Duals through the struct without being stripped by `convert`.

The `z` knot grid, `ell_prefactor`, and `limber_factor` are always
`Float64` — they're never differentiated wrt.
"""
mutable struct NumberCounts{T<:Real} <: AbstractComponents
    nz::Matrix{T}
    z::Vector{Float64}
    nz_norm::Matrix{T}
    bias::Matrix{T}
    nz_norms::Vector{T}
    Kernel::Matrix{T}
    ell_prefactor::Vector{Float64}
    limber_factor::Vector{Float64}
end

function NumberCounts(;
    nz::AbstractMatrix            = zeros(Float64, 1, 1),
    z::AbstractVector             = zeros(Float64, 1),
    bias::AbstractMatrix          = zeros(Float64, 1, 1),
    nz_norm::AbstractMatrix       = zeros(eltype(nz), 1, 1),
    nz_norms::AbstractVector      = _trapezoidal_norms(nz, z),
    Kernel::AbstractMatrix        = zeros(eltype(bias), 1, 1),
    ell_prefactor::AbstractVector = ones(size(Blast.full_ℓ_range, 1)),
    limber_factor::AbstractVector = ones(size(Blast.full_ℓ_range, 1)),
)
    T = promote_type(eltype(nz), eltype(bias), eltype(nz_norm),
                     eltype(nz_norms), eltype(Kernel))
    NumberCounts{T}(
        convert(Matrix{T},       nz),
        convert(Vector{Float64}, z),
        convert(Matrix{T},       nz_norm),
        convert(Matrix{T},       bias),
        convert(Vector{T},       nz_norms),
        convert(Matrix{T},       Kernel),
        convert(Vector{Float64}, ell_prefactor),
        convert(Vector{Float64}, limber_factor),
    )
end
has_nz(::NumberCounts) = true

@doc raw"""
    CosmicShear <: AbstractComponents

Weak-lensing shear component with kernel

```math
K_i^{\gamma}(\chi) = \frac{3 H_0^2 \Omega_m}{2 c^2}\, \chi (1+z)
\int_z^{\infty} \mathrm{d}z'\, \frac{H(z')}{c}\, \hat n_i(z')
\frac{\chi(z')-\chi}{\chi(z')}.
```
"""
mutable struct CosmicShear{T<:Real} <: AbstractComponents
    nz::Matrix{T}
    z::Vector{Float64}
    nz_norm::Matrix{T}
    nz_norms::Vector{T}
    Kernel::Matrix{T}
    ell_prefactor::Vector{Float64}
    limber_factor::Vector{Float64}
end

function CosmicShear(;
    nz::AbstractMatrix            = zeros(Float64, 1, 1),
    z::AbstractVector             = zeros(Float64, 1),
    nz_norm::AbstractMatrix       = zeros(eltype(nz), 1, 1),
    nz_norms::AbstractVector      = _trapezoidal_norms(nz, z),
    Kernel::AbstractMatrix        = zeros(eltype(nz), 1, 1),
    ell_prefactor::AbstractVector = sqrt.(factorial_frac.(Blast.full_ℓ_range)),
    limber_factor::AbstractVector = (Blast.full_ℓ_range .+ 0.5) .^ (-2),
)
    T = promote_type(eltype(nz), eltype(nz_norm), eltype(nz_norms), eltype(Kernel))
    CosmicShear{T}(
        convert(Matrix{T},       nz),
        convert(Vector{Float64}, z),
        convert(Matrix{T},       nz_norm),
        convert(Vector{T},       nz_norms),
        convert(Matrix{T},       Kernel),
        convert(Vector{Float64}, ell_prefactor),
        convert(Vector{Float64}, limber_factor),
    )
end
has_nz(::CosmicShear) = true

@doc raw"""
    CMBLensing <: AbstractComponents

CMB lensing convergence component with source plane at recombination:

```math
K^{\kappa_{\mathrm{CMB}}}(\chi) = \frac{3 H_0^2 \Omega_m}{2 c^2}\, \chi (1+z)
\left(1 - \frac{\chi}{\chi_*}\right).
```
"""
mutable struct CMBLensing{T<:Real} <: AbstractComponents
    Kernel::Matrix{T}
    ell_prefactor::Vector{Float64}
    limber_factor::Vector{Float64}
end

function CMBLensing(;
    Kernel::AbstractMatrix        = zeros(Float64, 1, 1),
    ell_prefactor::AbstractVector = Blast.full_ℓ_range .* (Blast.full_ℓ_range .+ 1),
    limber_factor::AbstractVector = (Blast.full_ℓ_range .+ 0.5) .^ (-2),
)
    T = eltype(Kernel)
    CMBLensing{T}(
        convert(Matrix{T},       Kernel),
        convert(Vector{Float64}, ell_prefactor),
        convert(Vector{Float64}, limber_factor),
    )
end

@doc raw"""
    RedshiftSpaceDistortions <: AbstractComponents

Redshift-space distortion component with kernel

```math
K_i^{\mathrm{RSD}}(z) = \frac{H(z)}{c}\, f(z)\, \hat n_i(z).
```
"""
mutable struct RedshiftSpaceDistortions{T<:Real} <: AbstractComponents
    nz::Matrix{T}
    z::Vector{Float64}
    nz_norm::Matrix{T}
    growth_rate::Vector{T}
    nz_norms::Vector{T}
    Kernel::Matrix{T}
    ell_prefactor::Vector{Float64}
    limber_factor::Vector{Float64}
end

function RedshiftSpaceDistortions(;
    nz::AbstractMatrix            = zeros(Float64, 1, 1),
    z::AbstractVector             = zeros(Float64, 1),
    nz_norm::AbstractMatrix       = zeros(eltype(nz), 1, 1),
    growth_rate::AbstractVector   = zeros(eltype(nz), 1),
    nz_norms::AbstractVector      = _trapezoidal_norms(nz, z),
    Kernel::AbstractMatrix        = zeros(eltype(nz), 1, 1),
    ell_prefactor::AbstractVector = ones(size(Blast.full_ℓ_range, 1)),
    limber_factor::AbstractVector = ones(size(Blast.full_ℓ_range, 1)),
)
    T = promote_type(eltype(nz), eltype(nz_norm), eltype(growth_rate),
                     eltype(nz_norms), eltype(Kernel))
    RedshiftSpaceDistortions{T}(
        convert(Matrix{T},       nz),
        convert(Vector{Float64}, z),
        convert(Matrix{T},       nz_norm),
        convert(Vector{T},       growth_rate),
        convert(Vector{T},       nz_norms),
        convert(Matrix{T},       Kernel),
        convert(Vector{Float64}, ell_prefactor),
        convert(Vector{Float64}, limber_factor),
    )
end
has_nz(::RedshiftSpaceDistortions) = true

@doc raw"""
    MagnificationBias <: AbstractComponents

Magnification-bias component weighted by the source-slope factor $5s-2$.

Implemented as the shear-like lensing kernel times the extra slope factor:

```math
K_i^{\mu}(\chi) = \frac{3 H_0^2 \Omega_m}{2 c^2}\, \chi (1+z)
\int_z^{\infty} \mathrm{d}z'\, \frac{H(z')}{c}\, \hat n_i(z')
\frac{\chi(z')-\chi}{\chi(z')}\,\left[5 s_i(z')-2\right].
```
"""
mutable struct MagnificationBias{T<:Real} <: AbstractComponents
    nz::Matrix{T}
    z::Vector{Float64}
    nz_norm::Matrix{T}
    s::Matrix{T}
    nz_norms::Vector{T}
    Kernel::Matrix{T}
    ell_prefactor::Vector{Float64}
    limber_factor::Vector{Float64}
end

function MagnificationBias(;
    nz::AbstractMatrix            = zeros(Float64, 1, 1),
    z::AbstractVector             = zeros(Float64, 1),
    nz_norm::AbstractMatrix       = zeros(eltype(nz), 1, 1),
    s::AbstractMatrix             = zeros(eltype(nz), 1, 1),
    nz_norms::AbstractVector      = _trapezoidal_norms(nz, z),
    Kernel::AbstractMatrix        = zeros(eltype(nz), 1, 1),
    ell_prefactor::AbstractVector = Blast.full_ℓ_range .* (Blast.full_ℓ_range .+ 1),
    limber_factor::AbstractVector = (Blast.full_ℓ_range .+ 0.5) .^ (-2),
)
    T = promote_type(eltype(nz), eltype(nz_norm), eltype(s),
                     eltype(nz_norms), eltype(Kernel))
    MagnificationBias{T}(
        convert(Matrix{T},       nz),
        convert(Vector{Float64}, z),
        convert(Matrix{T},       nz_norm),
        convert(Matrix{T},       s),
        convert(Vector{T},       nz_norms),
        convert(Matrix{T},       Kernel),
        convert(Vector{Float64}, ell_prefactor),
        convert(Vector{Float64}, limber_factor),
    )
end
has_nz(::MagnificationBias) = true

@doc raw"""
    IntrinsicAlignment <: AbstractComponents

Intrinsic-alignment component with NLA-style amplitude model:

Users can provide a generic amplitude table `A_IA` (same bin/redshift shape as
`nz_norm`) to be used directly in the kernel. If `A_IA` is not provided with the
required shape, BLAST falls back to the standard NLA model controlled by `A`.

```math
K_i^{\mathrm{IA}}(z) = \frac{H(z)}{c}\, A_{\mathrm{IA},i}(z)\, \hat n_i(z),
\qquad
A_{\mathrm{IA}}(z) = - A\, C_1\, \frac{\Omega_m}{D(z)}.
```
"""
mutable struct IntrinsicAlignment{T<:Real} <: AbstractComponents
    nz::Matrix{T}
    z::Vector{Float64}
    nz_norm::Matrix{T}
    A::T                 # Standard NLA amplitude (differentiable)
    C1::Float64          # NLA coupling constant; fixed
    A_IA::Matrix{T}
    nz_norms::Vector{T}
    Kernel::Matrix{T}
    ell_prefactor::Vector{Float64}
    limber_factor::Vector{Float64}
end

function IntrinsicAlignment(;
    nz::AbstractMatrix            = zeros(Float64, 1, 1),
    z::AbstractVector             = zeros(Float64, 1),
    nz_norm::AbstractMatrix       = zeros(eltype(nz), 1, 1),
    A::Real                       = 1.72,
    C1::Real                      = 0.0134,
    A_IA::AbstractMatrix          = zeros(eltype(nz), 1, 1),
    nz_norms::AbstractVector      = _trapezoidal_norms(nz, z),
    Kernel::AbstractMatrix        = zeros(eltype(nz), 1, 1),
    ell_prefactor::AbstractVector = sqrt.(factorial_frac.(Blast.full_ℓ_range)),
    limber_factor::AbstractVector = (Blast.full_ℓ_range .+ 0.5) .^ (-2),
)
    T = promote_type(eltype(nz), eltype(nz_norm), typeof(A), eltype(A_IA),
                     eltype(nz_norms), eltype(Kernel))
    IntrinsicAlignment{T}(
        convert(Matrix{T},       nz),
        convert(Vector{Float64}, z),
        convert(Matrix{T},       nz_norm),
        convert(T,               A),
        convert(Float64,         C1),
        convert(Matrix{T},       A_IA),
        convert(Vector{T},       nz_norms),
        convert(Matrix{T},       Kernel),
        convert(Vector{Float64}, ell_prefactor),
        convert(Vector{Float64}, limber_factor),
    )
end
has_nz(::IntrinsicAlignment) = true

@doc raw"""
    IntegratedSachsWolfe <: AbstractComponents

Integrated Sachs-Wolfe component sourced by the time variation of the
gravitational potential:

```math
K^{\mathrm{ISW}}(z) = \frac{3 T_{\mathrm{CMB}} H_0^2 \Omega_m}{c^3}\, H(z)\, \left[1-f(z)\right].
```
"""
mutable struct IntegratedSachsWolfe{T<:Real} <: AbstractComponents
    growth_rate::Vector{T}
    Kernel::Matrix{T}
    ell_prefactor::Vector{Float64}
    limber_factor::Vector{Float64}
end

function IntegratedSachsWolfe(;
    growth_rate::AbstractVector   = zeros(Float64, 1),
    Kernel::AbstractMatrix        = zeros(eltype(growth_rate), 1, 1),
    ell_prefactor::AbstractVector = ones(size(Blast.full_ℓ_range, 1)),
    limber_factor::AbstractVector = ones(size(Blast.full_ℓ_range, 1)),
)
    T = promote_type(eltype(growth_rate), eltype(Kernel))
    IntegratedSachsWolfe{T}(
        convert(Vector{T},       growth_rate),
        convert(Matrix{T},       Kernel),
        convert(Vector{Float64}, ell_prefactor),
        convert(Vector{Float64}, limber_factor),
    )
end

@doc raw"""
    PrimordialNonGaussianity <: AbstractComponents

Scale-dependent bias contribution induced by local-type primordial
non-Gaussianity:

```math
K_i^{\mathrm{PNG}}(z) = \frac{H(z)}{c}\, f_{\mathrm{NL}}\, b_{\Phi,i}(z)\, \hat n_i(z),
\qquad
b_{\Phi,i}(z) = 2\,\delta_c\,\left[b_i(z)-p\right].
```
"""
mutable struct PrimordialNonGaussianity{T<:Real} <: AbstractComponents
    nz::Matrix{T}
    z::Vector{Float64}
    nz_norm::Matrix{T}
    bias::Matrix{T}
    f_NL::T               # local-type fNL amplitude (differentiable)
    p::Float64            # tracer-type parameter; fixed
    nz_norms::Vector{T}
    Kernel::Matrix{T}
    ell_prefactor::Vector{Float64}
    limber_factor::Vector{Float64}
end

function PrimordialNonGaussianity(;
    nz::AbstractMatrix            = zeros(Float64, 1, 1),
    z::AbstractVector             = zeros(Float64, 1),
    nz_norm::AbstractMatrix       = zeros(eltype(nz), 1, 1),
    bias::AbstractMatrix          = zeros(eltype(nz), 1, 1),
    f_NL::Real                    = 0.0,
    p::Real                       = 0.0,
    nz_norms::AbstractVector      = _trapezoidal_norms(nz, z),
    Kernel::AbstractMatrix        = zeros(eltype(bias), 1, 1),
    ell_prefactor::AbstractVector = ones(size(Blast.full_ℓ_range, 1)),
    limber_factor::AbstractVector = ones(size(Blast.full_ℓ_range, 1)),
)
    T = promote_type(eltype(nz), eltype(nz_norm), eltype(bias), typeof(f_NL),
                     eltype(nz_norms), eltype(Kernel))
    PrimordialNonGaussianity{T}(
        convert(Matrix{T},       nz),
        convert(Vector{Float64}, z),
        convert(Matrix{T},       nz_norm),
        convert(Matrix{T},       bias),
        convert(T,               f_NL),
        convert(Float64,         p),
        convert(Vector{T},       nz_norms),
        convert(Matrix{T},       Kernel),
        convert(Vector{Float64}, ell_prefactor),
        convert(Vector{Float64}, limber_factor),
    )
end
has_nz(::PrimordialNonGaussianity) = true

"""
    GalaxyClustering <: AbstractCosmologicalProbes

Galaxy-clustering probe container.

It combines the density term `δ` with optional redshift-space distortions
`RSD`, magnification bias `μ`, and primordial non-Gaussianity `PNG`.
"""
@kwdef mutable struct GalaxyClustering <: AbstractCosmologicalProbes
    δ::NumberCounts
    RSD::Union{Nothing, RedshiftSpaceDistortions} = nothing
    μ::Union{Nothing, MagnificationBias} = nothing
    PNG::Union{Nothing, PrimordialNonGaussianity} = nothing
end

"""
    WeakLensing <: AbstractCosmologicalProbes

Weak-lensing probe container combining shear `γ` and optional intrinsic
alignment `IA`.
"""
@kwdef mutable struct WeakLensing <: AbstractCosmologicalProbes
    γ::CosmicShear
    IA::Union{Nothing, IntrinsicAlignment} = nothing
end

"""
    CMB <: AbstractCosmologicalProbes

CMB probe container combining lensing convergence `κ` and optional
Integrated Sachs-Wolfe contribution `ISW`.
"""
@kwdef mutable struct CMB <: AbstractCosmologicalProbes
    κ::CMBLensing
    ISW::Union{Nothing, IntegratedSachsWolfe} = nothing
end

@doc raw"""
    compute_kernel!(Component::NumberCounts, bg::Background)

Compute the galaxy number-counts kernel on the background redshift grid:

```math
K_i^{\delta}(z) = \frac{H(z)}{c}\, b_i(z)\, \hat n_i(z).
```
"""
function compute_kernel!(Component::NumberCounts, bg::Background) 
    check_and_normalize!(Component, bg.z)
    Component.Kernel = _number_counts_kernel(Component.bias, bg.H, Component.nz_norm)
end

# ----------------------------------------------------------------------------
# NumberCounts kernel: pure broadcast.
#
#   K[b, i] = bias[b, i] · H[i] / c · nz_norm[b, i]
#
# Extracted so AD backends can differentiate the kernel math on its own.
# Pure broadcast → no custom rrule needed; ForwardDiff / Zygote / Mooncake
# all handle `.*`-style ops natively.
# ----------------------------------------------------------------------------
function _number_counts_kernel(bias::AbstractMatrix, H::AbstractVector,
                               nz_norm::AbstractMatrix)
    return @. bias * (H' / C_LIGHT) * nz_norm
end

@doc raw"""
    compute_kernel!(Component::CosmicShear, bg::Background)

Compute the weak-lensing shear kernel using Simpson integration on the $\chi$
grid:

```math
K_i^{\gamma}(\chi) = \frac{3 H_0^2 \Omega_m}{2 c^2}\, \chi (1+z)
\int_z^{\infty} \mathrm{d}z'\, \frac{H(z')}{c}\, \hat n_i(z')
\frac{\chi(z')-\chi}{\chi(z')}.
```
"""
function compute_kernel!(Component::CosmicShear, bg::Background)
    check_and_normalize!(Component, bg.z)

    H0 = get_H0(bg)
    Ωm = get_Ωm(bg)
    prefac = 1.5 * H0^2 * Ωm / C_LIGHT^2

    simpson_matrix = simpson_weights_matrix(length(bg.z))
    Δχ = (bg.χ[end] - bg.χ[1]) / (length(bg.χ) - 1)

    Component.Kernel = _cosmic_shear_kernel_tullio(
        bg.H, bg.χ, bg.z, Component.nz_norm, simpson_matrix, Δχ, prefac)
end

# ----------------------------------------------------------------------------
# CosmicShear kernel contraction
#
#   K[b, i] = Σ_p (Δχ · prefac / c) · H[p] · S[i, p] · χ[i] · (1 + z[i]) ·
#                   n̂[b, p] · (χ[p] - χ[i]) / χ[p]
#
# Extracted for AD: a hand-written rrule lives in `src/chainrules.jl` and is
# registered as a Mooncake primitive via `@from_chainrules` in MooncakeExt.
# Mooncake's auto-differentiation of this Tullio works as of today, but the
# custom rrule pins semantics, insulates against upstream Tullio/Mooncake
# version drift, and gives predictable adjoint performance.
# ----------------------------------------------------------------------------
function _cosmic_shear_kernel_tullio(H::AbstractVector, χ::AbstractVector,
                                     z::AbstractVector, nz_norm::AbstractMatrix,
                                     simpson_matrix::AbstractMatrix,
                                     Δχ::Number, prefac::Number)
    @tullio K[b, i] := Δχ * (H[p] / C_LIGHT) * simpson_matrix[i, p] * prefac *
                       χ[i] * (1.0 + z[i]) * nz_norm[b, p] *
                       (χ[p] - χ[i]) / χ[p]
    return K
end

@doc raw"""
    compute_kernel!(Component::CMBLensing, bg::Background)

Compute the CMB lensing convergence kernel on the background grid:

```math
K^{\kappa_{\mathrm{CMB}}}(\chi) = \frac{3 H_0^2 \Omega_m}{2 c^2}\, \chi (1+z)
\left(1 - \frac{\chi}{\chi_*}\right),
\qquad \chi_* = \chi(z_*),\ z_* \simeq 1090.
```
"""
function compute_kernel!(Component::CMBLensing, bg::Background) 
    H0 = get_H0(bg)
    Ωm = get_Ωm(bg)
    prefac = 1.5 * H0^2 * Ωm / C_LIGHT^2

    χ_CMB = compute_χ(1090, bg.cosmo, order=120)
    Component.Kernel = _cmb_lensing_kernel(bg.χ, bg.z, χ_CMB, prefac)
end

# ----------------------------------------------------------------------------
# CMBLensing kernel: pure broadcast, source plane at recombination.
#
#   K[1, i] = prefac · χ[i] · (1 + z[i]) · (1 - χ[i] / χ_CMB)
#
# `χ_CMB` is treated as a plain Number input — differentiating wrt it
# propagates the χ_CMB dependence, which tracks H0/Ωm sensitivity through
# compute_χ(1090, cosmo) when the caller hooks that in.
# Output is shaped (1, nχ) for consistency with multi-bin kernels.
# ----------------------------------------------------------------------------
function _cmb_lensing_kernel(χ::AbstractVector, z::AbstractVector,
                             χ_CMB::Number, prefac::Number)
    kernel = @. prefac * χ * (1 + z) * (1 - χ / χ_CMB)
    return reshape(kernel, 1, length(kernel))
end

@doc raw"""
    compute_kernel!(Component::RedshiftSpaceDistortions, bg::Background)

Compute the redshift-space-distortion kernel from the linear growth rate and
normalized redshift distribution:

```math
K_i^{\mathrm{RSD}}(z) = \frac{H(z)}{c}\, f(z)\, \hat n_i(z).
```
"""
function compute_kernel!(Component::RedshiftSpaceDistortions, bg::Background) 
    check_and_normalize!(Component, bg.z)
    Component.Kernel = _rsd_kernel(bg.f, bg.H, Component.nz_norm)
end

# ----------------------------------------------------------------------------
# RedshiftSpaceDistortions kernel: pure broadcast.
#
#   K[b, i] = f[i] · H[i] / c · nz_norm[b, i]
#
# Pure broadcast → AD-native; no custom rrule required.
# ----------------------------------------------------------------------------
function _rsd_kernel(f::AbstractVector, H::AbstractVector,
                     nz_norm::AbstractMatrix)
    return @. f' * (H' / C_LIGHT) * nz_norm
end

@doc raw"""
    compute_kernel!(Component::MagnificationBias, bg::Background)

Compute the magnification-bias kernel with source-slope factor $5s-2$:

```math
K_i^{\mu}(\chi) = \frac{3 H_0^2 \Omega_m}{2 c^2}\, \chi (1+z)
\int_z^{\infty} \mathrm{d}z'\, \frac{H(z')}{c}\, \hat n_i(z')
\frac{\chi(z')-\chi}{\chi(z')}\, \left[5 s_i(z') - 2\right].
```
"""
function compute_kernel!(Component::MagnificationBias, bg::Background)
    check_and_normalize!(Component, bg.z)

    H0 = get_H0(bg)
    Ωm = get_Ωm(bg)
    prefac = 1.5 * H0^2 * Ωm / C_LIGHT^2

    simpson_matrix = simpson_weights_matrix(length(bg.z))
    Δχ = (bg.χ[end] - bg.χ[1]) / (length(bg.χ) - 1)

    Component.Kernel = _magnification_bias_kernel_tullio(
        bg.H, bg.χ, bg.z, Component.nz_norm, Component.s,
        simpson_matrix, Δχ, prefac)
end

# ----------------------------------------------------------------------------
# MagnificationBias kernel contraction (CosmicShear + (5·s - 2) factor).
#
#   K[b, i] = Σ_p (Δχ · prefac / c) · H[p] · S[i, p] · χ[i] · (1 + z[i]) ·
#                   n̂[b, p] · (χ[p] - χ[i]) / χ[p] · (5·s[b, p] - 2)
#
# Like `_cosmic_shear_kernel_tullio`, extracted so a hand-written rrule can
# be registered as a Mooncake primitive (see `src/chainrules.jl`).
# ----------------------------------------------------------------------------
function _magnification_bias_kernel_tullio(H::AbstractVector, χ::AbstractVector,
                                           z::AbstractVector, nz_norm::AbstractMatrix,
                                           s::AbstractMatrix,
                                           simpson_matrix::AbstractMatrix,
                                           Δχ::Number, prefac::Number)
    @tullio K[b, i] := Δχ * (H[p] / C_LIGHT) * simpson_matrix[i, p] * prefac *
                       χ[i] * (1.0 + z[i]) * nz_norm[b, p] *
                       (χ[p] - χ[i]) / χ[p] * (5.0 * s[b, p] - 2)
    return K
end

@doc raw"""
    compute_kernel!(Component::IntrinsicAlignment, bg::Background)

Compute the intrinsic-alignment kernel using the NLA amplitude model:

```math
K_i^{\mathrm{IA}}(z) = \frac{H(z)}{c}\, A_{\mathrm{IA},i}(z)\, \hat n_i(z),
\qquad
A_{\mathrm{IA}}(z) = - A\, C_1\, \frac{\Omega_m}{D(z)}.
```
"""
function compute_kernel!(Component::IntrinsicAlignment, bg::Background) 
    check_and_normalize!(Component, bg.z)

    # Dispatch: user-supplied A_IA matrix (shape matches nz_norm × bg.z) vs
    # NLA-model branch. The NLA branch composes A_IA construction and the
    # kernel multiplication into one pure function so AD can flow wrt
    # (A, C1, Ωm, bg.D, bg.H, nz_norm) without mutating Component.A_IA.
    user_A_IA = size(Component.A_IA) == (size(Component.nz_norm, 1), length(bg.z))

    if user_A_IA
        Component.Kernel = _ia_kernel(Component.A_IA, bg.H, Component.nz_norm)
    else
        Component.Kernel = _ia_kernel_nla(Component.A, Component.C1,
                                          get_Ωm(bg), bg.D,
                                          bg.H, Component.nz_norm)
    end
end

# ----------------------------------------------------------------------------
# IntrinsicAlignment kernel: pure broadcast.
#
#   K[b, i] = A_IA[b, i] · H[i] / c · nz_norm[b, i]
#
# Used directly when the user provides a full A_IA matrix.
# ----------------------------------------------------------------------------
function _ia_kernel(A_IA::AbstractMatrix, H::AbstractVector,
                    nz_norm::AbstractMatrix)
    return @. A_IA * (H' / C_LIGHT) * nz_norm
end

# ----------------------------------------------------------------------------
# IntrinsicAlignment kernel with inlined NLA amplitude model.
#
#   A_NLA(z) = -A · C1 · Ωm / D(z)                     (vector, length n_z)
#   K[b, i]  = A_NLA[i] · H[i] / c · nz_norm[b, i]
#
# Composes NLA_model's math and _ia_kernel into a single pure function so
# AD can differentiate wrt each of A, C1, Ωm, D, H, nz_norm independently
# without the A_IA struct-field mutation that the multi-step path required.
#
# The `repeat(nla_vals', n_bins, 1)` matrix materialization is skipped —
# `nla_vals'` broadcasts directly against the (n_bins, n_z) arrays, giving
# the same result with one fewer allocation.
# ----------------------------------------------------------------------------
function _ia_kernel_nla(A::Number, C1::Number, Ωm::Number,
                        D::AbstractVector, H::AbstractVector,
                        nz_norm::AbstractMatrix)
    nla_vals = @. -A * C1 * Ωm / D           # vector, length n_z
    return @. nla_vals' * (H' / C_LIGHT) * nz_norm
end

@doc raw"""
    compute_kernel!(Component::IntegratedSachsWolfe, bg::Background)

Compute the ISW kernel from the background expansion and growth history:

```math
K^{\mathrm{ISW}}(z) = \frac{3 T_{\mathrm{CMB}} H_0^2 \Omega_m}{c^3}\, H(z)\, \left[1-f(z)\right].
```
"""
function compute_kernel!(Component::IntegratedSachsWolfe, bg::Background) 
    H0 = get_H0(bg)
    Ωm = get_Ωm(bg)
    T_CMB = 2.726
    prefac = 3T_CMB * H0^2 * Ωm / C_LIGHT^3
    Component.Kernel = _isw_kernel(bg.H, bg.f, prefac)
end

# ----------------------------------------------------------------------------
# IntegratedSachsWolfe kernel: pure broadcast.
#
#   K[1, i] = prefac · H[i] · (1 - f[i])
#
# Output is shaped (1, nχ) for pipeline consistency.
# ----------------------------------------------------------------------------
function _isw_kernel(H::AbstractVector, f::AbstractVector, prefac::Number)
    kernel = @. prefac * H * (1 - f)
    return reshape(kernel, 1, length(kernel))
end

@doc raw"""
    compute_kernel!(Component::PrimordialNonGaussianity, bg::Background)

Compute the primordial non-Gaussianity scale-dependent-bias kernel:

```math
K_i^{\mathrm{PNG}}(z) = \frac{H(z)}{c}\, f_{\mathrm{NL}}\, b_{\Phi,i}(z)\, \hat n_i(z),
\qquad
b_{\Phi,i}(z) = 2\,\delta_c\,\left[b_i(z)-p\right].
```
"""
function compute_kernel!(Component::PrimordialNonGaussianity, bg::Background) 
    check_and_normalize!(Component, bg.z)
    Component.Kernel = _png_kernel(bg.H, Component.f_NL, Component.bias,
                                    Component.p, Component.nz_norm)
end

# ----------------------------------------------------------------------------
# PrimordialNonGaussianity kernel: pure broadcast with bΦ helper inlined.
#
#   b_phi[b, i] = 2 · δ_c · (bias[b, i] - p)
#   K[b, i]     = H[i] / c · f_NL · b_phi[b, i] · nz_norm[b, i]
#
# Takes bias and p directly (instead of a precomputed b_phi) so AD can flow
# wrt either bias, f_NL, or p in isolation.
# ----------------------------------------------------------------------------
function _png_kernel(H::AbstractVector, f_NL::Number,
                     bias::AbstractMatrix, p::Number,
                     nz_norm::AbstractMatrix)
    b_phi = bΦ(bias, p)
    return @. (H' / C_LIGHT) * f_NL * b_phi * nz_norm
end

function compute_kernel!(Component::Nothing, bg::Background) 
    return nothing
end

"""
    evaluate_components!(probe, bg)

Compute and store all kernels for the components of `probe` on the `Background` grid.
"""
# ────────────────────────────────────────────────────────────────────────────
# _promote_eltype (Phase B)
#
# When `bg::Background{U}` carries a wider eltype than a component was
# constructed with (U = Dual from ForwardDiff-ing cosmology parameters,
# component T = Float64), `compute_kernel!` would try to write `Matrix{Dual}`
# into `Matrix{Float64}` — silently stripping the Duals. This helper widens
# the component's T parameter by constructing a fresh struct with `convert`-ed
# fields *before* `compute_kernel!` runs.
#
# When T ≥ U (the common case — Float64 bg, Float64 components), returns the
# original struct unchanged (compile-time short-circuit).
# ────────────────────────────────────────────────────────────────────────────
_promote_eltype(::Nothing, ::Type) = nothing

# Table driving the 8 `_promote_eltype` methods below. Each entry gives a
# component's struct name and its field order, annotated by how each field
# should be handled when constructing a widened-eltype copy:
#
#   :T_mat    → convert(Matrix{W}, C.field)   (T-parametric matrix)
#   :T_vec    → convert(Vector{W}, C.field)   (T-parametric vector)
#   :T_scalar → convert(W,         C.field)   (T-parametric scalar)
#   :keep     → C.field                       (pass-through: Float64 arrays,
#                                              fixed scalars like C1 and p)
#
# The @eval loop below expands one method per row. Field order must match
# the positional constructor of each struct — that is the spec.
const _PROMOTE_COMPONENTS = (
    (:NumberCounts, (
        (:nz, :T_mat), (:z, :keep), (:nz_norm, :T_mat), (:bias, :T_mat),
        (:nz_norms, :T_vec),
        (:Kernel, :T_mat), (:ell_prefactor, :keep), (:limber_factor, :keep),
    )),
    (:CosmicShear, (
        (:nz, :T_mat), (:z, :keep), (:nz_norm, :T_mat),
        (:nz_norms, :T_vec),
        (:Kernel, :T_mat),
        (:ell_prefactor, :keep), (:limber_factor, :keep),
    )),
    (:CMBLensing, (
        (:Kernel, :T_mat), (:ell_prefactor, :keep), (:limber_factor, :keep),
    )),
    (:RedshiftSpaceDistortions, (
        (:nz, :T_mat), (:z, :keep), (:nz_norm, :T_mat), (:growth_rate, :T_vec),
        (:nz_norms, :T_vec),
        (:Kernel, :T_mat), (:ell_prefactor, :keep), (:limber_factor, :keep),
    )),
    (:MagnificationBias, (
        (:nz, :T_mat), (:z, :keep), (:nz_norm, :T_mat), (:s, :T_mat),
        (:nz_norms, :T_vec),
        (:Kernel, :T_mat), (:ell_prefactor, :keep), (:limber_factor, :keep),
    )),
    (:IntrinsicAlignment, (
        (:nz, :T_mat), (:z, :keep), (:nz_norm, :T_mat), (:A, :T_scalar),
        (:C1, :keep), (:A_IA, :T_mat),
        (:nz_norms, :T_vec),
        (:Kernel, :T_mat),
        (:ell_prefactor, :keep), (:limber_factor, :keep),
    )),
    (:IntegratedSachsWolfe, (
        (:growth_rate, :T_vec), (:Kernel, :T_mat),
        (:ell_prefactor, :keep), (:limber_factor, :keep),
    )),
    (:PrimordialNonGaussianity, (
        (:nz, :T_mat), (:z, :keep), (:nz_norm, :T_mat), (:bias, :T_mat),
        (:f_NL, :T_scalar), (:p, :keep),
        (:nz_norms, :T_vec),
        (:Kernel, :T_mat),
        (:ell_prefactor, :keep), (:limber_factor, :keep),
    )),
)

for (_struct, _fields) in _PROMOTE_COMPONENTS
    _args = [_kind === :T_mat    ? :(convert(Matrix{W}, C.$_field)) :
             _kind === :T_vec    ? :(convert(Vector{W}, C.$_field)) :
             _kind === :T_scalar ? :(convert(W,         C.$_field)) :
                                    :(C.$_field)
             for (_field, _kind) in _fields]
    @eval function _promote_eltype(C::$_struct{T}, ::Type{U}) where {T, U}
        W = promote_type(T, U)
        W == T && return C
        $_struct{W}($(_args...))
    end
end

function evaluate_components!(GC::GalaxyClustering, bg::Background{U}) where {U}
    GC.δ   = _promote_eltype(GC.δ, U)
    GC.RSD = _promote_eltype(GC.RSD, U)
    GC.μ   = _promote_eltype(GC.μ, U)
    GC.PNG = _promote_eltype(GC.PNG, U)
    compute_kernel!(GC.δ, bg)
    compute_kernel!(GC.RSD, bg)
    compute_kernel!(GC.μ, bg)
    compute_kernel!(GC.PNG, bg)
end

function evaluate_components!(WL::WeakLensing, bg::Background{U}) where {U}
    WL.γ  = _promote_eltype(WL.γ,  U)
    WL.IA = _promote_eltype(WL.IA, U)
    compute_kernel!(WL.γ, bg)
    compute_kernel!(WL.IA, bg)
end

function evaluate_components!(cmb::CMB, bg::Background{U}) where {U}
    cmb.κ   = _promote_eltype(cmb.κ,   U)
    cmb.ISW = _promote_eltype(cmb.ISW, U)
    compute_kernel!(cmb.κ, bg)
    compute_kernel!(cmb.ISW, bg)
end

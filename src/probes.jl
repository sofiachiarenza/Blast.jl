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
    prepare_nz_matrix(nz::AbstractMatrix, z::AbstractVector, z_grid::AbstractVector)

Interpolate and normalize an n(z) matrix bin-by-bin onto the calculation grid.
"""
function prepare_nz_matrix(nz::AbstractMatrix, z::AbstractVector, z_grid::AbstractVector)
    n_bins = size(nz, 1)
    z_lo, z_hi = first(z_grid), last(z_grid)

    # Per-bin normalization via Integrals.jl + QuadGKJL. SciMLSensitivity
    # supplies the rrule, unlike raw QuadGK.quadgk whose adaptive step
    # selection is primal-driven and gives wrong Dual-mode gradients.
    # Parameter `p` must be a plain Array for clean pullbacks; fixed knots
    # `z` are captured from the closure.
    # sensealg pins the Mooncake-native VJP path (IntegralsMooncakeExt),
    # avoiding the default ZygoteVJP which assumes scalar `p` and calls
    # `only(Δ)` on a vector tangent.
    sensealg = Integrals.ReCallVJP{Integrals.MooncakeVJP}(Integrals.MooncakeVJP())
    norms = map(1:n_bins) do b
        nz_row = nz[b, :]
        prob = IntegralProblem((x, p) -> akima_interpolation(p, z, x),
                               (z_lo, z_hi), nz_row)
        nz_norm_val = solve(prob, QuadGKJL(); sensealg=sensealg).u
        nz_norm_val > 0 || throw(ArgumentError(
            "n(z) for bin $b integrates to zero or negative ($nz_norm_val); " *
            "check that your redshift distribution is non-negative and overlaps z_grid."))
        nz_norm_val
    end

    # Vectorized Akima interpolation: ACE's matrix dispatch expects
    # (length(z), n_bins) layout and returns (length(z_grid), n_bins).
    interpolated = akima_interpolation(permutedims(nz), z, z_grid)
    # Normalize each bin, then return to (n_bins, length(z_grid)) layout.
    return permutedims(interpolated ./ reshape(norms, 1, :))
end

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
            Component.nz_norm = prepare_nz_matrix(Component.nz, Component.z, z_grid)
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
    Ωm = get_Ωm(bg.cosmo)
    return @. - A * C1 * Ωm / bg.D
end

@doc raw"""
    NumberCounts <: AbstractComponents

Galaxy number-counts density component.

This component stores the redshift distribution and linear bias entering

```math
K_i^{\delta}(z) = \frac{H(z)}{c}\, b_i(z)\, \hat n_i(z).
```
"""
@kwdef mutable struct NumberCounts <: AbstractComponents
    nz::Matrix{Float64} = zeros(1, 1)
    z::Vector{Float64} = zeros(1)
    nz_norm::Matrix{Float64} = zeros(1, 1)
    bias::Matrix{Float64} = zeros(1, 1)
    Kernel::Matrix{Float64} = zeros(1, 1)
    ell_prefactor::Vector{Float64} = ones(size(Blast.full_ℓ_range, 1))
    limber_factor::Vector{Float64} = ones(size(Blast.full_ℓ_range, 1))
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
@kwdef mutable struct CosmicShear <: AbstractComponents
    nz::Matrix{Float64} = zeros(1, 1)
    z::Vector{Float64} = zeros(1)
    nz_norm::Matrix{Float64} = zeros(1, 1)
    Kernel::Matrix{Float64} = zeros(1, 1)
    ell_prefactor::Vector{Float64} = @. sqrt(factorial_frac(Blast.full_ℓ_range))
    limber_factor::Vector{Float64} = (Blast.full_ℓ_range .+ 0.5) .^ (-2)
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
@kwdef mutable struct CMBLensing <: AbstractComponents
    Kernel::Matrix{Float64} = zeros(1, 1)
    ell_prefactor::Vector{Float64} = @. Blast.full_ℓ_range * (Blast.full_ℓ_range + 1)
    limber_factor::Vector{Float64} = (Blast.full_ℓ_range .+ 0.5) .^ (-2)
end

@doc raw"""
    RedshiftSpaceDistortions <: AbstractComponents

Redshift-space distortion component with kernel

```math
K_i^{\mathrm{RSD}}(z) = \frac{H(z)}{c}\, f(z)\, \hat n_i(z).
```
"""
@kwdef mutable struct RedshiftSpaceDistortions <: AbstractComponents
    nz::Matrix{Float64} = zeros(1, 1)
    z::Vector{Float64} = zeros(1)
    nz_norm::Matrix{Float64} = zeros(1, 1)
    growth_rate::Vector{Float64} = zeros(1)
    Kernel::Matrix{Float64} = zeros(1, 1)
    ell_prefactor::Vector{Float64} = ones(size(Blast.full_ℓ_range, 1))
    limber_factor::Vector{Float64} = ones(size(Blast.full_ℓ_range, 1))
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
@kwdef mutable struct MagnificationBias <: AbstractComponents
    nz::Matrix{Float64} = zeros(1, 1)
    z::Vector{Float64} = zeros(1)
    nz_norm::Matrix{Float64} = zeros(1, 1)
    s::Matrix{Float64} = zeros(1, 1)
    Kernel::Matrix{Float64} = zeros(1, 1)
    ell_prefactor::Vector{Float64} = @. Blast.full_ℓ_range * (Blast.full_ℓ_range + 1)
    limber_factor::Vector{Float64} = (Blast.full_ℓ_range .+ 0.5) .^ (-2)
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
@kwdef mutable struct IntrinsicAlignment <: AbstractComponents
    nz::Matrix{Float64} = zeros(1, 1)
    z::Vector{Float64} = zeros(1)
    nz_norm::Matrix{Float64} = zeros(1, 1)
    A::Float64 = 1.72     # Standard NLA amplitude
    C1::Float64 = 0.0134  # NLA coupling constant; almost always fixed
    A_IA::Matrix{Float64} = zeros(1, 1)
    Kernel::Matrix{Float64} = zeros(1, 1)
    ell_prefactor::Vector{Float64} = @. sqrt(factorial_frac(Blast.full_ℓ_range))
    limber_factor::Vector{Float64} = (Blast.full_ℓ_range .+ 0.5) .^ (-2)
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
@kwdef mutable struct IntegratedSachsWolfe <: AbstractComponents
    growth_rate::Vector{Float64} = zeros(1)
    Kernel::Matrix{Float64} = zeros(1, 1)
    ell_prefactor::Vector{Float64} = ones(size(Blast.full_ℓ_range, 1))
    limber_factor::Vector{Float64} = ones(size(Blast.full_ℓ_range, 1))
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
@kwdef mutable struct PrimordialNonGaussianity <: AbstractComponents
    nz::Matrix{Float64} = zeros(1, 1)
    z::Vector{Float64} = zeros(1)
    nz_norm::Matrix{Float64} = zeros(1, 1)
    bias::Matrix{Float64} = zeros(1, 1)
    f_NL::Float64 = 0.0
    p::Float64 = 0.0
    Kernel::Matrix{Float64} = zeros(1, 1)
    ell_prefactor::Vector{Float64} = ones(size(Blast.full_ℓ_range, 1))
    limber_factor::Vector{Float64} = ones(size(Blast.full_ℓ_range, 1))
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
    RSD::Union{RedshiftSpaceDistortions, Nothing} = nothing
    μ::Union{MagnificationBias, Nothing} = nothing
    PNG::Union{PrimordialNonGaussianity, Nothing} = nothing
end

"""
    WeakLensing <: AbstractCosmologicalProbes

Weak-lensing probe container combining shear `γ` and optional intrinsic
alignment `IA`.
"""
@kwdef mutable struct WeakLensing <: AbstractCosmologicalProbes
    γ::CosmicShear
    IA::Union{IntrinsicAlignment, Nothing} = nothing
end

"""
    CMB <: AbstractCosmologicalProbes

CMB probe container combining lensing convergence `κ` and optional
Integrated Sachs-Wolfe contribution `ISW`.
"""
@kwdef mutable struct CMB <: AbstractCosmologicalProbes
    κ::CMBLensing
    ISW::Union{IntegratedSachsWolfe, Nothing} = nothing
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

    H0 = get_H0(bg.cosmo)
    Ωm = get_Ωm(bg.cosmo)
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
    H0 = get_H0(bg.cosmo)
    Ωm = get_Ωm(bg.cosmo)
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

    H0 = get_H0(bg.cosmo)
    Ωm = get_Ωm(bg.cosmo)
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
                                          get_Ωm(bg.cosmo), bg.D,
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
    H0 = get_H0(bg.cosmo)
    Ωm = get_Ωm(bg.cosmo)
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
function evaluate_components!(GC::GalaxyClustering, bg::Background)
    compute_kernel!(GC.δ, bg)
    compute_kernel!(GC.RSD, bg)
    compute_kernel!(GC.μ, bg)
    compute_kernel!(GC.PNG, bg)
end

function evaluate_components!(WL::WeakLensing, bg::Background)
    compute_kernel!(WL.γ, bg)
    compute_kernel!(WL.IA, bg)
end

function evaluate_components!(cmb::CMB, bg::Background)
    compute_kernel!(cmb.κ, bg)
    compute_kernel!(cmb.ISW, bg)
end



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

function compute_w!(w::Nothing, c::PowerSpectrum)
    return nothing
end

@kwdef mutable struct w_2_00_ϕTT <: ProjectedMatterDensityComponent
    T̃ = T_tildes.T_2_00
    w::AbstractArray{<:Any, 3} = zeros(1, 1, 1)
end

function compute_w!(w::w_2_00_ϕTT, c::PowerSpectrum)
    w.w = w_ell_tullio(c.cϕTT.coefs, w.T̃)
end

@kwdef mutable struct w_minus2_00_ϕTT <: ProjectedMatterDensityComponent
    T̃ = T_tildes.T_minus2_00
    w::AbstractArray{<:Any, 3} = zeros(1, 1, 1)
end

function compute_w!(w::w_minus2_00_ϕTT, c::PowerSpectrum)
    w.w = w_ell_tullio(c.cϕTT.coefs, w.T̃)
end

@kwdef mutable struct w_0_00_ϕTT <: ProjectedMatterDensityComponent
    T̃ = T_tildes.T_0_00
    w::AbstractArray{<:Any, 3} = zeros(1, 1, 1)
end

function compute_w!(w::w_0_00_ϕTT, c::PowerSpectrum)
    w.w = w_ell_tullio(c.cϕTT.coefs, w.T̃)
end

@kwdef mutable struct w_0_02_ϕTT <: ProjectedMatterDensityComponent
    T̃ = T_tildes.T_0_02
    w::AbstractArray{<:Any, 3} = zeros(1, 1, 1)
end

function compute_w!(w::w_0_02_ϕTT, c::PowerSpectrum)
    w.w = w_ell_tullio(c.cϕTT.coefs, w.T̃)
end

@kwdef mutable struct w_0_20_ϕTT <: ProjectedMatterDensityComponent
    T̃ = T_tildes.T_0_20
    w::AbstractArray{<:Any, 3} = zeros(1, 1, 1)
end

function compute_w!(w::w_0_20_ϕTT, c::PowerSpectrum)
    w.w = w_ell_tullio(c.cϕTT.coefs, w.T̃)
end

@kwdef mutable struct w_2_02_ϕTT <: ProjectedMatterDensityComponent
    T̃ = T_tildes.T_2_02
    w::AbstractArray{<:Any, 3} = zeros(1, 1, 1)
end

function compute_w!(w::w_2_02_ϕTT, c::PowerSpectrum)
    w.w = w_ell_tullio(c.cϕTT.coefs, w.T̃)
end

@kwdef mutable struct w_2_20_ϕTT <: ProjectedMatterDensityComponent
    T̃ = T_tildes.T_2_20
    w::AbstractArray{<:Any, 3} = zeros(1, 1, 1)
end

function compute_w!(w::w_2_20_ϕTT, c::PowerSpectrum)
    w.w = w_ell_tullio(c.cϕTT.coefs, w.T̃)
end

@kwdef mutable struct w_2_22_ϕTT <: ProjectedMatterDensityComponent
    T̃ = T_tildes.T_2_22
    w::AbstractArray{<:Any, 3} = zeros(1, 1, 1)
end

function compute_w!(w::w_2_22_ϕTT, c::PowerSpectrum)
    w.w = w_ell_tullio(c.cϕTT.coefs, w.T̃)
end

@kwdef mutable struct w_2_00_ϕT <: ProjectedMatterDensityComponent
    T̃ = T_tildes.T_2_00
    w::AbstractArray{<:Any, 3} = zeros(1, 1, 1)
end

function compute_w!(w::w_2_00_ϕT, c::PowerSpectrum)
    w.w = w_ell_tullio(c.cϕT.coefs, w.T̃)
end

@kwdef mutable struct w_2_00_ϕT_R1 <: ProjectedMatterDensityComponent
    T̃ = T_tildes.T_2_00
    w::AbstractArray{<:Any, 3} = zeros(1, 1, 1)
end

function compute_w!(w::w_2_00_ϕT_R1, c::PowerSpectrum)
    coefs_R1 = c.cϕT.coefs[:,:, end]
    w.w = w_ell_tullio(coefs_R1, w.T̃)
end

@kwdef mutable struct w_0_00_ϕT <: ProjectedMatterDensityComponent
    T̃ = T_tildes.T_0_00
    w::AbstractArray{<:Any, 3} = zeros(1, 1, 1)
end

function compute_w!(w::w_0_00_ϕT, c::PowerSpectrum)
    w.w = w_ell_tullio(c.cϕT.coefs, w.T̃)
end

@kwdef mutable struct w_0_00_ϕT_R1 <: ProjectedMatterDensityComponent
    T̃ = T_tildes.T_0_00
    w::AbstractArray{<:Any, 3} = zeros(1, 1, 1)
end

function compute_w!(w::w_0_00_ϕT_R1, c::PowerSpectrum)
    coefs_R1 = c.cϕT.coefs[:,:, end]
    w.w = w_ell_tullio(coefs_R1, w.T̃)
end

@kwdef mutable struct w_2_02_ϕT <: ProjectedMatterDensityComponent
    T̃ = T_tildes.T_2_02
    w::AbstractArray{<:Any, 3} = zeros(1, 1, 1)
end

function compute_w!(w::w_2_02_ϕT, c::PowerSpectrum)
    w.w = w_ell_tullio(c.cϕT.coefs, w.T̃)
end

@kwdef mutable struct w_2_02_ϕT_R1 <: ProjectedMatterDensityComponent
    T̃ = T_tildes.T_2_02
    w::AbstractArray{<:Any, 3} = zeros(1, 1, 1)
end

function compute_w!(w::w_2_02_ϕT_R1, c::PowerSpectrum)
    coefs_R1 = c.cϕT.coefs[:,:, end]
    w.w = w_ell_tullio(coefs_R1, w.T̃)
end

@kwdef mutable struct w_2_20_ϕT <: ProjectedMatterDensityComponent
    T̃ = T_tildes.T_2_20
    w::AbstractArray{<:Any, 3} = zeros(1, 1, 1)
end

function compute_w!(w::w_2_20_ϕT, c::PowerSpectrum)
    w.w = w_ell_tullio(c.cϕT.coefs, w.T̃)
end

@kwdef mutable struct w_2_20_ϕT_R1 <: ProjectedMatterDensityComponent
    T̃ = T_tildes.T_2_20
    w::AbstractArray{<:Any, 3} = zeros(1, 1, 1)
end

function compute_w!(w::w_2_20_ϕT_R1, c::PowerSpectrum)
    coefs_R1 = c.cϕT.coefs[:,:,end]
    w.w = w_ell_tullio(coefs_R1, w.T̃)
end

@kwdef mutable struct w_2_00_ϕ <: ProjectedMatterDensityComponent
    T̃ = T_tildes.T_2_00
    w::AbstractArray{<:Any, 3} = zeros(1, 1, 1)
end

function compute_w!(w::w_2_00_ϕ, c::PowerSpectrum)
    w.w = w_ell_tullio(c.cϕ.coefs, w.T̃)
end

"""
    ProjectedMatterDensity

Container holding all projected matter density components required for all the active observables and components.

Each field corresponds to a specific kernel contribution (e.g. `w_2_00_ϕTT`,
`w_0_02_ϕT_R1`). Fields may be set to `nothing` if the corresponding contribution
is not required.

The container is populated during setup and filled by calling `compute_w!`
with a `PowerSpectrum` object.
"""
@kwdef mutable struct ProjectedMatterDensity <: AbstractProjectedMatterDensity
    w_2_00_ϕTT::Union{w_2_00_ϕTT, Nothing} = nothing
    w_minus2_00_ϕTT::Union{w_minus2_00_ϕTT, Nothing} = nothing
    w_0_00_ϕTT::Union{w_0_00_ϕTT, Nothing} = nothing
    w_0_02_ϕTT::Union{w_0_02_ϕTT, Nothing} = nothing
    w_0_20_ϕTT::Union{w_0_20_ϕTT, Nothing} = nothing
    w_2_02_ϕTT::Union{w_2_02_ϕTT, Nothing} = nothing
    w_2_20_ϕTT::Union{w_2_20_ϕTT, Nothing} = nothing
    w_2_22_ϕTT::Union{w_2_22_ϕTT, Nothing} = nothing
    w_2_00_ϕT::Union{w_2_00_ϕT, Nothing} = nothing
    w_2_00_ϕT_R1::Union{w_2_00_ϕT_R1, Nothing} = nothing
    w_0_00_ϕT::Union{w_0_00_ϕT, Nothing} = nothing
    w_0_00_ϕT_R1::Union{w_0_00_ϕT_R1, Nothing} = nothing
    w_2_02_ϕT::Union{w_2_02_ϕT, Nothing} = nothing
    w_2_02_ϕT_R1::Union{w_2_02_ϕT_R1, Nothing} = nothing
    w_2_20_ϕT::Union{w_2_20_ϕT, Nothing} = nothing
    w_2_20_ϕT_R1::Union{w_2_20_ϕT_R1, Nothing} = nothing
    w_2_00_ϕ::Union{w_2_00_ϕ, Nothing} = nothing
end 

"""
    compute_w!(W::ProjectedMatterDensity, c::PowerSpectrum)

Compute all active projected matter density components.
"""
function compute_w!(W::ProjectedMatterDensity, c::PowerSpectrum)
    compute_w!(W.w_2_00_ϕTT, c)
    compute_w!(W.w_minus2_00_ϕTT, c)
    compute_w!(W.w_0_00_ϕTT, c)
    compute_w!(W.w_0_02_ϕTT, c)
    compute_w!(W.w_0_20_ϕTT, c)
    compute_w!(W.w_2_02_ϕTT, c)
    compute_w!(W.w_2_20_ϕTT, c)
    compute_w!(W.w_2_22_ϕTT, c)
    compute_w!(W.w_2_00_ϕT, c)
    compute_w!(W.w_2_00_ϕT_R1, c)
    compute_w!(W.w_0_00_ϕT, c)
    compute_w!(W.w_0_00_ϕT_R1, c)
    compute_w!(W.w_2_02_ϕT, c)
    compute_w!(W.w_2_02_ϕT_R1, c)
    compute_w!(W.w_2_20_ϕT, c)
    compute_w!(W.w_2_20_ϕT_R1, c)
    compute_w!(W.w_2_00_ϕ, c)
end
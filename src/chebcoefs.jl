# FUNCTIONS FOR THE DECOMPOSITION OF THE POWER SPECTRUM ON THE BASIS OF THE CHEBYSHEV POLYNOMIALS
using FFTW
using ChainRulesCore
using ForwardDiff: Dual, value, partials, Partials, tagtype

"""
    ChebyshevPlan{ND, P, T}

A plan for computing the Chebyshev coefficients of a function evaluated at Chebyshev nodes.
"""
struct ChebyshevPlan{ND, P, T}
    fft_plan::P
    K::NTuple{ND, Int}
    nodes::NTuple{ND, Vector{T}}
    dim::NTuple{ND, Int}
end

"""
    chebpoints(K::Int, x_min::Number, x_max::Number)

Generate `K+1` Chebyshev nodes of the second kind (extrema) mapped to `[x_min, x_max]`.
This replaces `FastChebInterp.chebpoints` directly to preserve mathematical behavior natively.
"""
function chebpoints(K::Int, x_min::Number, x_max::Number)
    T = promote_type(typeof(x_min), typeof(x_max), Float64)
    k = 0:K
    # Cosine points in [-1, 1], descending from 1 to -1 exactly as FastChebInterp does
    nodes_std = cos.(π .* k ./ K)
    nodes = T(x_min) .+ T(0.5) .* (nodes_std .+ T(1.0)) .* T(x_max - x_min)
    return nodes
end

"""
    prepare_chebyshev_plan(x_mins, x_maxs, Ks; size_nd=nothing, dim=1)

Precomputes the Chebyshev nodes and the FFT plan required to compute coefficients.
K is the polynomial degree (K+1 nodes). For N-dimensional inputs, specify the `size_nd`
tuple and the target dimension `dim`.
"""
function prepare_chebyshev_plan(x_mins::NTuple{N, Any}, x_maxs::NTuple{N, Any}, Ks::NTuple{N, Int}; size_nd::Union{Tuple, Nothing}=nothing, dim::NTuple{N, Int}=ntuple(i->i, N)) where {N}
    T = promote_type(typeof.(x_mins)..., typeof.(x_maxs)..., Float64)
    nodes = ntuple(i -> chebpoints(Ks[i], x_mins[i], x_maxs[i]), N)

    if size_nd === nothing
        # If N > 1, we can only infer size_nd if dim is standard (1, 2, ..., N)
        if N > 1
            @assert dim == ntuple(i->i, N) "size_nd must be specified if dim is not (1, ..., N)"
        end
        size_nd = ntuple(i -> Ks[i] + 1, N)
    end

    for i in 1:N
        @assert size_nd[dim[i]] == Ks[i] + 1 "Size along target dimension $(dim[i]) must be Ks[$i]+1"
    end
    fft_plan = FFTW.plan_r2r(zeros(T, size_nd...), FFTW.REDFT00, dim; flags=FFTW.PATIENT, timelimit=Inf)

    return ChebyshevPlan(fft_plan, Ks, nodes, dim)
end

function prepare_chebyshev_plan(x_min::Number, x_max::Number, K::Int; size_nd::Union{Tuple, Nothing}=nothing, dim::Int=1)
    return prepare_chebyshev_plan((x_min,), (x_max,), (K,); size_nd=size_nd, dim=(dim,))
end

"""
    chebyshev_polynomials(x_grid, x_min, x_max, K)

Computes the matrix of Chebyshev polynomials evaluated on `x_grid`, mapped to `[-1, 1]` from `[x_min, x_max]`.
"""
function chebyshev_polynomials(x_grid::AbstractVector, x_min::Number, x_max::Number, K::Int)
    T = promote_type(eltype(x_grid), typeof(x_min), typeof(x_max), Float64)
    n = length(x_grid)
    map_to_domain(val) = T(2.0) * (T(val) - T(x_min)) / T(x_max - x_min) - T(1.0)
    z = map_to_domain.(x_grid)

    T_mat = zeros(T, n, K + 1)
    T_mat[:, 1] .= T(1.0)
    if K > 0
        T_mat[:, 2] .= z
    end
    for k in 3:(K + 1)
        T_mat[:, k] .= T(2.0) .* z .* T_mat[:, k-1] .- T_mat[:, k-2]
    end
    return T_mat
end

# Helper to apply scaling to raw FFT coefficients
function _scale_chebyshev_coeffs!(c, ND, plan_dim, plan_K)
    for i in 1:ND
        d = plan_dim[i]
        K_i = plan_K[i]
        n_d = size(c, d)
        
        # Scale endpoints by 1/(2K) and interiors by 1/K
        selectdim(c, d, 1) ./= (2 * K_i)
        selectdim(c, d, n_d) ./= (2 * K_i)
        if n_d > 2
            selectdim(c, d, 2:(n_d-1)) ./= K_i
        end
    end
    return c
end

"""
    chebyshev_decomposition(plan, f_vals)

Computes the Chebyshev coefficients for a function evaluated at the Chebyshev nodes.
Supports batched inputs (ranks higher than the plan rank) and ForwardDiff.Dual.
"""
function chebyshev_decomposition(plan::ChebyshevPlan{ND, P, T}, f_vals::AbstractArray) where {ND, P, T}
    # Rank for which the plan was created
    PR = length(size(plan.fft_plan))
    
    # Case 1: Rank matches exactly what the plan expects (spatial grid dimensions)
    if ndims(f_vals) == PR
        return _chebyshev_decomposition_single(plan, f_vals)
    end

    # Case 2: Batched input (Rank > PR)
    # Dimensions after PR are treated as batch dimensions.
    grid_size = size(f_vals)[1:PR]
    batch_size = size(f_vals)[PR+1:end]
    f_reshaped = reshape(f_vals, grid_size..., :)
    
    # First batch to get result size and type
    dummy_c = _chebyshev_decomposition_single(plan, copy(selectdim(f_reshaped, PR+1, 1)))
    c_reshaped = similar(f_reshaped, eltype(dummy_c), size(dummy_c)..., size(f_reshaped, PR+1))
    
    # Process each batch
    for i in 1:size(f_reshaped, PR+1)
        # We MUST use copy() to ensure contiguous memory for FFTW plan application
        f_slice = copy(selectdim(f_reshaped, PR+1, i))
        selectdim(c_reshaped, PR+1, i) .= _chebyshev_decomposition_single(plan, f_slice)
    end
    
    return reshape(c_reshaped, size(c_reshaped)[1:PR]..., batch_size...)
end

# Internal implementation for a single block (Rank == PR)
function _chebyshev_decomposition_single(plan::ChebyshevPlan{ND, P, T}, f_vals::AbstractArray{T}) where {ND, P, T}
    # f_vals rank must match plan rank (checked by caller)
    c = plan.fft_plan * f_vals
    return _scale_chebyshev_coeffs!(c, ND, plan.dim, plan.K)
end

function _chebyshev_decomposition_single(plan::ChebyshevPlan{ND, P, T}, f_vals::AbstractArray{<:Dual}) where {ND, P, T}
    vals = value.(f_vals)
    c_raw_val = plan.fft_plan * vals
    
    dual_type = eltype(f_vals)
    Tag = tagtype(dual_type)
    num_partials = length(partials(first(f_vals)))

    c_raw_partials = map(1:num_partials) do p
        p_vals = [partials(x)[p] for x in f_vals]
        plan.fft_plan * p_vals
    end

    c = map(CartesianIndices(c_raw_val)) do idx
        parts = Partials(ntuple(p -> c_raw_partials[p][idx], num_partials))
        Dual{Tag}(c_raw_val[idx], parts)
    end

    return _scale_chebyshev_coeffs!(c, ND, plan.dim, plan.K)
end

# AD rrule for Chebyshev decomposition
function ChainRulesCore.rrule(::typeof(chebyshev_decomposition), plan::ChebyshevPlan{ND, P, T}, f_vals::AbstractArray{T}) where {ND, P, T}
    c = chebyshev_decomposition(plan, f_vals)
    function chebyshev_decomposition_pullback(Δc_raw)
        Δf_vals = chebyshev_decomposition(plan, unthunk(Δc_raw))
        return NoTangent(), NoTangent(), Δf_vals
    end
    return c, chebyshev_decomposition_pullback
end

"""
    AbstractCoeff

Abstract supertype for Chebyshev-coefficient containers in BLAST.

Concrete subtypes of `AbstractCoeff` store Chebyshev expansion coefficients.

All coefficient types are constructed uniformly via and `build_coeff`.
"""
abstract type AbstractCoeff end

"""
    cϕTT <: AbstractCoeff

Container for Chebyshev coefficients of the unequal time power spectrum.

This type stores a three-dimensional array of the Chebyshev coefficients, where the decomposition is performed on the wakenumbers k.

Fields:
- `coefs::AbstractArray{<:Any,3}`: Chebyshev coefficients array of size (nk, nχ, nR).
"""
@kwdef mutable struct cϕTT <: AbstractCoeff
    coefs::AbstractArray{<:Any, 3} 
end

"""
    cϕT <: AbstractCoeff

Container for Chebyshev coefficients of the power spectrum built with only one transfer function (needed when PNGs are active).

This type stores a three-dimensional array of the Chebyshev coefficients, where the decomposition is performed on the wakenumbers k.

Fields:
- `coefs::AbstractArray{<:Any,3}`: Chebyshev coefficients array of size (nk, nχ, nR).
"""
@kwdef mutable struct cϕT <: AbstractCoeff
    coefs::AbstractArray{<:Any, 3}
end

"""
    cϕ <: AbstractCoeff

Container for Chebyshev coefficients of the primordial power spectrum (needed when PNGs are active).

This type stores a one-dimensional array of the Chebyshev coefficients, where the decomposition is performed on the wakenumbers k.

Fields:
- `coefs::AbstractArray{<:Any,3}`: Chebyshev coefficients array of size (nk).
"""
@kwdef mutable struct cϕ <: AbstractCoeff
    coefs::AbstractArray{<:Any, 1}
end

#Constructors
"""
    make_coeff(::Type{T}, c) where {T<:AbstractCoeff}

Construct a Chebyshev-coefficient container of type `T` from a coefficient array.

This is a lightweight constructor wrapper that standardizes the creation of
`AbstractCoeff` subtypes.

Arguments:
- `T<:AbstractCoeff`: The concrete coefficient type to construct (either cϕTT, cϕT, cϕ).
- `c`: Array of Chebyshev coefficients.

Returns:
- An instance of type `T` containing the provided coefficients.
"""
make_coeff(::Type{T}, c) where {T<:AbstractCoeff} = T(c)
"""
    chebinterp_native(c::AbstractVector, x_grid, x_min, x_max)

Evaluates a 1D Chebyshev expansion on a grid.
"""
function chebinterp_native(c::AbstractVector, x_grid, x_min, x_max)
    T = chebyshev_polynomials(x_grid, x_min, x_max, length(c) - 1)
    return T * c
end

"""
    build_coeff(::Type{T}, vals::AbstractArray, plan::Union{ChebyshevPlan,Nothing}) where {T<:AbstractCoeff}

Compute Chebyshev coefficients from input values and wrap them in a coefficient container.

This function applies an FFT-based Chebyshev transform to `vals` using the provided
`ChebyshevPlan`, then constructs a coefficient object of type `T`. If `plan` is `nothing`,
no computation is performed and `nothing` is returned.

Arguments:
- `T<:AbstractCoeff`: The concrete coefficient type to construct.
- `vals::AbstractArray`: Input array to be decomposed on the Chebyshev polynomials (i.e. the power spectrum).
- `plan::Union{ChebyshevPlan,Nothing}`: ChebyshevPlan used to compute coefficients.

Returns:
- An instance of type `T` containing the Chebyshev coefficients, or `nothing` if
  `plan` is `nothing` (i.e. if the coefficient is not active.).
"""
@inline function build_coeff(::Type{T}, vals::AbstractArray, plan::Union{ChebyshevPlan, Nothing}) where {T<:AbstractCoeff}
    c = chebyshev_decomposition(plan, vals)
    return make_coeff(T, c)
end

@inline function build_coeff(::Type{T}, vals::AbstractArray, plan::Nothing) where {T<:AbstractCoeff}
    return nothing
end

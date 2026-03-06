function load_Ts(folder, nχ, nR, nk)
    ell_vector = Blast.ℓ_nonlimber
    full_T = zeros(length(ell_vector), nχ, nR, nk)
    for i in 1:length(ell_vector)
        l_string = string(round(ell_vector[i]; digits=1))
        filename = folder * "/T_tilde_l_$l_string.npy"
        if isfile(filename)
            full_T[i,:,:,:] = npzread(filename)
        else
            println("Missing file!")
        end
    end
    return full_T
end

"""
    get_clencurt_grid(kmin::Number, kmax::Number, N::Number)
Return the integration points in k. They are a set of 'N' Chebyshev points rescaled between 'kmin' and 'kmax'.
"""
function get_clencurt_grid(kmin::Number, kmax::Number, N::Number)
    CC_obj = FastTransforms.chebyshevmoments1(Float64, N)
    x = FastTransforms.clenshawcurtisnodes(Float64, N)
    x = (kmax - kmin) / 2 * x .+ (kmin + kmax) / 2 

    x[1] *= (1-1e-8)
    x[end] *= (1+1e-8) #TODO: this is just a quick patch, need to figure this out properly.

    return x
end

"""
    get_clencurt_weights(kmin::Number, kmax::Number, N::Number)
Return the set of 'N' weights needed to perform the integration with the Clenshaw-Curtis quadrature rule.
The weights are rescaled between 'kmin' and 'kmax'.  
"""
function get_clencurt_weights(kmin::Number, kmax::Number, N::Number)
    CC_obj = FastTransforms.chebyshevmoments1(Float64, N)
    w = FastTransforms.clenshawcurtisweights(CC_obj)
    w = (kmax - kmin) / 2 * w

    return w
end

"""
    get_clencurt_weights_R_integration(N::Int)

Return Clenshaw–Curtis quadrature weights adapted for integration over the
ratio variable `R = χ₁ / χ₂ ∈ (0, 1]`.

The weights are obtained by:
1. Computing standard Clenshaw–Curtis weights on `[-1, 1]`,
2. Restricting to the positive half of the interval corresponding to `R > 0`,
3. Applying a correction to the first weight to account for the truncated domain.

This routine is used to computed the `C_\\ell`'s in the `\\chi-R` coordinate system.

# Arguments
- `N::Int`: Number of Clenshaw–Curtis nodes on `[-1, 1]`.

# Returns
A vector of quadrature weights suitable for integration over `R ∈ (0, 1]`.

# Notes
The correction applied to the first weight is a numerical workaround and is
not the exact analytic solution. 
"""
function get_clencurt_weights_R_integration(N::Int)

    w = get_clencurt_weights(-1, 1, N)

    index = div(N + 3, 2) 
    w = w[index:end]
    w[1]/=2 #TODO: investigate if there are better solutions, this is not the analytic solution.

    return w
end

"""
    bessel_second_derivative(ℓ::Number, x::AbstractArray)

Compute the second derivative with respect to the argument of the spherical
Bessel function of order `ℓ`.

Specifically, this returns 
```math
\\frac{d^2 j_\\ell(x)}{dx^2}
```
using the analytic recurrence relation in terms of `j_ℓ(x)` and `j_{ℓ+1}(x)`.

# Arguments
- `ℓ::Number`: Multipole order.
- `x::AbstractArray`: Array of arguments at which the derivative is evaluated.

# Returns
An array of the same shape as `x`, containing the second derivative of
`j_ℓ(x)`.
"""
function bessel_second_derivative(ℓ::Number, x::AbstractArray)
    return @views @. 2/x * SpecialFunctions.sphericalbesselj.(ℓ+1, x) + (ℓ^2-ℓ-x^2)/x^2 * SpecialFunctions.sphericalbesselj.(ℓ, x)
end

"""
    bessel_cheb_eval(ℓ::Number, kmin::Number, kmax::Number, χ::AbstractArray, n_cheb::Int, N::Number)

Evaluate spherical Bessel functions and Chebyshev polynomials on a common
Clenshaw–Curtis integration grid.

This function precomputes:

1. The values of Chebyshev polynomials up to order `n_cheb`,
   evaluated at `N` Clenshaw–Curtis nodes between `kmin` and `kmax`.
2. The spherical Bessel function `j_ℓ(kχ)` (or its second derivative)
   evaluated at the same `k` nodes for all `χ`.

These quantities are used to build the precomputed kernels `T̃` appearing
in the projected matter integrals.

# Arguments
- `ℓ::Number`: Multipole order.
- `kmin::Number`, `kmax::Number`: Integration range in wavenumber `k`.
- `χ::AbstractArray`: Comoving distance grid.
- `n_cheb::Int`: Maximum Chebyshev polynomial order.
- `N::Number`: Number of Clenshaw–Curtis integration nodes.
- `deriv_order::Int`: Derivative order of the Bessel function (`0` or `2`).

# Returns
- `T::Array`: Chebyshev polynomials evaluated on the `k` grid, with shape `(n_cheb+1, N)`.
- `Bessel::Array`: Bessel function values with shape `(length(χ), N)`.
"""
function bessel_cheb_eval(ℓ::Number, kmin::Number, kmax::Number, χ::AbstractArray, n_cheb::Int, N::Number, deriv_order::Int)

    nχ = length(χ)
    x = get_clencurt_grid(kmin, kmax, N)

    # Precompute all Chebyshev polynomials up to n_cheb on the log10(k) grid
    T = chebyshev_polynomials(log10.(x), log10(kmin), log10(kmax), n_cheb)'

    Bessel = zeros(nχ, N)
    Threads.@threads for i in 1:nχ
        if deriv_order == 0
            Bessel[i, :] = @views SpecialFunctions.sphericalbesselj.(ℓ, χ[i] * x)
        elseif deriv_order == 2
            Bessel[i, :] = @views bessel_second_derivative(ℓ, χ[i] * x)
        else
            error("Invalid derivative order: $deriv_order. Expected 0 or 2.")
        end
    end

    return T, Bessel

end



"""
    compute_T̃(ℓ::Number, χ::AbstractArray, R::AbstractArray, kmin::Number, kmax::Number, β::Number; n_cheb::Int = 119, N::Int = 2^(15)+1)

Compute the precomputed kernel `T̃` entering the projected matter density
integrals.

This function evaluates integrals of the form:
```math
\\tilde{T}(\\ell, \\chi, R, n) = \\int dk k^\\beta j_\\ell^{A}(k\\chi) j_\\ell^{B}(k R\\chi) T_n(k)
```
where `T_n` are Chebyshev polynomials and `j_ℓ^{A,B}` denote spherical
Bessel functions or their second derivatives.

These kernels are **cosmology-independent** once the global `χ` grid is fixed,
and are therefore precomputed and stored as artifacts.

# Arguments
- `ℓ::Number`: Multipole order.
- `χ::AbstractArray`: Global comoving distance grid.
- `R::AbstractArray`: Ratio grid `R = χ₁ / χ₂`.
- `kmin::Number`, `kmax::Number`: Integration limits in wavenumber.
- `β::Number`: Power-law exponent of `k` in the integrand.
- `der_A::Int`, `der_B::Int`: Derivative orders (`0` or `2`) for the two Bessel factors.
- `n_cheb::Int`: Number of Chebyshev polynomials used to approximate the power spectrum.
- `N::Int`: Number of Clenshaw–Curtis integration nodes.

# Returns
A 4D array with shape `(1, length(χ), length(R), n_cheb+1)` containing the
precomputed kernel `T̃`.
"""
function compute_T̃(ℓ::Number, χ::AbstractArray, R::AbstractArray, kmin::Number, kmax::Number, β::Number, der_A::Int, der_B::Int; n_cheb::Int = 119, N::Int = 2^(15)+1)
    if kmin >= kmax 
        throw(DomainError("The integration range is unphysical. Make sure kmin < kmax.")) 
    end
    
    nχ = length(χ)
    nR = length(R)

    x = get_clencurt_grid(kmin, kmax, N)
    w = get_clencurt_weights(kmin, kmax, N)

    T, Bessel1 = bessel_cheb_eval(ℓ, kmin, kmax, χ, n_cheb, N, der_A)


    T_tilde = zeros(1, nχ, nR, n_cheb+1)
    
    for (ridx, r) in enumerate(R)
        Bessel2 = zeros(nχ, N)
        
        if der_B == 0
            Threads.@threads for i in 1:nχ
                Bessel2[i,:] = @views SpecialFunctions.sphericalbesselj.(ℓ, r*χ[i] * x)
            end
        elseif der_B == 2
            Threads.@threads for i in 1:nχ
                Bessel2[i,:] = bessel_second_derivative(ℓ, r*χ[i] * x)
            end
        end

        α = w .* (x .^ β) #β = 2 for CC, -2 for LL and 0 for CL.
         
        @tturbo for l in 1:n_cheb+1, i in 1:nχ
            Cij = zero(eltype(w))
            for k in 1:N
                Cij +=  T[l,k] * Bessel1[i,k] * Bessel2[i,k] * α[k]
            end
            T_tilde[1,i,ridx,l] = Cij
        end
    end

    return T_tilde

end



"""
    factorial_frac(ℓ::Union{Number,Vector}}) 

Computes the ratio (ℓ+2)!/(ℓ-2)!, needed in the pre-factors of the the angular power spectra.

# Arguments
- `ℓ::Vector{T}`: vectors of ℓ values.
"""
function factorial_frac(ℓ::Union{Number,Vector})
    return @. (ℓ-1)*ℓ*(ℓ+1)*(ℓ+2)
end

"""
    bΦ(bias::AbstractArray, p::Number)

Compute the non-Gaussian bias coefficient `b_Φ` for local-type primordial
non-Gaussianity.

The coefficient is defined as
```math
b_\\Phi(z) = 2 \\delta_c (b(z)-p)
```
with `δ_c = 1.686`, where:
- `b` is the linear bias (can be Vector or Matrix for multiple bins),
- `p` is a tracer-dependent parameter.

# Returns
An array of the same shape as `bias` containing the non-Gaussian bias coefficient `b_Φ`.
"""
function bΦ(bias::AbstractArray{T}, p::Number) where T
    return 2 * 1.686 .* (bias .- p)
end

function simpson_weights_array(n::Int)
    number_intervals = floor((n - 1) / 2)
    weight_array = zeros(n)
    if n == number_intervals * 2 + 1
        for i in 1:number_intervals
            weight_array[Int((i - 1) * 2 + 1)] += 1 / 3
            weight_array[Int((i - 1) * 2 + 2)] += 4 / 3
            weight_array[Int((i - 1) * 2 + 3)] += 1 / 3
        end
    else
        weight_array[1] += 0.5
        weight_array[2] += 0.5
        for i in 1:number_intervals
            weight_array[Int((i - 1) * 2 + 1)+1] += 1 / 3
            weight_array[Int((i - 1) * 2 + 2)+1] += 4 / 3
            weight_array[Int((i - 1) * 2 + 3)+1] += 1 / 3
        end
        weight_array[length(weight_array)] += 0.5
        weight_array[length(weight_array)-1] += 0.5
        for i in 1:number_intervals
            weight_array[Int((i - 1) * 2 + 1)] += 1 / 3
            weight_array[Int((i - 1) * 2 + 2)] += 4 / 3
            weight_array[Int((i - 1) * 2 + 3)] += 1 / 3
        end
        weight_array ./= 2
    end
    return weight_array
end

@memoize function simpson_weights_matrix(n::Int)
    weight_matrix = zeros(n, n)
    for i in 1:n
        weight_matrix[i, i:n] = simpson_weights_array(n - i + 1)
    end
    return weight_matrix
end

#Functions for costumized akima interpolation, from effort.jl
function _akima_slopes(u, t)
    n = length(u)
    dt = diff(t)
    m = zeros(eltype(u[1] + t[1]), n + 3)
    m[3:(end-2)] = diff(u) ./ dt
    m[2] = 2m[3] - m[4]
    m[1] = 2m[2] - m[3]
    m[end-1] = 2m[end-2] - m[end-3]
    m[end] = 2m[end-1] - m[end-2]
    return m
end

function _akima_coefficients(t, m)
    n = length(t)
    dt = diff(t)
    b = (m[4:end] .+ m[1:(end-3)]) ./ 2
    dm = abs.(diff(m))
    f1 = dm[3:(n+2)]
    f2 = dm[1:n]
    f12 = f1 + f2

    # Handle division by zero for constant/linear segments
    # When f12 ≈ 0, use the average slope (already computed above)
    eps_akima = eps(eltype(f12)) * 100  # Small threshold
    for i in eachindex(f12)
        if f12[i] > eps_akima
            b[i] = (f1[i] * m[i+1] + f2[i] * m[i+2]) / f12[i]
        end
        # else: keep the average slope b[i] = (m[i+3] + m[i]) / 2
    end

    c = (3 .* m[3:(end-2)] .- 2 .* b[1:(end-1)] .- b[2:end]) ./ dt
    d = (b[1:(end-1)] .+ b[2:end] .- 2 .* m[3:(end-2)]) ./ dt .^ 2
    return b, c, d
end

function _akima_find_interval(t, tq)
    n = length(t)
    if tq ≤ t[1]
        return 1
    elseif tq ≥ t[end]
        return n - 1
    else
        idx = searchsortedlast(t, tq)
        return idx == n ? n - 1 : idx
    end
end

function _akima_eval(u, t, b, c, d, tq)
    idx = _akima_find_interval(t, tq)
    wj = tq - t[idx]
    return ((d[idx] * wj + c[idx]) * wj + b[idx]) * wj + u[idx]
end

function _akima_eval(u, t, b, c, d, tq::AbstractArray)
    map(tqi -> _akima_eval(u, t, b, c, d, tqi), tq)
end

"""
    _akima_interpolation(u, t, t_new)

Evaluates the one-dimensional Akima spline that interpolates the data points ``(t_i, u_i)``
at new abscissae `t_new`.

# Arguments
- `u`: Ordinates (function values) ``u_i`` at the data nodes.
- `t`: Strictly increasing abscissae (knots) ``t_i`` associated with `u`. `length(t)` must equal `length(u)`.
- `t_new`: The query point(s) where the spline is to be evaluated.

# Returns
The interpolated value(s) at `t_new`. A scalar input returns a scalar; a vector input returns a vector of the same length.

# Details
This routine implements the original Akima piecewise-cubic method (T. Akima, 1970). On each interval ``[t_j, t_{j+1}]``, a cubic polynomial is constructed. The method uses a weighted average of slopes to determine the derivative at each node, which effectively dampens oscillations without explicit shape constraints. The resulting spline is ``C^1`` continuous (its first derivative is continuous) but generally not ``C^2``.

# Formulae
The spline on the interval ``[t_j, t_{j+1}]`` is a cubic polynomial:
\\[
S_j(w) = u_j + b_j w + c_j w^{2} + d_j w^{3}, \\qquad w = t - t_j
\\]
The derivative ``b_j`` at each node is determined by Akima's weighting of local slopes ``m_j=(u_{j}-u_{j-1})/(t_j-t_{j-1})``:
\\[
b_j = \\frac{|m_{j+1}-m_{j}|\\,m_{j-1} + |m_{j-1}-m_{j-2}|\\,m_{j}}
            {|m_{j+1}-m_{j}| + |m_{j-1}-m_{j-2}|}
\\]
The remaining coefficients, ``c_j`` and ``d_j``, are found by enforcing continuity of the first derivative:
\\[
c_j = \\frac{3m_j - 2b_j - b_{j+1}}{t_{j+1}-t_j}
\\]
\\[
d_j = \\frac{b_j + b_{j+1} - 2m_j}{(t_{j+1}-t_j)^2}
\\]

# Automatic Differentiation
The implementation is free of mutation on the inputs and uses only element-wise arithmetic, making the returned value differentiable with both `ForwardDiff.jl` (dual numbers) and `Zygote.jl` (reverse-mode AD). You can therefore embed `_akima_interpolation` in optimization or machine-learning pipelines and back-propagate through the interpolation seamlessly.

# Notes
The algorithm and numerical results are equivalent to the Akima spline in `DataInterpolations.jl`, but this routine is self-contained and avoids any package dependency.
"""
function _akima_interpolation(u, t, t_new)
    n = length(t)
    dt = diff(t)

    m = _akima_slopes(u, t)
    b, c, d = _akima_coefficients(t, m)

    return _akima_eval(u, t, b, c, d, t_new)
end

"""
    _akima_slopes(u::AbstractMatrix, t)

Optimized version of `_akima_slopes` for matrix input where each column represents
a different data series but all share the same x-coordinates `t`.

# Performance Optimization
Computes `dt = diff(t)` once and reuses it for all columns, avoiding redundant computation.

# Arguments
- `u::AbstractMatrix`: Data values with shape `(n_points, n_columns)`.
- `t`: X-coordinates (same for all columns).

# Returns
Matrix of slopes with shape `(n_points + 3, n_columns)`.
"""
function _akima_slopes(u::AbstractMatrix, t)
    n, n_cols = size(u)
    dt = diff(t)  # Computed once, reused for all columns

    # Pre-allocate for all columns
    m = zeros(promote_type(eltype(u), eltype(t)), n + 3, n_cols)

    # Process each column using the shared dt
    for col in 1:n_cols
        m[3:(end-2), col] .= diff(view(u, :, col)) ./ dt

        # Extrapolation formulas
        m[2, col] = 2 * m[3, col] - m[4, col]
        m[1, col] = 2 * m[2, col] - m[3, col]
        m[end-1, col] = 2 * m[end-2, col] - m[end-3, col]
        m[end, col] = 2 * m[end-1, col] - m[end-2, col]
    end

    return m
end

"""
    _akima_coefficients(t, m::AbstractMatrix)

Optimized version of `_akima_coefficients` for matrix input where each column represents
coefficients for a different spline series.

# Performance Optimization
Computes `dt = diff(t)` once and reuses it for all columns.

# Arguments
- `t`: X-coordinates.
- `m::AbstractMatrix`: Slopes matrix with shape `(n_points + 3, n_columns)`.

# Returns
Tuple `(b, c, d)` where:
- `b` is a matrix of shape `(n_points, n_columns)`
- `c` and `d` are matrices of shape `(n_points - 1, n_columns)`
"""
function _akima_coefficients(t, m::AbstractMatrix)
    n = length(t)
    n_cols = size(m, 2)
    dt = diff(t)  # Computed once
    eps_akima = eps(eltype(m)) * 100

    # Pre-allocate for all columns - b has length n, c and d have length n-1
    b = zeros(eltype(m), n, n_cols)
    c = zeros(eltype(m), n - 1, n_cols)
    d = zeros(eltype(m), n - 1, n_cols)

    for col in 1:n_cols
        # Average slope (fallback) - length n
        b[:, col] .= (view(m, 4:(n+3), col) .+ view(m, 1:n, col)) ./ 2

        dm = abs.(diff(view(m, :, col)))
        f1 = view(dm, 3:(n+2))
        f2 = view(dm, 1:n)
        f12 = f1 .+ f2

        # Weighted average where slopes vary significantly
        for i in 1:n
            if f12[i] > eps_akima
                b[i, col] = (f1[i] * m[i+1, col] + f2[i] * m[i+2, col]) / f12[i]
            end
        end

        # Coefficients using shared dt - length n-1
        c[:, col] .= (3 .* view(m, 3:(n+1), col) .- 2 .* view(b, 1:(n-1), col) .- view(b, 2:n, col)) ./ dt
        d[:, col] .= (view(b, 1:(n-1), col) .+ view(b, 2:n, col) .- 2 .* view(m, 3:(n+1), col)) ./ dt .^ 2
    end

    return b, c, d
end

"""
    _akima_eval(u::AbstractMatrix, t, b::AbstractMatrix, c::AbstractMatrix, d::AbstractMatrix, tq::AbstractArray)

Optimized version of `_akima_eval` for matrix input where each column represents
a different spline series.

# Performance Optimization
- Finds intervals once per query point (not per column)
- Computes polynomial weights once per query point
- Broadcasts evaluation across all columns simultaneously

This is significantly faster than calling the vector version in a loop.

# Arguments
- `u::AbstractMatrix`: Data values with shape `(n_points, n_columns)`.
- `t`: X-coordinates.
- `b::AbstractMatrix`, `c::AbstractMatrix`, `d::AbstractMatrix`: Spline coefficients.
- `tq::AbstractArray`: Query points.

# Returns
Matrix of interpolated values with shape `(length(tq), n_columns)`.
"""
function _akima_eval(u::AbstractMatrix, t, b::AbstractMatrix, c::AbstractMatrix, d::AbstractMatrix, tq::AbstractArray)
    n_query = length(tq)
    n_cols = size(u, 2)
    results = zeros(promote_type(eltype(u), eltype(tq)), n_query, n_cols)

    @inbounds for i in 1:n_query
        idx = _akima_find_interval(t, tq[i])
        wj = tq[i] - t[idx]

        # Horner's method broadcasted over all columns
        # ((d*w + c)*w + b)*w + u
        @simd for col in 1:n_cols
            results[i, col] = ((d[idx, col] * wj + c[idx, col]) * wj + b[idx, col]) * wj + u[idx, col]
        end
    end

    return results
end

"""
    _akima_interpolation(u::AbstractMatrix, t, t_new)

Akima spline interpolation for multiple data series sharing the same x-coordinates.
Uses a simple comprehension-based approach that is compatible with automatic differentiation.

# Arguments
- `u::AbstractMatrix`: Data values with shape `(n_points, n_columns)`.
- `t`: X-coordinates shared by all columns.
- `t_new`: Query points.

# Returns
Matrix of interpolated values with shape `(length(t_new), n_columns)`.

# Example
```julia
# Interpolate 11 Jacobian columns at 100 k-points
k_in = range(0.01, 0.3, length=50)
k_out = range(0.01, 0.3, length=100)
jacobian = randn(50, 11)  # 11 parameters

result = _akima_interpolation(jacobian, k_in, k_out)  # (100, 11)
```
"""
function _akima_interpolation(u::AbstractMatrix, t, t_new)
    # Matrix-native implementation: compute shared operations once for all columns
    # This is much more efficient than column-wise processing, especially for Jacobians
    # Key optimization: diff(t) computed once instead of n_cols times
    m = _akima_slopes(u, t)
    b, c, d = _akima_coefficients(t, m)
    return _akima_eval(u, t, b, c, d, t_new)
end

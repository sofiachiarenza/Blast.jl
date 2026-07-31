function load_Ts(folder, nχ, nR, nk)
    ell_vector = Blast.ℓ_nonlimber
    full_T = zeros(length(ell_vector), nχ, nR, nk)
    for i in 1:length(ell_vector)
        l_string = string(round(ell_vector[i]; digits=1))
        filename = folder * "/T_tilde_l_$l_string.npy"
        if isfile(filename)
            full_T[i,:,:,:] = npzread(filename)
        else
            error("Missing T_tilde file: $filename. Cannot load T̃ kernels.")
        end
    end
    return full_T
end

"""
    get_clencurt_grid(kmin::Number, kmax::Number, N::Number)
Return the integration points in k. They are a set of 'N' Chebyshev points rescaled between 'kmin' and 'kmax'.
"""
function get_clencurt_grid(kmin::Number, kmax::Number, N::Number)
    x = FastTransforms.clenshawcurtisnodes(Float64, Int(N))
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
function get_clencurt_weights(kmin::Number, kmax::Number, N::Int)
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

    # Always build the quadrature rule at the full (untruncated) nR=64 grid
    # size, regardless of N (which reflects the possibly-truncated Blast.R
    # via _R_KEEP_IDX) — R-truncation drops nodes, it never re-derives a
    # lower-order quadrature rule. See constants.jl's R_TRUNCATION_FRAC.
    N_full = 2 * length(_R_full64) + 1
    w = get_clencurt_weights(-1, 1, N_full)

    index = div(N_full + 3, 2)
    w = w[index:end]
    w[1]/=2 #TODO: investigate if there are better solutions, this is not the analytic solution.

    return w[_R_KEEP_IDX]
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


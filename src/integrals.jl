"""
    SimpsonWeightArray(n::Int; T=Float64)

Computes the weights for the Simpson quadrature rule for numerical integration based on the input number of points `n`.

# Arguments
- `n::Int`: The number of points (must be at least 2).
- `T`: (optional) The type of the output array. Defaults to `Float64`.

# Returns
- An array of length `n` with the weights of type `T` for the Simpson quadrature rule.
"""
function SimpsonWeightArray(n::Int; T=Float64)
    @assert n > 1 "You cannot integrate with only 1 sampling point."
    number_intervals = floor((n-1)/2)
    weight_array = zeros(n)
    if n == number_intervals*2+1
        for i in 1:number_intervals
            weight_array[Int((i-1)*2+1)] += 1/3
            weight_array[Int((i-1)*2+2)] += 4/3
            weight_array[Int((i-1)*2+3)] += 1/3
        end
    else
        weight_array[1] += 0.5
        weight_array[2] += 0.5
        for i in 1:number_intervals
            weight_array[Int((i-1)*2+1)+1] += 1/3
            weight_array[Int((i-1)*2+2)+1] += 4/3
            weight_array[Int((i-1)*2+3)+1] += 1/3
        end
        weight_array[length(weight_array)]   += 0.5
        weight_array[length(weight_array)-1] += 0.5
        for i in 1:number_intervals
            weight_array[Int((i-1)*2+1)] += 1/3
            weight_array[Int((i-1)*2+2)] += 4/3
            weight_array[Int((i-1)*2+3)] += 1/3
        end
        weight_array ./= 2
    end
    return T.(weight_array)
end

"""
    compute_Cℓ(w::AbstractArray{T, 3}, K::AbstractArray{T, 4}, χ::AbstractVector, R::AbstractVector)

Computes the Cℓ's by performing the two outer integrals in χ and R. The integration in χ is performed using the Simpson quadrature rule, while the integration in R is performed using the Clenshaw-Curtis quadrature rule.

# Arguments
- `w::AbstractArray{T, 3}`: A 3D array representing the projected matter densities, i.e. the inner integrals in k. The three axis are (ℓ, χ, R).
- `K::AbstractArray{T, 4}`: A 4D array representing the kernel function, with dimensions (i,j,χ,R). i and j are the tomographic bins.
- `χ::AbstractVector`: A 1D array containing the χ values.
- `R::AbstractVector`: A 1D array containing the R values.

# Returns
- A multi-dimensional array `Cℓ` with axis (ℓ, i, j) containing the angular power spectrum coefficients in every combination of tomographic bins.
"""
function compute_Cℓ(w::AbstractArray{T, 3}, K::AbstractArray{T, 4}, χ::AbstractVector, R::AbstractVector) where T
    nχ = length(χ)
    nR = length(R)
    @assert nχ == size(w, 2) "Dimension mismatch: the χ array passed doesn't correspond to the one used in the evaluation of the inner k integral."

    #Integration in χ is peformed using the Simpson quadrature rule
    Δχ = ((last(χ)-first(χ))/(nχ-1))
    pesi_χ = SimpsonWeightArray(nχ)

    #Integration in R is performed using the Clenshaw-Curtis quadrature rule
    CC_obj = FastTransforms.chebyshevmoments1(Float64, 2*nR+1)
    pesi_R = FastTransforms.clenshawcurtisweights(CC_obj)
    pesi_R = pesi_R[nR+2:end]
    pesi_R[1]/=2 #TODO: investigate if there are better solutions

    @tullio Cℓ[l,i,j] := χ[n]*K[i,j,n,m]*w[l,n,m]*pesi_χ[n]*pesi_R[m]*Δχ

    return Cℓ

end
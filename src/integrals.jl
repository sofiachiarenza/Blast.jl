#TODO: missing documentation of these new functions 

function make_grid(BackgroundQuantities::BackgroundQuantities, R::Vector{T}) where T
    return vec(BackgroundQuantities.χz_array * R')
end

function grid_interpolator(AbstractCosmologicalProbes::AbstractCosmologicalProbes, 
    BackgroundQuantities::BackgroundQuantities, grid::Vector{T}) where T

    kernel_interpolated = zeros(T, AbstractCosmologicalProbes.n_bins, length(BackgroundQuantities.χz_array))

    for b in 1: AbstractCosmologicalProbes.n_bins
        interp = AkimaInterpolation(AbstractCosmologicalProbes.Kernel[b,:], BackgroundQuantities.χz_array, extrapolate=true)
        kernel_interpolated[b, :] = interp.(grid)
    end

    return kernel_interpolated
end

function get_kernel_array(AbstractCosmologicalProbes::GalaxyKernel, 
    BackgroundQuantities::BackgroundQuantities, R::Vector{T}) where T

    nχ = length(BackgroundQuantities.χz_array)
    nR = length(R)
    n_bins = AbstractCosmologicalProbes.n_bins

    W_array = reshape(grid_interpolator(GalaxyKernel, BackgroundQuantities, make_grid(BackgroundQuantities.χz_array, R)), n_bins, nχ, nR)

    return W_array
end

function get_kernel_array(AbstractCosmologicalProbes::Union{ShearKernel, CMBLensingKernel}, 
    BackgroundQuantities::BackgroundQuantities, R::Vector{T}) where T

    nχ = length(BackgroundQuantities.χz_array)
    nR = length(R)
    n_bins = AbstractCosmologicalProbes.n_bins

    W_L = grid_interpolator(AbstractCosmologicalProbes, BackgroundQuantities, make_grid(BackgroundQuantities.χz_array, R))
    χ2_app = make_grid(BackgroundQuantities.χz_array, R) .^ 2
    W_array = reshape( W_L./χ2_app , n_bins, nχ, nR)

    return W_array
end

function combine_kernels(ProbeA::Union{GalaxyKernel, ShearKernel, CMBLensingKernel}, 
    ProbeB::Union{GalaxyKernel, ShearKernel, CMBLensingKernel},
    BackgroundQuantities::BackgroundQuantities, R::Vector{T}) where T

    W_A = get_kernel_array(ProbeA, BackgroundQuantities, R)
    W_A_r1 = W_A[:,:,end]

    W_B = get_kernel_array(ProbeB, BackgroundQuantities, R)
    W_B_r1 = W_B[:,:,end]

    @tullio K[i,j,c,r] := W_A_r1[i,c] * W_B[j,c,r] + W_A[i,c,r]*W_B_r1[j,c]

    return K
end



"""
    simpson_weight_array(n::Int; T=Float64)

Computes the weights for the Simpson quadrature rule for numerical integration based on the input number of points `n`.

# Arguments
- `n::Int`: The number of points (must be at least 2).
- `T`: (optional) The type of the output array. Defaults to `Float64`.

# Returns
- An array of length `n` with the weights of type `T` for the Simpson quadrature rule.
"""
function simpson_weight_array(n::Int; T=Float64)
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
function compute_Cℓ(w::AbstractArray{T, 3}, K::AbstractArray{T, 4}, BackgroundQuantities::BackgroundQuantities, R::AbstractVector) where T
    nχ = length(BackgroundQuantities.χz_array)
    nR = length(R)

    if nχ != size(w, 2)
        throw(DimensionMismatch("Dimension mismatch: the χ array passed doesn't correspond to the one used in the evaluation of the inner k integral. Expected $nχ, got $(size(w, 2))."))
    end

    if nR != size(w, 3)
        throw(DimensionMismatch("Dimension mismatch: the R array passed doesn't correspond to the one used in the evaluation of the inner k integral. Expected $nR, got $(size(w, 3))."))
    end

    #Integration in χ is performed using the Simpson quadrature rule
    Δχ = ((last(BackgroundQuantities.χz_array)-first(BackgroundQuantities.χz_array))/(nχ-1))
    w_χ = simpson_weight_array(nχ)

    #Integration in R is performed using the Clenshaw-Curtis quadrature rule
    #CC_obj = FastTransforms.chebyshevmoments1(Float64, 2*nR+1)
    #w_R = FastTransforms.clenshawcurtisweights(CC_obj)
    w_R = get_clencurt_weights(-1, 1, 2*nR+1 )
    w_R = w_R[nR+2:end]
    w_R[1]/=2 #TODO: investigate if there are better solutions, this is not the analytic solution.

    @tullio Cℓ[l,i,j] := BackgroundQuantities.χz_array[n]*K[i,j,n,m]*w[l,n,m]*w_χ[n]*w_R[m]*Δχ

    return Cℓ

end
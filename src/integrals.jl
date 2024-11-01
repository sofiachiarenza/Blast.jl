#TODO: missing documentation of these new functions 

function make_grid(BackgroundQuantities::BackgroundQuantities, R::Vector{T}) where T
    return vec(BackgroundQuantities.χz_array * R')
end

function grid_interpolator(AbstractCosmologicalProbes::Union{GalaxyKernel, ShearKernel, CMBLensingKernel}, 
    BackgroundQuantities::BackgroundQuantities, grid::Vector{T}) where T

    n_bins = size(AbstractCosmologicalProbes.Kernel, 1)

    kernel_interpolated = zeros(n_bins, length(grid))

    for b in 1:n_bins
        interp = AkimaInterpolation(AbstractCosmologicalProbes.Kernel[b,:], BackgroundQuantities.χz_array, extrapolate=true)
        kernel_interpolated[b, :] = interp.(grid)
    end

    return kernel_interpolated
end

function get_kernel_array(AbstractCosmologicalProbes::GalaxyKernel, 
    BackgroundQuantities::BackgroundQuantities, R::Vector{T}) where T

    nχ = size(AbstractCosmologicalProbes.Kernel, 2)
    nR = length(R)
    n_bins = size(AbstractCosmologicalProbes.Kernel, 1)

    W_array = reshape(grid_interpolator(AbstractCosmologicalProbes, BackgroundQuantities, make_grid(BackgroundQuantities, R)), n_bins, nχ, nR)

    return W_array
end

function get_kernel_array(AbstractCosmologicalProbes::Union{ShearKernel, CMBLensingKernel}, 
    BackgroundQuantities::BackgroundQuantities, R::Vector{T}) where T

    nχ = size(AbstractCosmologicalProbes.Kernel, 2)
    nR = length(R)
    n_bins = size(AbstractCosmologicalProbes.Kernel, 1)

    W_L = grid_interpolator(AbstractCosmologicalProbes, BackgroundQuantities, make_grid(BackgroundQuantities, R))

    χ2_app = zeros(n_bins, nχ*nR)
    for i in 1:n_bins
        χ2_app[i,:] = make_grid(BackgroundQuantities, R) .^ 2
    end
    
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

function factorial_frac(n::Vector{T}) where T
    return @. (n-1)*n*(n+1)*(n+2)
end

#TODO: this function is HORRIBLE, pelase come up with something better!!!
function get_ell_prefactor(ProbeA::Union{GalaxyKernel, ShearKernel, CMBLensingKernel}, 
    ProbeB::Union{GalaxyKernel, ShearKernel, CMBLensingKernel}, ℓ_list::Vector)

    if isa(ProbeA, GalaxyKernel) && isa(ProbeB, GalaxyKernel)
        prefactor = 2 / π * ones(length(ℓ_list))
    elseif isa(ProbeA, GalaxyKernel) && isa(ProbeB, ShearKernel)
        prefactor =  2 / π * sqrt.(factorial_frac(ℓ_list))
    elseif isa(ProbeA, ShearKernel) && isa(ProbeB, GalaxyKernel)
        prefactor =  2 / π * sqrt.(factorial_frac(ℓ_list))
    elseif isa(ProbeA, ShearKernel) && isa(ProbeB, ShearKernel)
        prefactor = 2 / π * factorial_frac(ℓ_list)
    elseif isa(ProbeA, CMBLensingKernel) && isa(ProbeB, ShearKernel)
        prefactor = 2 / π * ones(length(ℓ_list)) #TODO: figure out all cmb lensing prefactors!!!!!!
    elseif isa(ProbeA, CMBLensingKernel) && isa(ProbeB, CMBLensingKernel)
        prefactor = 2 / π * ones(length(ℓ_list))
    elseif isa(ProbeA, CMBLensingKernel) && isa(ProbeB, GalaxyKernel)
        prefactor = 2 / π * ones(length(ℓ_list))
    end

    return prefactor
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
function compute_Cℓ(w::AbstractArray{T, 3}, ProbeA::Union{GalaxyKernel, ShearKernel, CMBLensingKernel}, 
    ProbeB::Union{GalaxyKernel, ShearKernel, CMBLensingKernel}, BackgroundQuantities::BackgroundQuantities, R::AbstractVector, ℓ_list::AbstractArray{T,1} = Blast.ℓ) where T

    nχ = length(BackgroundQuantities.χz_array)
    nR = length(R)

    K = combine_kernels(ProbeA, ProbeB, BackgroundQuantities, R)

    #Integration in χ is performed using the Simpson quadrature rule
    Δχ = ((last(BackgroundQuantities.χz_array)-first(BackgroundQuantities.χz_array))/(nχ-1))
    w_χ = simpson_weight_array(nχ)

    #Integration in R is performed using the Clenshaw-Curtis quadrature rule
    w_R = get_clencurt_weights(-1, 1, 2*nR+1 )
    w_R = w_R[nR+2:end]
    w_R[1]/=2 #TODO: investigate if there are better solutions, this is not the analytic solution.

    ell_prefactor = get_ell_prefactor(ProbeA, ProbeB, ℓ_list)

    @tullio Cℓ[l,i,j] := ell_prefactor[l]*BackgroundQuantities.χz_array[n]*K[i,j,n,m]*w[l,n,m]*w_χ[n]*w_R[m]*Δχ

    return Cℓ 
end
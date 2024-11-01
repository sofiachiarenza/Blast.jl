"""
    make_grid(BackgroundQuantities::BackgroundQuantities, R::Vector{T}) where T

Constructs a grid by multiplying the `χz_array` from `BackgroundQuantities` with the vector `R`.

# Arguments
- `BackgroundQuantities::BackgroundQuantities`: An instance of the `BackgroundQuantities` type that contains the `χz_array`.
- `R::Vector{T}`: A vector of values to be used in the grid construction, where `T` can be any type.
"""
function make_grid(BackgroundQuantities::BackgroundQuantities, R::Vector{T}) where T
    return vec(BackgroundQuantities.χz_array * R')
end

"""
    grid_interpolator(AbstractCosmologicalProbes::Union{GalaxyKernel, ShearKernel, CMBLensingKernel}, 
        BackgroundQuantities::BackgroundQuantities, grid::Vector{T}) where T

Interpolates the kernel values for a given grid based on the specified cosmological probes. 
Returns a 2D array of interpolated kernel values, where rows correspond to the number of bins and columns correspond to the grid points.

# Arguments
- `AbstractCosmologicalProbes::Union{GalaxyKernel, ShearKernel, CMBLensingKernel}`: The kernel data to be interpolated.
- `BackgroundQuantities::BackgroundQuantities`: An instance of the `BackgroundQuantities` type that contains the `χz_array`.
- `grid::Vector{T}`: A vector of values where the interpolated kernel values will be evaluated.
"""
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

"""
    get_kernel_array(AbstractCosmologicalProbes::GalaxyKernel, 
        BackgroundQuantities::BackgroundQuantities, R::Vector{T}) where T

Obtains the kernel array for the `GalaxyKernel` probe, with dimensions (bins, nχ, nR).

# Arguments
- `AbstractCosmologicalProbes::GalaxyKernel`: An instance of the `GalaxyKernel` type.
- `BackgroundQuantities::BackgroundQuantities`: An instance of the `BackgroundQuantities` type.
- `R::Vector{T}`: A vector of values for which the kernel array is to be computed.
"""
function get_kernel_array(AbstractCosmologicalProbes::GalaxyKernel, 
    BackgroundQuantities::BackgroundQuantities, R::Vector{T}) where T

    n_bins = size(AbstractCosmologicalProbes.Kernel, 1)
    nχ = size(AbstractCosmologicalProbes.Kernel, 2)
    nR = length(R)
    
    W_array = reshape(grid_interpolator(AbstractCosmologicalProbes, BackgroundQuantities, make_grid(BackgroundQuantities, R)), n_bins, nχ, nR)

    return W_array
end

"""
    get_kernel_array(AbstractCosmologicalProbes::Union{ShearKernel, CMBLensingKernel}, 
        BackgroundQuantities::BackgroundQuantities, R::Vector{T}) where T

Obtains the kernel array for the `ShearKernel` or `CMBLensingKernel` probe, with dimensions (bins, nχ, nR)
The difference with respect to the galaxy case is that these kernels are divided by χ².

# Arguments
- `AbstractCosmologicalProbes::Union{ShearKernel, CMBLensingKernel}`: An instance of either `ShearKernel` or `CMBLensingKernel`.
- `BackgroundQuantities::BackgroundQuantities`: An instance of the `BackgroundQuantities` type.
- `R::Vector{T}`: A vector of values for which the kernel array is to be computed.
"""
function get_kernel_array(AbstractCosmologicalProbes::Union{ShearKernel, CMBLensingKernel}, 
    BackgroundQuantities::BackgroundQuantities, R::Vector{T}) where T

    n_bins = size(AbstractCosmologicalProbes.Kernel, 1)
    nχ = size(AbstractCosmologicalProbes.Kernel, 2)
    nR = length(R)
    
    W_L = grid_interpolator(AbstractCosmologicalProbes, BackgroundQuantities, make_grid(BackgroundQuantities, R))

    χ2_app = zeros(n_bins, nχ*nR)
    for i in 1:n_bins
        χ2_app[i,:] = make_grid(BackgroundQuantities, R) .^ 2
    end
    
    W_array = reshape( W_L./χ2_app , n_bins, nχ, nR)

    return W_array
end

"""
    combine_kernels(ProbeA::Union{GalaxyKernel, ShearKernel, CMBLensingKernel}, 
        ProbeB::Union{GalaxyKernel, ShearKernel, CMBLensingKernel}, 
        BackgroundQuantities::BackgroundQuantities, R::Vector{T}) where T

Combines the kernels from two different cosmological probes into a single array. 
This is needed to perform the integration in the χ-R coordinates.

# Arguments
- `ProbeA::Union{GalaxyKernel, ShearKernel, CMBLensingKernel}`: The first cosmological probe to combine.
- `ProbeB::Union{GalaxyKernel, ShearKernel, CMBLensingKernel}`: The second cosmological probe to combine.
- `BackgroundQuantities::BackgroundQuantities`: An instance of the `BackgroundQuantities` type.
- `R::Vector{T}`: A vector of values used in the combination.
"""
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
    factorial_frac(ℓ::Union{Number,Vector{T}}}) where T

Computes the ratio (ℓ+2)!/(ℓ-2)!, needed in the pre-factors of the the angular power spectra.

# Arguments
- `ℓ::Vector{T}`: vectors of ℓ values.
"""
function factorial_frac(ℓ::Union{Number,Vector{T}}) where T
    return @. (ℓ-1)*ℓ*(ℓ+1)*(ℓ+2)
end

"""
    get_ell_prefactor(ProbeA::Union{GalaxyKernel, ShearKernel, CMBLensingKernel}, 
        ProbeB::Union{GalaxyKernel, ShearKernel, CMBLensingKernel}, ℓ_list::Vector)

Calculates the prefactor for the angular power spectrum based on the types of the two probes.

# Arguments
- `ProbeA::Union{GalaxyKernel, ShearKernel, CMBLensingKernel}`: The first probe used to determine the prefactor.
- `ProbeB::Union{GalaxyKernel, ShearKernel, CMBLensingKernel}`: The second probe used to determine the prefactor.
- `ℓ_list::Vector`: A vector of angular multipole values.

# Returns
- A vector of prefactor values corresponding to the input `ℓ_list`.
"""
function get_ell_prefactor(ProbeA::Union{GalaxyKernel, ShearKernel, CMBLensingKernel}, 
    ProbeB::Union{GalaxyKernel, ShearKernel, CMBLensingKernel}, ℓ_list::Vector)

    #TODO: this function is HORRIBLE, pelase come up with something better!!!
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
    compute_Cℓ(w::AbstractArray{T, 3}, 
               ProbeA::Union{GalaxyKernel, ShearKernel, CMBLensingKernel}, 
               ProbeB::Union{GalaxyKernel, ShearKernel, CMBLensingKernel}, 
               BackgroundQuantities::BackgroundQuantities, 
               R::AbstractVector, 
               ℓ_list::AbstractArray{T,1} = Blast.ℓ) where T

Computes the angular power spectrum `Cℓ` by performing two outer integrals over `χ` and `R`. 
The Simpson quadrature rule is used for integration over `χ`, while Clenshaw-Curtis quadrature is used for `R`.

# Arguments
- `w::AbstractArray{T, 3}`: A 3D array representing the projected matter densities, containing the inner integrals over `k`. The array dimensions are (ℓ, χ, R).
- `ProbeA::Union{GalaxyKernel, ShearKernel, CMBLensingKernel}`: The first cosmological probe kernel type.
- `ProbeB::Union{GalaxyKernel, ShearKernel, CMBLensingKernel}`: The second cosmological probe kernel type.
- `BackgroundQuantities::BackgroundQuantities`: Contains background information, including `χz_array` for distances in the cosmology.
- `R::AbstractVector`: A 1D array representing the radial grid values for integration.
- `ℓ_list::AbstractArray{T,1}`: An optional 1D array of angular multipole values, defaulting to the global variable Blast.ℓ, the list of ℓ values used for the precomputed part.

# Returns
- A 3D array `Cℓ` with dimensions (ℓ, i, j), where `i` and `j` represent the tomographic bins. The array contains the computed angular power spectrum coefficients for each combination of `ℓ` values and tomographic bins.
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
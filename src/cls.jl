#TODO: re-think background handling
function make_grid(χ::Array{<:Any, 1}, R::Array{<:Any, 1}) 
    return vec(χ * R') 
end

#Change of coordinated to χ-R grid
function grid_interpolator(Probe::AbstractComponents, grid::Array{<:Any, 1}) 

    n_bins = size(Probe.Kernel, 1)
    kernel_interpolated = zeros(n_bins, length(grid))

    for b in 1:n_bins
        interp = AkimaInterpolation(Probe.Kernel[b,:], Blast.χ, extrapolation=ExtrapolationType.Extension)
        kernel_interpolated[b, :] = interp.(grid)
    end

    return kernel_interpolated
end

function get_kernel_array(Probe::Union{NumberCounts, RedshiftSpaceDistortions, PrimordialNonGaussianity, IntegratedSachsWolfe})

    n_bins = size(Probe.Kernel, 1)
    nχ = size(Blast.χ, 1)
    nR = size(Blast.R, 1)
    
    W_array = reshape(grid_interpolator(Probe, make_grid(Blast.χ, Blast.R)), n_bins, nχ, nR)

    return W_array
end

function get_kernel_array(Probe::Union{CosmicShear, IntrinsicAlignment, MagnificationBias, CMBLensing})

    n_bins = size(Probe.Kernel, 1)
    nχ = size(Blast.χ, 1)
    nR = size(Blast.R, 1)
    
    W_L = grid_interpolator(Probe, make_grid(Blast.χ, Blast.R))

    χ2_app = zeros(n_bins, nχ*nR)
    for i in 1:n_bins
        χ2_app[i,:] = make_grid(Blast.χ, Blast.R) .^ 2
    end
    
    W_array = reshape( W_L./χ2_app , n_bins, nχ, nR)

    return W_array
end

#Symmetryzation of the integral
function combine_kernels(ProbeA::AbstractComponents, ProbeB::AbstractComponents)

    W_A = get_kernel_array(ProbeA)
    W_A_r1 = W_A[:,:,end]

    W_B = get_kernel_array(ProbeB)
    W_B_r1 = W_B[:,:,end]

    @tullio K[i,j,c,r] := W_A_r1[i,c] * W_B[j,c,r] + W_A[i,c,r] * W_B_r1[j,c]

    return K
end

#Now the nightmare begins
function compute_Cℓ(Component1::AbstractComponents, Component2::AbstractComponents, w::ProjectedMatterDensityComponent) 

    nχ = size(Blast.χ, 1)
    nR = size(Blast.R, 1)

    K = combine_kernels(Component1, Component2)

    #Integration in χ is performed using the Simpson quadrature rule
    Δχ = ((last(Blast.χ)-first(Blast.χ))/(nχ-1))
    w_χ = simpson_weight_array(nχ)

    #Integration in R is performed using the Clenshaw-Curtis quadrature rule
    w_R = get_clencurt_weights_R_integration(2*nR+1)

    prefactor = Component1.ell_prefactor .* Component2.ell_prefactor
    pmj = w.w

    @tullio Cℓ[l,i,j] := prefactor[l]*Blast.χ[n]*K[i,j,n,m]*pmj[l,n,m]*w_χ[n]*w_R[m]*Δχ

    return Cℓ 
end

function compute_Cℓ(Component1::AbstractComponents, Component2::AbstractComponents, 
    w02::ProjectedMatterDensityComponent, w20::ProjectedMatterDensityComponent) 

    nχ = size(Blast.χ, 1)
    nR = size(Blast.R, 1)

    W_A = get_kernel_array(Component1)
    W_A_r1 = W_A[:,:,end]

    W_B = get_kernel_array(Component2)
    W_B_r1 = W_B[:,:,end]

    pmj02 = w02.w
    pmj20 = w20.w

    @tullio K[l,i,j,c,r] := W_A_r1[i,c] * W_B[j,c,r] * pmj02[l,c,r] + W_A[i,c,r] * W_B_r1[j,c] * pmj20[l,c,r]

    #Integration in χ is performed using the Simpson quadrature rule
    Δχ = ((last(Blast.χ)-first(Blast.χ))/(nχ-1))
    w_χ = simpson_weight_array(nχ)

    #Integration in R is performed using the Clenshaw-Curtis quadrature rule
    w_R = get_clencurt_weights_R_integration(2*nR+1)

    prefactor = Component1.ell_prefactor .* Component2.ell_prefactor

    @tullio Cℓ[l,i,j] := prefactor[l]*Blast.χ[n]*K[l,i,j,n,m]*w_χ[n]*w_R[m]*Δχ
end


function get_Cℓ(Component1::NumberCounts, Component2::NumberCounts, W::ProjectedMatterDensity)
    return compute_Cℓ(Component1, Component2, W.w_2_00_ϕTT)
end

function get_Cℓ(Component1::NumberCounts, Component2::RedshiftSpaceDistortions, W::ProjectedMatterDensity)
    return compute_Cℓ(Component1, Component2, W.w_2_02_ϕTT, W.w_2_20_ϕTT)
end

function get_Cℓ(Component1::RedshiftSpaceDistortions, Component2::NumberCounts, W::ProjectedMatterDensity)
    return compute_Cℓ(Component1, Component2, W.w_2_20_ϕTT, W.w_2_02_ϕTT)
end

function get_Cℓ(Component1::RedshiftSpaceDistortions, Component2::RedshiftSpaceDistortions, W::ProjectedMatterDensity)
    return compute_Cℓ(Component1, Component2, W.w_2_22_ϕTT)
end

function get_Cℓ(Component1::NumberCounts, Component2::MagnificationBias, W::ProjectedMatterDensity)
    return compute_Cℓ(Component1, Component2, W.w_0_00_ϕTT)
end

function get_Cℓ(Component1::MagnificationBias, Component2::NumberCounts, W::ProjectedMatterDensity)
    return compute_Cℓ(Component1, Component2, W.w_0_00_ϕTT)
end

function get_Cℓ(Component1::MagnificationBias, Component2::MagnificationBias, W::ProjectedMatterDensity)
    return compute_Cℓ(Component1, Component2, W.w_minus2_00_ϕTT)
end

function get_Cℓ(Component1::MagnificationBias, Component2::RedshiftSpaceDistortions, W::ProjectedMatterDensity)
    return compute_Cℓ(Component1, Component2, W.w_0_02_ϕTT, W.w_0_20_ϕTT)
end

function get_Cℓ(Component1::RedshiftSpaceDistortions, Component2::MagnificationBias, W::ProjectedMatterDensity)
    return compute_Cℓ(Component1, Component2, W.w_0_20_ϕTT, W.w_0_02_ϕTT)
end

function get_Cℓ(Component1::NumberCounts, Component2::PrimordialNonGaussianity, W::ProjectedMatterDensity)
    return compute_Cℓ(Component1, Component2, W.w_2_00_ϕT)
end

function get_Cℓ(Component1::PrimordialNonGaussianity, Component2::NumberCounts, W::ProjectedMatterDensity)
    return compute_Cℓ(Component1, Component2, W.w_2_00_ϕT)
end

function get_Cℓ(Component1::PrimordialNonGaussianity, Component2::RedshiftSpaceDistortions, W::ProjectedMatterDensity)
    return compute_Cℓ(Component1, Component2, W.w_2_02_ϕT, W.w_2_20_ϕT)
end

function get_Cℓ(Component1::RedshiftSpaceDistortions, Component2::PrimordialNonGaussianity, W::ProjectedMatterDensity)
    return compute_Cℓ(Component1, Component2, W.w_2_20_ϕT, W.w_2_02_ϕT)
end

function get_Cℓ(Component1::MagnificationBias, Component2::PrimordialNonGaussianity, W::ProjectedMatterDensity)
    return compute_Cℓ(Component1, Component2, W.w_0_00_ϕT)
end

function get_Cℓ(Component1::PrimordialNonGaussianity, Component2::MagnificationBias, W::ProjectedMatterDensity)
    return compute_Cℓ(Component1, Component2, W.w_0_00_ϕT)
end

function get_Cℓ(Component1::PrimordialNonGaussianity, Component2::PrimordialNonGaussianity, W::ProjectedMatterDensity)
    return compute_Cℓ(Component1, Component2, W.w_2_00_ϕ)
end

function get_Cℓ(G::GalaxyClustering, W::ProjectedMatterDensity)
    Cℓ_δδ = get_Cℓ(G.δ, G.δ, W)
    Cℓ_δRSD = get_Cℓ(G.δ, G.RSD, W)
    Cℓ_RSDδ = get_Cℓ(G.RSD, G.δ, W)
    Cℓ_RSDRSD = get_Cℓ(G.RSD, G.RSD, W)
    Cℓ_δμ = get_Cℓ(G.δ, G.μ, W)
    Cℓ_μδ = get_Cℓ(G.μ, G.δ, W)
    Cℓ_μμ = get_Cℓ(G.μ, G.μ, W)
    Cℓ_μRSD = get_Cℓ(G.μ, G.RSD, W)
    Cℓ_RSDμ = get_Cℓ(G.RSD, G.μ, W)
    Cℓ_δfNL = get_Cℓ(G.δ, G.PNG, W)
    Cℓ_fNLδ = get_Cℓ(G.PNG, G.δ, W)
    Cℓ_fNLRSD =get_Cℓ(G.PNG, G.RSD, W)
    Cℓ_RSDfNL =get_Cℓ(G.RSD, G.PNG, W)
    Cℓ_μfNL =get_Cℓ(G.μ, G.PNG, W)
    Cℓ_fNLμ =get_Cℓ(G.PNG, G.μ, W)
    Cℓ_fNLfNL =get_Cℓ(G.PNG, G.PNG, W)

    return @. Cℓ_δδ - Cℓ_δRSD - Cℓ_RSDδ + Cℓ_RSDRSD + Cℓ_δμ + Cℓ_μδ + Cℓ_μμ - Cℓ_μRSD - Cℓ_RSDμ + Cℓ_δfNL + Cℓ_fNLδ - Cℓ_fNLRSD - Cℓ_RSDfNL + Cℓ_μfNL + Cℓ_fNLμ + Cℓ_fNLfNL
end


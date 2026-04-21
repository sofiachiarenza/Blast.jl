"""
    make_grid(χ::Array{<:Any, 1}, R::Array{<:Any, 1}) 

Construct a flattened χ–R grid corresponding to all combinations of comoving distance `χ` and ratio `R`.
"""
function make_grid(χ::Array{<:Any, 1}, R::Array{<:Any, 1}) 
    return vec(χ * R') 
end

"""
    grid_interpolator(Probe, bg)

Interpolate a probe kernel from the χ grid onto a χ–R grid using Background interpolators.
"""
function grid_interpolator(Probe::AbstractComponents, bg::Background) 
    kernel = Probe.Kernel
    chi_grid = bg.χ
    grid_flat = make_grid(Blast.χ, Blast.R)
    return collect(akima_interpolation(collect(kernel'), chi_grid, grid_flat)')
end

"""
    get_kernel_array(Probe, bg)

Return the probe kernel evaluated on the χ–R grid.
"""
function get_kernel_array(Probe::Union{NumberCounts, RedshiftSpaceDistortions, PrimordialNonGaussianity, IntegratedSachsWolfe}, bg::Background)
    n_bins = size(Probe.Kernel, 1)
    nχ = length(bg.χ)
    nR = length(Blast.R)
    W_array = reshape(grid_interpolator(Probe, bg), n_bins, nχ, nR)
    return W_array
end

function get_kernel_array(Probe::Union{CosmicShear, IntrinsicAlignment, MagnificationBias, CMBLensing}, bg::Background)
    n_bins = size(Probe.Kernel, 1)
    nχ = length(bg.χ)
    nR = length(Blast.R)
    W_L = grid_interpolator(Probe, bg)
    χ2_app = zeros(n_bins, nχ*nR)
    grid_flat = make_grid(Blast.χ, Blast.R)
    for i in 1:n_bins
        χ2_app[i,:] = grid_flat .^ 2
    end
    W_array = reshape( W_L./χ2_app , n_bins, nχ, nR)
    return W_array
end

function _combine_kernels_tullio(W_A, W_B)
    W_A_r1 = W_A[:, :, end]
    W_B_r1 = W_B[:, :, end]
    @tullio K[idx_i, idx_j, idx_c, idx_r] := W_A_r1[idx_i, idx_c] * W_B[idx_j, idx_c, idx_r] + W_A[idx_i, idx_c, idx_r] * W_B_r1[idx_j, idx_c]
    return K
end

"""
    combine_kernels(ProbeA, ProbeB, bg)
"""
function combine_kernels(ProbeA::AbstractComponents, ProbeB::AbstractComponents, bg::Background)
    W_A = get_kernel_array(ProbeA, bg)
    W_B = get_kernel_array(ProbeB, bg)
    return _combine_kernels_tullio(W_A, W_B)
end

function _compute_Cℓ_tullio(K, pmj, w_χ, w_R, prefactor, Δχ, χ_grid)
    @tullio Cℓ[idx_l, idx_i, idx_j] := prefactor[idx_l] * χ_grid[idx_n] * K[idx_i, idx_j, idx_n, idx_m] * pmj[idx_l, idx_n, idx_m] * w_χ[idx_n] * w_R[idx_m] * Δχ
    return Cℓ
end

function _compute_Cℓ_rsd_tullio(W_A_r1, W_B, pmj02, W_A, W_B_r1, pmj20, w_χ, w_R, prefactor, Δχ, χ_grid)
    @tullio K[idx_l, idx_i, idx_j, idx_c, idx_r] := W_A_r1[idx_i, idx_c] * W_B[idx_j, idx_c, idx_r] * pmj02[idx_l, idx_c, idx_r] + W_A[idx_i, idx_c, idx_r] * W_B_r1[idx_j, idx_c] * pmj20[idx_l, idx_c, idx_r]
    @tullio Cℓ[idx_l, idx_i, idx_j] := prefactor[idx_l] * χ_grid[idx_n] * K[idx_l, idx_i, idx_j, idx_n, idx_m] * w_χ[idx_n] * w_R[idx_m] * Δχ
    return Cℓ
end

# Shared finalize step: given `full_Cℓ` stacked on `full_ℓ_range`, interpolate
# onto the requested ℓ grid via a Chebyshev decomposition of ℓ² · C_ℓ (the
# ℓ² weighting smooths the spectrum before decomposition; we divide it back
# out after chebinterp_native). Called by every `get_Cℓ(ℓ, ...)` mixing the
# non-Limber and Limber paths.
#
# `eltype(full_Cℓ)` is used for the zeros(...) initializer so ForwardDiff.Dual
# threads through the loop assignment; plain `zeros(...)` would silently strip
# Duals at the `Cℓ_final[:, i, j] = ...` store.
function _finalize_Cℓ(full_Cℓ, ℓ, nbins_A::Integer, nbins_B::Integer, P::Union{FFTPlans, Nothing})
    Cℓ_final = zeros(eltype(full_Cℓ), size(ℓ, 1), nbins_A, nbins_B)
    plan_interp = isnothing(P) ? prepare_chebyshev_plan(2.0, 2000.0, 100) : P.plan_ℓ
    for i in 1:nbins_A, j in 1:nbins_B
        c_coeffs = chebyshev_decomposition(plan_interp, reverse(full_Cℓ[:, i, j] .* (Blast.full_ℓ_range .^ 2.)))
        Cℓ_final[:, i, j] = chebinterp_native(c_coeffs, ℓ, 2.0, 2000.0) ./ (ℓ .^ 2)
    end
    return Cℓ_final
end

"""
    compute_Cℓ(Component1, Component2, w, bg)
"""
function compute_Cℓ(Component1::AbstractComponents, Component2::AbstractComponents, w::ProjectedMatterDensityComponent, bg::Background) 
    nχ = size(bg.χ, 1)
    nR = size(Blast.R, 1)
    K = combine_kernels(Component1, Component2, bg)
    Δχ = ((last(bg.χ)-first(bg.χ))/(nχ-1))
    w_χ = simpson_weights_array(nχ)
    w_R = get_clencurt_weights_R_integration(2*nR+1)
    prefactor = 2/π .* Component1.ell_prefactor[1:length(ℓ_nonlimber)] .* Component2.ell_prefactor[1:length(ℓ_nonlimber)]
    pmj = w.w
    χ_grid = bg.χ
    return _compute_Cℓ_tullio(K, pmj, w_χ, w_R, prefactor, Δχ, χ_grid)
end

"""
    compute_Cℓ(Component1, Component2, w02, w20, bg)
"""
function compute_Cℓ(Component1::AbstractComponents, Component2::AbstractComponents, 
    w02::ProjectedMatterDensityComponent, w20::ProjectedMatterDensityComponent, bg::Background) 
    nχ = size(bg.χ, 1)
    nR = size(Blast.R, 1)
    W_A = get_kernel_array(Component1, bg)
    W_A_r1 = W_A[:,:,end]
    W_B = get_kernel_array(Component2, bg)
    W_B_r1 = W_B[:,:,end]
    pmj02 = w02.w
    pmj20 = w20.w
    Δχ = ((last(bg.χ)-first(bg.χ))/(nχ-1))
    w_χ = simpson_weights_array(nχ)
    w_R = get_clencurt_weights_R_integration(2*nR+1)
    prefactor = 2/π .* Component1.ell_prefactor[1:length(ℓ_nonlimber)] .* Component2.ell_prefactor[1:length(ℓ_nonlimber)]
    χ_grid = bg.χ
    return _compute_Cℓ_rsd_tullio(W_A_r1, W_B, pmj02, W_A, W_B_r1, pmj20, w_χ, w_R, prefactor, Δχ, χ_grid)
end

"""
    get_Cℓ(...)
"""
function get_Cℓ(Component1::Union{AbstractComponents, Nothing}, Component2::Union{AbstractComponents, Nothing}, W::ProjectedMatterDensity, bg::Background)
    return 0.
end

## Galaxy clustering auto
function get_Cℓ(Component1::NumberCounts, Component2::NumberCounts, W::ProjectedMatterDensity, bg::Background)
    return compute_Cℓ(Component1, Component2, W.w_2_00_ϕTT, bg)
end

function get_Cℓ(Component1::NumberCounts, Component2::RedshiftSpaceDistortions, W::ProjectedMatterDensity, bg::Background)
    return compute_Cℓ(Component1, Component2, W.w_2_02_ϕTT, W.w_2_20_ϕTT, bg)
end

function get_Cℓ(Component1::RedshiftSpaceDistortions, Component2::NumberCounts, W::ProjectedMatterDensity, bg::Background)
    return compute_Cℓ(Component1, Component2, W.w_2_20_ϕTT, W.w_2_02_ϕTT, bg)
end

function get_Cℓ(Component1::RedshiftSpaceDistortions, Component2::RedshiftSpaceDistortions, W::ProjectedMatterDensity, bg::Background)
    return compute_Cℓ(Component1, Component2, W.w_2_22_ϕTT, bg)
end

function get_Cℓ(Component1::NumberCounts, Component2::MagnificationBias, W::ProjectedMatterDensity, bg::Background)
    return compute_Cℓ(Component1, Component2, W.w_0_00_ϕTT, bg)
end

function get_Cℓ(Component1::MagnificationBias, Component2::NumberCounts, W::ProjectedMatterDensity, bg::Background)
    return compute_Cℓ(Component1, Component2, W.w_0_00_ϕTT, bg)
end

function get_Cℓ(Component1::MagnificationBias, Component2::MagnificationBias, W::ProjectedMatterDensity, bg::Background)
    return compute_Cℓ(Component1, Component2, W.w_minus2_00_ϕTT, bg)
end

function get_Cℓ(Component1::MagnificationBias, Component2::RedshiftSpaceDistortions, W::ProjectedMatterDensity, bg::Background)
    return compute_Cℓ(Component1, Component2, W.w_0_02_ϕTT, W.w_0_20_ϕTT, bg)
end

function get_Cℓ(Component1::RedshiftSpaceDistortions, Component2::MagnificationBias, W::ProjectedMatterDensity, bg::Background)
    return compute_Cℓ(Component1, Component2, W.w_0_20_ϕTT, W.w_0_02_ϕTT, bg)
end

function get_Cℓ(Component1::NumberCounts, Component2::PrimordialNonGaussianity, W::ProjectedMatterDensity, bg::Background)
    return compute_Cℓ(Component1, Component2, W.w_2_00_ϕT_R1, bg)
end

function get_Cℓ(Component1::PrimordialNonGaussianity, Component2::NumberCounts, W::ProjectedMatterDensity, bg::Background)
    return compute_Cℓ(Component1, Component2, W.w_2_00_ϕT, bg)
end

function get_Cℓ(Component1::PrimordialNonGaussianity, Component2::RedshiftSpaceDistortions, W::ProjectedMatterDensity, bg::Background)
    return compute_Cℓ(Component1, Component2, W.w_2_02_ϕT, W.w_2_20_ϕT, bg)
end

function get_Cℓ(Component1::RedshiftSpaceDistortions, Component2::PrimordialNonGaussianity, W::ProjectedMatterDensity, bg::Background)
    return compute_Cℓ(Component1, Component2, W.w_2_20_ϕT_R1, W.w_2_02_ϕT_R1, bg)
end

function get_Cℓ(Component1::MagnificationBias, Component2::PrimordialNonGaussianity, W::ProjectedMatterDensity, bg::Background)
    return compute_Cℓ(Component1, Component2, W.w_0_00_ϕT_R1, bg)
end

function get_Cℓ(Component1::PrimordialNonGaussianity, Component2::MagnificationBias, W::ProjectedMatterDensity, bg::Background)
    return compute_Cℓ(Component1, Component2, W.w_0_00_ϕT, bg)
end

function get_Cℓ(Component1::PrimordialNonGaussianity, Component2::PrimordialNonGaussianity, W::ProjectedMatterDensity, bg::Background)
    return compute_Cℓ(Component1, Component2, W.w_2_00_ϕ, bg)
end

function get_Cℓ(ℓ::AbstractArray{<:Any,1}, G::GalaxyClustering, Pk::PowerSpectrum, W::ProjectedMatterDensity, bg::Background, P::Union{FFTPlans, Nothing}=nothing)
    Cℓ_δδ = get_Cℓ(G.δ, G.δ, W, bg)
    Cℓ_δRSD = get_Cℓ(G.δ, G.RSD, W, bg)
    Cℓ_RSDδ = get_Cℓ(G.RSD, G.δ, W, bg)
    Cℓ_RSDRSD = get_Cℓ(G.RSD, G.RSD, W, bg)
    Cℓ_δμ = get_Cℓ(G.δ, G.μ, W, bg)
    Cℓ_μδ = get_Cℓ(G.μ, G.δ, W, bg)
    Cℓ_μμ = get_Cℓ(G.μ, G.μ, W, bg)
    Cℓ_μRSD = get_Cℓ(G.μ, G.RSD, W, bg)
    Cℓ_RSDμ = get_Cℓ(G.RSD, G.μ, W, bg)
    Cℓ_δfNL = get_Cℓ(G.δ, G.PNG, W, bg)
    Cℓ_fNLδ = get_Cℓ(G.PNG, G.δ, W, bg)
    Cℓ_fNLRSD =get_Cℓ(G.PNG, G.RSD, W, bg)
    Cℓ_RSDfNL =get_Cℓ(G.RSD, G.PNG, W, bg)
    Cℓ_μfNL =get_Cℓ(G.μ, G.PNG, W, bg)
    Cℓ_fNLμ =get_Cℓ(G.PNG, G.μ, W, bg)
    Cℓ_fNLfNL =get_Cℓ(G.PNG, G.PNG, W, bg)

    Cℓ_nonlimber = @. Cℓ_δδ - Cℓ_δRSD - Cℓ_RSDδ + Cℓ_RSDRSD + Cℓ_δμ + Cℓ_μδ + Cℓ_μμ - Cℓ_μRSD - Cℓ_RSDμ + Cℓ_δfNL + Cℓ_fNLδ - Cℓ_fNLRSD - Cℓ_RSDfNL + Cℓ_μfNL + Cℓ_fNLμ + Cℓ_fNLfNL
    Cℓ_correction = get_limber_correction(G, Pk)
    Cℓ_limber = get_limber_Cℓ(G, Pk)
    full_Cℓ = cat(Cℓ_nonlimber .+ Cℓ_correction, Cℓ_limber; dims=1) 

    nbins = size(G.δ.Kernel, 1)
    return _finalize_Cℓ(full_Cℓ, ℓ, nbins, nbins, P)
end

## Weak lensing auto
function get_Cℓ(Component1::CosmicShear, Component2::CosmicShear, W::ProjectedMatterDensity, bg::Background)
    return compute_Cℓ(Component1, Component2, W.w_minus2_00_ϕTT, bg)
end

function get_Cℓ(Component1::CosmicShear, Component2::IntrinsicAlignment, W::ProjectedMatterDensity, bg::Background)
    return compute_Cℓ(Component1, Component2, W.w_minus2_00_ϕTT, bg)
end

function get_Cℓ(Component1::IntrinsicAlignment, Component2::CosmicShear, W::ProjectedMatterDensity, bg::Background)
    return compute_Cℓ(Component1, Component2, W.w_minus2_00_ϕTT, bg)
end

function get_Cℓ(Component1::IntrinsicAlignment, Component2::IntrinsicAlignment, W::ProjectedMatterDensity, bg::Background)
    return compute_Cℓ(Component1, Component2, W.w_minus2_00_ϕTT, bg)
end

function get_Cℓ(ℓ::AbstractArray{<:Any, 1}, L::WeakLensing, Pk::PowerSpectrum, W::ProjectedMatterDensity, bg::Background, P::Union{FFTPlans, Nothing}=nothing)
    Cℓ_γγ = get_Cℓ(L.γ, L.γ, W, bg)
    Cℓ_γI = get_Cℓ(L.γ, L.IA, W, bg)
    Cℓ_Iγ = get_Cℓ(L.IA, L.γ, W, bg)
    Cℓ_II = get_Cℓ(L.IA, L.IA, W, bg)

    Cℓ_nonlimber = @. Cℓ_γγ + Cℓ_γI + Cℓ_Iγ + Cℓ_II
    Cℓ_correction = get_limber_correction(L, Pk)
    Cℓ_limber = get_limber_Cℓ(L, Pk)
    full_Cℓ = cat(Cℓ_nonlimber .+ Cℓ_correction, Cℓ_limber; dims=1) 

    nbins = size(L.γ.Kernel, 1)
    return _finalize_Cℓ(full_Cℓ, ℓ, nbins, nbins, P)
end

## Cross clustering-lensing
function get_Cℓ(Component1::NumberCounts, Component2::CosmicShear, W::ProjectedMatterDensity, bg::Background)
    return compute_Cℓ(Component1, Component2, W.w_0_00_ϕTT, bg)
end

function get_Cℓ(Component1::NumberCounts, Component2::IntrinsicAlignment, W::ProjectedMatterDensity, bg::Background)
    return compute_Cℓ(Component1, Component2, W.w_0_00_ϕTT, bg)
end

function get_Cℓ(Component1::RedshiftSpaceDistortions, Component2::CosmicShear, W::ProjectedMatterDensity, bg::Background)
    return compute_Cℓ(Component1, Component2, W.w_0_20_ϕTT, W.w_0_02_ϕTT, bg)
end

function get_Cℓ(Component1::RedshiftSpaceDistortions, Component2::IntrinsicAlignment, W::ProjectedMatterDensity, bg::Background)
    return compute_Cℓ(Component1, Component2, W.w_0_20_ϕTT, W.w_0_02_ϕTT, bg)
end

function get_Cℓ(Component1::MagnificationBias, Component2::CosmicShear, W::ProjectedMatterDensity, bg::Background)
    return compute_Cℓ(Component1, Component2, W.w_minus2_00_ϕTT, bg)
end

function get_Cℓ(Component1::MagnificationBias, Component2::IntrinsicAlignment, W::ProjectedMatterDensity, bg::Background)
    return compute_Cℓ(Component1, Component2, W.w_minus2_00_ϕTT, bg)
end

function get_Cℓ(Component1::PrimordialNonGaussianity, Component2::CosmicShear, W::ProjectedMatterDensity, bg::Background)
    return compute_Cℓ(Component1, Component2, W.w_0_00_ϕT, bg)
end

function get_Cℓ(Component1::PrimordialNonGaussianity, Component2::IntrinsicAlignment, W::ProjectedMatterDensity, bg::Background)
    return compute_Cℓ(Component1, Component2, W.w_0_00_ϕT, bg)
end

function get_Cℓ(ℓ::AbstractArray{<:Any, 1}, G::GalaxyClustering, L::WeakLensing, Pk::PowerSpectrum, W::ProjectedMatterDensity, bg::Background, P::Union{FFTPlans, Nothing}=nothing)
    Cℓ_δγ = get_Cℓ(G.δ, L.γ, W, bg)
    Cℓ_δI = get_Cℓ(G.δ, L.IA, W, bg)
    Cℓ_RSDγ = get_Cℓ(G.RSD, L.γ, W, bg)
    Cℓ_RSDI = get_Cℓ(G.RSD, L.IA, W, bg)
    Cℓ_μγ = get_Cℓ(G.μ, L.γ, W, bg)
    Cℓ_μI = get_Cℓ(G.μ, L.IA, W, bg)
    Cℓ_fNLγ = get_Cℓ(G.PNG, L.γ, W, bg)
    Cℓ_fNLI = get_Cℓ(G.PNG, L.IA, W, bg)

    Cℓ_nonlimber = @. Cℓ_δγ + Cℓ_δI - Cℓ_RSDγ - Cℓ_RSDI + Cℓ_μγ + Cℓ_μI + Cℓ_fNLγ + Cℓ_fNLI
    Cℓ_correction = get_limber_correction(G, L, Pk)
    Cℓ_limber = get_limber_Cℓ(G, L, Pk)
    full_Cℓ = cat(Cℓ_nonlimber .+ Cℓ_correction, Cℓ_limber; dims=1) 

    nbins_A = size(G.δ.Kernel, 1)
    nbins_B = size(L.γ.Kernel, 1)
    return _finalize_Cℓ(full_Cℓ, ℓ, nbins_A, nbins_B, P)
end

function get_Cℓ(ℓ::AbstractArray{<:Any, 1}, L::WeakLensing, G::GalaxyClustering, Pk::PowerSpectrum, W::ProjectedMatterDensity, bg::Background, P::Union{FFTPlans, Nothing}=nothing)
    return get_Cℓ(ℓ, G, L, Pk, W, bg, P)
end

## Cross clustering - CMB 
function get_Cℓ(Component1::CMBLensing, Component2::NumberCounts, W::ProjectedMatterDensity, bg::Background)
    return compute_Cℓ(Component1, Component2, W.w_0_00_ϕTT, bg)
end

function get_Cℓ(Component1::CMBLensing, Component2::RedshiftSpaceDistortions, W::ProjectedMatterDensity, bg::Background)
    return compute_Cℓ(Component1, Component2, W.w_0_02_ϕTT, W.w_0_20_ϕTT, bg)
end

function get_Cℓ(Component1::CMBLensing, Component2::MagnificationBias, W::ProjectedMatterDensity, bg::Background)
    return compute_Cℓ(Component1, Component2, W.w_minus2_00_ϕTT, bg)
end

function get_Cℓ(Component1::CMBLensing, Component2::PrimordialNonGaussianity, W::ProjectedMatterDensity, bg::Background)
    return compute_Cℓ(Component1, Component2, W.w_0_00_ϕT_R1, bg)
end

function get_Cℓ(Component1::IntegratedSachsWolfe, Component2::NumberCounts, W::ProjectedMatterDensity, bg::Background)
    return compute_Cℓ(Component1, Component2, W.w_0_00_ϕTT, bg)
end

function get_Cℓ(Component1::IntegratedSachsWolfe, Component2::RedshiftSpaceDistortions, W::ProjectedMatterDensity, bg::Background)
    return compute_Cℓ(Component1, Component2, W.w_0_02_ϕTT, W.w_0_20_ϕTT, bg)
end

function get_Cℓ(Component1::IntegratedSachsWolfe, Component2::MagnificationBias, W::ProjectedMatterDensity, bg::Background)
    return compute_Cℓ(Component1, Component2, W.w_minus2_00_ϕTT, bg)
end

function get_Cℓ(Component1::IntegratedSachsWolfe, Component2::PrimordialNonGaussianity, W::ProjectedMatterDensity, bg::Background)
    return compute_Cℓ(Component1, Component2, W.w_0_00_ϕT_R1, bg)
end

function get_Cℓ(ℓ::AbstractArray{<:Any, 1}, K::CMB, G::GalaxyClustering, Pk::PowerSpectrum, W::ProjectedMatterDensity, bg::Background, P::Union{FFTPlans, Nothing}=nothing)
    Cℓ_κδ = get_Cℓ(K.κ, G.δ, W, bg)
    Cℓ_κRSD = get_Cℓ(K.κ, G.RSD, W, bg)
    Cℓ_κμ = get_Cℓ(K.κ, G.μ, W, bg)
    Cℓ_κfNL = get_Cℓ(K.κ, G.PNG, W, bg)
    Cℓ_Tδ = get_Cℓ(K.ISW, G.δ, W, bg)
    Cℓ_TRSD = get_Cℓ(K.ISW, G.RSD, W, bg)
    Cℓ_Tμ = get_Cℓ(K.ISW, G.μ, W, bg)
    Cℓ_TfNL = get_Cℓ(K.ISW, G.PNG, W, bg)

    Cℓ_nonlimber = @. Cℓ_κδ - Cℓ_κRSD + Cℓ_κμ + Cℓ_κfNL + Cℓ_Tδ - Cℓ_TRSD + Cℓ_Tμ + Cℓ_TfNL
    Cℓ_correction = get_limber_correction(K, G, Pk)
    Cℓ_limber = get_limber_Cℓ(K, G, Pk)
    full_Cℓ = cat(Cℓ_nonlimber .+ Cℓ_correction, Cℓ_limber; dims=1) 

    nbins_A = size(K.κ.Kernel, 1)
    nbins_B = size(G.δ.Kernel, 1)
    return _finalize_Cℓ(full_Cℓ, ℓ, nbins_A, nbins_B, P)
end

function get_Cℓ(ℓ::AbstractArray{<:Any, 1}, G::GalaxyClustering, K::CMB, Pk::PowerSpectrum, W::ProjectedMatterDensity, bg::Background, P::Union{FFTPlans, Nothing}=nothing)
    return get_Cℓ(ℓ, K, G, Pk, W, bg, P)
end

## Cross lensing - CMB
function get_Cℓ(Component1::CMBLensing, Component2::CosmicShear, W::ProjectedMatterDensity, bg::Background)
    return compute_Cℓ(Component1, Component2, W.w_minus2_00_ϕTT, bg)
end

function get_Cℓ(Component1::CMBLensing, Component2::IntrinsicAlignment, W::ProjectedMatterDensity, bg::Background)
    return compute_Cℓ(Component1, Component2, W.w_minus2_00_ϕTT, bg)
end

function get_Cℓ(Component1::IntegratedSachsWolfe, Component2::CosmicShear, W::ProjectedMatterDensity, bg::Background)
    return compute_Cℓ(Component1, Component2, W.w_minus2_00_ϕTT, bg)
end

function get_Cℓ(Component1::IntegratedSachsWolfe, Component2::IntrinsicAlignment, W::ProjectedMatterDensity, bg::Background)
    return compute_Cℓ(Component1, Component2, W.w_minus2_00_ϕTT, bg)
end

function get_Cℓ(ℓ::AbstractArray{<:Any, 1}, K::CMB, L::WeakLensing, Pk::PowerSpectrum, W::ProjectedMatterDensity, bg::Background, P::Union{FFTPlans, Nothing}=nothing)
    Cℓ_κγ = get_Cℓ(K.κ, L.γ, W, bg)
    Cℓ_κI = get_Cℓ(K.κ, L.IA, W, bg)
    Cℓ_Tγ = get_Cℓ(K.ISW, L.γ, W, bg)
    Cℓ_TI = get_Cℓ(K.ISW, L.IA, W, bg)

    Cℓ_nonlimber = @. Cℓ_κγ + Cℓ_κI + Cℓ_Tγ + Cℓ_TI
    Cℓ_correction = get_limber_correction(K, L, Pk)
    Cℓ_limber = get_limber_Cℓ(K, L, Pk)
    full_Cℓ = cat(Cℓ_nonlimber .+ Cℓ_correction, Cℓ_limber; dims=1) 

    nbins_A = size(K.κ.Kernel, 1)
    nbins_B = size(L.γ.Kernel, 1)
    return _finalize_Cℓ(full_Cℓ, ℓ, nbins_A, nbins_B, P)
end

function get_Cℓ(ℓ::AbstractArray{<:Any, 1}, S::WeakLensing, K::CMB, Pk::PowerSpectrum, W::ProjectedMatterDensity, bg::Background, P::Union{FFTPlans, Nothing}=nothing)
    return get_Cℓ(ℓ, K, S, Pk, W, bg, P)
end

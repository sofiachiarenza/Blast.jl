"""
    FFTPlans

Container for pre-allocated FFTW plans used to decompose the power spectrum on the Chebyshev basis.

These plans are reused across power spectrum evaluations to avoid repeated FFTW planning overhead.

# Fields
- `plan_ϕTT`: FFT plan for the unequal-time power spectrum `P(k, χ₁, χ₂)`.
- `plan_ϕT`: Optional FFT plan for the power spectrum `P(k, χ₁)` built with a single transfer function.
- `plan_ϕ`: Optional FFT plan for the primordial power spectrum.
"""
@kwdef mutable struct FFTPlans 
    plan_ϕTT::FFTW.r2rFFTWPlan 
    plan_ϕT::Union{FFTW.r2rFFTWPlan, Nothing} = nothing
    plan_ϕ::Union{FFTW.r2rFFTWPlan, Nothing} = nothing
end

"""
    SetUp(probes...)

Construct the projected matter density containers and FFT plans required
for the given set of cosmological probes.

The returned objects depend on which physical effects are active in the probes:
- Galaxy clustering: density, RSD, magnification, PNG
- Weak lensing: shear and intrinsic alignment
- CMB: lensing convergence and ISW

Only the required projected matter components are instantiated; all others
are set to `nothing` avoiding useless overhead.

# Arguments
- `probes...`: Any combination of `GalaxyClustering`, `WeakLensing`, and `CMB`.

# Returns
- `W::ProjectedMatterDensity`: Container holding all required  inner k integrals.
- `P::FFTPlans`: Pre-allocated FFT plans used for Chebyshev decomposition.

# Notes
- The ordering of probes does not matter.
- Conceptually, this function computes and stores everything that is only needed once.
"""
function SetUp(G::GalaxyClustering)
    plan_ϕTT = plan_fft(randn(size(Blast.k_cheb, 1), size(Blast.χ, 1), size(Blast.R, 1)),1)
    plan_ϕT = nothing
    plan_ϕ = nothing

    w_δ = w_2_00_ϕTT()
    w_μ_B = nothing
    w_μ_A = nothing
    w_μxRSD_A = nothing
    w_μxRSD_B = nothing
    w_RSD_A = nothing
    w_RSD_B = nothing
    w_RSD_C = nothing
    w_PNG_A = nothing
    w_PNG_B = nothing
    w_μxPNG_A = nothing
    w_μxPNG_B = nothing
    w_RSDxPNG_A = nothing
    w_RSDxPNG_C = nothing
    w_RSDxPNG_B = nothing
    w_RSDxPNG_D = nothing
    w_PNG_C = nothing

    if !isnothing(G.RSD)
        w_RSD_A = w_2_02_ϕTT()
        w_RSD_B = w_2_20_ϕTT()
        w_RSD_C = w_2_22_ϕTT()
    end

    if !isnothing(G.μ)
        w_μ_A = w_0_00_ϕTT()
        w_μ_B = w_minus2_00_ϕTT()
    end

    if !isnothing(G.PNG)
        plan_ϕT = plan_ϕTT
        plan_ϕ = plan_fft(randn(size(Blast.k_cheb, 1)),1)
        w_PNG_A = w_2_00_ϕT()
        w_PNG_B = w_2_00_ϕT_R1()
        w_PNG_C = w_2_00_ϕ()
    end

    if !isnothing(G.RSD) && !isnothing(G.μ)
        w_μxRSD_A = w_0_02_ϕTT()
        w_μxRSD_B = w_0_20_ϕTT()
    end

    if !isnothing(G.RSD) && !isnothing(G.PNG)
        w_RSDxPNG_A = w_2_02_ϕT()
        w_RSDxPNG_B = w_2_20_ϕT()
        w_RSDxPNG_C = w_2_02_ϕT_R1()
        w_RSDxPNG_D = w_2_20_ϕT_R1()
    end

    if !isnothing(G.μ) && !isnothing(G.PNG)
        w_μxPNG_A = w_0_00_ϕT()
        w_μxPNG_B = w_0_00_ϕT_R1()
    end
    
    Plans = FFTPlans(plan_ϕTT, plan_ϕT, plan_ϕ)
    W = ProjectedMatterDensity(w_δ, w_μ_B, w_μ_A, w_μxRSD_A, w_μxRSD_B, w_RSD_A, w_RSD_B, w_RSD_C, w_PNG_A, w_PNG_B, 
                                w_μxPNG_A, w_μxPNG_B, w_RSDxPNG_A, w_RSDxPNG_C, w_RSDxPNG_B, w_RSDxPNG_D, w_PNG_C)
    return W, Plans
end

function SetUp(L::WeakLensing)
    
    plan_ϕTT = plan_fft(randn(size(Blast.k_cheb, 1), size(Blast.χ, 1), size(Blast.R, 1)),1)
    plan_ϕT = nothing
    plan_ϕ = nothing

    w_γ = w_minus2_00_ϕTT()    
    
    Plans = FFTPlans(plan_ϕTT, plan_ϕT, plan_ϕ)
    W = ProjectedMatterDensity(nothing, w_γ, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing)
    return W, Plans
end

function SetUp(G::GalaxyClustering, L::WeakLensing)
    
    plan_ϕTT = plan_fft(randn(size(Blast.k_cheb, 1), size(Blast.χ, 1), size(Blast.R, 1)),1)
    plan_ϕT = nothing
    plan_ϕ = nothing

    w_δ = w_2_00_ϕTT()
    w_μ_B = w_minus2_00_ϕTT()
    w_μ_A = w_0_00_ϕTT()
    w_μxRSD_A = nothing
    w_μxRSD_B = nothing
    w_RSD_A = nothing
    w_RSD_B = nothing
    w_RSD_C = nothing
    w_PNG_A = nothing
    w_PNG_B = nothing
    w_μxPNG_A = nothing
    w_μxPNG_B = nothing
    w_RSDxPNG_A = nothing
    w_RSDxPNG_C = nothing
    w_RSDxPNG_B = nothing
    w_RSDxPNG_D = nothing
    w_PNG_C = nothing

    if !isnothing(G.RSD)
        w_RSD_A = w_2_02_ϕTT()
        w_RSD_B = w_2_20_ϕTT()
        w_RSD_C = w_2_22_ϕTT()
    end

    if !isnothing(G.PNG)
        plan_ϕT = plan_ϕTT
        plan_ϕ = plan_fft(randn(size(Blast.k_cheb, 1)),1)
        w_PNG_A = w_2_00_ϕT()
        w_PNG_B = w_2_00_ϕT_R1()
        w_PNG_C = w_2_00_ϕ()
    end

    if !isnothing(G.RSD) && (!isnothing(G.μ) || !isnothing(L.γ))
        w_μxRSD_A = w_0_02_ϕTT()
        w_μxRSD_B = w_0_20_ϕTT()
    end

    if !isnothing(G.RSD) && !isnothing(G.PNG)
        w_RSDxPNG_A = w_2_02_ϕT()
        w_RSDxPNG_B = w_2_20_ϕT()
        w_RSDxPNG_C = w_2_02_ϕT_R1()
        w_RSDxPNG_D = w_2_20_ϕT_R1()
    end

    if (!isnothing(G.μ) || !isnothing(L.γ)) && !isnothing(G.PNG)
        w_μxPNG_A = w_0_00_ϕT()
        w_μxPNG_B = w_0_00_ϕT_R1()
    end
    
    Plans = FFTPlans(plan_ϕTT, plan_ϕT, plan_ϕ)
    W = ProjectedMatterDensity(w_δ, w_μ_B, w_μ_A, w_μxRSD_A, w_μxRSD_B, w_RSD_A, w_RSD_B, w_RSD_C, w_PNG_A, w_PNG_B, 
                                w_μxPNG_A, w_μxPNG_B, w_RSDxPNG_A, w_RSDxPNG_C, w_RSDxPNG_B, w_RSDxPNG_D, w_PNG_C)
    return W, Plans
end

function SetUp(L::WeakLensing, G::GalaxyClustering)
    return SetUp(G, L)
end

function SetUp(G::GalaxyClustering, C::CMB)
    
    plan_ϕTT = plan_fft(randn(size(Blast.k_cheb, 1), size(Blast.χ, 1), size(Blast.R, 1)),1)
    plan_ϕT = nothing
    plan_ϕ = nothing

    w_δ = w_2_00_ϕTT()
    w_μ_B = nothing
    w_μ_A = w_0_00_ϕTT()
    w_μxRSD_A = nothing
    w_μxRSD_B = nothing
    w_RSD_A = nothing
    w_RSD_B = nothing
    w_RSD_C = nothing
    w_PNG_A = nothing
    w_PNG_B = nothing
    w_μxPNG_A = nothing
    w_μxPNG_B = nothing
    w_RSDxPNG_A = nothing
    w_RSDxPNG_C = nothing
    w_RSDxPNG_B = nothing
    w_RSDxPNG_D = nothing
    w_PNG_C = nothing

    if !isnothing(G.RSD)
        w_RSD_A = w_2_02_ϕTT()
        w_RSD_B = w_2_20_ϕTT()
        w_RSD_C = w_2_22_ϕTT()
    end

    if !isnothing(G.μ)
        w_μ_B = w_minus2_00_ϕTT()
    end

    if !isnothing(G.PNG)
        plan_ϕT = plan_ϕTT
        plan_ϕ = plan_fft(randn(size(Blast.k_cheb, 1)),1)
        w_PNG_A = w_2_00_ϕT()
        w_PNG_B = w_2_00_ϕT_R1()
        w_PNG_C = w_2_00_ϕ()
    end

    if !isnothing(G.RSD) && (!isnothing(G.μ) || !isnothing(C.κ))
        w_μxRSD_A = w_0_02_ϕTT()
        w_μxRSD_B = w_0_20_ϕTT()
    end

    if !isnothing(G.RSD) && !isnothing(G.PNG)
        w_RSDxPNG_A = w_2_02_ϕT()
        w_RSDxPNG_B = w_2_20_ϕT()
        w_RSDxPNG_C = w_2_02_ϕT_R1()
        w_RSDxPNG_D = w_2_20_ϕT_R1()
    end

    if (!isnothing(G.μ) || !isnothing(C.κ)) && !isnothing(G.PNG)
        w_μxPNG_A = w_0_00_ϕT()
        w_μxPNG_B = w_0_00_ϕT_R1()
    end
    
    Plans = FFTPlans(plan_ϕTT, plan_ϕT, plan_ϕ)
    W = ProjectedMatterDensity(w_δ, w_μ_B, w_μ_A, w_μxRSD_A, w_μxRSD_B, w_RSD_A, w_RSD_B, w_RSD_C, w_PNG_A, w_PNG_B, 
                                w_μxPNG_A, w_μxPNG_B, w_RSDxPNG_A, w_RSDxPNG_C, w_RSDxPNG_B, w_RSDxPNG_D, w_PNG_C)
    return W, Plans
end

function SetUp(C::CMB, G::GalaxyClustering)
    return SetUp(G, C)
end

function SetUp(L::WeakLensing, C::CMB)
    
    plan_ϕTT = plan_fft(randn(size(Blast.k_cheb, 1), size(Blast.χ, 1), size(Blast.R, 1)),1)
    plan_ϕT = nothing
    plan_ϕ = nothing

    w_γ = w_minus2_00_ϕTT()    
    
    Plans = FFTPlans(plan_ϕTT, plan_ϕT, plan_ϕ)
    W = ProjectedMatterDensity(nothing, w_γ, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing)
    return W, Plans
end

function SetUp(C::CMB, L::WeakLensing)
    return SetUp(L, C)
end

function SetUp(G::GalaxyClustering, L::WeakLensing, C::CMB)
    
    plan_ϕTT = plan_fft(randn(size(Blast.k_cheb, 1), size(Blast.χ, 1), size(Blast.R, 1)),1)
    plan_ϕT = nothing
    plan_ϕ = nothing

    w_δ = w_2_00_ϕTT()
    w_μ_B = w_minus2_00_ϕTT()
    w_μ_A = w_0_00_ϕTT()
    w_μxRSD_A = nothing
    w_μxRSD_B = nothing
    w_RSD_A = nothing
    w_RSD_B = nothing
    w_RSD_C = nothing
    w_PNG_A = nothing
    w_PNG_B = nothing
    w_μxPNG_A = nothing
    w_μxPNG_B = nothing
    w_RSDxPNG_A = nothing
    w_RSDxPNG_C = nothing
    w_RSDxPNG_B = nothing
    w_RSDxPNG_D = nothing
    w_PNG_C = nothing

    if !isnothing(G.RSD)
        w_RSD_A = w_2_02_ϕTT()
        w_RSD_B = w_2_20_ϕTT()
        w_RSD_C = w_2_22_ϕTT()
    end

    if !isnothing(G.μ)
        w_μ_A = w_0_00_ϕTT()
        w_μ_B = w_minus2_00_ϕTT()
    end

    if !isnothing(G.PNG)
        plan_ϕT = plan_ϕTT
        plan_ϕ = plan_fft(randn(size(Blast.k_cheb, 1)),1)
        w_PNG_A = w_2_00_ϕT()
        w_PNG_B = w_2_00_ϕT_R1()
        w_PNG_C = w_2_00_ϕ()
    end

    if !isnothing(G.RSD) && (!isnothing(G.μ) || !isnothing(L.γ))
        w_μxRSD_A = w_0_02_ϕTT()
        w_μxRSD_B = w_0_20_ϕTT()
    end

    if !isnothing(G.RSD) && !isnothing(G.PNG)
        w_RSDxPNG_A = w_2_02_ϕT()
        w_RSDxPNG_B = w_2_20_ϕT()
        w_RSDxPNG_C = w_2_02_ϕT_R1()
        w_RSDxPNG_D = w_2_20_ϕT_R1()
    end

    if (!isnothing(G.μ) || !isnothing(L.γ)) && !isnothing(G.PNG)
        w_μxPNG_A = w_0_00_ϕT()
        w_μxPNG_B = w_0_00_ϕT_R1()
    end
    
    Plans = FFTPlans(plan_ϕTT, plan_ϕT, plan_ϕ)
    W = ProjectedMatterDensity(w_δ, w_μ_B, w_μ_A, w_μxRSD_A, w_μxRSD_B, w_RSD_A, w_RSD_B, w_RSD_C, w_PNG_A, w_PNG_B, 
                                w_μxPNG_A, w_μxPNG_B, w_RSDxPNG_A, w_RSDxPNG_C, w_RSDxPNG_B, w_RSDxPNG_D, w_PNG_C)
    return W, Plans
end

function SetUp(G::GalaxyClustering, C::CMB, L::WeakLensing)
    return SetUp(G,L,C)
end

function SetUp(L::WeakLensing, G::GalaxyClustering, C::CMB)
    return SetUp(G,L,C)
end

function SetUp(L::WeakLensing, C::CMB, G::GalaxyClustering)
    return SetUp(G,L,C)
end

function SetUp(C::CMB, G::GalaxyClustering, L::WeakLensing)
    return SetUp(G,L,C)
end

function SetUp(C::CMB, L::WeakLensing, G::GalaxyClustering)
    return SetUp(G,L,C)
end

"""
    PowerSpectrum

A mutable struct acting as a container for all processed spectral data required for the angular power spectrum measurements.

This struct caches both the Chebyshev coefficients used in the high-precision 
non-Limber integration and the pre-interpolated power spectra used for 
the Limber approximation.

# Fields
- `cϕTT`, `cϕT`, `cϕ`: Structs storing all the necessary Chebyshev coefficients.
- `ΔP_limber`: A 2D array of the non-linear correction to the matter power spectrum.  
- `Pδ_limber`: A 2D array of the total matter power spectrum.
"""
@kwdef mutable struct PowerSpectrum 
    cϕTT::cϕTT
    cϕT::Union{cϕT, Nothing}
    cϕ::Union{cϕ, Nothing}
    ΔP_limber::AbstractArray{<:Any, 2} 
    Pδ_limber::AbstractArray{<:Any, 2} 
end

"""
    get_PΦ(k::AbstractArray, cosmo::AbstractCosmology)

Computes the primordial power spectrum \$P_\\Phi(k)\$.

The implementation follows the standard inflationary prediction:
```math
P_\\Phi(k) = \\frac{9}{25} \\frac{2\\pi^2 A_s}{k^3} \\left( \\frac{k}{k_{pivot}}
```
where the pivot scale is fixed at `0.05 Mpc^{-1}``.

# Returns
- A 1D array of the primordial potential power spectrum at the requested scales. 
"""
function get_PΦ(k::AbstractArray{<:Any,1}, cosmo::AbstractCosmology)
    return @. 9/25 * 2 * π^2 * cosmo.As / (k^3) * (k/0.05)^(cosmo.ns - 1.)
end

""" get_Tm(pk::AbstractArray{<:Any, 2}, k::AbstractArray, cosmo::AbstractCosmology)

Extracts the matter transfer function T_m(k,z) from a matter power spectrum.

The transfer function is defined as the square root of the ratio between the processed matter power spectrum and the primordial potential power spectrum:
```math
T_m(k, z) = \\sqrt{\\frac{P(k, z)}{P_\\Phi(k)}}
```
"""
function get_Tm(pk::AbstractArray{<:Any,2}, k::AbstractArray{<:Any, 1}, cosmo::AbstractCosmology)
    prim_pk = get_PΦ(k , cosmo)
    return sqrt.(pk ./ reshape(prim_pk, 1, :))
end

"""
    transform_to_R_frame(matrix::AbstractArray{<:Any,2}, bg::BackgroundQuantities, grid::AbstractCosmologicalGrid)

Transforms a 2D array (typically a transfer function) from standard coordinates `(χ, k)` 
to a ratio-based coordinate system `(χ, Rχ, k)`, where `0<R<1`. This coordinate change is a specialized optimization for the non-Limber integration.

# Returns
- A 3D array of shape to `(nχ, nR, nk)`.
"""
function transform_to_R_frame(matrix::AbstractArray{<:Any,2}, bg::BackgroundQuantities, grid::AbstractCosmologicalGrid)
    new_χs = make_grid(bg, R)
    x = resample_redshifts(bg, grid, new_χs)
    interp = Blast._akima_interpolation(matrix, Blast.z_lin, x)  
    return reshape( interp,  size(bg.χz_array, 1), size(R,1), size(interp,2))
end

#TODO: this function will interface with mapse once ready.
"""
    prepare_pk_workspace(P::FFTPlans, 
                         pk::AbstractArray{<:Any, 2}, pk_limber_lin::AbstractArray{<:Any, 2}, 
                         pk_limber_nonlin::AbstractArray{<:Any, 2}, 
                         cosmo::AbstractCosmology, bg::BackgroundQuantities, grid::AbstractCosmologicalGrid)

Assembles the `PowerSpectrum` struct by processing primordial, linear, and 
non-linear power spectra into the format needed for `C_ℓ` calculation.

The function performs the following operations:
1. Decomposes the matter power spectrum into the primordial power spectrum `P_Φ`
   and the transfer function `T_m(k, χ)`, which are used to define the unequal time power spectrum `P(k, χ1, χ2)`.
2. Perform change of coordinates to use a more optimal integration basis.
3. Computes 3D Chebyshev coefficients for non-Limber integration.
4. Pre-interpolates the Limber-limit power spectra `P(k = \\frac{(ℓ+1/2)}{χ}, χ)`` on the correct grid.

# Parameters:
- `P`: Pre-allocated `FFTPlans` for Chebyshev transforms.
- `pk`: Linear matter power spectrum for the non-Limber integration.
- `pk_limber_lin/nonlin`: Linear and non linear matter power spectra for the Limber integration.
"""
function prepare_pk_workspace(P::FFTPlans, pk::AbstractArray{<:Any, 2}, pk_limber_lin::AbstractArray{<:Any, 2}, pk_limber_nonlin::AbstractArray{<:Any, 2}, cosmo::AbstractCosmology, bg::BackgroundQuantities, grid::AbstractCosmologicalGrid)
    #Treating the non-limber power spectrum
    P_ϕ = get_PΦ(10 .^ Blast.k_cheb, cosmo)
    transfer_func = get_Tm(pk, 10 .^ Blast.k_cheb, cosmo)
    
    transfer_func_χR = transform_to_R_frame(transfer_func, bg, grid)
    transfer_func_χ1 = transfer_func_χR[:,end,:]

    @tullio P_ϕTT[k, i, j] := P_ϕ[k] * transfer_func_χ1[i,k] * transfer_func_χR[i, j, k] 
    @tullio P_ϕT[k, i, j] := P_ϕ[k]* transfer_func_χR[i, j, k]

    c1 = build_coeff(cϕTT, P_ϕTT, P.plan_ϕTT)
    c2 = build_coeff(cϕT,  P_ϕT,  P.plan_ϕT)   
    c3 = build_coeff(cϕ,   P_ϕ,   P.plan_ϕ)    

    lb, ub = [minimum(Blast.z_cheb),minimum(Blast.k_limber)], [maximum(Blast.z_cheb), maximum(Blast.k_limber)] # lower and upper bounds of the domain, respectively

    P_ϕ_limber = get_PΦ(10 .^ Blast.k_limber, cosmo)'
    T_m_limber_lin = get_Tm(pk_limber_lin, 10 .^ Blast.k_limber, cosmo)
    P_limber_lin = P_ϕ_limber .* T_m_limber_lin .^ 2.
    limber_pk_linear = chebinterp(log10.(P_limber_lin), lb, ub)
    
    T_m_limber_nonlin = get_Tm(pk_limber_nonlin, 10 .^ Blast.k_limber, cosmo)
    P_limber_nonlin = P_ϕ_limber .* T_m_limber_nonlin .^ 2.
    limber_pk_nonlinear = chebinterp(log10.(P_limber_nonlin), lb, ub)

    #TODO: handle this better, having this here sucks
    z_of_χ = AkimaInterpolation(grid.z_range, bg.χz_array, extrapolation=ExtrapolationType.Extension)

    ℓ_grid = reshape(Blast.full_ℓ_range .+ 0.5, :, 1)   
    χ_grid = reshape(Blast.χ, 1, :)          
    k = ℓ_grid ./ χ_grid               
    ΔP_limber = [10. ^ limber_pk_nonlinear(SVector(z_of_χ.(χ)[j], log10(k[i, j]))) .- 10. ^ limber_pk_linear(SVector(z_of_χ.(χ)[j], log10(k[i, j]))) for i in 1:size(ℓ_grid, 1), j in 1:size(Blast.χ, 1)]
    Pδ_limber = [10. ^ limber_pk_nonlinear(SVector(z_of_χ.(χ)[j], log10(k[i, j]))) for i in 1:size(ℓ_grid, 1), j in 1:size(Blast.χ, 1)]

    return PowerSpectrum(c1, c2, c3,  ΔP_limber, Pδ_limber)
end
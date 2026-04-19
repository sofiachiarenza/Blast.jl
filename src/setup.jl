"""
    FFTPlans

Container for pre-allocated FFTW plans used to decompose the power spectrum on the Chebyshev basis.

These plans are reused across power spectrum evaluations to avoid repeated FFTW planning overhead.

# Fields
- `plan_ϕTT`: FFT plan for the unequal-time power spectrum `P(k, χ₁, χ₂)`.
- `plan_ϕT`: Optional FFT plan for the power spectrum `P(k, χ₁)` built with a single transfer function.
- `plan_ϕ`: Optional FFT plan for the primordial power spectrum.
- `plan_limber`: FFT plan for the Limber power spectrum.
- `T_k_limber`: Precomputed k-basis polynomials for the Limber grid.
- `plan_ℓ`: Plan for the final C_ℓ interpolation.
"""
@kwdef mutable struct FFTPlans 
    plan_ϕTT::ChebyshevPlan 
    plan_ϕT::Union{ChebyshevPlan, Nothing} = nothing
    plan_ϕ::Union{ChebyshevPlan, Nothing} = nothing
    plan_limber::ChebyshevPlan
    T_k_limber::Array{Float64, 3}
    plan_ℓ::ChebyshevPlan
end

function _setup_limber_plan()
    lk_min, lk_max = minimum(k_limber), maximum(k_limber)
    z_min, z_max = minimum(z_cheb), maximum(z_cheb)
    K_k, K_z = length(k_limber) - 1, length(z_cheb) - 1

    plan_limber = prepare_chebyshev_plan((lk_min, z_min), (lk_max, z_max), (K_k, K_z))
    
    T_k = get_limber_k_polynomials(plan_limber, Blast.full_ℓ_range, Blast.χ; is_log_k=true)
    
    plan_ℓ = prepare_chebyshev_plan(2, 2000, 100)
    
    return plan_limber, T_k, plan_ℓ
end

"""
    SetUp(probes...)

Construct the projected matter density containers and FFT plans required
for the given set of cosmological probes.

# Returns
- `W::ProjectedMatterDensity`: Container holding all required inner k integrals.
- `P::FFTPlans`: Pre-allocated FFT plans used for Chebyshev decomposition.
"""
function SetUp(G::GalaxyClustering)
    plan_ϕTT = prepare_chebyshev_plan(minimum(k_cheb), maximum(k_cheb), length(k_cheb) - 1; size_nd=(size(Blast.k_cheb, 1), size(Blast.χ, 1), size(Blast.R, 1)), dim=1)
    plan_ϕT = nothing
    plan_ϕ = nothing
    plan_limber, T_k, plan_ℓ = _setup_limber_plan()

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
        plan_ϕ = prepare_chebyshev_plan(minimum(k_cheb), maximum(k_cheb), length(k_cheb) - 1; size_nd=(size(Blast.k_cheb, 1),), dim=1)
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
    
    Plans = FFTPlans(plan_ϕTT, plan_ϕT, plan_ϕ, plan_limber, T_k, plan_ℓ)
    W = ProjectedMatterDensity(w_δ, w_μ_B, w_μ_A, w_μxRSD_A, w_μxRSD_B, w_RSD_A, w_RSD_B, w_RSD_C, w_PNG_A, w_PNG_B, 
                                w_μxPNG_A, w_μxPNG_B, w_RSDxPNG_A, w_RSDxPNG_C, w_RSDxPNG_B, w_RSDxPNG_D, w_PNG_C)
    return W, Plans
end

function SetUp(L::WeakLensing)
    
    plan_ϕTT = prepare_chebyshev_plan(minimum(k_cheb), maximum(k_cheb), length(k_cheb) - 1; size_nd=(size(Blast.k_cheb, 1), size(Blast.χ, 1), size(Blast.R, 1)), dim=1)
    plan_ϕT = nothing
    plan_ϕ = nothing
    plan_limber, T_k, plan_ℓ = _setup_limber_plan()

    w_γ = w_minus2_00_ϕTT()    
    
    Plans = FFTPlans(plan_ϕTT, plan_ϕT, plan_ϕ, plan_limber, T_k, plan_ℓ)
    W = ProjectedMatterDensity(nothing, w_γ, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing)
    return W, Plans
end

function SetUp(G::GalaxyClustering, L::WeakLensing)
    
    plan_ϕTT = prepare_chebyshev_plan(minimum(k_cheb), maximum(k_cheb), length(k_cheb) - 1; size_nd=(size(Blast.k_cheb, 1), size(Blast.χ, 1), size(Blast.R, 1)), dim=1)
    plan_ϕT = nothing
    plan_ϕ = nothing
    plan_limber, T_k, plan_ℓ = _setup_limber_plan()

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
        plan_ϕ = prepare_chebyshev_plan(minimum(k_cheb), maximum(k_cheb), length(k_cheb) - 1; size_nd=(size(Blast.k_cheb, 1),), dim=1)
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
    
    Plans = FFTPlans(plan_ϕTT, plan_ϕT, plan_ϕ, plan_limber, T_k, plan_ℓ)
    W = ProjectedMatterDensity(w_δ, w_μ_B, w_μ_A, w_μxRSD_A, w_μxRSD_B, w_RSD_A, w_RSD_B, w_RSD_C, w_PNG_A, w_PNG_B, 
                                w_μxPNG_A, w_μxPNG_B, w_RSDxPNG_A, w_RSDxPNG_C, w_RSDxPNG_B, w_RSDxPNG_D, w_PNG_C)
    return W, Plans
end

function SetUp(L::WeakLensing, G::GalaxyClustering)
    return SetUp(G, L)
end

function SetUp(G::GalaxyClustering, C::CMB)
    
    plan_ϕTT = prepare_chebyshev_plan(minimum(k_cheb), maximum(k_cheb), length(k_cheb) - 1; size_nd=(size(Blast.k_cheb, 1), size(Blast.χ, 1), size(Blast.R, 1)), dim=1)
    plan_ϕT = nothing
    plan_ϕ = nothing
    plan_limber, T_k, plan_ℓ = _setup_limber_plan()

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
        plan_ϕ = prepare_chebyshev_plan(minimum(k_cheb), maximum(k_cheb), length(k_cheb) - 1; size_nd=(size(Blast.k_cheb, 1),), dim=1)
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
    
    Plans = FFTPlans(plan_ϕTT, plan_ϕT, plan_ϕ, plan_limber, T_k, plan_ℓ)
    W = ProjectedMatterDensity(w_δ, w_μ_B, w_μ_A, w_μxRSD_A, w_μxRSD_B, w_RSD_A, w_RSD_B, w_RSD_C, w_PNG_A, w_PNG_B, 
                                w_μxPNG_A, w_μxPNG_B, w_RSDxPNG_A, w_RSDxPNG_C, w_RSDxPNG_B, w_RSDxPNG_D, w_PNG_C)
    return W, Plans
end

function SetUp(C::CMB, G::GalaxyClustering)
    return SetUp(G, C)
end

function SetUp(L::WeakLensing, C::CMB)
    
    plan_ϕTT = prepare_chebyshev_plan(minimum(k_cheb), maximum(k_cheb), length(k_cheb) - 1; size_nd=(size(Blast.k_cheb, 1), size(Blast.χ, 1), size(Blast.R, 1)), dim=1)
    plan_ϕT = nothing
    plan_ϕ = nothing
    plan_limber, T_k, plan_ℓ = _setup_limber_plan()

    w_γ = w_minus2_00_ϕTT()    
    
    Plans = FFTPlans(plan_ϕTT, plan_ϕT, plan_ϕ, plan_limber, T_k, plan_ℓ)
    W = ProjectedMatterDensity(nothing, w_γ, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing)
    return W, Plans
end

function SetUp(C::CMB, L::WeakLensing)
    return SetUp(L, C)
end

function SetUp(G::GalaxyClustering, L::WeakLensing, C::CMB)
    
    plan_ϕTT = prepare_chebyshev_plan(minimum(k_cheb), maximum(k_cheb), length(k_cheb) - 1; size_nd=(size(Blast.k_cheb, 1), size(Blast.χ, 1), size(Blast.R, 1)), dim=1)
    plan_ϕT = nothing
    plan_ϕ = nothing
    plan_limber, T_k, plan_ℓ = _setup_limber_plan()

    w_δ = w_2_00_ϕTT()
    # w_minus2_00_ϕTT is always needed: γ×γ (L auto), κ×κ (C auto), γ×κ (L×C cross)
    w_μ_B = w_minus2_00_ϕTT()
    # w_0_00_ϕTT is always needed: δ×γ, δ×IA (G×L cross), δ×κ (G×C cross)
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
        plan_ϕ = prepare_chebyshev_plan(minimum(k_cheb), maximum(k_cheb), length(k_cheb) - 1; size_nd=(size(Blast.k_cheb, 1),), dim=1)
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
    
    Plans = FFTPlans(plan_ϕTT, plan_ϕT, plan_ϕ, plan_limber, T_k, plan_ℓ)
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
    cϕTT::cϕTT{Float64}
    cϕT::Union{cϕT{Float64}, Nothing}
    cϕ::Union{cϕ{Float64}, Nothing}
    ΔP_limber::Matrix{Float64}
    Pδ_limber::Matrix{Float64}
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
    As = get_As(cosmo)
    ns = get_ns(cosmo)
    return @. 9/25 * 2 * π^2 * As / (k^3) * (k/0.05)^(ns - 1.)
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
    transform_to_R_frame(matrix::AbstractArray{<:Any,2}, bg::Background)

Transforms a 2D array from standard coordinates `(z, k)` to a ratio-based 
coordinate system `(χ, Rχ, k)` using pre-built Background interpolators.
"""
function transform_to_R_frame(matrix::AbstractArray{<:Any,2}, bg::Background)
    new_χs = make_grid(Blast.χ, Blast.R)
    x = bg.z_of_χ.(new_χs)
    interp = akima_interpolation(matrix, Blast.z_lin, x)
    return reshape(interp, length(Blast.χ), length(Blast.R), size(interp, 2))
end

#TODO: this function will interface with mapse once ready.
"""
    prepare_pk_workspace(P::FFTPlans, 
                         pk::AbstractArray{<:Any, 2}, pk_limber_lin::AbstractArray{<:Any, 2}, 
                         pk_limber_nonlin::AbstractArray{<:Any, 2}, 
                         bg::Background)
"""
function prepare_pk_workspace(P::FFTPlans, pk::AbstractArray{<:Any, 2}, pk_limber_lin::AbstractArray{<:Any, 2}, pk_limber_nonlin::AbstractArray{<:Any, 2}, bg::Background)
    #Treating the non-limber power spectrum
    P_ϕ = get_PΦ(10 .^ Blast.k_cheb, bg.cosmo)
    transfer_func = get_Tm(pk, 10 .^ Blast.k_cheb, bg.cosmo)
    
    transfer_func_χR = transform_to_R_frame(transfer_func, bg)
    transfer_func_χ1 = transfer_func_χR[:,end,:]

    P_ϕTT = _p_phi_TT_tullio(P_ϕ, transfer_func_χ1, transfer_func_χR)
    P_ϕT  = _p_phi_T_tullio(P_ϕ, transfer_func_χR)

    c1 = build_coeff(cϕTT, P_ϕTT, P.plan_ϕTT)
    c2 = build_coeff(cϕT,  P_ϕT,  P.plan_ϕT)   
    c3 = build_coeff(cϕ,   P_ϕ,   P.plan_ϕ)    


    P_ϕ_limber = get_PΦ(10 .^ Blast.k_limber, bg.cosmo)'
    
    T_m_limber_lin = get_Tm(pk_limber_lin, 10 .^ Blast.k_limber, bg.cosmo)
    P_limber_lin = P_ϕ_limber .* T_m_limber_lin .^ 2.
    # permutedims materializes the transpose into a plain Matrix, avoiding
    # an Adjoint wrapper across the chebyshev_decomposition rrule boundary
    # (Mooncake's upstream rule in AbstractCosmologicalEmulators doesn't
    # support Adjoint tangents; plain Matrix works).
    c_lin = chebyshev_decomposition(P.plan_limber, permutedims(log10.(P_limber_lin)))

    T_m_limber_nonlin = get_Tm(pk_limber_nonlin, 10 .^ Blast.k_limber, bg.cosmo)
    P_limber_nonlin = P_ϕ_limber .* T_m_limber_nonlin .^ 2.
    c_nonlin = chebyshev_decomposition(P.plan_limber, permutedims(log10.(P_limber_nonlin)))

    T_z_limber = get_limber_coords_polynomials(P.plan_limber, bg.z)

    P_lin_grid = 10.0 .^ limber_eval(c_lin, T_z_limber, P.T_k_limber)
    P_nonlin_grid = 10.0 .^ limber_eval(c_nonlin, T_z_limber, P.T_k_limber)

    ΔP_limber = P_nonlin_grid .- P_lin_grid
    Pδ_limber = P_nonlin_grid

    return PowerSpectrum(c1, c2, c3, ΔP_limber, Pδ_limber)
end

# ----------------------------------------------------------------------------
# prepare_pk_workspace tensor products.  Both are pure outer-products (no
# contraction), so these are really broadcasts dressed as @tullio.  We extract
# them as named helpers and register rrules + Mooncake primitives so the AD
# contract is explicit and stable under future Tullio/Mooncake changes.
#
#   P_ϕTT[k,i,j] = P_ϕ[k] · T_χ1[i,k] · T_χR[i,j,k]
#   P_ϕT[k,i,j]  = P_ϕ[k]              · T_χR[i,j,k]
# ----------------------------------------------------------------------------
function _p_phi_TT_tullio(P_ϕ::AbstractVector, T_χ1::AbstractMatrix,
                          T_χR::AbstractArray{<:Any, 3})
    @tullio out[k, i, j] := P_ϕ[k] * T_χ1[i, k] * T_χR[i, j, k]
    return out
end

function _p_phi_T_tullio(P_ϕ::AbstractVector, T_χR::AbstractArray{<:Any, 3})
    @tullio out[k, i, j] := P_ϕ[k] * T_χR[i, j, k]
    return out
end
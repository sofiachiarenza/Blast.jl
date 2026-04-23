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
# Parametrized on each plan's concrete type so every field access infers to
# a concrete plan, not the abstract `ChebyshevPlan` UnionAll. Without this,
# `prepare_pk_workspace(::FFTPlans, …)` could not infer its return to
# `PowerSpectrum{T}` — only to the bare `PowerSpectrum` UnionAll — which
# poisoned `f_full`'s inference all the way to the top.
@kwdef mutable struct FFTPlans{PϕTT<:ChebyshevPlan,
                               PϕT<:Union{ChebyshevPlan, Nothing},
                               Pϕ<:Union{ChebyshevPlan, Nothing},
                               PL<:ChebyshevPlan,
                               Pℓ<:ChebyshevPlan}
    plan_ϕTT::PϕTT
    plan_ϕT::PϕT = nothing
    plan_ϕ::Pϕ = nothing
    plan_limber::PL
    T_k_limber::Array{Float64, 3}
    plan_ℓ::Pℓ
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

# Unequal-time (ϕTT) FFT plan over (k, χ, R). Identical for every SetUp.
function _setup_plan_ϕTT()
    return prepare_chebyshev_plan(minimum(k_cheb), maximum(k_cheb), length(k_cheb) - 1;
                                  size_nd=(size(Blast.k_cheb, 1), size(Blast.χ, 1), size(Blast.R, 1)),
                                  dim=1)
end

# Primordial (ϕ) FFT plan over k only. Needed when PNG is active.
function _setup_plan_ϕ()
    return prepare_chebyshev_plan(minimum(k_cheb), maximum(k_cheb), length(k_cheb) - 1;
                                  size_nd=(size(Blast.k_cheb, 1),),
                                  dim=1)
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
    plan_ϕTT = _setup_plan_ϕTT()
    plan_limber, T_k, plan_ℓ = _setup_limber_plan()

    # Each `have_*` folds to a compile-time constant once G's concrete type
    # params are known, so every ternary below is const-propagated away and
    # the returned `ProjectedMatterDensity{...}` has fully concrete type
    # parameters (never `Union{T, Nothing}` nor UnionAll).
    have_RSD = !isnothing(G.RSD)
    have_μ   = !isnothing(G.μ)
    have_PNG = !isnothing(G.PNG)

    plan_ϕT = have_PNG ? plan_ϕTT : nothing
    plan_ϕ  = have_PNG ? _setup_plan_ϕ() : nothing

    W = ProjectedMatterDensity(
        w_2_00_ϕTT      = w_2_00_ϕTT(),  # δ (always)
        w_minus2_00_ϕTT = have_μ ? w_minus2_00_ϕTT() : nothing,
        w_0_00_ϕTT      = have_μ ? w_0_00_ϕTT() : nothing,
        w_0_02_ϕTT      = (have_RSD && have_μ) ? w_0_02_ϕTT() : nothing,
        w_0_20_ϕTT      = (have_RSD && have_μ) ? w_0_20_ϕTT() : nothing,
        w_2_02_ϕTT      = have_RSD ? w_2_02_ϕTT() : nothing,
        w_2_20_ϕTT      = have_RSD ? w_2_20_ϕTT() : nothing,
        w_2_22_ϕTT      = have_RSD ? w_2_22_ϕTT() : nothing,
        w_2_00_ϕT       = have_PNG ? w_2_00_ϕT() : nothing,
        w_2_00_ϕT_R1    = have_PNG ? w_2_00_ϕT_R1() : nothing,
        w_0_00_ϕT       = (have_μ && have_PNG) ? w_0_00_ϕT() : nothing,
        w_0_00_ϕT_R1    = (have_μ && have_PNG) ? w_0_00_ϕT_R1() : nothing,
        w_2_02_ϕT       = (have_RSD && have_PNG) ? w_2_02_ϕT() : nothing,
        w_2_02_ϕT_R1    = (have_RSD && have_PNG) ? w_2_02_ϕT_R1() : nothing,
        w_2_20_ϕT       = (have_RSD && have_PNG) ? w_2_20_ϕT() : nothing,
        w_2_20_ϕT_R1    = (have_RSD && have_PNG) ? w_2_20_ϕT_R1() : nothing,
        w_2_00_ϕ        = have_PNG ? w_2_00_ϕ() : nothing,
    )

    Plans = FFTPlans(plan_ϕTT, plan_ϕT, plan_ϕ, plan_limber, T_k, plan_ℓ)
    return W, Plans
end

function SetUp(L::WeakLensing)
    plan_ϕTT = _setup_plan_ϕTT()
    plan_limber, T_k, plan_ℓ = _setup_limber_plan()

    W = ProjectedMatterDensity(
        w_minus2_00_ϕTT = w_minus2_00_ϕTT(),  # γ (always under L)
    )

    Plans = FFTPlans(plan_ϕTT, nothing, nothing, plan_limber, T_k, plan_ℓ)
    return W, Plans
end

function SetUp(G::GalaxyClustering, L::WeakLensing)
    plan_ϕTT = _setup_plan_ϕTT()
    plan_limber, T_k, plan_ℓ = _setup_limber_plan()

    have_RSD = !isnothing(G.RSD)
    have_μ   = !isnothing(G.μ)
    have_PNG = !isnothing(G.PNG)
    have_γ   = !isnothing(L.γ)  # always true — L.γ is required

    plan_ϕT = have_PNG ? plan_ϕTT : nothing
    plan_ϕ  = have_PNG ? _setup_plan_ϕ() : nothing

    W = ProjectedMatterDensity(
        w_2_00_ϕTT      = w_2_00_ϕTT(),            # δ
        w_minus2_00_ϕTT = w_minus2_00_ϕTT(),       # γ
        w_0_00_ϕTT      = w_0_00_ϕTT(),            # δ×γ cross
        w_0_02_ϕTT      = (have_RSD && (have_μ || have_γ)) ? w_0_02_ϕTT() : nothing,
        w_0_20_ϕTT      = (have_RSD && (have_μ || have_γ)) ? w_0_20_ϕTT() : nothing,
        w_2_02_ϕTT      = have_RSD ? w_2_02_ϕTT() : nothing,
        w_2_20_ϕTT      = have_RSD ? w_2_20_ϕTT() : nothing,
        w_2_22_ϕTT      = have_RSD ? w_2_22_ϕTT() : nothing,
        w_2_00_ϕT       = have_PNG ? w_2_00_ϕT() : nothing,
        w_2_00_ϕT_R1    = have_PNG ? w_2_00_ϕT_R1() : nothing,
        w_0_00_ϕT       = ((have_μ || have_γ) && have_PNG) ? w_0_00_ϕT() : nothing,
        w_0_00_ϕT_R1    = ((have_μ || have_γ) && have_PNG) ? w_0_00_ϕT_R1() : nothing,
        w_2_02_ϕT       = (have_RSD && have_PNG) ? w_2_02_ϕT() : nothing,
        w_2_02_ϕT_R1    = (have_RSD && have_PNG) ? w_2_02_ϕT_R1() : nothing,
        w_2_20_ϕT       = (have_RSD && have_PNG) ? w_2_20_ϕT() : nothing,
        w_2_20_ϕT_R1    = (have_RSD && have_PNG) ? w_2_20_ϕT_R1() : nothing,
        w_2_00_ϕ        = have_PNG ? w_2_00_ϕ() : nothing,
    )

    Plans = FFTPlans(plan_ϕTT, plan_ϕT, plan_ϕ, plan_limber, T_k, plan_ℓ)
    return W, Plans
end

function SetUp(L::WeakLensing, G::GalaxyClustering)
    return SetUp(G, L)
end

function SetUp(G::GalaxyClustering, C::CMB)
    plan_ϕTT = _setup_plan_ϕTT()
    plan_limber, T_k, plan_ℓ = _setup_limber_plan()

    have_RSD = !isnothing(G.RSD)
    have_μ   = !isnothing(G.μ)
    have_PNG = !isnothing(G.PNG)
    have_κ   = !isnothing(C.κ)  # always true — C.κ is required

    plan_ϕT = have_PNG ? plan_ϕTT : nothing
    plan_ϕ  = have_PNG ? _setup_plan_ϕ() : nothing

    W = ProjectedMatterDensity(
        w_2_00_ϕTT      = w_2_00_ϕTT(),     # δ
        w_minus2_00_ϕTT = have_μ ? w_minus2_00_ϕTT() : nothing,
        w_0_00_ϕTT      = w_0_00_ϕTT(),     # κ×δ cross (always under G×C)
        w_0_02_ϕTT      = (have_RSD && (have_μ || have_κ)) ? w_0_02_ϕTT() : nothing,
        w_0_20_ϕTT      = (have_RSD && (have_μ || have_κ)) ? w_0_20_ϕTT() : nothing,
        w_2_02_ϕTT      = have_RSD ? w_2_02_ϕTT() : nothing,
        w_2_20_ϕTT      = have_RSD ? w_2_20_ϕTT() : nothing,
        w_2_22_ϕTT      = have_RSD ? w_2_22_ϕTT() : nothing,
        w_2_00_ϕT       = have_PNG ? w_2_00_ϕT() : nothing,
        w_2_00_ϕT_R1    = have_PNG ? w_2_00_ϕT_R1() : nothing,
        w_0_00_ϕT       = ((have_μ || have_κ) && have_PNG) ? w_0_00_ϕT() : nothing,
        w_0_00_ϕT_R1    = ((have_μ || have_κ) && have_PNG) ? w_0_00_ϕT_R1() : nothing,
        w_2_02_ϕT       = (have_RSD && have_PNG) ? w_2_02_ϕT() : nothing,
        w_2_02_ϕT_R1    = (have_RSD && have_PNG) ? w_2_02_ϕT_R1() : nothing,
        w_2_20_ϕT       = (have_RSD && have_PNG) ? w_2_20_ϕT() : nothing,
        w_2_20_ϕT_R1    = (have_RSD && have_PNG) ? w_2_20_ϕT_R1() : nothing,
        w_2_00_ϕ        = have_PNG ? w_2_00_ϕ() : nothing,
    )

    Plans = FFTPlans(plan_ϕTT, plan_ϕT, plan_ϕ, plan_limber, T_k, plan_ℓ)
    return W, Plans
end

function SetUp(C::CMB, G::GalaxyClustering)
    return SetUp(G, C)
end

function SetUp(L::WeakLensing, C::CMB)
    plan_ϕTT = _setup_plan_ϕTT()
    plan_limber, T_k, plan_ℓ = _setup_limber_plan()

    W = ProjectedMatterDensity(
        w_minus2_00_ϕTT = w_minus2_00_ϕTT(),  # γ (always); κ×γ cross uses this
    )

    Plans = FFTPlans(plan_ϕTT, nothing, nothing, plan_limber, T_k, plan_ℓ)
    return W, Plans
end

function SetUp(C::CMB, L::WeakLensing)
    return SetUp(L, C)
end

function SetUp(G::GalaxyClustering, L::WeakLensing, C::CMB)
    plan_ϕTT = _setup_plan_ϕTT()
    plan_limber, T_k, plan_ℓ = _setup_limber_plan()

    have_RSD = !isnothing(G.RSD)
    have_μ   = !isnothing(G.μ)
    have_PNG = !isnothing(G.PNG)
    have_γ   = !isnothing(L.γ)  # always true — L.γ required
    have_κ   = !isnothing(C.κ)  # always true — C.κ required

    plan_ϕT = have_PNG ? plan_ϕTT : nothing
    plan_ϕ  = have_PNG ? _setup_plan_ϕ() : nothing

    W = ProjectedMatterDensity(
        w_2_00_ϕTT      = w_2_00_ϕTT(),             # δ
        w_minus2_00_ϕTT = w_minus2_00_ϕTT(),        # γ×γ, γ×κ
        w_0_00_ϕTT      = w_0_00_ϕTT(),             # δ×γ, δ×IA, δ×κ
        w_0_02_ϕTT      = (have_RSD && (have_μ || have_γ)) ? w_0_02_ϕTT() : nothing,
        w_0_20_ϕTT      = (have_RSD && (have_μ || have_γ)) ? w_0_20_ϕTT() : nothing,
        w_2_02_ϕTT      = have_RSD ? w_2_02_ϕTT() : nothing,
        w_2_20_ϕTT      = have_RSD ? w_2_20_ϕTT() : nothing,
        w_2_22_ϕTT      = have_RSD ? w_2_22_ϕTT() : nothing,
        w_2_00_ϕT       = have_PNG ? w_2_00_ϕT() : nothing,
        w_2_00_ϕT_R1    = have_PNG ? w_2_00_ϕT_R1() : nothing,
        w_0_00_ϕT       = ((have_μ || have_γ) && have_PNG) ? w_0_00_ϕT() : nothing,
        w_0_00_ϕT_R1    = ((have_μ || have_γ) && have_PNG) ? w_0_00_ϕT_R1() : nothing,
        w_2_02_ϕT       = (have_RSD && have_PNG) ? w_2_02_ϕT() : nothing,
        w_2_02_ϕT_R1    = (have_RSD && have_PNG) ? w_2_02_ϕT_R1() : nothing,
        w_2_20_ϕT       = (have_RSD && have_PNG) ? w_2_20_ϕT() : nothing,
        w_2_20_ϕT_R1    = (have_RSD && have_PNG) ? w_2_20_ϕT_R1() : nothing,
        w_2_00_ϕ        = have_PNG ? w_2_00_ϕ() : nothing,
    )

    Plans = FFTPlans(plan_ϕTT, plan_ϕT, plan_ϕ, plan_limber, T_k, plan_ℓ)
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
mutable struct PowerSpectrum{T<:Real}
    cϕTT::cϕTT{T}
    cϕT::Union{cϕT{T}, Nothing}
    cϕ::Union{cϕ{T}, Nothing}
    ΔP_limber::Matrix{T}
    Pδ_limber::Matrix{T}
end

# Outer constructor: infer T from cϕTT (which is always present), promote
# ΔP_limber/Pδ_limber to match. cϕT/cϕ may be Nothing — they're typed via
# the struct's Union{..., Nothing}.
function PowerSpectrum(cϕTT::cϕTT{T},
                       cϕT::Union{cϕT{T}, Nothing},
                       cϕ::Union{cϕ{T}, Nothing},
                       ΔP_limber::AbstractMatrix,
                       Pδ_limber::AbstractMatrix) where {T<:Real}
    PowerSpectrum{T}(cϕTT, cϕT, cϕ,
                     convert(Matrix{T}, ΔP_limber),
                     convert(Matrix{T}, Pδ_limber))
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

# Background-dispatched overloads: not a convenience, a type-stability path.
# The cosmo-only methods above return `::Any` because `w0waCDMCosmology`
# declares its parameter fields as the abstract `::Number` (upstream, in
# AbstractCosmologicalEmulators), so `get_As(cosmo)` / `get_ns(cosmo)` come
# back as `Any` and poison everything downstream (→ `prepare_pk_workspace`
# → `PowerSpectrum` UnionAll → `f_full` inferred as `Any`).
#
# Routing through `Background{T}` uses the typed `get_As(bg)` / `get_ns(bg)`
# from cosmo.jl, whose `convert(T, …)::T` pins the return at the concrete
# eltype `T`. Fixing it at the `w0waCDMCosmology` layer is an upstream
# change and is deliberately not attempted here.
function get_PΦ(k::AbstractArray{<:Any,1}, bg::Background{T}) where {T}
    As = get_As(bg)
    ns = get_ns(bg)
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

# Background-dispatched overload — see note on get_PΦ above.
function get_Tm(pk::AbstractArray{<:Any,2}, k::AbstractArray{<:Any, 1},
                bg::Background)
    prim_pk = get_PΦ(k, bg)
    return sqrt.(pk ./ reshape(prim_pk, 1, :))
end

"""
    transform_to_R_frame(matrix::AbstractArray{<:Any,2}, bg::Background)

Transforms a 2D array from standard coordinates `(z, k)` to a ratio-based 
coordinate system `(χ, Rχ, k)` using pre-built Background interpolators.
"""
function transform_to_R_frame(matrix::AbstractArray{<:Any,2}, bg::Background)
    new_χs = make_grid(Blast.χ, Blast.R)
    # Pass the whole χ vector to the z_of_χ akima closure in ONE call rather
    # than scalar-broadcasting (`bg.z_of_χ.(new_χs)`). The dotted form
    # rebuilt the akima coefficient arrays (b, c, d of length ~1000) on every
    # single one of the 6144 scalar queries, costing ~100 ms and ~1 GiB per
    # Background evaluation. The vector form shares a single setup across
    # queries via `_akima_eval(..., tq::AbstractArray)` — ~1 ms, ~200 KiB.
    x = bg.z_of_χ(new_χs)
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
    P_ϕ = get_PΦ(10 .^ Blast.k_cheb, bg)
    transfer_func = get_Tm(pk, 10 .^ Blast.k_cheb, bg)
    
    transfer_func_χR = transform_to_R_frame(transfer_func, bg)
    transfer_func_χ1 = transfer_func_χR[:,end,:]

    P_ϕTT = _p_phi_TT_tullio(P_ϕ, transfer_func_χ1, transfer_func_χR)
    P_ϕT  = _p_phi_T_tullio(P_ϕ, transfer_func_χR)

    c1 = build_coeff(cϕTT, P_ϕTT, P.plan_ϕTT)
    c2 = build_coeff(cϕT,  P_ϕT,  P.plan_ϕT)   
    c3 = build_coeff(cϕ,   P_ϕ,   P.plan_ϕ)    


    P_ϕ_limber = get_PΦ(10 .^ Blast.k_limber, bg)'
    
    T_m_limber_lin = get_Tm(pk_limber_lin, 10 .^ Blast.k_limber, bg)
    P_limber_lin = P_ϕ_limber .* T_m_limber_lin .^ 2.
    # permutedims materializes the transpose into a plain Matrix, avoiding
    # an Adjoint wrapper across the chebyshev_decomposition rrule boundary
    # (Mooncake's upstream rule in AbstractCosmologicalEmulators doesn't
    # support Adjoint tangents; plain Matrix works).
    c_lin = chebyshev_decomposition(P.plan_limber, permutedims(log10.(P_limber_lin)))

    T_m_limber_nonlin = get_Tm(pk_limber_nonlin, 10 .^ Blast.k_limber, bg)
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
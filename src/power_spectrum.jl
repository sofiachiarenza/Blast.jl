@kwdef mutable struct PowerSpectrum 
    PΦ::AbstractArray{<:Any, 1} = zeros(1)
    Tm::AbstractArray{<:Any, 2} = zeros(1,1)
    unequal_time::AbstractArray{<:Any, 3} = zeros(1,1,1)
    cϕTT::cϕTT = cϕTT()
    cϕT::Union{cϕT, NullCoeff} = NullCoeff()
    cϕ::Union{cϕ, NullCoeff} = NullCoeff()
    ΔP_limber::AbstractArray{<:Any, 2} = zeros(1,1)
    Pδ_limber::AbstractArray{<:Any, 2} = zeros(1,1)
end

function get_PΦ(k::AbstractArray{<:Any,1}, cosmo::AbstractCosmology)
    return @. 9/25 * 2 * π^2 * cosmo.As / (k^3) * (k/0.05)^(cosmo.ns - 1.)
end

function get_Tm(pk::AbstractArray{<:Any,2}, cosmo::AbstractCosmology)
    prim_pk = get_PΦ( 10 .^ Blast.k_cheb , cosmo)
    return sqrt.(pk ./ reshape(prim_pk, 1, :))
end

function to_χR_frame(matrix::AbstractArray{<:Any,2}, plan::FFTW.r2rFFTWPlan, bg::BackgroundQuantities, grid::AbstractCosmologicalGrid)
    coefs = fast_chebcoefs(matrix, plan)
    new_χs = make_grid(bg, R)
    x = resample_redshifts(bg, grid, new_χs)
    chebyshevs = chebyshev_polynomials(x, z_cheb)
    pk_interp = coefs' * chebyshevs  
    return reshape( pk_interp, size(k_cheb,1),  size(bg.χz_array, 1), size(R,1) )
end

function PowerSpectrumSetUp(pk::AbstractArray{<:Any, 2}, pk_limber_lin::AbstractArray{<:Any, 2}, pk_limber_nonlin::AbstractArray{<:Any, 2}, f_NL::Number, cosmo::AbstractCosmology, bg::BackgroundQuantities, grid::AbstractCosmologicalGrid)
    primordial_pk = get_PΦ(10 .^ Blast.k_cheb, cosmo)
    transfer_func = get_Tm(pk, cosmo)
    plan_z = plan_fft(log10.(transfer_func), 1)
    transfer_func_χR = 10 .^ to_χR_frame(log10.(transfer_func), plan_z, bg, grid)

    transfer_func_χ1 = transfer_func_χR[:,:,end]
    @tullio P_ϕTT[k, i, j] := primordial_pk[k] * transfer_func_χ1[k,i] * transfer_func_χR[k, i, j] 
    @tullio P_ϕT[k, i, j] := primordial_pk[k]* transfer_func_χR[k, i, j]

    plan = Blast.plan_fft(P_ϕTT,1)

    c1 = cϕTT()
    c = Blast.fast_chebcoefs(P_ϕTT, plan)
    c1.coefs = permutedims(c, (2,3,1))

    if f_NL == 0
        c2 = NullCoeff()
        c3 = NullCoeff()
    else
        c2 = cϕT()
        c = Blast.fast_chebcoefs(P_ϕT, plan)
        c2.coefs = permutedims(c, (2,3,1))

        c3 = cϕ()
        plan_primordial = Blast.plan_fft(primordial_pk,1) 
        c3.coefs = Blast.fast_chebcoefs(primordial_pk, plan_primordial)
    end

    lb, ub = [minimum(Blast.z_cheb),minimum(Blast.k_limber)], [maximum(Blast.z_cheb), maximum(Blast.k_limber)] # lower and upper bounds of the domain, respectively
    limber_pk_linear = chebinterp(log10.(pk_limber_lin), lb, ub)
    limber_pk_nonlinear = chebinterp(log10.(pk_limber_nonlin), lb, ub)

    z_of_χ = AkimaInterpolation(grid.z_range, bg.χz_array, extrapolation=ExtrapolationType.Extension)

    ℓ_grid = reshape(Blast.full_ℓ_range .+ 0.5, :, 1)   
    χ_grid = reshape(Blast.χ, 1, :)          
    k = ℓ_grid ./ χ_grid               
    ΔP_limber = [10. ^ limber_pk_nonlinear(SVector(z_of_χ.(χ)[j], log10(k[i, j]))) .- 10. ^ limber_pk_linear(SVector(z_of_χ.(χ)[j], log10(k[i, j]))) for i in 1:size(ℓ_grid, 1), j in 1:size(Blast.χ, 1)]
    Pδ_limber = [10. ^ limber_pk_nonlinear(SVector(z_of_χ.(χ)[j], log10(k[i, j]))) for i in 1:size(ℓ_grid, 1), j in 1:size(Blast.χ, 1)]

    PK = PowerSpectrum(primordial_pk, transfer_func, P_ϕTT, c1, c2, c3,  ΔP_limber, Pδ_limber)
end
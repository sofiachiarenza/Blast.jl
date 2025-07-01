function Cℓ_limber(pk, ℓ, bg::BackgroundQuantities, ProbeA::Blast.GalaxyKernel, ProbeB::Blast.GalaxyKernel )
    χ = bg.χz_array
    n = length(χ)
    
    F = 1

    Δχ = ((χ[n]-χ[1])/(n-1))
    pesi = Blast.simpson_weight_array(n)

    pk_over_chi = pk ./ (χ .^ 2)

    KA = ProbeA.Kernel
    KB = ProbeB.Kernel

    @tullio Cℓ[i,j] := Δχ*pk_over_chi[m]*KA[i,m]*KB[j,m]*pesi[m]*F
    return Cℓ
end

function Cℓ_limber(pk, ℓ, bg::BackgroundQuantities, ProbeA::Blast.GalaxyKernel, ProbeB::Blast.ShearKernel )
    χ = bg.χz_array
    n = length(χ)
    
    F = @. sqrt(Blast.factorial_frac(ℓ)) / (ℓ+0.5)^2

    Δχ = ((χ[n]-χ[1])/(n-1))
    pesi = Blast.simpson_weight_array(n)

    pk_over_chi = pk ./ (χ .^ 2)

    KA = ProbeA.Kernel
    KB = ProbeB.Kernel

    @tullio Cℓ[i,j] := Δχ*pk_over_chi[m]*KA[i,m]*KB[j,m]*pesi[m]*F
    return Cℓ
end

function Cℓ_limber(pk, ℓ, bg::BackgroundQuantities, ProbeA::Blast.ShearKernel, ProbeB::Blast.ShearKernel )
    χ = bg.χz_array
    n = length(χ)
    
    F = @. Blast.factorial_frac(ℓ) / (ℓ+0.5)^4

    Δχ = ((χ[n]-χ[1])/(n-1))
    pesi = Blast.simpson_weight_array(n)

    pk_over_chi = pk ./ (χ .^ 2)

    KA = ProbeA.Kernel
    KB = ProbeB.Kernel

    @tullio Cℓ[i,j] := Δχ*pk_over_chi[m]*KA[i,m]*KB[j,m]*pesi[m]*F
    return Cℓ
end

function Cℓ_limber(pk, ℓ, bg::BackgroundQuantities, ProbeA::Blast.CMBLensingKernel, ProbeB::Blast.ShearKernel )
    χ = bg.χz_array
    n = length(χ)
    
    F = @. sqrt(Blast.factorial_frac(ℓ)) * ℓ * (ℓ+1) / (ℓ+0.5)^4

    Δχ = ((χ[n]-χ[1])/(n-1))
    pesi = Blast.simpson_weight_array(n)

    pk_over_chi = pk ./ (χ .^ 2)

    KA = reshape(ProbeA.Kernel, 1, n)
    KB = ProbeB.Kernel

    @tullio Cℓ[i,j] := Δχ*pk_over_chi[m]*KA[i,m]*KB[j,m]*pesi[m]*F
    return Cℓ
end

function Cℓ_limber(pk, ℓ, bg::BackgroundQuantities, ProbeA::Blast.CMBLensingKernel, ProbeB::Blast.GalaxyKernel )
    χ = bg.χz_array
    n = length(χ)
    bins = size(ProbeB.Kernel)[1]
    F = @. ℓ * (ℓ + 1) / (ℓ+0.5)^2

    Δχ = ((χ[n]-χ[1])/(n-1))
    pesi = Blast.simpson_weight_array(n)

    pk_over_chi = pk ./ (χ .^ 2)

    KA = reshape(ProbeA.Kernel, 1, length(χ))
    KB = ProbeB.Kernel

    @tullio Cℓ[i,j] := Δχ*pk_over_chi[m]*KA[i,m]*KB[j,m]*pesi[m]*F
    return Cℓ
end

function Cℓ_limber(pk, ℓ, bg::BackgroundQuantities, ProbeA::Blast.CMBLensingKernel, ProbeB::Blast.RSDKernel )
    χ = bg.χz_array
    n = length(χ)
    bins = size(ProbeB.Kernel)[1]
    F = @. ℓ * (ℓ + 1) / (ℓ+0.5)^2

    Δχ = ((χ[n]-χ[1])/(n-1))
    pesi = Blast.simpson_weight_array(n)

    pk_over_chi = pk ./ (χ .^ 2)

    KA = reshape(ProbeA.Kernel, 1, length(χ))
    KB = reshape(limber_rsd_kernel(ℓ, bg, ProbeB), bins, n)

    @tullio Cℓ[i,j] := Δχ*pk_over_chi[m]*KA[i,m]*KB[j,m]*pesi[m]*F
    return Cℓ
end

function Cℓ_limber(pk, ℓ, bg::BackgroundQuantities, ProbeA::Blast.CMBLensingKernel, ProbeB::Blast.LensMagKernel )
    χ = bg.χz_array
    n = length(χ)
    F = @. (ℓ * (ℓ + 1))^2 / (ℓ+0.5)^4

    Δχ = ((χ[n]-χ[1])/(n-1))
    pesi = Blast.simpson_weight_array(n)

    pk_over_chi = pk ./ (χ .^ 2)

    KA = reshape(ProbeA.Kernel, 1, length(χ))
    KB = ProbeB.Kernel

    @tullio Cℓ[i,j] := Δχ*pk_over_chi[m]*KA[i,m]*KB[j,m]*pesi[m]*F
    return Cℓ
end

function Cℓ_limber(pk, ℓ, bg::BackgroundQuantities, ProbeA::Blast.GalaxyKernel, ProbeB::Blast.RSDKernel )
    χ = bg.χz_array
    n = length(χ)
    F = 1

    Δχ = ((χ[n]-χ[1])/(n-1))
    pesi = Blast.simpson_weight_array(n)

    pk_over_chi = pk ./ (χ .^ 2)

    KA = ProbeA.Kernel
    KB = reshape(limber_rsd_kernel(ℓ, bg, ProbeB), size(ProbeB.Kernel)[1], n)

    @tullio Cℓ[i,j] := Δχ*pk_over_chi[m]*KA[i,m]*KB[j,m]*pesi[m]*F
    return Cℓ
end

function Cℓ_limber(pk, ℓ, bg::BackgroundQuantities, ProbeA::Blast.RSDKernel, ProbeB::Blast.RSDKernel )
    χ = bg.χz_array
    n = length(χ)
    bins = size(ProbeA.Kernel)[1]
    F = 1

    Δχ = ((χ[n]-χ[1])/(n-1))
    pesi = Blast.simpson_weight_array(n)

    pk_over_chi = pk ./ (χ .^ 2)

    KA = reshape(limber_rsd_kernel(ℓ, bg, ProbeA), bins, n)
    KB = reshape(limber_rsd_kernel(ℓ, bg, ProbeB), bins, n)

    @tullio Cℓ[i,j] := Δχ*pk_over_chi[m]*KA[i,m]*KB[j,m]*pesi[m]*F
    return Cℓ
end

function Cℓ_limber(pk, ℓ, bg::BackgroundQuantities, ProbeA::Blast.GalaxyKernel, ProbeB::Blast.LensMagKernel )
    χ = bg.χz_array
    n = length(χ)
    
    F = @. ℓ * (ℓ + 1) / (ℓ+0.5)^2

    Δχ = ((χ[n]-χ[1])/(n-1))
    pesi = Blast.simpson_weight_array(n)

    pk_over_chi = pk ./ (χ .^ 2)

    KA = ProbeA.Kernel
    KB = ProbeB.Kernel

    @tullio Cℓ[i,j] := Δχ*pk_over_chi[m]*KA[i,m]*KB[j,m]*pesi[m]*F
    return Cℓ
end

function Cℓ_limber(pk, ℓ, bg::BackgroundQuantities, ProbeA::Blast.LensMagKernel, ProbeB::Blast.LensMagKernel )
    χ = bg.χz_array
    n = length(χ)
    
    F = @. (ℓ * (ℓ + 1))^2 / (ℓ+0.5)^4

    Δχ = ((χ[n]-χ[1])/(n-1))
    pesi = Blast.simpson_weight_array(n)

    pk_over_chi = pk ./ (χ .^ 2)

    KA = ProbeA.Kernel
    KB = ProbeB.Kernel

    @tullio Cℓ[i,j] := Δχ*pk_over_chi[m]*KA[i,m]*KB[j,m]*pesi[m]*F
    return Cℓ
end

function Cℓ_limber(pk, ℓ, bg::BackgroundQuantities, ProbeA::Blast.RSDKernel, ProbeB::Blast.LensMagKernel )
    χ = bg.χz_array
    n = length(χ)
    bins = size(ProbeA.Kernel)[1]
    F = @. ℓ * (ℓ + 1) / (ℓ+0.5)^2

    Δχ = ((χ[n]-χ[1])/(n-1))
    pesi = Blast.simpson_weight_array(n)

    pk_over_chi = pk ./ (χ .^ 2)

    KA = reshape(limber_rsd_kernel(ℓ, bg, ProbeA), bins, n)
    KB = ProbeB.Kernel

    @tullio Cℓ[i,j] := Δχ*pk_over_chi[m]*KA[i,m]*KB[j,m]*pesi[m]*F
    return Cℓ
end;

function Cℓ_limber(pk, ℓ, bg::BackgroundQuantities, ProbeA::Blast.ShearKernel, ProbeB::Blast.LensMagKernel )
    χ = bg.χz_array
    n = length(χ)
    bins = size(ProbeA.Kernel)[1]
    F = @. sqrt(factorial_frac(ℓ)) / (ℓ+0.5)^4 * ℓ * (ℓ + 1) 

    Δχ = ((χ[n]-χ[1])/(n-1))
    pesi = Blast.simpson_weight_array(n)

    pk_over_chi = pk ./ (χ .^ 2)

    KA = ProbeA.Kernel
    KB = ProbeB.Kernel

    @tullio Cℓ[i,j] := Δχ*pk_over_chi[m]*KA[i,m]*KB[j,m]*pesi[m]*F
    return Cℓ
end

function Cℓ_limber(pk, ℓ, bg::BackgroundQuantities, ProbeA::Blast.ShearKernel, ProbeB::Blast.RSDKernel )
    χ = bg.χz_array
    n = length(χ)
    F = @. sqrt(factorial_frac(ℓ)) / (ℓ+0.5)^2

    Δχ = ((χ[n]-χ[1])/(n-1))
    pesi = Blast.simpson_weight_array(n)

    pk_over_chi = pk ./ (χ .^ 2)

    KA = ProbeA.Kernel
    KB = reshape(limber_rsd_kernel(ℓ, bg, ProbeB), size(ProbeB.Kernel)[1], n)

    @tullio Cℓ[i,j] := Δχ*pk_over_chi[m]*KA[i,m]*KB[j,m]*pesi[m]*F
    return Cℓ
end


#Limber Inclusion in the new code. All the above will be deprecated. 
function get_limber_kernel(Component::AbstractComponents)
    total_prefactor = Component.ell_prefactor .* Component.limber_factor
    kernel = Component.Kernel'
    kernel = reshape(kernel, 1, size(kernel, 1), size(kernel, 2))  # (1, 200, nbins)
    prefactor = reshape(total_prefactor, :, 1, 1)  # (101, 1, 1)
    return prefactor .* kernel  # Result: (101, 200, nbins)
end

function get_limber_kernel(Component::Nothing)
    return 0.
end

function get_limber_kernel(G::GalaxyClustering)
    #TODO: in the correction i'm currently excluding RSDs. The limber implementation still doesn't work, but also the contribution at the scales of interest is null basically.
    return get_limber_kernel(G.δ) .+ get_limber_kernel(G.μ) # .+ get_limber_kernel(G.PNG) .+ get_limber_kernel(G.RSD)
end

function get_limber_kernel(L::WeakLensing)
    return get_limber_kernel(L.γ) .+ get_limber_kernel(L.IA)
end

function get_limber_kernel(C::CMB)
    return get_limber_kernel(C.κ) .+ get_limber_kernel(C.ISW)
end

function get_limber_correction(Probe::Union{GalaxyClustering, WeakLensing, CMB}, pk::PowerSpectrum)
    ΔP_over_χ2 = pk.ΔP_limber[1:size(Blast.ℓ_nonlimber, 1), :] ./ reshape(Blast.χ, 1, :) .^ 2

    n = size(Blast.χ, 1)
    Δχ = ((χ[n]-χ[1])/(n-1))
    weights = Blast.simpson_weight_array(n)

    K = get_limber_kernel(Probe)[1:size(Blast.ℓ_nonlimber, 1), :, :]

    @tullio Cℓ[l,i,j] := ΔP_over_χ2[l,m]*K[l,m,i]*K[l,m,j]*weights[m]*Δχ
end

function get_limber_correction(ProbeA::Union{GalaxyClustering, WeakLensing, CMB}, ProbeB::Union{GalaxyClustering, WeakLensing, CMB}, pk::PowerSpectrum)
    ΔP_over_χ2 = pk.ΔP_limber[1:size(Blast.ℓ_nonlimber, 1), :] ./ reshape(Blast.χ, 1, :) .^ 2

    n = size(Blast.χ, 1)
    Δχ = ((χ[n]-χ[1])/(n-1))
    weights = Blast.simpson_weight_array(n)

    KA = get_limber_kernel(ProbeA)[1:size(Blast.ℓ_nonlimber, 1), :, :]
    KB = get_limber_kernel(ProbeB)[1:size(Blast.ℓ_nonlimber, 1), :, :]

    @tullio Cℓ[l,i,j] := ΔP_over_χ2[l,m]*KA[l,m,i]*KB[l,m,j]*weights[m]*Δχ
end

#this is to compute the Cℓ's for ℓ>215
function get_limber_Cℓ(Probe::Union{GalaxyClustering, WeakLensing, CMB}, pk::PowerSpectrum)
    Pδ_over_χ2 = pk.Pδ_limber[size(Blast.ℓ_nonlimber, 1)+1:end, :] ./ reshape(Blast.χ, 1, :) .^ 2

    n = size(Blast.χ, 1)
    Δχ = ((χ[n]-χ[1])/(n-1))
    weights = Blast.simpson_weight_array(n)

    K = get_limber_kernel(Probe)[size(Blast.ℓ_nonlimber, 1)+1:end, :, :]

    @tullio Cℓ[l,i,j] := Pδ_over_χ2[l,m]*K[l,m,i]*K[l,m,j]*weights[m]*Δχ
end

function get_limber_Cℓ(ProbeA::Union{GalaxyClustering, WeakLensing, CMB}, ProbeB::Union{GalaxyClustering, WeakLensing, CMB}, pk::PowerSpectrum)
    Pδ_over_χ2 = pk.Pδ_limber[size(Blast.ℓ_nonlimber, 1)+1:end, :] ./ reshape(Blast.χ, 1, :) .^ 2

    n = size(Blast.χ, 1)
    Δχ = ((χ[n]-χ[1])/(n-1))
    weights = Blast.simpson_weight_array(n)

    KA = get_limber_kernel(ProbeA)[size(Blast.ℓ_nonlimber, 1)+1:end, :, :]
    KB = get_limber_kernel(ProbeB)[size(Blast.ℓ_nonlimber, 1)+1:end, :, :]

    @tullio Cℓ[l,i,j] := Pδ_over_χ2[l,m]*KA[l,m,i]*KB[l,m,j]*weights[m]*Δχ
end



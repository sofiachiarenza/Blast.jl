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

function get_limber_kernel(G::GalaxyClustering)
    #TODO: in the correction i'm currently excluding RSDs. The limber implementation still doesn't work, but also the contribution at the scales of interest is null basically.
    return G.δ.Kernel  .* reshape(G.δ.ell_prefactor, :,1,1) .* reshape(G.δ.limber_factor, :,1,1)  .+ G.μ.Kernel  .* reshape(G.μ.ell_prefactor, :,1,1) .* reshape(G.μ.limber_factor, :,1,1) #.+ G.PNG.Kernel  .* reshape(G.PNG.ell_prefactor, :,1,1) .* reshape(G.PNG.limber_factor, :,1,1)
end

function get_limber_kernel(L::WeakLensing)
    return L.γ.Kernel .* reshape(L.γ.ell_prefactor, :, 1, 1) .* reshape(L.γ.limber_factor, :, 1, 1) #.+ L.IA.Kernel .* reshape(L.IA.ell_prefactor, :, 1, 1) .* reshape(L.IA.limber_factor, :, 1, 1)
end

function get_limber_kernel(C::CMB)
    return C.κ.Kernel .* reshape(C.κ.ell_prefactor, :, 1, 1) .* reshape(C.κ.limber_factor, :, 1, 1) .+ C.ISW.Kernel .* reshape(C.ISW.ell_prefactor, :, 1, 1) .* reshape(C.ISW.limber_factor, :, 1, 1)
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



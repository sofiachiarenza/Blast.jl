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


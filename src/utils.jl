function load_Ts(folder, nχ, nR, nk)
    ell_vector = Blast.ℓ
    full_T = zeros(length(ell_vector), nχ, nR, nk)
    for i in 1:length(ell_vector)
        l_string = string(round(ell_vector[i]; digits=1))
        filename = folder * "/T_tilde_l_$l_string.npy"
        if isfile(filename)
            full_T[i,:,:,:] = npzread(filename)
        else
            println("Missing file!")
        end
    end
    return full_T
end

function stoopid_2D_interpolator(x::AbstractArray{T,1}, y::AbstractArray{T,1}, f::AbstractArray{T,2}, logx::Bool, logy::Bool, logf::Bool ) where T
    if logx
        x_resc = LinRange(log10(first(x)),log10(last(x)), length(x))
    else
        x_resc = LinRange(first(x),last(x), length(x))
    end 

    if logy
        y_resc = LinRange(log10(first(y)),log10(last(y)), length(y))
    else
        y_resc = LinRange(first(y),last(y), length(y))
    end

    if logf 
        Interpolator = Interpolations.interpolate(log10.(f),BSpline(Cubic(Line(OnGrid()))))
        Interpolator = scale(Interpolator, (x_resc, y_resc))
        Interpolator = Interpolations.extrapolate(Interpolator, Line());
    else
        Interpolator = Interpolations.interpolate(f,BSpline(Cubic(Line(OnGrid()))))
        Interpolator = scale(Interpolator, (x_resc, y_resc))
        Interpolator = Interpolations.extrapolate(Interpolator, Line());
    end

    return Interpolator
end

function stoopid_unequal_time_pk(interpolator, k::AbstractArray{T,1}, z1::Number, z2::Number) where T
    #TODO: this only works if the interpolator is created with false true true
    return @. sqrt(10^interpolator(z1,log10(k)) * 10^interpolator(z2,log10(k)))
end

function limber_rsd_kernel(ℓ::Number, bg::BackgroundQuantities, RSDK::Blast.RSDKernel)
    χ = bg.χz_array
    nbins = size(RSDK.Kernel)[1]
    rds_kernels = zeros( nbins, length(χ) )

    for b in 1:nbins
        kernel_interp = AkimaInterpolation(RSDK.Kernel[b, :], χ, extrapolation=ExtrapolationType.Extension)
        piece1 = @. (2*ℓ^2 + 2*ℓ - 1) / ((2*ℓ - 1)*(2*ℓ + 3)) * RSDK.Kernel[b, :]
        piece2 = @. (ℓ - 1)*ℓ / ((2*ℓ - 1) * sqrt(2*ℓ - 3)*(2*ℓ + 1)) * kernel_interp.((2*ℓ - 3)/(2*ℓ + 1) * χ)
        piece3 = @. (ℓ + 1)*(ℓ + 2) / ((2*ℓ + 3) * sqrt((2*ℓ + 1)*(2*ℓ + 5))) * kernel_interp.((2*ℓ + 5)/(2*ℓ + 1) * χ)
        rds_kernels[b, :] .= piece1 .- piece2 .- piece3
    end

    return rds_kernels
end
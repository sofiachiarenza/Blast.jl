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
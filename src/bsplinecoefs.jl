"""
    build_coeff(C, values, plan::CubicBSplinePlan)

Fit source values with ACE's cubic B-spline collocation solve and wrap the
coefficients in the existing source-family container `C`.

The first axis is the spectral sampling axis. Remaining axes are independent
source series, e.g. `(k, χ, R)` for `cϕTT` and `cϕT`. Collocation values are
not B-spline coefficients: a solve is required even when the sample and
evaluation sites coincide. The result keeps the input shape and scalar type
promotion supplied by ACE, including ForwardDiff dual numbers.

This method does not normalize sources or select tensor banks. A ratio backend
must supply consistently normalized values and matching source-weighted banks.
"""
function build_coeff(::Type{C}, values::AbstractArray,
                     plan::AbstractCosmologicalEmulators.CubicBSplinePlan) where {C<:AbstractCoeff}
    coefficients = if values isa AbstractVector
        AbstractCosmologicalEmulators.bspline_coefficients(plan, values)
    else
        series = reshape(values, size(values, 1), :)
        fitted = AbstractCosmologicalEmulators.bspline_coefficients(plan, series)
        reshape(fitted, size(values))
    end
    return make_coeff(C, coefficients)
end

function _geometry_signature(ell,chi,r)
    io=IOBuffer();write(io,"blast-fixed-axes-v1")
    for a in (ell,chi,r)
        write(io,htol(UInt64(length(a))))
        write(io,reinterpret(UInt8,htol.(reinterpret(UInt64,collect(Float64,a)))))
    end
    w=ltoh.(reinterpret(UInt64,SHA.sha256(take!(io))))
    return (w[1],w[2],w[3],w[4])
end
const PRODUCTION_GEOMETRY_SIGNATURE=_geometry_signature(ℓ_nonlimber,χ,R)

"""
    RadialQuadraturePlan(chi, weights, suffix; R=Blast.R, R_weights=nothing, ell=Blast.ℓ_nonlimber)

Fixed physical-distance quadrature data. `weights` and `suffix` include the
physical integration measure in Mpc. `suffix[i,j]` multiplies source value j
in the integral from `chi[i]` to the right endpoint; partial-panel rows need
not be triangular. By default R uses the unchanged production rule. Custom
R nodes require explicit weights and an exact final R=1 source evaluation;
that endpoint may have zero integration weight. Weights are never normalized.

Construction copies inputs and builds geometry caches once, outside AD. Plan
arrays are read-only by contract and may be shared between inference chains.
This constructor neither loads nor resamples a tensor bank. `ell` declares the
bank multipoles for geometry compatibility; it does not change radial weights
or select a final multipole interpolator.
"""
struct RadialQuadraturePlan
    chi::Vector{Float64}
    weights::Vector{Float64}
    suffix::Matrix{Float64}
    R::Vector{Float64}
    R_weights::Vector{Float64}
    chiR_flat::Vector{Float64}
    inv_chiR2_row::Matrix{Float64}
    inv_chi2_row::Matrix{Float64}
    signature::NTuple{4,UInt64}
    legacy::Bool
    geometry::NTuple{4,UInt64}
end

function RadialQuadraturePlan(chi::AbstractVector, weights::AbstractVector,
                              suffix::AbstractMatrix; legacy=false,
                              R::AbstractVector=Blast.R, R_weights=nothing,
                              ell::AbstractVector=Blast.ℓ_nonlimber)
    x, w, s = collect(Float64,chi), collect(Float64,weights), Matrix{Float64}(suffix)
    ls=collect(Float64,ell)
    !isempty(ls) && all(isfinite,ls) && all(>(0),ls) ||
        throw(ArgumentError("bank multipoles must be finite and positive"))
    n=length(x)
    n>=2 && all(isfinite,x) && first(x)>0 && all(>(0),diff(x)) ||
        throw(ArgumentError("chi nodes must be finite, positive and strictly increasing"))
    length(w)==n && size(s)==(n,n) || throw(DimensionMismatch("radial weights/suffix must match chi"))
    all(isfinite,w) && all(isfinite,s) || throw(ArgumentError("radial weights must be finite"))
    if legacy
        x==collect(LinRange(first(x),last(x),n)) || throw(ArgumentError("legacy quadrature needs uniform chi"))
        step=(last(x)-first(x))/(n-1)
        w==step.*simpson_weights_array(n) && s==step.*simpson_weights_matrix(n) ||
            throw(ArgumentError("legacy weights must equal the existing Simpson rule"))
    end
    r=collect(Float64,R)
    length(r)>=2 && all(isfinite,r) && first(r)>0 && last(r)==1 && all(>(0),diff(r)) ||
        throw(ArgumentError("R nodes must be finite, positive, strictly increasing and end at exactly 1"))
    wr=if isnothing(R_weights)
        r==Blast.R || throw(ArgumentError("custom R nodes require explicit R_weights"))
        get_clencurt_weights_R_integration(2length(r)+1)
    else
        collect(Float64,R_weights)
    end
    length(wr)==length(r) || throw(DimensionMismatch("R weights must match R nodes"))
    all(isfinite,wr) && all(>=(0),wr) && any(>(0),wr) ||
        throw(ArgumentError("R weights must be finite, nonnegative and not all zero"))
    xr=vec(x*transpose(r))
    io=IOBuffer();write(io,"blast-radial-plan-v1")
    for a in (x,w,s,r,wr)
        write(io,htol(UInt64(length(a))))
        write(io,reinterpret(UInt8,htol.(reinterpret(UInt64,vec(a)))))
    end
    words=ltoh.(reinterpret(UInt64,SHA.sha256(take!(io))))
    return RadialQuadraturePlan(x,w,s,r,wr,xr,reshape(inv.(xr.^2),1,:),
        reshape(inv.(x.^2),1,:),(words[1],words[2],words[3],words[4]),legacy,
        _geometry_signature(ls,x,r))
end

const LEGACY_RADIAL_PLAN = let
    step=(last(χ)-first(χ))/(length(χ)-1)
    RadialQuadraturePlan(χ,step.*simpson_weights_array(length(χ)),
        step.*simpson_weights_matrix(length(χ));legacy=true)
end

_default_radial_plan(chi) = chi==Blast.χ ? LEGACY_RADIAL_PLAN : nothing

function _background_geometry(chi_grid, plan::Union{Nothing,RadialQuadraturePlan})
    grid=isnothing(chi_grid) ? (isnothing(plan) ? Blast.χ : plan.chi) : chi_grid
    if isnothing(plan)
        plan=_default_radial_plan(grid)
    else
        grid==plan.chi || throw(ArgumentError("Background chi_grid differs from its radial plan"))
    end
    return grid,plan
end

function _kernel_quadrature(chi, plan::Nothing)
    n=length(chi)
    chi==collect(LinRange(first(chi),last(chi),n)) || throw(ArgumentError(
        "nonuniform Background kernels require an explicit RadialQuadraturePlan"))
    return simpson_weights_matrix(n),(last(chi)-first(chi))/(n-1)
end
function _kernel_quadrature(chi, plan::RadialQuadraturePlan)
    if plan.legacy
        return simpson_weights_matrix(length(chi)),(last(chi)-first(chi))/(length(chi)-1)
    end
    return plan.suffix,1.0
end

function _require_radial_background(bg,plan::RadialQuadraturePlan)
    length(bg.χ)==length(plan.chi) || throw(DimensionMismatch("Background and radial plan chi lengths differ"))
    bg.χ==plan.chi || throw(ArgumentError("Background chi coordinates differ from the radial plan"))
    !isnothing(bg.radial_plan) && bg.radial_plan.signature==plan.signature ||
        throw(ArgumentError("Background and spectrum use different radial weights/suffix policies"))
    return nothing
end

function _radial_outer(plan::RadialQuadraturePlan)
    if plan.legacy
        return simpson_weights_array(length(plan.chi)),(last(plan.chi)-first(plan.chi))/(length(plan.chi)-1)
    end
    return plan.weights,1.0
end

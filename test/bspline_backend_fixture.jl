# Manufactured PRODUCTION-GEOMETRY direct spline banks (22 ell x 128 chi x 64 R x count) with analytic,
# nonconstant tensors and a private dataset id. They prove plumbing and AD, never physical accuracy.
using Blast
using Blast.AbstractCosmologicalEmulators: CubicBSplinePlan, bspline_basis, knot_vector, basis_matrix

function manufactured_grid(n, power)
    logsites = log10(5e-5) .+ (log10(16.0) - log10(5e-5)) .* range(0, 1; length=n) .^ power
    sites = 10 .^ logsites; sites[1] = 5e-5; sites[end] = 16.0
    logsites = log10.(sites)
    knots = knot_vector(bspline_basis(CubicBSplinePlan(logsites, logsites; extrapolation=:throw)))
    return Blast.SplineGrid(sites, logsites, knots)
end

function manufactured_direct_bank(; n=9, power=1.25, offset=0.0,
        axes=(ell=Blast.ℓ_nonlimber,chi=Blast.χ,R=Blast.R))
    grid = manufactured_grid(n, power); sites, logsites, knots = grid.sites, grid.logsites, grid.knots
    nℓ, nχ, nR = length(axes.ell), length(axes.chi), length(axes.R)
    tensor(g) = [(1 + 0.1g + offset) * exp(-0.02i) * (1 + 0.3sin(0.05j + 0.1g)) * (0.5 + 0.01k) * cos(0.3c + 0.07g)
                 for i in 1:nℓ, j in 1:nχ, k in 1:nR, c in 1:n]
    tensors = NamedTuple{Blast._BANK_GEOMETRIES}(Tuple(tensor(g) for g in 1:8))
    bank = Blast.TensorBank(tensors, logsites; basis=:bspline_cubic_notaknot, dataset_id="private:manufactured-direct-$n-$offset",axes)
    return Blast.DirectSplineBank(bank, grid)
end

# Manufactured production-geometry RATIO bank: 8 TT + 4 T + 1 PHI deliberately DIFFERENT nonconstant arrays and
# fixed positive analytic q (generator-style phi, smooth Pmm). Algebra/AD evidence only, never physics.
function manufactured_ratio_bank(; n=9, power=1.35, offset=0.0)
    grid = manufactured_grid(n, power)
    nℓ, nχ, nR = length(Blast.ℓ_nonlimber), length(Blast.χ), length(Blast.R)
    tensor(g) = [(1 + 0.07g + offset) * exp(-0.015i) * (1 + 0.25sin(0.04j + 0.13g)) * (0.4 + 0.012k) * sin(0.21c + 0.05g + 0.3)
                 for i in 1:nℓ, j in 1:nχ, k in 1:nR, c in 1:n]
    TT = NamedTuple{Blast._BANK_GEOMETRIES}(Tuple(tensor(g) for g in 1:8))
    T = NamedTuple{(:T_2_00, :T_0_00, :T_2_02, :T_2_20)}(Tuple(tensor(10 + g) for g in 1:4))
    PHI = tensor(20)
    pmm = [fixture_pmm(k, 0.0) for k in grid.sites]
    q = Blast._normalizers(pmm, Blast._fiducial_phi(grid.sites, 3.044, 0.9649))
    sig = Blast.bank_signature_digest(:bspline_cubic_notaknot, "private:manufactured-ratio-$n-$offset", grid.logsites)
    return Blast.RatioSplineBank(TT, T, PHI, grid, q, sig)
end

# Smooth positive analytic Pmm, as in the characterization workload.
fixture_pmm(k, z) = 2.0e4 * (k / 0.2)^0.96 / (1 + (k / 0.2)^2)^2 / (1 + z)^2
fixture_pk(grid) = [fixture_pmm(k, z) for z in grid.z, k in grid.k]

# Literal contraction w[i,j,k] = sum_l c[l,(j),(k)] T[i,j,k,l] (explicit loops; independent of Tullio/adapters).
function literal_w(c, T)
    nℓ, nχ, nR, n = size(T)
    w = zeros(promote_type(eltype(c), eltype(T)), nℓ, nχ, nR)
    for k in 1:nR, j in 1:nχ, i in 1:nℓ
        s = zero(eltype(w))
        for l in 1:n
            s += (c isa AbstractVector ? c[l] : c isa AbstractMatrix ? c[l, j] : c[l, j, k]) * T[i, j, k, l]
        end
        w[i, j, k] = s
    end
    return w
end

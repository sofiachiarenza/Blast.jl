# Unfolded (χ₁, χ₂) check of the δ×PNG transfer-leg ordering (issue #38).
#
#   julia --project=. -t auto validation/png_external/unfolded_check.jl
#
# The δ×PNG cross-spectrum carries one transfer function, on the δ leg:
#
#   C^{δP}_ij = ∫∫ dχ₁ dχ₂ W^δ_i(χ₁) W^P_j(χ₂) D(χ₁) F_ℓ(χ₁, χ₂),
#   F_ℓ(χ₁, χ₂) = ∫ dk k² P(k) j_ℓ(kχ₁) j_ℓ(kχ₂)   (symmetric).
#
# The reference integrates this over the full square, with no folding, so it
# makes no assumption about which projected tensor goes on which triangle.
#
# Blast folds onto (χ, R) and feeds pre-contracted tensors to its Cℓ helpers.
# From prepare_pk_workspace / _PROJECTED_MATTER_COMPONENTS:
#
#   w_ϕT[ℓ, χ, R]    ↔ T(k, Rχ)   (P_ϕT[k,i,j] = P_ϕ[k] T_χR[i,j,k])
#   w_ϕT_R1[ℓ, χ, R] ↔ T(k, χ)    (coefs[:, :, end], i.e. R = 1)
#   w_ϕTT[ℓ, χ, R]   ↔ T(k, χ) T(k, Rχ)
#
# With linear growth T(k, χ) = D(χ) T₀(k), those are exactly D(Rχ) F, D(χ) F
# and D(χ) D(Rχ) F evaluated at (χ, Rχ). The script builds them on Blast's
# production χ and R grids and runs Blast's own fused contraction helpers with
# Blast's own quadrature weights. Prefactors are set to 1 on both sides.
#
# Three comparisons:
#   control  δ×δ with w_ϕTT (symmetric, no ordering question). This measures
#            the folded-quadrature floor on these kernels.
#   old      pre-fix call sites: δ×PNG with w_ϕT_R1, PNG×δ with w_ϕT, each
#            through the one-tensor helper (same tensor on both triangles).
#   fixed    current call sites (issue #38): two-tensor helper with
#            (w_ϕT_R1, w_ϕT) for δ×PNG and (w_ϕT, w_ϕT_R1) for PNG×δ.
# The same three are repeated with D ≡ 1, where old and fixed must coincide
# because the PNG tensors become symmetric. test/test_png_transfer_ordering.jl
# runs the same comparison through the production get_Cℓ call sites.

using Blast, SpecialFunctions, Printf, LinearAlgebra

# ----------------------------------------------------------------------------
# Toy linear model. Shapes are illustrative; nothing below depends on them
# being realistic, only on F being symmetric and D evolving.
# ----------------------------------------------------------------------------
const ℓs = [2, 10, 30]

growth(χ) = exp(-χ / 4000)                                    # D(χ)
pk_cross(k) = inv(k) / (1 + (k / 0.02)^2) * exp(-(k / 0.1)^2)  # ~ k² T₀ P_Φ shape

gauss(χ, μ, σ) = exp(-(χ - μ)^2 / (2σ^2))
const bins = ((1500.0, 200.0), (3000.0, 300.0))               # unequal bins
W_δ(i, χ)   = (1 + χ / 3000) * gauss(χ, bins[i]...)            # b(χ) n_i(χ)
W_PNG(i, χ) = (0.5 + χ / 1500) * gauss(χ, bins[i]...)          # f_NL b_Φ(χ) n_i(χ)
const nbins = length(bins)

# k quadrature defining F. Only consistency matters: the reference and the
# folded tensors evaluate the same discrete sum at different (χ₁, χ₂).
const k = collect(2e-4:2e-4:0.4)
const wk = fill(step(2e-4:2e-4:0.4), length(k)) .* k .^ 2 .* pk_cross.(k)

bessel_row(ℓ, χ) = sphericalbesselj.(ℓ, k .* χ)

# ----------------------------------------------------------------------------
# Unfolded reference: 2-D Simpson on a fine uniform grid covering all kernels.
# ----------------------------------------------------------------------------
const n_ref = 1041
const χ_ref = collect(range(300.0, 5500.0; length=n_ref))
const s_ref = let h = step(range(300.0, 5500.0; length=n_ref))
    s = [isodd(a) ? 2.0 : 4.0 for a in 1:n_ref]
    s[1] = s[end] = 1.0
    s .* h ./ 3
end

function F_unfolded(ℓ)
    J = Matrix{Float64}(undef, n_ref, length(k))
    Threads.@threads for a in 1:n_ref
        J[a, :] = bessel_row(ℓ, χ_ref[a])
    end
    return J * Diagonal(wk) * J'
end

# ∫∫ dχ₁ dχ₂ A_i(χ₁) B_j(χ₂) T(χ₁, χ₂) F(χ₁, χ₂)
function unfolded(F, A, B, transfer)
    C = zeros(nbins, nbins)
    for i in 1:nbins, j in 1:nbins
        acc = 0.0
        for b in 1:n_ref, a in 1:n_ref
            χ₁, χ₂ = χ_ref[a], χ_ref[b]
            acc += s_ref[a] * s_ref[b] * A(i, χ₁) * B(j, χ₂) * transfer(χ₁, χ₂) * F[a, b]
        end
        C[i, j] = acc
    end
    return C
end

# ----------------------------------------------------------------------------
# Folded Blast side on the production grids.
# ----------------------------------------------------------------------------
const χ = Blast.χ
const R = Blast.R
const nχ, nR = length(χ), length(R)
const Δχ = (last(χ) - first(χ)) / (nχ - 1)
const w_χ = Blast.simpson_weights_array(nχ)
const w_R = Blast.get_clencurt_weights_R_integration(2nR + 1)
const prefactor = ones(length(ℓs))

# F(χ_n, R_m χ_n) for every ℓ: (nℓ, nχ, nR)
function F_folded()
    F = zeros(length(ℓs), nχ, nR)
    for (iℓ, ℓ) in enumerate(ℓs)
        Threads.@threads for n in 1:nχ
            v = wk .* bessel_row(ℓ, χ[n])
            for m in 1:nR
                F[iℓ, n, m] = dot(bessel_row(ℓ, R[m] * χ[n]), v)
            end
        end
    end
    return F
end

# Blast's get_kernel_array layout: W[b, n, m] = W_b(R_m χ_n)
kernel_array(W) = [W(b, R[m] * χ[n]) for b in 1:nbins, n in 1:nχ, m in 1:nR]

one_tensor(WA, WB, pmj) =
    Blast._compute_Cℓ_fused_tullio(WA, WB, pmj, w_χ, w_R, prefactor, Δχ, χ)
two_tensor(WA, WB, pmj_tri1, pmj_tri2) =
    Blast._compute_Cℓ_rsd_fused_tullio(WA, WB, pmj_tri1, pmj_tri2, w_χ, w_R, prefactor, Δχ, χ)

# ----------------------------------------------------------------------------
# Run
# ----------------------------------------------------------------------------
relerr(x, ref) = abs(x - ref) / abs(ref)
maxnorm(x, ref) = maximum(abs, x .- ref) / maximum(abs, ref)

function run_case(label, D, Ffold, Frefs)
    Dχ  = reshape(D.(χ), 1, nχ, 1)
    DRχ = reshape([D(R[m] * χ[n]) for n in 1:nχ, m in 1:nR], 1, nχ, nR)
    w_ϕTT, w_ϕT, w_ϕT_R1 = Ffold .* Dχ .* DRχ, Ffold .* DRχ, Ffold .* Dχ

    Wδ, WP = kernel_array(W_δ), kernel_array(W_PNG)
    blast_dd = one_tensor(Wδ, Wδ, w_ϕTT)
    old   = one_tensor(Wδ, WP, w_ϕT_R1) .+ one_tensor(WP, Wδ, w_ϕT)
    fixed = two_tensor(Wδ, WP, w_ϕT_R1, w_ϕT) .+ two_tensor(WP, Wδ, w_ϕT, w_ϕT_R1)

    println("\n== $label ==")
    @printf("%4s %5s %13s %11s %11s %11s\n", "ℓ", "(i,j)", "ref δP+Pδ", "control", "old", "fixed")
    worst = (control=0.0, old=0.0, fixed=0.0)
    for (iℓ, ℓ) in enumerate(ℓs)
        F = Frefs[iℓ]
        ref_dd = unfolded(F, W_δ, W_δ, (a, b) -> D(a) * D(b))
        # T on the δ leg: χ₁ for δ×PNG, χ₂ for PNG×δ
        ref_png = unfolded(F, W_δ, W_PNG, (a, b) -> D(a)) .+
                  unfolded(F, W_PNG, W_δ, (a, b) -> D(b))
        for i in 1:nbins, j in 1:nbins
            @printf("%4d %5s %13.5e %11.2e %11.2e %11.2e\n", ℓ, "($i,$j)", ref_png[i, j],
                    relerr(blast_dd[iℓ, i, j], ref_dd[i, j]),
                    relerr(old[iℓ, i, j], ref_png[i, j]),
                    relerr(fixed[iℓ, i, j], ref_png[i, j]))
        end
        worst = (control = max(worst.control, maxnorm(blast_dd[iℓ, :, :], ref_dd)),
                 old     = max(worst.old,     maxnorm(old[iℓ, :, :], ref_png)),
                 fixed   = max(worst.fixed,   maxnorm(fixed[iℓ, :, :], ref_png)))
    end
    @printf("max|Δ|/max|ref| over bins, worst ℓ:  control %.2e   old %.2e   fixed %.2e\n",
            worst.control, worst.old, worst.fixed)
    return worst
end

println("threads: ", Threads.nthreads())
t = @elapsed begin
    Ffold = F_folded()
    Frefs = [F_unfolded(ℓ) for ℓ in ℓs]
end
@printf("F tensors built in %.1f s\n", t)

evolving = run_case("D(χ) = exp(-χ/4000)", growth, Ffold, Frefs)
frozen   = run_case("D ≡ 1 (no transfer evolution)", _ -> 1.0, Ffold, Frefs)

println("\nVerdict:")
floor_ = 10 * evolving.control
for (name, err) in (("old", evolving.old), ("fixed", evolving.fixed))
    println("  ", rpad(name, 9), err <= floor_ ? "matches" : "DOES NOT match",
            " the unfolded reference (", @sprintf("%.2e", err),
            " vs 10× control floor ", @sprintf("%.2e", floor_), ")")
end

using LinearAlgebra

isdefined(@__MODULE__, :get_test_bg) || include("test_fixtures.jl")

# Mixed-PNG transfer-leg ordering (issue #38).
#
# A X×PNG pair carries one transfer function, on the non-PNG leg X. Folding
# the (χ₁, χ₂) integral onto (χ, R) therefore needs ϕT_R1 (T at χ) on the
# triangle where X sits at χ and ϕT (T at Rχ) on the other. These tests fill
# the projected-matter workspace with synthetic tensors built from a symmetric
# F(χ₁, χ₂), a growth factor D(χ) and an RSD derivative weight h(χ) on the RSD
# leg, run the production get_Cℓ paths on real components, and compare with
# an unfolded integral over the full (χ₁, χ₂) square in which D and h are put
# on the right leg explicitly. Nothing in the reference knows about the
# folding.

# Spherical Bessel j₂, with a series below x = 0.1 where the closed form cancels.
_j2(x) = x < 0.1 ? x^2 / 15 - x^4 / 210 + x^6 / 7560 :
                   (3 / x^2 - 1) * sin(x) / x - 3 * cos(x) / x^2

function _simpson_weights(x::AbstractRange)
    n = length(x)
    isodd(n) || throw(ArgumentError("Simpson needs an odd number of nodes"))
    s = [isodd(a) ? 2.0 : 4.0 for a in 1:n]
    s[1] = s[end] = 1.0
    return s .* step(x) ./ 3
end

@testset "Mixed-PNG transfer-leg ordering" begin
    bg = get_test_bg()
    χ, R = Blast.χ, Blast.R
    nχ, nR, nℓ = length(χ), length(R), length(Blast.ℓ_nonlimber)

    # ---- toy linear model ---------------------------------------------------
    D(x) = exp(-x / 4000)       # T(k, χ) = D(χ) T₀(k)
    h(x) = 1 + x / 5000         # RSD derivative weight, on the RSD leg
    ks = 2e-4:2e-4:0.4
    wk = step(ks) .* ks .^ 2 .* (inv.(ks) ./ (1 .+ (ks ./ 0.02) .^ 2) .* exp.(-(ks ./ 0.1) .^ 2))
    bessel(xs) = [_j2(k * x) for x in xs, k in ks]

    # F(χ_n, R_m χ_n) on Blast's folded grid
    J_χ = bessel(χ)
    F_fold = zeros(nχ, nR)
    for n in 1:nχ
        F_fold[n, :] = bessel(R .* χ[n]) * (wk .* J_χ[n, :])
    end
    χ_at = repeat(χ, 1, nR)         # position of triangle-1's A leg
    Rχ_at = χ .* R'                 # position of triangle-1's B leg

    lift(A) = repeat(reshape(A, 1, nχ, nR), nℓ, 1, 1)
    wrap(T, A) = T{Float64}(zeros(1, 1, 1, 1), lift(A))
    # ϕT ↔ T(k, Rχ), ϕT_R1 ↔ T(k, χ); _02 has the derivative at Rχ, _20 at χ.
    toy = (
        w_2_00_ϕT    = wrap(Blast.w_2_00_ϕT,    F_fold .* D.(Rχ_at)),
        w_2_00_ϕT_R1 = wrap(Blast.w_2_00_ϕT_R1, F_fold .* D.(χ_at)),
        w_0_00_ϕT    = wrap(Blast.w_0_00_ϕT,    F_fold .* D.(Rχ_at)),
        w_0_00_ϕT_R1 = wrap(Blast.w_0_00_ϕT_R1, F_fold .* D.(χ_at)),
        w_2_02_ϕT    = wrap(Blast.w_2_02_ϕT,    F_fold .* h.(Rχ_at) .* D.(Rχ_at)),
        w_2_02_ϕT_R1 = wrap(Blast.w_2_02_ϕT_R1, F_fold .* h.(Rχ_at) .* D.(χ_at)),
        w_2_20_ϕT    = wrap(Blast.w_2_20_ϕT,    F_fold .* h.(χ_at) .* D.(Rχ_at)),
        w_2_20_ϕT_R1 = wrap(Blast.w_2_20_ϕT_R1, F_fold .* h.(χ_at) .* D.(χ_at)),
    )
    W = Blast.ProjectedMatterDensity(; toy...)

    # ---- components: two unequal bins ---------------------------------------
    z = collect(0.0:0.02:3.6)
    n_bg = length(bg.z)
    gaussian_bins(centres) = reduce(vcat, (permutedims(exp.(-0.5 .* ((z .- c) ./ σ) .^ 2))
                                           for (c, σ) in centres))
    nz_g = gaussian_bins(((0.5, 0.10), (1.2, 0.15)))
    nz_s = gaussian_bins(((0.8, 0.20), (1.4, 0.25)))
    bias = vcat(permutedims(1.2 .+ 0.1 .* bg.z), permutedims(1.6 .+ 0.1 .* bg.z))

    G = Blast.GalaxyClustering(
        δ   = Blast.NumberCounts(nz=nz_g, z=z, bias=bias),
        RSD = Blast.RedshiftSpaceDistortions(nz=nz_g, z=z, growth_rate=bg.f),
        μ   = Blast.MagnificationBias(nz=nz_g, z=z, s=fill(0.3, 2, n_bg)),
        PNG = Blast.PrimordialNonGaussianity(nz=nz_g, z=z, bias=bias, f_NL=1.0, p=1.0),
    )
    L = Blast.WeakLensing(γ  = Blast.CosmicShear(nz=nz_s, z=z),
                          IA = Blast.IntrinsicAlignment(nz=nz_s, z=z))
    C = Blast.CMB(κ=Blast.CMBLensing(), ISW=Blast.IntegratedSachsWolfe())
    Blast.evaluate_components!(G, bg)
    Blast.evaluate_components!(L, bg)
    Blast.evaluate_components!(C, bg)

    # ---- unfolded reference -------------------------------------------------
    χ_ref = range(1.0, last(χ); length=2801)
    s_ref = _simpson_weights(χ_ref)
    J_ref = bessel(χ_ref)
    F_ref = J_ref * Diagonal(wk) * J_ref'
    lensing_like = Union{Blast.CosmicShear, Blast.IntrinsicAlignment,
                         Blast.MagnificationBias, Blast.CMBLensing}

    # Same interpolant get_kernel_array uses, evaluated off the folded grid.
    function kernel_ref(P)
        K = collect(Blast.akima_interpolation(collect(P.Kernel'), bg.χ, collect(χ_ref))')
        P isa lensing_like && (K ./= χ_ref' .^ 2)
        return K
    end

    # pref_ℓ ∫∫ dχ₁ dχ₂ A_i(χ₁) B_j(χ₂) F(χ₁, χ₂) × D on the non-PNG leg
    # × h on the RSD leg (if any).
    function unfolded(A, B)
        x = collect(χ_ref)
        weight(P) = (P isa Blast.PrimordialNonGaussianity ? ones(length(x)) : D.(x)) .*
                    (P isa Blast.RedshiftSpaceDistortions ? h.(x) : ones(length(x)))
        KA = kernel_ref(A) .* (s_ref .* weight(A))'
        KB = kernel_ref(B) .* (s_ref .* weight(B))'
        pair = KA * F_ref * KB'
        return reshape(Blast._pair_prefactor(A, B), :, 1, 1) .* reshape(pair, 1, size(pair)...)
    end
    maxnorm(x, ref) = maximum(abs, x .- ref) / maximum(abs, ref)

    # The pre-fix call sites miss the reference by 1.7e-2 to 9e-2 on every
    # pair here. The folded quadrature itself (checked with D ≡ 1, where there
    # is no ordering question) agrees to ~1e-5, except for the broad μ, γ
    # and κ kernels that reach down to χ → 0, where its floor is ~4e-4.
    # Each tolerance sits between that floor and the pre-fix error.
    tol = 2e-4
    tol_broad = 2e-3
    broad = Union{Blast.MagnificationBias, Blast.CosmicShear, Blast.CMBLensing}
    tolerance(A, B) = (A isa broad || B isa broad) ? tol_broad : tol

    @testset "component methods: $(nameof(typeof(A))) × $(nameof(typeof(B)))" for (A, B) in (
        (G.δ, G.PNG), (G.PNG, G.δ),
        (G.RSD, G.PNG), (G.PNG, G.RSD),
        (G.μ, G.PNG), (G.PNG, G.μ),
        (G.PNG, L.γ), (G.PNG, L.IA),
        (C.κ, G.PNG), (C.ISW, G.PNG),
    )
        @test maxnorm(Blast.get_Cℓ(A, B, W, bg), unfolded(A, B)) < tolerance(A, B)
    end

    @testset "reverse pairs are index transposes" begin
        for (A, B) in ((G.δ, G.PNG), (G.RSD, G.PNG), (G.μ, G.PNG))
            AB = Blast.get_Cℓ(A, B, W, bg)
            BA = Blast.get_Cℓ(B, A, W, bg)
            @test BA ≈ permutedims(AB, (1, 3, 2)) rtol = 1e-12
        end
    end

    @testset "the old one-tensor contraction is detected" begin
        # Same tensor on both triangles: what main used before the fix.
        old = Blast.compute_Cℓ(G.δ, G.PNG, W.w_2_00_ϕT_R1, bg)
        @test maxnorm(old, unfolded(G.δ, G.PNG)) > 100 * tol
    end

    @testset "cached GG paths" begin
        integ = Blast._prepare_nonlimber_integration(bg)
        K(P) = Blast.get_kernel_array(P, bg)
        δfNL, fNLδ, fNLRSD, RSDfNL, μfNL, fNLμ, _ = Blast._gc_nonlimber_png_pairs(
            G, W, integ, K(G.δ), K(G.RSD), K(G.μ), K(G.PNG))
        @test maxnorm(δfNL, unfolded(G.δ, G.PNG)) < tol
        @test maxnorm(fNLδ, unfolded(G.PNG, G.δ)) < tol
        @test maxnorm(RSDfNL, unfolded(G.RSD, G.PNG)) < tol
        @test maxnorm(fNLRSD, unfolded(G.PNG, G.RSD)) < tol
        @test maxnorm(μfNL, unfolded(G.μ, G.PNG)) < tol_broad
        @test maxnorm(fNLμ, unfolded(G.PNG, G.μ)) < tol_broad

    end

    @testset "cached GL path" begin
        # GC×WL with only δ + PNG on the galaxy side. Every non-PNG tensor is
        # zero, so the PNG-on minus PNG-off difference isolates Cℓ_fNLγ +
        # Cℓ_fNLI; the Limber parts cancel and _finalize_Cℓ_parts is linear.
        G_gl = Blast.GalaxyClustering(δ=G.δ, PNG=G.PNG)
        _, Plans = Blast.SetUp(G_gl, L)
        pk = ones(length(Blast.z_lin), length(Blast.k_cheb))
        pk_limber = ones(length(Blast.z_cheb), length(Blast.k_limber))
        PS = Blast.prepare_pk_workspace(Plans, pk, pk_limber, pk_limber, bg)

        zero_ϕTT = wrap(Blast.w_0_00_ϕTT, zeros(nχ, nR))
        W_on = Blast.ProjectedMatterDensity(w_0_00_ϕTT=zero_ϕTT,
                                            w_0_00_ϕT=toy.w_0_00_ϕT, w_0_00_ϕT_R1=toy.w_0_00_ϕT_R1)
        W_off = Blast.ProjectedMatterDensity(w_0_00_ϕTT=zero_ϕTT,
                                             w_0_00_ϕT=wrap(Blast.w_0_00_ϕT, zeros(nχ, nR)),
                                             w_0_00_ϕT_R1=wrap(Blast.w_0_00_ϕT_R1, zeros(nχ, nR)))

        ℓ = Blast.full_ℓ_range
        on = Blast.get_Cℓ(ℓ, G_gl, L, PS, W_on, bg, Plans)
        off = Blast.get_Cℓ(ℓ, G_gl, L, PS, W_off, bg, Plans)

        n_limber = length(Blast.full_ℓ_range) - nℓ
        expected = Blast._finalize_Cℓ_parts(unfolded(G.PNG, L.γ) .+ unfolded(G.PNG, L.IA), 0.0,
                                            zeros(n_limber, 2, 2), ℓ, 2, 2, Plans)
        @test maxnorm(on .- off, expected) < tol_broad
    end
end

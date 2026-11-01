# Phase 6 mandatory canary: non-square (2 ell, 3 chi, 5 R, 7 coefficients) ratio contraction with deliberately
# different, nonseparable TT/T/PHI tensors. Independent literal sums for all 17 outputs; the legacy shared-TT
# fusion would be wrong for ratio banks and is shown to be detectably wrong. Selector: `ratio-canary`.
using Test, Blast, LinearAlgebra
using AbstractCosmologicalEmulators, FastChebInterp      # same imports as the full suite
isdefined(@__MODULE__, :literal_w) || include("bspline_backend_fixture.jl")

const CNL, CNX, CNR, CN = 2, 3, 5, 7
canary_array(seed) = [sin(0.7i + 0.3j + 0.1seed) * cos(1.1k - 0.2seed) + 0.1c * i * j + 0.01seed * k * c + seed for i in 1:CNL, j in 1:CNX, k in 1:CNR, c in 1:CN]
canary_coef3(seed) = [cos(0.3l + 0.5j - 0.2k + seed) + 0.05l * j + 0.02seed for l in 1:CN, j in 1:CNX, k in 1:CNR]
const CANARY_SIG = Blast.BankSignature(:bspline_cubic_notaknot, CN, (0x1, 0x2, 0x3, 0x4))

function canary_ratio_bank()
    TT = NamedTuple{Blast._BANK_GEOMETRIES}(Tuple(canary_array(10g) for g in 1:8))
    T = NamedTuple{(:T_2_00, :T_0_00, :T_2_02, :T_2_20)}(Tuple(canary_array(100 + 10g) for g in 1:4))
    grid = Blast.SplineGrid(collect(range(1.0, 2.0; length=CN)), collect(range(0.0, 0.3; length=CN)), zeros(CN + 4))
    q = Blast.FiducialNormalizers(ones(CN), ones(CN), ones(CN))
    return Blast.RatioSplineBank(TT, T, canary_array(200), grid, q, CANARY_SIG)
end

canary_W(bank; only=nothing) = Blast.ProjectedMatterDensity(;
    (name => getfield(Blast, name)(bank) for (name, _, _, _) in Blast._PROJECTED_MATTER_COMPONENTS
     if isnothing(only) || name in only)...)
canary_pk(seed) = Blast.PowerSpectrum(Blast.cϕTT(canary_coef3(seed)), Blast.cϕT(canary_coef3(seed + 1)),
                                      Blast.cϕ([cos(0.9l + seed) + 0.3 for l in 1:CN]), rand(3, 4), rand(3, 4), CANARY_SIG)
canary_tensor(bank, tfield, src) = src === :ϕTT ? getfield(bank.TT, tfield) : src === :ϕT ? bank.T[tfield] : bank.PHI

function canary_literal(bank, Pk, name, tfield, src, r1)
    c = src === :ϕTT ? Pk.cϕTT.coefs : src === :ϕT ? Pk.cϕT.coefs : Pk.cϕ.coefs
    return literal_w(r1 ? c[:, :, end] : c, canary_tensor(bank, tfield, src))      # R1 slices the SOURCE coefficients
end

@testset "Ratio canary: 17 outputs, shared-tensor guard, poisoned buffers, changed inputs" begin
    bank = canary_ratio_bank()
    W = canary_W(bank)
    @test length(fieldnames(typeof(W))) == 17 && all(!isnothing, getfield.(Ref(W), fieldnames(typeof(W))))
    # Components hold their ACTUAL family array; R1 holds the T-family tensor (alias of the ϕT component).
    for (name, tfield, src, r1) in Blast._PROJECTED_MATTER_COMPONENTS
        @test getfield(W, name).T̃ === canary_tensor(bank, tfield, src)
    end
    @test W.w_2_00_ϕT_R1.T̃ === W.w_2_00_ϕT.T̃ !== W.w_2_00_ϕTT.T̃
    @test !Blast._shares_tensor(W.w_2_00_ϕTT, W.w_2_00_ϕT, W.w_2_00_ϕT_R1, W.w_2_00_ϕ)
    @test !Blast._shares_tensor(W.w_0_00_ϕTT, W.w_0_00_ϕT, W.w_0_00_ϕT_R1)

    ws = Blast.allocate_compute_w(W, canary_pk(0.0))
    for seed in (0.0, 0.7, 0.0)                                           # changed inputs, then the first again
        Pk = canary_pk(seed)
        Wf = Blast.compute_w(W, Pk)
        for f in fieldnames(typeof(ws)); fill!(getfield(ws, f).w, NaN); end      # poisoned reused buffers
        Blast.compute_w!(ws, Pk)
        for (name, tfield, src, r1) in Blast._PROJECTED_MATTER_COMPONENTS
            ref = canary_literal(bank, Pk, name, tfield, src, r1)
            @test getfield(Wf, name).w ≈ ref rtol=1e-12
            @test all(isfinite, getfield(ws, name).w)
            @test getfield(ws, name).w ≈ ref rtol=1e-12
        end
    end

    # The legacy fused scan reads the TT tensor for every family: demonstrably wrong for ratio banks.
    Pk = canary_pk(0.0)
    w_tt, w_t, w_t1, w_phi = (zeros(CNL, CNX, CNR) for _ in 1:4)
    Blast._w_fused4_inplace!(w_tt, w_t, w_t1, w_phi, Pk.cϕTT.coefs, Pk.cϕT.coefs, Pk.cϕT.coefs[:, :, end], Pk.cϕ.coefs, bank.TT.T_2_00)
    correct_t = literal_w(Pk.cϕT.coefs, bank.T.T_2_00)
    @test w_tt ≈ literal_w(Pk.cϕTT.coefs, bank.TT.T_2_00) rtol=1e-12          # TT output is right ...
    @test norm(w_t - correct_t) > 0.1 * norm(correct_t)                       # ... T, R1 and PHI are not
    @test norm(w_t1 - literal_w(Pk.cϕT.coefs[:, :, end], bank.T.T_2_00)) > 0.1 * norm(w_t1)
    @test norm(w_phi - literal_w(Pk.cϕ.coefs, bank.PHI)) > 0.1 * norm(w_phi)

    # Absent fields / single outputs: groups that are incomplete use the single scans.
    for only in ((:w_2_00_ϕT,), (:w_2_00_ϕTT, :w_2_00_ϕT), (:w_0_00_ϕTT, :w_0_00_ϕT_R1), (:w_2_22_ϕTT, :w_2_20_ϕT, :w_2_20_ϕT_R1))
        Ws = canary_W(bank; only)
        wsub = Blast.allocate_compute_w(Ws, Pk); foreach(f -> isnothing(getfield(wsub, f)) || fill!(getfield(wsub, f).w, NaN), fieldnames(typeof(wsub)))
        Blast.compute_w!(wsub, Pk)
        for (name, tfield, src, r1) in Blast._PROJECTED_MATTER_COMPONENTS
            name in only || continue
            @test getfield(wsub, name).w ≈ canary_literal(bank, Pk, name, tfield, src, r1) rtol=1e-12
        end
    end

    # Direct/Chebyshev-style SHARED tensors still take the fused path and are exact.
    shared = Blast.TensorBank(bank.TT.T_2_00, bank.TT.T_0_00, bank.TT.T_minus2_00, bank.TT.T_0_02, bank.TT.T_0_20,
                              bank.TT.T_2_02, bank.TT.T_2_20, bank.TT.T_2_22, collect(range(0.0, 0.3; length=CN)), CANARY_SIG)
    Wd = canary_W(shared)
    @test Blast._shares_tensor(Wd.w_2_00_ϕTT, Wd.w_2_00_ϕT, Wd.w_2_00_ϕT_R1, Wd.w_2_00_ϕ)
    @test Blast._shares_tensor(Wd.w_2_02_ϕTT, Wd.w_2_02_ϕT, Wd.w_2_02_ϕT_R1)
    wsd = Blast.allocate_compute_w(Wd, Pk); foreach(f -> fill!(getfield(wsd, f).w, NaN), fieldnames(typeof(wsd)))
    Blast.compute_w!(wsd, Pk)
    for (name, tfield, src, r1) in Blast._PROJECTED_MATTER_COMPONENTS
        c = src === :ϕTT ? Pk.cϕTT.coefs : src === :ϕT ? Pk.cϕT.coefs : Pk.cϕ.coefs
        @test getfield(wsd, name).w ≈ literal_w(r1 ? c[:, :, end] : c, getfield(bank.TT, tfield)) rtol=1e-12
    end
    # Missing T-family / PHI geometries are rejected, never replaced by dummies.
    @test_throws ArgumentError Blast._ratio_tensor(bank, Val(:ϕT), :T_2_22)
    @test_throws ArgumentError Blast._ratio_tensor(bank, Val(:ϕ), :T_0_00)
end

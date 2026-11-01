# PowerSpectrum constructor promotion (mixed-scalar public AD cases). Focused selector: `promotion`.
using Test, Blast, ForwardDiff, LinearAlgebra
using AbstractCosmologicalEmulators, FastChebInterp      # same imports as the full suite
isdefined(@__MODULE__, :build_characterization_state) || include("characterization_workload.jl")
isdefined(@__MODULE__, :get_test_cosmo) || include("test_fixtures.jl")

const PromoDual = ForwardDiff.Dual{Nothing,Float64,2}
mkdual(v, p1, p2) = PromoDual(v, ForwardDiff.Partials{2,Float64}((p1, p2)))
dual_like(A) = reshape([mkdual(a, 0.3 * i, -0.1 * i) for (i, a) in enumerate(A)], size(A))

@testset "PowerSpectrum promotion: constructor boundary" begin
    nk = length(Blast.k_cheb)
    sig = Blast.LEGACY_BANK_SIGNATURE
    a, b, c = randn(nk, 2, 2), randn(nk, 2, 2), randn(nk)
    cTT, cT, cP = Blast.cϕTT(a), Blast.cϕT(b), Blast.cϕ(c)
    L1, L2 = rand(3, 4), rand(3, 4)

    # All Float64: same objects, no conversion; the 5-argument legacy form still tags the legacy signature.
    ps = Blast.PowerSpectrum(cTT, cT, cP, L1, L2, sig)
    @test ps isa Blast.PowerSpectrum{Float64} && ps.cϕTT === cTT && ps.cϕT === cT && ps.cϕ === cP
    @test ps.ΔP_limber === L1 && ps.Pδ_limber === L2
    Blast.PowerSpectrum(cTT, cT, cP, L1, L2, sig)                       # warm
    @test (@allocated Blast.PowerSpectrum(cTT, cT, cP, L1, L2, sig)) < 512     # no hidden copies
    @test Blast.PowerSpectrum(cTT, nothing, nothing, L1, L2).signature == sig
    @test_throws ArgumentError Blast.PowerSpectrum(Blast.cϕTT(randn(5, 2, 2)), nothing, nothing, L1, L2)   # not legacy sites

    # Fixed Float64 phi with Dual TT/T (Float64 background, Dual Pmm): promoted, zero partials, same values.
    dTT, dT = Blast.cϕTT(dual_like(a)), Blast.cϕT(dual_like(b))
    ps = Blast.PowerSpectrum(dTT, dT, cP, dual_like(L1), dual_like(L2), sig)
    @test ps isa Blast.PowerSpectrum{PromoDual}
    @test ps.cϕTT === dTT && ps.cϕT === dT                               # already the right T: identity kept
    @test eltype(ps.cϕ.coefs) == PromoDual && ForwardDiff.value.(ps.cϕ.coefs) == c
    @test all(x -> all(iszero, ForwardDiff.partials(x)), ps.cϕ.coefs)    # constant family: zero derivative, not fake values
    # Inactive families stay Nothing.
    ps = Blast.PowerSpectrum(dTT, nothing, nothing, dual_like(L1), dual_like(L2), sig)
    @test isnothing(ps.cϕT) && isnothing(ps.cϕ) && ps isa Blast.PowerSpectrum{PromoDual}

    # Limber-only Dual: constant Float64 non-Limber coefficients are promoted, the derivative is NOT stripped.
    dL1, dL2 = dual_like(L1), dual_like(L2)
    ps = Blast.PowerSpectrum(cTT, cT, cP, dL1, dL2, sig)
    @test ps isa Blast.PowerSpectrum{PromoDual}
    @test ps.ΔP_limber === dL1 && ps.Pδ_limber === dL2 && any(x -> any(!iszero, ForwardDiff.partials(x)), ps.Pδ_limber)
    for (f, orig) in ((ps.cϕTT, a), (ps.cϕT, b), (ps.cϕ, c))
        @test eltype(f.coefs) == PromoDual && ForwardDiff.value.(f.coefs) == orig
        @test all(x -> all(iszero, ForwardDiff.partials(x)), f.coefs)
    end
    @test_throws ArgumentError Blast.PowerSpectrum(Blast.cϕTT(complex.(a)), nothing, nothing, L1, L2, sig)
end

@testset "PowerSpectrum promotion: public prepare_pk_workspace (fixed Float64 background)" begin
    state = build_characterization_state()           # 17-component PNG workload, Chebyshev backend
    (; bg, plans, pmm) = state
    ref = Blast.prepare_pk_workspace(plans, pmm.pk_nonlimber, pmm.pk_limber_linear, pmm.pk_limber_nonlinear, bg)
    @test ref isa Blast.PowerSpectrum{Float64} && !isnothing(ref.cϕ)

    # Dual non-Limber Pmm only: previously a MethodError (cϕ stayed Float64).
    pk = [mkdual(v, i == 1 ? v : 0.0, 0.0) for (i, v) in enumerate(vec(pmm.pk_nonlimber))]
    pk = reshape(pk, size(pmm.pk_nonlimber))
    ps = Blast.prepare_pk_workspace(plans, pk, pmm.pk_limber_linear, pmm.pk_limber_nonlinear, bg)
    @test ps isa Blast.PowerSpectrum{PromoDual}
    @test ForwardDiff.value.(ps.cϕ.coefs) == ref.cϕ.coefs                              # same primal values
    @test all(x -> all(iszero, ForwardDiff.partials(x)), ps.cϕ.coefs)                  # constant primordial family
    @test any(x -> any(!iszero, ForwardDiff.partials(x)), ps.cϕTT.coefs)               # differentiated family keeps its derivative
    @test ForwardDiff.value.(ps.Pδ_limber) == ref.Pδ_limber

    # Limber-only Dual: fixed non-Limber Float64 Pmm; Limber matrices carry the derivative.
    ln = reshape([mkdual(v, v, 0.0) for v in vec(pmm.pk_limber_nonlinear)], size(pmm.pk_limber_nonlinear))
    ps = Blast.prepare_pk_workspace(plans, pmm.pk_nonlimber, pmm.pk_limber_linear, ln, bg)
    @test ps isa Blast.PowerSpectrum{PromoDual}
    @test any(x -> any(!iszero, ForwardDiff.partials(x)), ps.Pδ_limber)                # not stripped
    @test ForwardDiff.value.(ps.cϕTT.coefs) == ref.cϕTT.coefs
    @test all(x -> all(iszero, ForwardDiff.partials(x)), ps.cϕTT.coefs)
    # The functional contraction then runs with the promoted Dual coefficients.
    W = Blast.compute_w(state.weights_template, ps)
    @test W.w_2_00_ϕTT isa Blast.w_2_00_ϕTT{PromoDual}
end

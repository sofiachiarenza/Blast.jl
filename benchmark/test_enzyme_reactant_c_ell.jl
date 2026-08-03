using Blast
using Enzyme
using ForwardDiff
using Reactant
using Random

function main()
    Reactant.set_default_backend("cpu")
    Random.seed!(12345)
    nl, nχ, nR, nA, nB = 4, 5, 3, 2, 3
    WA = randn(nA, nχ, nR)
    WB = randn(nB, nχ, nR)
    pmj = randn(nl, nχ, nR)
    pref = randn(nl)
    wχ = rand(nχ)
    wR = rand(nR)
    χ = rand(nχ)
    Δχ = 0.7

    loss(p) = sum(Blast.reactant_nonlimber_c_ell(WA, WB, p, pref, wχ, wR, χ, Δχ))
    reference = ForwardDiff.gradient(loss, pmj)

    loss_dynamic(p, WA_, WB_, pref_, wχ_, wR_, χ_) = sum(
        Blast.reactant_nonlimber_c_ell(WA_, WB_, p, pref_, wχ_, wR_, χ_, Δχ),
    )
    grad_dynamic(p, WA_, WB_, pref_, wχ_, wR_, χ_) = Enzyme.gradient(
        Reverse,
        loss_dynamic,
        p,
        Const(WA_), Const(WB_), Const(pref_), Const(wχ_), Const(wR_), Const(χ_),
    )[1]
    WA_r, WB_r = Reactant.to_rarray(WA), Reactant.to_rarray(WB)
    pmj_r, pref_r = Reactant.to_rarray(pmj), Reactant.to_rarray(pref)
    wχ_r, wR_r, χ_r = map(Reactant.to_rarray, (wχ, wR, χ))
    compiled = Reactant.@compile sync=true grad_dynamic(pmj_r, WA_r, WB_r, pref_r, wχ_r, wR_r, χ_r)
    result = compiled(pmj_r, WA_r, WB_r, pref_r, wχ_r, wR_r, χ_r)
    Reactant.synchronize(result)
    println("Reactant=", pkgversion(Reactant))
    println("gradient_shape=", size(reference))
    println("gradient_error=", maximum(abs, Array(result) .- reference))
    println("gradient_norm=", maximum(abs, reference))
end

main()

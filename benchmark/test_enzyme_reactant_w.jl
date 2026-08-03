using Blast
using Enzyme
using ForwardDiff
using Reactant
using Random

function main()
    Reactant.set_default_backend("cpu")
    Random.seed!(97531)
    nl, ni, nj, nk = 3, 2, 4, 3
    cTT = randn(nl, nj, nk)
    cT = randn(nl, nj, nk)
    cϕ = randn(nl, nj)
    Ts = ntuple(_ -> randn(nl, ni, nj, nk), 8)

    loss(c) = begin
        r = Blast.reactant_compute_w(c, cT, cϕ, Ts...)
        sum(sum, r)
    end
    reference = ForwardDiff.gradient(loss, cTT)

    function grad_fun(c, cT_, cϕ_, T1, T2, T3, T4, T5, T6, T7, T8)
        loss_dynamic(c_, cT_, cϕ_, T1_, T2_, T3_, T4_, T5_, T6_, T7_, T8_) = begin
            r = Blast.reactant_compute_w(c_, cT_, cϕ_, T1_, T2_, T3_, T4_, T5_, T6_, T7_, T8_)
            sum(sum, r)
        end
        Enzyme.gradient(
            Reverse, loss_dynamic, c,
            Const(cT_), Const(cϕ_), Const(T1), Const(T2), Const(T3), Const(T4),
            Const(T5), Const(T6), Const(T7), Const(T8),
        )[1]
    end
    args = (
        Reactant.to_rarray(cTT), Reactant.to_rarray(cT), Reactant.to_rarray(cϕ),
        map(Reactant.to_rarray, Ts)...,
    )
    compiled = Reactant.@compile sync=true grad_fun(args...)
    result = compiled(args...)
    Reactant.synchronize(result)
    println("Reactant=", pkgversion(Reactant))
    println("gradient_shape=", size(reference))
    println("gradient_error=", maximum(abs, Array(result) .- reference))
    println("gradient_norm=", maximum(abs, reference))
end

main()

using AbstractCosmologicalEmulators
using Blast
using Enzyme
using ForwardDiff
using Reactant

const ROOT = "/home/marcobonici/Desktop/work/CosmologicalEmulators/blast_play/for_marco_copy"
include(joinpath(ROOT, "performance", "workload.jl"))

function main()
    Reactant.set_default_backend("cpu")
    workload = build_workload(ROOT, 8)
    plan = prepare_chebyshev_plan(
        minimum(Blast.k_cheb), maximum(Blast.k_cheb), length(Blast.k_cheb) - 1;
        size_nd=(length(Blast.k_cheb),), dim=1,
    )
    transform = Blast.reactant_chebyshev_matrix(plan)
    pk = workload.pk

    loss(p, M) = sum(Blast.reactant_prepare_nonlimber_spectrum(p, workload.bg, M)[1])
    reference = ForwardDiff.gradient(p -> loss(p, transform), pk)
    loss_enzyme(p, M) = Enzyme.gradient(Reverse, loss, p, Const(M))[1]
    grad_fun(p, M) = loss_enzyme(p, M)

    pk_r = Reactant.to_rarray(pk)
    M_r = Reactant.to_rarray(transform)
    compiled = Reactant.@compile sync=true grad_fun(pk_r, M_r)
    result = compiled(pk_r, M_r)
    Reactant.synchronize(result)
    println("Reactant=", pkgversion(Reactant))
    println("pk_shape=", size(pk), " gradient_shape=", size(reference))
    println("gradient_error=", maximum(abs, Array(result) .- reference))
    println("gradient_norm=", maximum(abs, reference))
end

main()

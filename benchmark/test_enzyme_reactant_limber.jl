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
    plan = workload.plans.plan_limber
    Kk, Kz = plan.K
    Mk = Blast.reactant_chebyshev_matrix(prepare_chebyshev_plan(
        minimum(plan.nodes[1]), maximum(plan.nodes[1]), Kk;
        size_nd=(Kk + 1,), dim=1,
    ))
    Mz = Blast.reactant_chebyshev_matrix(prepare_chebyshev_plan(
        minimum(plan.nodes[2]), maximum(plan.nodes[2]), Kz;
        size_nd=(Kz + 1,), dim=1,
    ))
    Tz = Blast.get_limber_coords_polynomials(plan, workload.bg.z)
    Tk = workload.plans.T_k_limber
    pk_nonlin = workload.pk_limber_nonlin
    loss(pk) = sum(Blast.reactant_prepare_limber(pk, pk_nonlin, workload.bg, Mk, Mz, Tz, Tk)[2])
    reference = ForwardDiff.gradient(loss, workload.pk_limber_lin)

    loss_dynamic(pk, pknl, Mk_, Mz_, Tz_, Tk_) = sum(
        Blast.reactant_prepare_limber(pk, pknl, workload.bg, Mk_, Mz_, Tz_, Tk_)[2],
    )
    grad_fun(pk, pknl, Mk_, Mz_, Tz_, Tk_) = Enzyme.gradient(
        Reverse, loss_dynamic, pk,
        Const(pknl), Const(Mk_), Const(Mz_), Const(Tz_), Const(Tk_),
    )[1]
    args = (
        Reactant.to_rarray(workload.pk_limber_lin),
        Reactant.to_rarray(pk_nonlin), Reactant.to_rarray(Mk), Reactant.to_rarray(Mz),
        Reactant.to_rarray(Tz), Reactant.to_rarray(Tk),
    )
    compiled = Reactant.@compile sync=true grad_fun(args...)
    result = compiled(args...)
    Reactant.synchronize(result)
    println("Reactant=", pkgversion(Reactant))
    println("pk_shape=", size(workload.pk_limber_lin))
    println("gradient_error=", maximum(abs, Array(result) .- reference))
    println("gradient_norm=", maximum(abs, reference))
end

main()

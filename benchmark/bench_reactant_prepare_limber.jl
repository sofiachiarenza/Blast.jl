using AbstractCosmologicalEmulators
using BenchmarkTools
using Blast
using Reactant

const ROOT = "/home/marcobonici/Desktop/work/CosmologicalEmulators/blast_play/for_marco_copy"
include(joinpath(ROOT, "performance", "workload.jl"))

function main()
    Reactant.set_default_backend("cpu")
    workload = build_workload(ROOT, 8)
    reference = prepare_spectrum(workload)
    plan = workload.plans.plan_limber
    Kk, Kz = plan.K
    plan_k = prepare_chebyshev_plan(
        minimum(plan.nodes[1]), maximum(plan.nodes[1]), Kk;
        size_nd=(Kk + 1,), dim=1,
    )
    plan_z = prepare_chebyshev_plan(
        minimum(plan.nodes[2]), maximum(plan.nodes[2]), Kz;
        size_nd=(Kz + 1,), dim=1,
    )
    transform_k = Blast.reactant_chebyshev_matrix(plan_k)
    transform_z = Blast.reactant_chebyshev_matrix(plan_z)
    T_z = Blast.get_limber_coords_polynomials(plan, workload.bg.z)
    T_k = workload.plans.T_k_limber

    pk_lin1 = 1.02 .* workload.pk_limber_lin .+ 0.001
    pk_nonlin1 = 0.97 .* workload.pk_limber_nonlin .+ 0.002
    args = (
        Reactant.to_rarray(workload.pk_limber_lin),
        Reactant.to_rarray(workload.pk_limber_nonlin),
        Reactant.to_rarray(transform_k), Reactant.to_rarray(transform_z),
        Reactant.to_rarray(T_z), Reactant.to_rarray(T_k),
    )
    args1 = (
        Reactant.to_rarray(pk_lin1), Reactant.to_rarray(pk_nonlin1), args[3:end]...
    )
    prepare_fun = let bg=workload.bg
        (pk_lin, pk_nonlin, Mk, Mz, Tz, Tk) -> Blast.reactant_prepare_limber(
            pk_lin, pk_nonlin, bg, Mk, Mz, Tz, Tk,
        )
    end
    compiled = Reactant.compile(prepare_fun, args; sync=true)
    result = map(Array, compiled(args...))
    result1 = map(Array, compiled(args1...))
    ordinary1 = Blast.prepare_pk_workspace(
        workload.plans, workload.pk, pk_lin1, pk_nonlin1, workload.bg,
    )
    println("Reactant=", pkgversion(Reactant))
    println("plan_K=", plan.K, " T_k_shape=", size(T_k))
    println("delta_parity_error=", maximum(abs, result[1] .- reference.ΔP_limber))
    println("Pdelta_parity_error=", maximum(abs, result[2] .- reference.Pδ_limber))
    println("delta_new_input_error=", maximum(abs, result1[1] .- ordinary1.ΔP_limber))
    println("Pdelta_new_input_error=", maximum(abs, result1[2] .- ordinary1.Pδ_limber))
    println("delta_compiled_delta=", maximum(abs, result1[1] .- result[1]))
    println("delta_reference_delta=", maximum(abs, ordinary1.ΔP_limber .- reference.ΔP_limber))
    println("Pdelta_compiled_delta=", maximum(abs, result1[2] .- result[2]))
    println("Pdelta_reference_delta=", maximum(abs, ordinary1.Pδ_limber .- reference.Pδ_limber))
    host_call = () -> Blast.prepare_pk_workspace(
        workload.plans, workload.pk, workload.pk_limber_lin,
        workload.pk_limber_nonlin, workload.bg,
    )
    reactant_call = () -> compiled(args...)
    println("host_ms=", minimum(@benchmark $host_call() samples=3 evals=1).time / 1e6)
    println("reactant_ms=", minimum(@benchmark $reactant_call() samples=5 evals=1).time / 1e6)
end

main()

using Blast
using Reactant
using NPZ

const ROOT = "/home/marcobonici/Desktop/work/CosmologicalEmulators/blast_play/for_marco_copy"
include(joinpath(ROOT, "performance", "workload.jl"))

function rss_kib()
    for line in eachline("/proc/self/status")
        startswith(line, "VmRSS:") && return parse(Int, split(line)[2])
    end
    error("VmRSS missing")
end

function invoke_discard(plan, args)
    result = plan(args...)
    foreach(Reactant.synchronize, values(result))
    return nothing
end

function main()
    Reactant.set_default_backend("cpu")
    workload = build_workload(ROOT, 1)
    plan = Blast.build_reactant_full_plan(workload)
    args = map(Reactant.to_rarray,
        (workload.pk, workload.pk_limber_lin, workload.pk_limber_nonlin))
    invoke_discard(plan, args)
    GC.gc(true)
    println("rss_after_build_gc_kib=", rss_kib())
    for i in 1:20
        invoke_discard(plan, args)
        i % 5 == 0 && println("rss_no_gc_iteration_$(i)_kib=", rss_kib())
    end
    GC.gc(true)
    println("rss_after_no_gc_sweep_kib=", rss_kib())
    for i in 1:20
        invoke_discard(plan, args)
        GC.gc(true)
        i % 5 == 0 && println("rss_with_gc_iteration_$(i)_kib=", rss_kib())
    end
end

main()

using Blast
using DelimitedFiles
using NPZ

get(ENV, "BLAST_REGENERATE_FULL_FIXTURES", "0") == "1" ||
    error("Set BLAST_REGENERATE_FULL_FIXTURES=1 explicitly")

const ROOT = "/home/marcobonici/Desktop/work/CosmologicalEmulators/blast_play/for_marco_copy"
const OUT = joinpath(@__DIR__, "reference_fixtures", "reactant_full")
include(joinpath(ROOT, "performance", "workload.jl"))

function main()
    workload = build_workload(ROOT, 1)
    pk0, lin0, non0 = workload.pk, workload.pk_limber_lin, workload.pk_limber_nonlin
    pk1 = pk0 .* reshape(1 .+ 1e-5 .* sin.(1:size(pk0, 1)), :, 1)
    lin1 = lin0 .* reshape(1 .+ 1e-5 .* sin.(1:size(lin0, 1)), :, 1)
    non1 = non0 .* reshape(1 .+ 1e-5 .* cos.(1:size(non0, 1)), :, 1)
    cases = (
        ("baseline", pk0, lin0, non0),
        ("nonlimber", pk1, lin0, non0),
        ("limber_linear", pk0, lin1, non0),
        ("limber_nonlinear", pk0, lin0, non1),
    )
    mkpath(OUT)
    metadata = String[]
    for (label, pk, lin, non) in cases
        result = full_power_to_cls(merge(workload,
            (; pk, pk_limber_lin=lin, pk_limber_nonlin=non)); pass_plans=true)
        for name in (:cl_gg, :cl_gs, :cl_ss)
            value = getproperty(result, name)
            writedlm(joinpath(OUT, "$(label)_$(name).txt"), vec(value))
            push!(metadata, "$(label)_$(name) size=$(join(size(value), ','))")
        end
    end
    open(joinpath(OUT, "metadata.txt"), "w") do io
        println(io, "generator=ordinary full_power_to_cls")
        println(io, "perturbation_amplitude=1e-5")
        foreach(line -> println(io, line), metadata)
    end
    println("wrote fixtures to $OUT")
end

main()

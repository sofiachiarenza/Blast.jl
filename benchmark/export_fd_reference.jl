# Create a NEW FD-reference fixture from a known-passed objective-check log.
# Does not re-run FD.   julia --project=benchmark benchmark/export_fd_reference.jl \
#     LOG OUTPUT.toml   (env: BLAST_REALISTIC_DATA_ROOT, BLAST_REALISTIC_WITH_PNG=1,
#     BLAST_FD_ORIGIN="baseline=e2b47d6 ace=0.11.0 ...")
using Blast, NPZ, ForwardDiff
include(joinpath(@__DIR__, "..", "test", "characterization_workload.jl"))
include(joinpath(@__DIR__, "realistic_workload.jl"))
include(joinpath(@__DIR__, "fd_reference_fixture.jl"))
length(ARGS) == 2 || error("usage: export_fd_reference.jl LOG OUTPUT.toml")
root = ENV["BLAST_REALISTIC_DATA_ROOT"]
workload = build_realistic_workload(root, 1)   # context only (PNG flag, active components)
context = fd_fixture_context(root, [1.0, 1.0, 1.0], [1.07, 0.93, 1.11];
                             with_png=workload.with_png, active_components=length(workload.active_fields))
origin = Dict{String,Any}("note" => get(ENV, "BLAST_FD_ORIGIN", ""), "log" => abspath(ARGS[1]))
write_fd_reference(ARGS[2], ARGS[1], context; origin)
println("FD_REFERENCE_WRITTEN ", ARGS[2], " sha256=", file_sha256(ARGS[2]))

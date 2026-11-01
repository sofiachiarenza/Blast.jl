# Saved independent ForwardDiff reference for the realistic 3-Pmm-amplitude
# objective, produced once by the baseline code (zero-partial lifted background)
# and used as a regression oracle for later runs. Fail closed: any stale,
# missing, non-finite or mismatching entry is an error; there is no fallback and
# no regeneration here. Select with BLAST_BENCH_FD_REFERENCE=<file.toml>;
# unset means the FD reference is computed as before.
using TOML, SHA

const FD_FIXTURE_FORMAT = "blast-fd-reference-v1"
const FD_DATASET_FILES = ("lsst_nz.npz", "ell_n5k.npz", "pk_nonlimber_161.npz",
                          "pk_limber_lin.npz", "pk_limber_nonlin.npz")

file_sha256(path) = open(io -> bytes2hex(sha256(io)), path)
floats_sha256(v) = bytes2hex(sha256(reinterpret(UInt8, htol.(collect(Float64, v)))))

"""
Context a reference is only valid for: inputs x0/x1, PNG/17-component workload, R
count, SHA-256 of the five data files, exact axis fingerprints, and SHA-256 of
the objective/workload harness code.
"""
function fd_fixture_context(root, x0, x1; with_png, active_components, bench_dir=@__DIR__)
    return Dict{String,Any}(
        "x0" => collect(Float64, x0), "x1" => collect(Float64, x1),
        "with_png" => with_png, "active_components" => active_components,
        "n_R" => length(Blast.R),
        "datasets" => Dict{String,Any}(f => file_sha256(joinpath(root, "data", f)) for f in FD_DATASET_FILES),
        "axes" => Dict{String,Any}("ell_nonlimber" => floats_sha256(Blast.ℓ_nonlimber), "chi" => floats_sha256(Blast.χ),
                                   "R" => floats_sha256(Blast.R), "k_cheb" => floats_sha256(Blast.k_cheb)),
        "harness" => Dict{String,Any}("realistic_workload.jl" => file_sha256(joinpath(bench_dir, "realistic_workload.jl")),
                                      "characterization_workload.jl" => file_sha256(joinpath(bench_dir, "..", "test", "characterization_workload.jl"))))
end

_parse_vec(s) = parse.(Float64, strip.(split(s, ",")))

"""
Parse the two accepted records from a known-passed objective-check log, verifying
EXIT=0 and exactly one record each for x0 and x1; values are taken verbatim.
"""
function parse_fd_log(text)
    occursin(r"^EXIT=0$"m, text) || error("log has no EXIT=0")
    occursin("workload built (with_png=true, active=17)", text) || error("log does not show the 17-component PNG workload")
    records = Dict{String,Any}()
    for m in eachmatch(r"FD point (x0|x1): primal=(\S+) reference=(\S+) gradient=\[([^\]]*)\]", text)
        haskey(records, m[1]) && error("duplicate record for $(m[1])")
        records[m[1]] = Dict{String,Any}("primal_fd" => parse(Float64, m[2]), "float_reference" => parse(Float64, m[3]),
                                         "gradient" => _parse_vec(m[4]))
    end
    Set(keys(records)) == Set(["x0", "x1"]) || error("expected exactly the x0 and x1 records; got $(sort(collect(keys(records))))")
    return records
end

"Write a NEW fixture atomically from a passed log; refuses to overwrite."
function write_fd_reference(path, logpath, context; origin)
    isfile(path) && error("refusing to overwrite $path")
    text = read(logpath, String)
    records = parse_fd_log(text)
    for (k, r) in records
        all(isfinite, r["gradient"]) && isfinite(r["primal_fd"]) && isfinite(r["float_reference"]) || error("non-finite record $k")
        length(r["gradient"]) == 3 || error("gradient of $k is not length 3")
    end
    fixture = Dict{String,Any}("format" => FD_FIXTURE_FORMAT,
                               "origin" => merge(Dict{String,Any}(origin), Dict{String,Any}("log_sha256" => file_sha256(logpath))),
                               "context" => context, "point" => records)
    tmp = path * ".partial." * string(getpid())
    open(io -> TOML.print(io, fixture), tmp, "w")
    mv(tmp, path)
    return path
end

function _finite_vec(v, n, what)
    v isa AbstractVector && length(v) == n && all(x -> x isa Float64 && isfinite(x), v) || error("fixture $what is not $n finite Float64 values")
    return Vector{Float64}(v)
end

"Load and validate; returns (; x0, x1, points, origin, sha256)."
function load_fd_reference(path, context)
    isfile(path) || error("FD reference fixture missing: $path")
    fixture = TOML.parsefile(path)
    get(fixture, "format", nothing) == FD_FIXTURE_FORMAT || error("fixture format is not $FD_FIXTURE_FORMAT")
    for key in ("origin", "context", "point")
        haskey(fixture, key) && fixture[key] isa AbstractDict || error("fixture lacks section $key")
    end
    stored = fixture["context"]
    mismatches = String[]
    for key in sort(collect(union(keys(stored), keys(context))))
        if !haskey(stored, key) || !haskey(context, key)
            push!(mismatches, key)
        elseif stored[key] isa AbstractDict
            for sub in sort(collect(union(keys(stored[key]), keys(context[key]))))
                get(stored[key], sub, nothing) == get(context[key], sub, missing) || push!(mismatches, "$key.$sub")
            end
        elseif stored[key] != context[key]
            push!(mismatches, key)
        end
    end
    isempty(mismatches) || error("stale/mismatching FD reference fixture; differing context: $(join(mismatches, ", "))")
    points = Dict{String,Any}()
    for k in ("x0", "x1")
        haskey(fixture["point"], k) || error("fixture lacks point $k")
        r = fixture["point"][k]
        points[k] = (; primal_fd=r["primal_fd"], float_reference=r["float_reference"], gradient=_finite_vec(r["gradient"], 3, "$k gradient"))
        all(x -> x isa Float64 && isfinite(x), (points[k].primal_fd, points[k].float_reference)) || error("fixture $k primal is not finite")
    end
    return (; x0=_finite_vec(stored["x0"], 3, "x0"), x1=_finite_vec(stored["x1"], 3, "x1"), points,
            origin=Dict{String,Any}(fixture["origin"]), sha256=file_sha256(path))
end

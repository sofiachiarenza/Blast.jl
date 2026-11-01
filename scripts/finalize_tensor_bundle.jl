# Finalize a portable B-spline tensor bundle manifest for already-merged campaign banks.
#
#   julia --project=. scripts/finalize_tensor_bundle.jl --campaign DIR --bankroot DIR --counts 33[,49,...] \
#         [--manifest-name bundle_N033.toml] [--runtime-ace-rev REV] [--generator-ace-version ID]
#   julia --project=. scripts/finalize_tensor_bundle.jl --help
#
# Reads the campaign's real artifacts (inputs/manifest.json, inputs/*.txt, bank_mapping.json and
# the per-count merge manifests `<bankroot>/N<count>/manifest.json` written by merge.py, whose
# paths are relative to the bank ROOT). It verifies every selected payload (SHA-256, byte count,
# NPY header) and the input/generator identities, adds ONLY small metadata NPY files under
# `<bankroot>/meta/`, validates the staged manifest with the public loader (both representations
# of every selected count, one at a time), and publishes the manifest LAST, atomically. Tensor
# payloads are never copied, moved or rewritten; existing metadata is reused only on an exact
# hash match, otherwise the run fails. Status is always "candidate"; coverage is exactly the
# listed counts (no claim about other counts, no science certification).
module BundleFinalizer

using Blast, JSON, SHA, TOML
using Blast.AbstractCosmologicalEmulators: CubicSpline, CubicBSplinePlan, bspline_basis, knot_vector

fin_error(msg) = throw(ArgumentError("finalizer: " * msg))

sha_file(path) = bytes2hex(open(sha256, path))
parse_floats(path) = [parse(Float64, strip(l)) for l in eachline(path) if !isempty(strip(l))]
function parse_table(path)
    rows = [parse.(Float64, split(strip(l))) for l in eachline(path) if !isempty(strip(l))]
    all(length(r) == 2 for r in rows) || fin_error("$path: expected two columns")
    return permutedims(reduce(hcat, rows))
end

# Minimal NPY v1 writer: little-endian Float64, Fortran order (what the loader requires).
function npy_write(path, A::AbstractArray{Float64})
    shape = join(size(A), ", ") * (ndims(A) == 1 ? "," : "")
    header = "{'descr': '<f8', 'fortran_order': True, 'shape': ($shape), }"
    header *= " "^(mod(-(10 + length(header) + 1), 64)) * "\n"
    open(path, "w") do io
        write(io, b"\x93NUMPY", UInt8(1), UInt8(0), htol(UInt16(length(header))), header)
        write(io, htol.(reinterpret(UInt64, vec(A))))       # explicit little-endian, bit-exact
    end
    return path
end
npy_header_bytes(path) = open(io -> (skip(io, 8); 10 + Int(ltoh(read(io, UInt16)))), path)

function file_ref(root, rel)
    path = joinpath(root, rel)
    return Dict{String,Any}("file" => rel, "sha256" => sha_file(path), "bytes" => filesize(path),
                            "header_bytes" => npy_header_bytes(path))
end

# Independent restatement of the generator's normalizer (generator.jl: natural ln/ln CubicSpline,
# As = exp(ln10_10_As)*1e-10), used only for the saved witness. Not the loader's helper.
function witness_q(table, sites, ln10_10_As, n_s)
    spline = CubicSpline(log.(table[:, 2]), log.(table[:, 1]))
    pmm = exp.(spline(log.(sites)))
    As = exp(ln10_10_As) * 1e-10
    phi = @. (9 / 25) * 2 * pi^2 * As / sites^3 * (sites / 0.05)^(n_s - 1)
    return hcat(pmm, sqrt.(phi .* pmm), phi)
end

# A directory below `root` that exists as a real directory or is created here; every component is
# checked (no `..`, no empty/absolute/malformed parts, no symbolic links), never resolved away.
function safe_dir!(root, rel)
    Blast._check_relpath(rel)
    path = root
    for segment in split(rel, '/')
        path = joinpath(path, segment)
        islink(path) && fin_error("`$rel`: component `$segment` is a symbolic link")
        if ispath(path)
            isdir(path) || fin_error("`$rel`: component `$segment` is not a directory")
        else
            mkdir(path)
        end
    end
    return path
end

# Read a campaign file by relative name with the loader's path policy (no escapes, no symlinks).
read_path(root_real, rel) = Blast._resolve(root_real, Blast._check_relpath(rel))

# Remove only files this run created, only inside the declared root, never through a link.
function cleanup_owned!(owned, root)
    for p in owned
        (islink(p) || ispath(p)) || continue
        isdir(dirname(p)) && startswith(realpath(dirname(p)) * "/", root * "/") && rm(p; force = true)
    end
end

"Structural checks of one merge manifest (merge.py layout); returns (identity, listed banks)."
function check_merge_manifest(mm, n, shape_tail, inputs_sha)
    mm["count"] == n || fin_error("count $n: merge manifest is for count $(mm["count"])")
    mm["shape"] == [shape_tail..., n] || fin_error("count $n: merge manifest shape $(mm["shape"])")
    lowercase(String(mm["dtype"])) == "float64" || fin_error("count $n: dtype $(mm["dtype"])")
    mm["axis_order"] == ["ell", "chi", "R", "bspline_coefficient"] || fin_error("count $n: unexpected axis order")
    mm["inputs_manifest_sha256"] == inputs_sha || fin_error("count $n: merged from different campaign inputs")
    listed = [(String(b["family"]), String(b["geometry"])) for b in mm["banks"]]
    (length(listed) == length(ALL_BANKS) && Set(listed) == Set(ALL_BANKS)) ||
        fin_error("count $n: merge manifest must list exactly the 21 bank families")
    return (String(mm["generator_sha256"]), String(mm["inputs_manifest_sha256"]), Int(mm["quadrature_N"])), listed
end

const ALL_BANKS = (Blast._DIRECT_BANKS..., Blast._RATIO_BANKS...)

"Component -> source family names of the campaign's bank_mapping.json, checked against the loader's table."
function source_map_rows(mapping_path)
    mapping = JSON.parsefile(mapping_path)
    rows = Dict{String,Any}[]
    for (name, tfield, source, r1) in Blast._PROJECTED_MATTER_COMPONENTS
        component = Blast._ascii_component(name)
        haskey(mapping, component) || fin_error("bank_mapping.json lacks `$component`")
        entry = mapping[component]
        src = Blast._SOURCE_NAME[source]
        campaign_source = src == "T" ? (r1 ? "T1" : "TR") : src
        entry["source"] == campaign_source || fin_error("bank_mapping `$component`: source `$(entry["source"])` != `$campaign_source`")
        entry["direct"]["family"] == "direct" && entry["direct"]["geometry"] == String(tfield) || fin_error("bank_mapping `$component`: wrong direct bank")
        entry["ratio"]["family"] == Blast._RATIO_FAMILY[src] && entry["ratio"]["geometry"] == String(tfield) || fin_error("bank_mapping `$component`: wrong ratio bank")
        push!(rows, Dict{String,Any}("component" => component, "source" => src, "orientation" => r1 ? "R1" : src,
                                     "family_direct" => ["direct", String(tfield)], "family_ratio" => [Blast._RATIO_FAMILY[src], String(tfield)]))
    end
    length(mapping) == length(rows) || fin_error("bank_mapping.json has components other than the 17 projected ones")
    return rows
end

function runtime_ace_rev_from_manifest()
    path = joinpath(dirname(Base.active_project()), "Manifest.toml")
    isfile(path) || return nothing
    for entry in get(get(TOML.parsefile(path), "deps", Dict()), "AbstractCosmologicalEmulators", [])
        haskey(entry, "repo-rev") && return String(entry["repo-rev"])
    end
    return nothing
end

"""
    finalize(; campaign, bankroot, counts, manifest_name = "bundle.toml", status = "candidate",
             axes = Blast.production_axes(), expected_k_range = Blast.PRODUCTION_K_RANGE,
             runtime_ace_rev = nothing, generator_ace_version = nothing, log = println)

Returns the manifest path. Throws `ArgumentError` (publishing nothing, removing its own staging
files) on any inconsistency.
"""
function finalize(; campaign, bankroot, counts, manifest_name = "bundle.toml", status = "candidate",
                  axes = Blast.production_axes(), expected_k_range = Blast.PRODUCTION_K_RANGE,
                  runtime_ace_rev = nothing, generator_ace_version = nothing, log = println)
    status == "candidate" || fin_error("only status \"candidate\" can be produced here")
    (!isempty(manifest_name) && basename(manifest_name) == manifest_name && endswith(manifest_name, ".toml")) ||
        fin_error("--manifest-name must be a plain *.toml file name")
    counts = sort(unique(collect(Int, counts)))
    isempty(counts) && fin_error("no counts requested")
    isdir(bankroot) || fin_error("bank root not found")
    root = realpath(bankroot)
    campaign_real = realpath(campaign)
    isdir(joinpath(campaign_real, "inputs")) || fin_error("campaign inputs/ not found")
    inputs = realpath(joinpath(campaign_real, "inputs"))
    owned = String[]
    tag = string(getpid(), ".", bytes2hex(rand(UInt8, 4)))
    try
        # ---- campaign identity: every listed input file, and the manifest's own hash ----
        im = JSON.parsefile(read_path(inputs, "manifest.json"))
        declared = im["input_sha256"]
        for (name, sha) in declared                   # names are validated BEFORE any read or hash
            path = read_path(inputs, String(name))
            sha_file(path) == sha || fin_error("campaign input `$name` fails its checksum")
        end
        # Every input this run reads must be declared (and was therefore hash-verified above).
        needed = vcat(["ell.txt", "chi.txt", "R.txt", "fiducial_pk.txt"],
                      ["$(k)_$(r)_$(n).txt" for k in ("sites", "knots"), r in ("direct", "ratio"), n in counts][:])
        for name in needed
            haskey(declared, name) || fin_error("campaign input `$name` is not declared in inputs/manifest.json input_sha256")
        end
        inputs_sha = sha_file(read_path(inputs, "manifest.json"))
        kmin, kmax = Float64(im["kmin"]), Float64(im["kmax"])
        (abs(kmin / expected_k_range[1] - 1) ≤ 4eps() && abs(kmax / expected_k_range[2] - 1) ≤ 4eps()) ||
            fin_error("inputs k range ($kmin, $kmax) is not the trusted range $expected_k_range")
        params = im["fiducial_parameters"]
        ln10_10_As = Float64(params["ln10^{10}A_s"]); n_s = Float64(params["n_s"])
        ell, chi, R = (parse_floats(read_path(inputs, n * ".txt")) for n in ("ell", "chi", "R"))
        (ell == axes.ell && chi == axes.chi && R == axes.R) || fin_error("campaign axes differ from the expected axes")
        table = parse_table(read_path(inputs, "fiducial_pk.txt"))
        source_rows = source_map_rows(read_path(campaign_real, "bank_mapping.json"))
        gen_ace = generator_ace_version
        if isnothing(gen_ace)
            vp = joinpath(campaign_real, "vendor_provenance.json")
            (islink(vp) || !isfile(vp)) && fin_error("pass --generator-ace-version (no regular vendor_provenance.json)")
            gen_ace = "vendor_head:" * String(JSON.parsefile(vp)["ACE_local_head"])
        end
        rt_rev = isnothing(runtime_ace_rev) ? runtime_ace_rev_from_manifest() : runtime_ace_rev
        isnothing(rt_rev) && fin_error("pass --runtime-ace-rev (not found in the active Manifest)")

        # ---- verify every selected count's merged banks BEFORE touching anything ----
        shape_tail = (length(axes.ell), length(axes.chi), length(axes.R))
        identities = Set{Tuple{String,String,Int}}()
        verified = Dict{Int,Vector{Dict{String,Any}}}()
        for n in counts
            dir = joinpath(root, "N$(lpad(n, 3, '0'))")
            (isfile(joinpath(dir, "manifest.json")) && !islink(dir)) || fin_error("count $n has no completed merge manifest")
            mm = JSON.parsefile(Blast._resolve(root, "N$(lpad(n, 3, '0'))/manifest.json"))
            identity, _ = check_merge_manifest(mm, n, shape_tail, inputs_sha)
            push!(identities, identity)
            entries = Dict{String,Any}[]
            for b in mm["banks"]
                family, geometry = String(b["family"]), String(b["geometry"])
                rel = "N$(lpad(n, 3, '0'))/$family/$geometry.npy"
                b["path"] == rel || fin_error("count $n: bank path `$(b["path"])` is not `$rel` (relative to the bank root)")
                path = Blast._resolve(root, rel)
                filesize(path) == b["bytes"] || fin_error("count $n: `$rel` has $(filesize(path)) bytes, manifest $(b["bytes"])")
                sha_file(path) == b["sha256"] || fin_error("count $n: `$rel` fails its SHA-256")
                header = Blast._npy_header(path)
                (header.descr == "<f8" && header.fortran && header.shape == (shape_tail..., n)) ||
                    fin_error("count $n: `$rel` is not Fortran-order Float64 $((shape_tail..., n))")
                filesize(path) == header.header_bytes + 8 * prod(header.shape) || fin_error("count $n: `$rel` is truncated")
                push!(entries, Dict{String,Any}("family" => family, "geometry" => geometry,
                      "file_ref" => Dict{String,Any}("file" => rel, "sha256" => b["sha256"], "bytes" => Int(b["bytes"]),
                                                     "header_bytes" => header.header_bytes)))
            end
            verified[n] = entries
        end
        length(identities) == 1 || fin_error("selected counts disagree on generator/inputs/quadrature identity")
        generator_sha, _, quadrature_N = only(identities)
        log("verified merged banks for counts $counts (21 payloads each)")

        # ---- small metadata under <bankroot>/meta (exact-hash reuse only) ----
        meta_dir = safe_dir!(root, "meta")                # every parent checked before any write
        function publish_meta!(name, A)
            rel = "meta/$name.npy"; target = joinpath(meta_dir, name * ".npy")
            islink(target) && fin_error("`$rel` is a symbolic link")
            tmp = target * ".partial." * tag
            (islink(tmp) || ispath(tmp)) && fin_error("staging path `$tmp` already exists")
            push!(owned, tmp)
            npy_write(tmp, A)
            if ispath(target)
                sha_file(target) == sha_file(tmp) || fin_error("existing `$rel` differs from the expected content; refusing to overwrite")
                rm(tmp)
            else
                mv(tmp, target)
            end
            return file_ref(root, rel)
        end
        axes_refs = Dict{String,Any}(n => publish_meta!("axes_$n", Float64.(v)) for (n, v) in (("ell", ell), ("chi", chi), ("R", R)))
        fiducial_ref = publish_meta!("fiducial_pmm", table)
        evidence = String[]
        count_entries = Dict{String,Any}[]
        for n in counts, rep in ("direct", "ratio")
            sites = parse_floats(read_path(inputs, "sites_$(rep)_$n.txt"))
            knots = parse_floats(read_path(inputs, "knots_$(rep)_$n.txt"))
            logsites = log10.(sites)
            length(sites) == n && length(knots) == n + 4 || fin_error("count $n/$rep: wrong sites/knots length")
            knot_vector(bspline_basis(CubicBSplinePlan(logsites, logsites; extrapolation = :throw))) == knots ||
                fin_error("count $n/$rep: saved knots differ from the runtime ACE not-a-knot basis")
            push!(evidence, "$rep:$n")
            tag_n = "N$(lpad(n, 3, '0'))_$rep"
            entry = Dict{String,Any}("count" => n, "representation" => rep,
                "sites_phys" => publish_meta!("sites_phys_$tag_n", sites),
                "sites_log10" => publish_meta!("sites_log10_$tag_n", logsites),
                "knots" => publish_meta!("knots_$tag_n", knots),
                "bank" => filter(e -> (rep == "direct") == (e["family"] == "direct"), verified[n]))
            rep == "ratio" && (entry["witness"] = publish_meta!("witness_$tag_n", witness_q(table, sites, ln10_10_As, n_s)))
            push!(count_entries, entry)
        end

        manifest = Dict{String,Any}(
            "schema" => Blast.BUNDLE_SCHEMA, "status" => status,
            "convention" => Dict{String,Any}(
                "coordinate" => Blast.BUNDLE_CONVENTION.coordinate, "basis" => Blast.BUNDLE_CONVENTION.basis,
                "coefficient" => Blast.BUNDLE_CONVENTION.coefficient, "layout" => Blast.BUNDLE_CONVENTION.layout,
                "endpoint" => Blast.BUNDLE_CONVENTION.endpoint, "endpoint_epsilon" => Blast.BUNDLE_ENDPOINT_EPSILON,
                "direct_source" => Blast.BUNDLE_CONVENTION.direct_source, "ratio_source" => Blast.BUNDLE_CONVENTION.ratio_source),
            "provenance" => Dict{String,Any}("generator_sha256" => generator_sha, "inputs_manifest_sha256" => inputs_sha,
                "quadrature_N" => quadrature_N, "generator_ace_version" => gen_ace, "runtime_ace_rev" => String(rt_rev),
                "basis_compat_evidence" => "saved knots equal the runtime ACE not-a-knot basis exactly for designs " * join(evidence, ",")),
            "range" => Dict{String,Any}("kmin" => kmin, "kmax" => kmax),
            "axes" => axes_refs,
            "fiducial" => Dict{String,Any}("table" => fiducial_ref, "rows" => size(table, 1), "ln10_10_As" => ln10_10_As,
                "n_s" => n_s, "pivot" => Blast.BUNDLE_PIVOT,
                "pmm_interpolator" => Blast.BUNDLE_FIDUCIAL.pmm_interpolator, "phi_formula" => Blast.BUNDLE_FIDUCIAL.phi_formula,
                "q_TT" => Blast.BUNDLE_FIDUCIAL.q_TT, "q_T" => Blast.BUNDLE_FIDUCIAL.q_T, "q_PHI" => Blast.BUNDLE_FIDUCIAL.q_PHI),
            "count" => count_entries, "source_map" => source_rows)   # no timestamp: byte-deterministic

        staged = joinpath(root, manifest_name * ".partial." * tag)
        (islink(staged) || ispath(staged)) && fin_error("staging path `$staged` already exists")
        islink(joinpath(root, manifest_name)) && fin_error("`$manifest_name` is a symbolic link")
        push!(owned, staged)
        open(io -> TOML.print(io, manifest; sorted = true), staged, "w")
        for n in counts, rep in (:direct, :ratio)       # one selection at a time; released before the next
            Blast.load_tensor_bank(staged; count = n, representation = rep, axes, expected_k_range)
            GC.gc()
        end
        log("staged manifest validated by the public loader (direct and ratio, counts $counts)")

        target = joinpath(root, manifest_name)
        if ispath(target)
            read(target) == read(staged) || fin_error("`$manifest_name` already exists with different content; refusing to overwrite")
            rm(staged)
            log("already published, byte-identical: $target")
        else
            mv(staged, target)
            log("published $target")
        end
        log("COVERAGE counts=$counts status=$status sha256=$(sha_file(target))")
        return target
    catch
        cleanup_owned!(owned, root)
        rethrow()
    end
end

const USAGE = """
usage: finalize_tensor_bundle.jl --campaign DIR --bankroot DIR --counts N[,N...] [options]
  --campaign DIR            campaign root (inputs/, bank_mapping.json, vendor_provenance.json)
  --bankroot DIR            merged banks root (N033/, ... with merge manifests); metadata is added under meta/
  --counts 33,49            counts to include (each must be completely merged)
  --manifest-name NAME.toml manifest file written LAST in the bank root (default bundle.toml)
  --runtime-ace-rev REV     runtime ACE revision recorded (default: repo-rev in the active Manifest.toml)
  --generator-ace-version V generator ACE identity (default: vendor_head:<vendor_provenance.json ACE_local_head>)
  --help                    show this message
Status is always "candidate"; existing metadata/manifests are never overwritten.
"""

function main(args)
    if any(a -> a in ("-h", "--help"), args)
        print(USAGE); return nothing
    end
    opts = Dict{String,String}()
    i = 1
    while i ≤ length(args)
        key = args[i]
        startswith(key, "--") && key[3:end] in ("campaign", "bankroot", "counts", "manifest-name", "runtime-ace-rev", "generator-ace-version") ||
            fin_error("unknown option `$key`\n" * USAGE)
        i < length(args) || fin_error("option `$key` needs a value")
        opts[key[3:end]] = args[i + 1]; i += 2
    end
    for required in ("campaign", "bankroot", "counts")
        haskey(opts, required) || fin_error("missing --$required\n" * USAGE)
    end
    counts = try
        parse.(Int, split(opts["counts"], ','))
    catch
        fin_error("--counts must be comma-separated integers")
    end
    return finalize(; campaign = opts["campaign"], bankroot = opts["bankroot"], counts,
                    manifest_name = get(opts, "manifest-name", "bundle.toml"),
                    runtime_ace_rev = get(opts, "runtime-ace-rev", nothing),
                    generator_ace_version = get(opts, "generator-ace-version", nothing))
end

end # module BundleFinalizer

if abspath(PROGRAM_FILE) == @__FILE__
    BundleFinalizer.main(ARGS)
end

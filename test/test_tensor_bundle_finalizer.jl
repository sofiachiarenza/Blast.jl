# Phase 4B: finalizer tests on a fake campaign built from synthetic metadata but the REAL
# JSON layouts (inputs/manifest.json, bank_mapping.json, merge.py per-count manifests with
# paths relative to the bank root). No loader monkeypatching; tiny arrays, temporary files.
using Test, Blast, JSON, SHA, TOML, LinearAlgebra
using AbstractCosmologicalEmulators, FastChebInterp      # same imports as the full suite
isdefined(@__MODULE__, :write_synthetic_bundle) || include("tensor_bundle_fixture.jl")
isdefined(@__MODULE__, :BundleFinalizer) || include(joinpath(@__DIR__, "..", "scripts", "finalize_tensor_bundle.jl"))

const FIN_ARGS = (axes = SYN_AXES, expected_k_range = (SYN_KMIN, SYN_KMAX), runtime_ace_rev = "synthetic-rev")

function write_json(path, obj); mkpath(dirname(path)); open(io -> JSON.print(io, obj, 2), path, "w"); end
write_lines(path, v) = (mkpath(dirname(path)); open(io -> foreach(x -> println(io, repr(x)), v), path, "w"))

"Fake campaign + merged banks for `counts`, in the real campaign file layouts."
function write_fake_campaign(dir; counts=(9,), generator="generator-A")
    campaign = joinpath(dir, "campaign"); inputs = joinpath(campaign, "inputs"); banks = joinpath(dir, "banks")
    mkpath(inputs); mkpath(banks)
    for (name, v) in pairs(SYN_AXES); write_lines(joinpath(inputs, "$name.txt"), v); end
    table = syn_table()
    open(io -> foreach(r -> println(io, repr(r[1]), "\t", repr(r[2])), eachrow(table)), joinpath(inputs, "fiducial_pk.txt"), "w")
    for n in (9, 11), kind in ("direct", "ratio")
        sites = syn_sites(n, kind); x = log10.(sites)
        write_lines(joinpath(inputs, "sites_$(kind)_$n.txt"), sites)
        write_lines(joinpath(inputs, "knots_$(kind)_$n.txt"),
                    knot_vector(bspline_basis(CubicBSplinePlan(x, x; extrapolation=:throw))))
    end
    input_hashes = Dict(name => syn_sha(joinpath(inputs, name)) for name in readdir(inputs))
    write_json(joinpath(inputs, "manifest.json"), Dict("basis" => "synthetic", "kmin" => SYN_KMIN, "kmax" => SYN_KMAX,
        "counts" => collect(counts), "dimensions" => [3, 4, 5], "dtype" => "Float64",
        "fiducial_parameters" => Dict("ln10^{10}A_s" => SYN_LN10AS, "n_s" => SYN_NS, "h" => 0.67),
        "input_sha256" => input_hashes))
    mapping = Dict(Blast._ascii_component(name) => Dict("direct" => Dict("family" => "direct", "geometry" => string(t)),
        "ratio" => Dict("family" => Blast._RATIO_FAMILY[Blast._SOURCE_NAME[src]], "geometry" => string(t)),
        "source" => (Blast._SOURCE_NAME[src] == "T" ? (r1 ? "T1" : "TR") : Blast._SOURCE_NAME[src]),
        "coefficient_source_orientation" => "x") for (name, t, src, r1) in Blast._PROJECTED_MATTER_COMPONENTS)
    write_json(joinpath(campaign, "bank_mapping.json"), mapping)
    write_json(joinpath(campaign, "vendor_provenance.json"), Dict("ACE_local_head" => "abc123"))
    inputs_sha = syn_sha(joinpath(inputs, "manifest.json"))
    for n in counts
        entries = Dict{String,Any}[]
        for (family, geometry) in vcat(collect(SYN_DIRECT), collect(SYN_RATIO))
            rel = "N$(lpad(n, 3, '0'))/$family/$geometry.npy"
            syn_write_npy(joinpath(banks, rel), syn_tensor(family, geometry, n))
            push!(entries, Dict("family" => family, "geometry" => geometry, "path" => rel,
                                "sha256" => syn_sha(joinpath(banks, rel)), "bytes" => filesize(joinpath(banks, rel))))
        end
        write_json(joinpath(banks, "N$(lpad(n, 3, '0'))", "manifest.json"), Dict(
            "count" => n, "shape" => [3, 4, 5, n], "dtype" => "float64",
            "axis_order" => ["ell", "chi", "R", "bspline_coefficient"], "generator_sha256" => syn_sha_text(generator),
            "inputs_manifest_sha256" => inputs_sha, "quadrature_N" => 131073, "banks" => entries))
    end
    return (; campaign, banks, inputs)
end

fin(c; counts=(9,), name="bundle_synthetic.toml", kw...) =
    BundleFinalizer.finalize(; campaign=c.campaign, bankroot=c.banks, counts, manifest_name=name, log=_ -> nothing, FIN_ARGS..., kw...)
frejects(c; why=r"", kw...) = try; fin(c; kw...); false; catch e
    e isa ArgumentError && occursin(why, e.msg) || (@warn "unexpected failure" e; false)
end
leftovers(c) = [f for (r, _, fs) in walkdir(c.banks) for f in fs if occursin(".partial.", f)]
snapshot(c) = Dict(joinpath(r, f) => (syn_sha(joinpath(r, f)), mtime(joinpath(r, f))) for (r, _, fs) in walkdir(c.banks) for f in fs)
function with_campaign(f; kw...)
    dir = mktempdir(); c = write_fake_campaign(dir; kw...); f(c, dir)
end
const FIN_LOAD = (axes = SYN_AXES, expected_k_range = (SYN_KMIN, SYN_KMAX))

@testset "Finalizer: positive path, portable load, coverage, idempotence" begin
    with_campaign() do c, dir
        before = snapshot(c)
        manifest = fin(c)
        @test isfile(manifest) && isempty(leftovers(c))
        after = snapshot(c)
        for (path, state) in before; @test after[path] == state; end            # every payload/merge manifest untouched (hash+mtime)
        m = TOML.parsefile(manifest)
        @test m["status"] == "candidate" && sort([e["count"] for e in m["count"]]) == [9, 9]   # exactly the requested count (direct+ratio)
        direct = Blast.load_tensor_bank(manifest; count=9, representation=:direct, FIN_LOAD...)
        ratio = Blast.load_tensor_bank(manifest; count=9, representation=:ratio, FIN_LOAD...)
        @test direct.tensors.T_2_00 == syn_tensor("direct", "T_2_00", 9)
        @test ratio.TT.T_0_20 == syn_tensor("ratio_TT", "T_0_20", 9) && ratio.PHI == syn_tensor("ratio_PHI", "T_2_00", 9)
        pmm, qT, phi = syn_q(syn_table(), ratio.grid.sites)
        @test ratio.q.qPHI == phi && isapprox(ratio.q.qTT, pmm; rtol=8eps())
        @test_throws ArgumentError Blast.load_tensor_bank(manifest; count=11, representation=:direct, FIN_LOAD...)   # no coverage claim for others
        # Relocation as a unit keeps the signature.
        moved = joinpath(mktempdir(), "banks_moved"); cp(c.banks, moved)
        relocated = Blast.load_tensor_bank(joinpath(moved, "bundle_synthetic.toml"); count=9, representation=:ratio, FIN_LOAD...)
        @test Blast.bank_signature(relocated) == Blast.bank_signature(ratio)
        # Exact idempotence: identical rerun leaves everything byte- and mtime-identical.
        state = snapshot(c)
        @test fin(c) == manifest
        @test snapshot(c) == state && isempty(leftovers(c))
        @test m["provenance"]["generator_ace_version"] == "vendor_head:abc123"
        @test occursin("direct:9,ratio:9", m["provenance"]["basis_compat_evidence"])
    end
end

@testset "Finalizer: several counts" begin
    with_campaign(counts=(9, 11)) do c, dir
        manifest = fin(c; counts=(9, 11))
        m = TOML.parsefile(manifest)
        @test sort([(e["count"], e["representation"]) for e in m["count"]]) == [(9, "direct"), (9, "ratio"), (11, "direct"), (11, "ratio")]
        @test Blast.load_tensor_bank(manifest; count=11, representation=:ratio, FIN_LOAD...) isa Blast.RatioSplineBank
    end
end

@testset "Finalizer: rejections publish nothing, overwrite nothing, leave no staging files" begin
    function rejects_cleanly(mutate; why=r"", kw...)
        with_campaign() do c, dir
            mutate(c)
            before = snapshot(c)
            ok = frejects(c; why, kw...)
            after = snapshot(c)
            ok && !isfile(joinpath(c.banks, "bundle_synthetic.toml")) && isempty(leftovers(c)) &&
                all(p -> haskey(after, p) && (!startswith(p, joinpath(c.banks, "N")) || after[p] == before[p]), keys(before))
        end
    end
    @test rejects_cleanly(c -> nothing; counts=(9, 11), why=r"count 11 has no completed merge manifest")     # missing count
    @test rejects_cleanly(c -> nothing; status="certified", why=r"candidate")
    @test rejects_cleanly(c -> nothing; manifest_name="../x.toml", why=r"plain")
    @test rejects_cleanly(c -> nothing; counts=Int[], why=r"no counts")
    @test rejects_cleanly(c -> nothing; expected_k_range=(2 * SYN_KMIN, SYN_KMAX), why=r"trusted range")
    @test rejects_cleanly(c -> nothing; axes=(; SYN_AXES..., chi=[1.0, 2.5, 4.0, 8.0]), why=r"axes differ")
    # campaign input corruption
    @test rejects_cleanly(c -> (p = joinpath(c.inputs, "sites_ratio_9.txt"); write(p, read(p, String) * "\n# tampered\n")); why=r"fails its checksum")
    @test rejects_cleanly(c -> rm(joinpath(c.inputs, "knots_direct_9.txt")); why=r"missing file")
    # bank payload corruption: checksum, truncation (short data), wrong header
    bankfile(c, rel) = joinpath(c.banks, rel)
    @test rejects_cleanly(c -> (p = bankfile(c, "N009/direct/T_2_00.npy"); b = read(p); b[end] ⊻= 0x01; write(p, b)); why=r"fails its SHA-256")
    @test rejects_cleanly(c -> (p = bankfile(c, "N009/ratio_TT/T_0_00.npy"); write(p, read(p)[1:end-8])); why=r"bytes")
    @test rejects_cleanly(c -> rm(bankfile(c, "N009/ratio_PHI/T_2_00.npy")); why=r"missing file")
    # merge manifest inconsistencies
    function edit_merge(c, f)
        p = joinpath(c.banks, "N009", "manifest.json"); m = JSON.parsefile(p); f(m); write_json(p, m)
    end
    @test rejects_cleanly(c -> edit_merge(c, m -> m["shape"] = [3, 4, 5, 10]); why=r"shape")
    @test rejects_cleanly(c -> edit_merge(c, m -> m["inputs_manifest_sha256"] = "0"^64); why=r"different campaign inputs")
    @test rejects_cleanly(c -> edit_merge(c, m -> pop!(m["banks"])); why=r"exactly the 21")
    @test rejects_cleanly(c -> edit_merge(c, m -> m["banks"][1]["path"] = "N009/direct/other.npy"); why=r"is not `N009")
    @test rejects_cleanly(c -> edit_merge(c, m -> m["banks"][1]["path"] = "direct/T_2_00.npy"); why=r"relative to the bank root")   # old N-less path
    @test rejects_cleanly(c -> edit_merge(c, m -> m["dtype"] = "float32"); why=r"dtype")
    @test rejects_cleanly(c -> edit_merge(c, m -> m["banks"][2]["sha256"] = "0"^64); why=r"SHA-256")
    # bank_mapping disagreeing with the loader's component table
    @test rejects_cleanly(c -> (p = joinpath(c.campaign, "bank_mapping.json"); m = JSON.parsefile(p);
                               m["w_2_00_phiTT"]["ratio"]["family"] = "ratio_T"; write_json(p, m)); why=r"wrong ratio bank")
    # counts merged from different generators
    with_campaign(counts=(9, 11)) do c, dir
        p = joinpath(c.banks, "N011", "manifest.json"); m = JSON.parsefile(p); m["generator_sha256"] = syn_sha_text("generator-B"); write_json(p, m)
        @test frejects(c; counts=(9, 11), why=r"disagree on generator")
        @test !isfile(joinpath(c.banks, "bundle_synthetic.toml")) && isempty(leftovers(c))
    end
end

@testset "Finalizer: existing metadata and manifests are never overwritten" begin
    with_campaign() do c, dir
        manifest = fin(c)
        # A different existing manifest of the same name: refuse, bytes and mtime unchanged.
        before = (read(manifest), mtime(manifest))
        write(manifest, read(manifest, String) * "\n# edited\n"); edited = (read(manifest), mtime(manifest))
        @test frejects(c; why=r"already exists with different content")
        @test (read(manifest), mtime(manifest)) == edited && isempty(leftovers(c))
        write(manifest, before[1])
        # Existing identical metadata is reused without rewriting.
        axes_file = joinpath(c.banks, "meta", "axes_chi.npy"); state = (syn_sha(axes_file), mtime(axes_file))
        @test fin(c; name="bundle_second.toml") == joinpath(realpath(c.banks), "bundle_second.toml")
        @test (syn_sha(axes_file), mtime(axes_file)) == state
        # Existing metadata with different content: fail, file untouched.
        write(axes_file, vcat(read(axes_file)[1:end-8], zeros(UInt8, 8)))
        tampered = (syn_sha(axes_file), mtime(axes_file))
        @test frejects(c; name="bundle_third.toml", why=r"refusing to overwrite")
        @test (syn_sha(axes_file), mtime(axes_file)) == tampered
        @test !isfile(joinpath(c.banks, "bundle_third.toml")) && isempty(leftovers(c))
    end
end

@testset "Finalizer: command line contract" begin
    @test BundleFinalizer.main(["--help"]) === nothing
    @test_throws ArgumentError BundleFinalizer.main(["--bogus", "1"])
    @test_throws ArgumentError BundleFinalizer.main(["--campaign", "x"])                 # missing required options
    @test_throws ArgumentError BundleFinalizer.main(["--campaign", "x", "--bankroot", "y", "--counts", "a,b"])
    @test_throws ArgumentError BundleFinalizer.main(["--campaign"])                       # option without value
end

# ---- publication-path safety: nothing may be written, read-through or removed outside the declared root ----
dir_state(d) = Dict(relpath(joinpath(r, f), d) => (read(joinpath(r, f)), mtime(joinpath(r, f)))
                    for (r, _, fs) in walkdir(d) for f in fs)
function try_symlink(target, link)
    try; symlink(target, link); return true; catch err; err isa Base.IOError || rethrow(); return false; end
end
no_partials(dirs...) = all(isempty([f for (r, _, fs) in walkdir(d) for f in fs if occursin(".partial.", f)]) for d in dirs)

@testset "Finalizer: symbolic links and traversal never escape the declared root" begin
    # meta/ is a symlink to a directory OUTSIDE the bank root.
    with_campaign() do c, dir
        outside = mkpath(joinpath(dir, "outside")); write(joinpath(outside, "keep.txt"), "keep")
        if try_symlink(outside, joinpath(c.banks, "meta"))
            before = dir_state(outside); banks_before = snapshot(c)
            @test frejects(c; why=r"symbolic link")
            @test dir_state(outside) == before && no_partials(outside, c.banks)      # outside names, bytes, mtimes unchanged
            @test !isfile(joinpath(c.banks, "bundle_synthetic.toml"))
            @test all(p -> snapshot(c)[p] == banks_before[p], keys(banks_before))
        end
    end
    # The manifest target is a symlink to a file outside.
    with_campaign() do c, dir
        outside = mkpath(joinpath(dir, "outside")); write(joinpath(outside, "target.toml"), "outside manifest")
        if try_symlink(joinpath(outside, "target.toml"), joinpath(c.banks, "bundle_synthetic.toml"))
            before = dir_state(outside)
            @test frejects(c; why=r"symbolic link")
            @test dir_state(outside) == before && no_partials(outside, c.banks)
        end
    end
    # A meta file inside a real meta/ directory that is a symlink to a file outside.
    with_campaign() do c, dir
        mkpath(joinpath(c.banks, "meta")); outside = mkpath(joinpath(dir, "outside")); write(joinpath(outside, "axes.npy"), "outside")
        if try_symlink(joinpath(outside, "axes.npy"), joinpath(c.banks, "meta", "axes_ell.npy"))
            before = dir_state(outside)
            @test frejects(c; why=r"symbolic link")
            @test dir_state(outside) == before && no_partials(outside, c.banks)
        end
    end
    # inputs/manifest.json names that traverse, are absolute, or are symlinks are rejected BEFORE any read/hash.
    function with_input_name(name, expect; setup=nothing)
        with_campaign() do c, dir
            outside = mkpath(joinpath(dir, "outside")); write(joinpath(outside, "secret.txt"), "secret\n")
            isnothing(setup) || setup(c, outside) === false && return
            p = joinpath(c.inputs, "manifest.json"); m = JSON.parsefile(p)
            m["input_sha256"][name] = syn_sha(joinpath(outside, "secret.txt")); write_json(p, m)
            before = dir_state(outside)
            @test frejects(c; why=expect)
            @test dir_state(outside) == before && no_partials(outside, c.banks, c.inputs) && !isfile(joinpath(c.banks, "bundle_synthetic.toml"))
        end
    end
    with_input_name("../../outside/secret.txt", r"escapes or is not normalized")
    with_input_name("/etc/passwd", r"absolute")
    with_input_name("sub\\x.txt", r"malformed")
    with_input_name("link.txt", r"symbolic link"; setup=(c, outside) -> try_symlink(joinpath(outside, "secret.txt"), joinpath(c.inputs, "link.txt")))
    # A declared input file replaced by a symlink (same bytes) is rejected, not followed.
    with_campaign() do c, dir
        outside = mkpath(joinpath(dir, "outside")); copy = joinpath(outside, "ell.txt"); cp(joinpath(c.inputs, "ell.txt"), copy)
        rm(joinpath(c.inputs, "ell.txt"))
        if try_symlink(copy, joinpath(c.inputs, "ell.txt"))
            before = dir_state(outside)
            @test frejects(c; why=r"symbolic link")
            @test dir_state(outside) == before && no_partials(outside, c.banks)
        end
    end
    # Every input the run reads must be declared in input_sha256.
    with_campaign() do c, dir
        p = joinpath(c.inputs, "manifest.json"); m = JSON.parsefile(p); delete!(m["input_sha256"], "knots_ratio_9.txt"); write_json(p, m)
        @test frejects(c; why=r"not declared")
        @test no_partials(c.banks) && !isfile(joinpath(c.banks, "bundle_synthetic.toml"))
    end
end

@testset "Finalizer: the real N033 merge manifest layout fits the parser (metadata only)" begin
    mm = JSON.parsefile(joinpath(@__DIR__, "real_merge_manifest_N033_shape.json"))
    identity, listed = BundleFinalizer.check_merge_manifest(mm, 33, (22, 128, 64), mm["inputs_manifest_sha256"])
    @test identity == (mm["generator_sha256"], mm["inputs_manifest_sha256"], 131073)
    @test length(listed) == 21 && Set(listed) == Set(BundleFinalizer.ALL_BANKS)
    @test all(b -> b["path"] == "N033/$(b["family"])/$(b["geometry"]).npy", mm["banks"])    # paths relative to the bank root
    @test all(b -> b["bytes"] > 0 && occursin(r"^[0-9a-f]{64}$", b["sha256"]), mm["banks"])
    # Wrong expectations are still rejected against the same real layout.
    @test_throws ArgumentError BundleFinalizer.check_merge_manifest(mm, 49, (22, 128, 64), mm["inputs_manifest_sha256"])
    @test_throws ArgumentError BundleFinalizer.check_merge_manifest(mm, 33, (22, 128, 52), mm["inputs_manifest_sha256"])
    @test_throws ArgumentError BundleFinalizer.check_merge_manifest(mm, 33, (22, 128, 64), "0"^64)
    @test BundleFinalizer.npy_header_bytes isa Function
end

@testset "Finalizer: NPY metadata writer is explicitly little-endian and bit exact" begin
    path = joinpath(mktempdir(), "v.npy"); v = [1.0, -0.0, 2.5e-300, Inf, 123456.789]
    BundleFinalizer.npy_write(path, v)
    bytes = read(path); hb = BundleFinalizer.npy_header_bytes(path)
    @test bytes[1:6] == b"\x93NUMPY" && hb % 64 == 0 && length(bytes) == hb + 8 * length(v)
    @test ltoh.(reinterpret(UInt64, bytes[hb+1:end])) == reinterpret(UInt64, v)
    header = Blast._npy_header(path)
    @test header.descr == "<f8" && header.fortran && header.shape == (5,) && header.header_bytes == hb
end

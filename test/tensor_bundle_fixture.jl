# Tiny public-shaped synthetic tensor bundles, written into temporary directories.
# Index-coded tensors (value encodes ell, chi, R, coefficient, and family) make any
# permutation or family swap visible. Algebra/API evidence only, not bank science.
using Blast, TOML, SHA
using Blast.AbstractCosmologicalEmulators: CubicSpline, CubicBSplinePlan, bspline_basis, knot_vector

const SYN_AXES = (ell = [9.0, 5.0, 2.0], chi = [1.0, 2.5, 4.0, 7.0], R = [0.2, 0.5, 0.8, 0.95, 1.0])
const SYN_KMIN, SYN_KMAX = 1.0e-3, 10.0
const SYN_LN10AS, SYN_NS = 3.044, 0.9649
const SYN_ALL_COUNTS = (9, 11)
const SYN_DIRECT = Tuple(("direct", string(g)) for g in Blast._BANK_GEOMETRIES)
const SYN_RATIO = Blast._RATIO_BANKS

syn_sha(path) = bytes2hex(sha256(read(path)))
syn_sha_text(text) = bytes2hex(sha256(text))

function syn_write_npy(path, A::AbstractArray; descr="<f8", fortran=true)
    T = descr == "<f8" ? Float64 : Float32
    shape = join(size(A), ", ") * (ndims(A) == 1 ? "," : "")
    header = "{'descr': '$descr', 'fortran_order': $(fortran ? "True" : "False"), 'shape': ($shape), }"
    total = 10 + length(header) + 1
    header *= " "^(mod(-total, 64)) * "\n"
    mkpath(dirname(path))
    open(path, "w") do io
        write(io, b"\x93NUMPY", UInt8(1), UInt8(0), UInt16(length(header)), header)
        write(io, fortran ? T.(A) : T.(permutedims(A, ndims(A):-1:1)))
    end
    return path
end

# Exact NPY header byte count read from the file itself (10 + header length for NPY v1).
syn_header_bytes(path) = open(io -> (skip(io, 8); 10 + Int(read(io, UInt16))), path)

function syn_ref(root, rel; shape=nothing)
    path = joinpath(root, rel)
    return Dict{String,Any}("file" => rel, "sha256" => syn_sha(path), "bytes" => filesize(path),
                            "header_bytes" => syn_header_bytes(path))
end

syn_sites(n, representation; kmin=SYN_KMIN, kmax=SYN_KMAX) = begin
    p = representation == "direct" ? 1.3 : 1.7
    x = log10(kmin) .+ (log10(kmax) - log10(kmin)) .* range(0, 1; length=n) .^ p
    phys = 10 .^ x; phys[1] = kmin; phys[end] = kmax
    phys
end

# Independent (non-loader) fiducial functions for witnesses.
syn_table(kmin=SYN_KMIN, kmax=SYN_KMAX) = (k = 10 .^ range(log10(kmin), log10(kmax); length=64); k[1] = kmin; k[end] = kmax;
               hcat(k, @. 2.0e4 * k^0.96 / (1 + (k / 0.2)^2)^2))
function syn_q(table, sites)
    pmm = exp.(CubicSpline(log.(table[:, 2]), log.(table[:, 1]))(log.(sites)))
    As = exp(SYN_LN10AS) * 1e-10
    phi = @. (9 / 25) * 2 * pi^2 * As / sites^3 * (sites / 0.05)^(SYN_NS - 1)
    return pmm, sqrt.(phi .* pmm), phi
end

syn_code(family, geometry) = 1.0e5 * (findfirst(==((family, geometry)), vcat(collect(SYN_DIRECT), collect(SYN_RATIO))))
syn_tensor(family, geometry, n; dims=(3,4,5)) = [1000.0a + 100.0b + 10.0c + d + syn_code(family, geometry)
                                   for a in 1:dims[1], b in 1:dims[2], c in 1:dims[3], d in 1:n]

function write_synthetic_bundle(root; counts=SYN_ALL_COUNTS, witness=true, kmin=SYN_KMIN, kmax=SYN_KMAX,
                                grid_axes=SYN_AXES)
    mkpath(root)
    dims=map(length,values(grid_axes))
    for (name, v) in pairs(grid_axes); syn_write_npy(joinpath(root, "meta", "axes_$name.npy"), v); end
    table = syn_table(kmin, kmax); syn_write_npy(joinpath(root, "meta", "fiducial_pmm.npy"), table)
    entries = Dict{String,Any}[]
    for n in counts, rep in ("direct", "ratio")
        sites = syn_sites(n, rep; kmin, kmax); logsites = log10.(sites)
        knots = knot_vector(bspline_basis(CubicBSplinePlan(logsites, logsites; extrapolation=:throw)))
        tag = "N$(lpad(n, 3, '0'))_$rep"
        syn_write_npy(joinpath(root, "meta", "sites_phys_$tag.npy"), sites)
        syn_write_npy(joinpath(root, "meta", "sites_log10_$tag.npy"), logsites)
        syn_write_npy(joinpath(root, "meta", "knots_$tag.npy"), knots)
        entry = Dict{String,Any}("count" => n, "representation" => rep,
            "sites_phys" => syn_ref(root, "meta/sites_phys_$tag.npy"),
            "sites_log10" => syn_ref(root, "meta/sites_log10_$tag.npy"),
            "knots" => syn_ref(root, "meta/knots_$tag.npy"))
        if witness && rep == "ratio"
            q = syn_q(table, sites)
            syn_write_npy(joinpath(root, "meta", "witness_$tag.npy"), hcat(q...))
            entry["witness"] = syn_ref(root, "meta/witness_$tag.npy")
        end
        entry["bank"] = Dict{String,Any}[]
        for (family, geometry) in (rep == "direct" ? SYN_DIRECT : SYN_RATIO)
            rel = "N$(lpad(n, 3, '0'))/$family/$geometry.npy"
            syn_write_npy(joinpath(root, rel), syn_tensor(family, geometry, n;dims))
            push!(entry["bank"], Dict{String,Any}("family" => family, "geometry" => geometry,
                                                    "file_ref" => syn_ref(root, rel; shape=(dims...,n))))
        end
        push!(entries, entry)
    end
    src(s) = s == :ϕTT ? "TT" : s == :ϕT ? "T" : "PHI"
    fam = Dict("TT" => "ratio_TT", "T" => "ratio_T", "PHI" => "ratio_PHI")
    source_map = [Dict{String,Any}("component" => Blast._ascii_component(name), "source" => src(source),
                                   "orientation" => r1 ? "R1" : src(source),
                                   "family_direct" => ["direct", string(tfield)],
                                   "family_ratio" => [fam[src(source)], string(tfield)])
                  for (name, tfield, source, r1) in Blast._PROJECTED_MATTER_COMPONENTS]
    manifest = Dict{String,Any}(
        "schema" => Blast.BUNDLE_SCHEMA, "status" => "candidate", "created_utc" => "2026-01-01T00:00:00",
        "convention" => Dict{String,Any}(
            "coordinate" => Blast.BUNDLE_CONVENTION.coordinate, "basis" => Blast.BUNDLE_CONVENTION.basis,
            "coefficient" => Blast.BUNDLE_CONVENTION.coefficient, "layout" => Blast.BUNDLE_CONVENTION.layout,
            "endpoint" => Blast.BUNDLE_CONVENTION.endpoint, "endpoint_epsilon" => Blast.BUNDLE_ENDPOINT_EPSILON,
            "direct_source" => Blast.BUNDLE_CONVENTION.direct_source, "ratio_source" => Blast.BUNDLE_CONVENTION.ratio_source),
        "provenance" => Dict{String,Any}("generator_sha256" => syn_sha_text("generator"), "inputs_manifest_sha256" => syn_sha_text("inputs"),
            "quadrature_N" => 131073, "generator_ace_version" => "0.12.0", "runtime_ace_rev" => "synthetic",
            "basis_compat_evidence" => "synthetic fixture"),
        "range" => Dict{String,Any}("kmin" => kmin, "kmax" => kmax),
        "axes" => Dict{String,Any}(string(name) => syn_ref(root, "meta/axes_$name.npy") for name in keys(grid_axes)),
        "fiducial" => Dict{String,Any}("table" => syn_ref(root, "meta/fiducial_pmm.npy"), "rows" => size(table, 1),
            "ln10_10_As" => SYN_LN10AS, "n_s" => SYN_NS, "pivot" => 0.05,
            "pmm_interpolator" => Blast.BUNDLE_FIDUCIAL.pmm_interpolator, "phi_formula" => Blast.BUNDLE_FIDUCIAL.phi_formula,
            "q_TT" => Blast.BUNDLE_FIDUCIAL.q_TT, "q_T" => Blast.BUNDLE_FIDUCIAL.q_T, "q_PHI" => Blast.BUNDLE_FIDUCIAL.q_PHI),
        "count" => entries, "source_map" => source_map)
    save_synthetic_manifest(root, manifest)
    return manifest
end

save_synthetic_manifest(root, manifest) = open(io -> TOML.print(io, manifest), joinpath(root, "bundle.toml"), "w")
load_synthetic_manifest(root) = TOML.parsefile(joinpath(root, "bundle.toml"))
syn_bank(manifest, n, rep, family, geometry) = only(b for e in manifest["count"] if e["count"] == n && e["representation"] == rep
                                                     for b in e["bank"] if b["family"] == family && b["geometry"] == geometry)
syn_entry(manifest, n, rep) = only(e for e in manifest["count"] if e["count"] == n && e["representation"] == rep)
syn_refresh!(root, ref) = (p = joinpath(root, ref["file"]); ref["sha256"] = syn_sha(p); ref["bytes"] = filesize(p); ref["header_bytes"] = syn_header_bytes(p); ref)
syn_load(root; n=9, rep=:direct, axes=SYN_AXES, expected_k_range=(SYN_KMIN, SYN_KMAX)) =
    Blast.load_tensor_bank(joinpath(root, "bundle.toml"); count=n, representation=rep, axes=axes, expected_k_range=expected_k_range)

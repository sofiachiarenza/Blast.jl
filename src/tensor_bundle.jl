# Cold loader for portable cubic B-spline tensor bundles (phase 4A).
#
# A bundle is a directory holding `bundle.toml`, small `meta/*.npy` payloads and
# merged tensor payloads `N<count>/<family>/<geometry>.npy`. Paths in the
# manifest are relative to the bundle root, so the directory moves as a unit.
# Nothing here runs in differentiated code: all validation, hashing and the
# fixed normalizers are computed once, cold.

const BUNDLE_SCHEMA = "blast.tensor_bank_bundle/1"
const _BUNDLE_SIGNATURE_FORMAT = "blast-tensor-bank-signature-v2"
const _BUNDLE_STATUSES = ("candidate", "certified")
const _BUNDLE_REPRESENTATIONS = (:direct, :ratio)
const _RATIO_T_GEOMETRIES = ("T_2_00", "T_0_00", "T_2_02", "T_2_20")
const _DIRECT_BANKS = Tuple(("direct", String(g)) for g in _BANK_GEOMETRIES)
const _RATIO_BANKS = (Tuple(("ratio_TT", String(g)) for g in _BANK_GEOMETRIES)...,
                      Tuple(("ratio_T", g) for g in _RATIO_T_GEOMETRIES)...,
                      ("ratio_PHI", "T_2_00"))

# Supported semantic contract: machine IDs plus explicit numeric parameters. Each value
# is what THIS runtime implements; a manifest claiming anything else is rejected, never
# merely hashed. (generator: ln/ln natural fiducial spline, As = exp(l)*1e-10, raw q.)
const BUNDLE_CONVENTION = (
    coordinate = "log10_k_per_mpc",
    basis = "ace_cubic_bspline_notaknot_degree3_extrapolation_throw",
    coefficient = "bspline_coefficients_collocation_solve",
    layout = "float64_ell_chi_r_coef_fortran",
    endpoint = "inward_quadrature_nodes_relative_epsilon",
    direct_source = "unweighted_shared_bank",
    ratio_source = "fixed_q_divide_full_physical_source_before_solve",
)
const BUNDLE_ENDPOINT_EPSILON = 1.0e-8
const BUNDLE_FIDUCIAL = (
    pmm_interpolator = "ln_ln_natural_cubic_spline_ace",
    phi_formula = "generator_phi_v1_9_25_2pi2_As_over_k3_k_over_pivot_pow_ns_minus_1_As_exp_l_times_1e_minus_10",
    q_TT = "pmm_fid",
    q_T = "sqrt_phi_fid_times_pmm_fid",
    q_PHI = "phi_fid",
)
const BUNDLE_PIVOT = 0.05
"Trusted production k range [1/Mpc]: the legacy `k_cheb` endpoints (constants.jl)."
const PRODUCTION_K_RANGE = (5.0e-5, 16.0)
const _MAX_COUNT = 1 << 16
const _MAX_TABLE_ROWS = 1 << 24
const _HEX64 = r"^[0-9a-f]{64}$"

const _TOP_KEYS = ("schema", "status", "created_utc", "convention", "provenance", "range", "axes", "fiducial", "count", "source_map")
const _CONVENTION_KEYS = (keys(BUNDLE_CONVENTION)..., :endpoint_epsilon)
const _PROVENANCE_KEYS = ("generator_sha256", "inputs_manifest_sha256", "quadrature_N", "generator_ace_version",
                    "runtime_ace_rev", "basis_compat_evidence")
const _FIDUCIAL_KEYS = ("table", "rows", "ln10_10_As", "n_s", "pivot", "pmm_interpolator", "phi_formula", "q_TT", "q_T", "q_PHI")
const _COUNT_KEYS = ("count", "representation", "sites_phys", "sites_log10", "knots", "witness", "bank")

"Production (ℓ, χ, R) axes the bundle must match exactly (tests pass explicit synthetic axes)."
production_axes() = (ell = ℓ_nonlimber, chi = χ, R = R)

_bundle_error(message) = throw(ArgumentError("tensor bundle: " * message))

function _only_keys(table, allowed, where)
    for key in keys(table)
        String(key) in String.(collect(allowed)) || _bundle_error("unsupported key `$key` in $where")
    end
end

function _req_int(table, key; min=1, max=typemax(Int))
    haskey(table, key) && table[key] isa Integer && !(table[key] isa Bool) || _bundle_error("missing or non-integer `$key`")
    min ≤ table[key] ≤ max || _bundle_error("`$key` = $(table[key]) is out of range [$min, $max]")
    return Int(table[key])
end

function _req(table, key, ::Type{T}) where {T}
    haskey(table, key) && table[key] isa T || _bundle_error("missing or mistyped `$key`")
    return table[key]
end

# ── Paths and files ──────────────────────────────────────────────────────────
function _check_relpath(rel::AbstractString)
    (isempty(rel) || startswith(rel, "/") || isabspath(rel) || occursin('\\', rel) || occursin('\0', rel)) &&
        _bundle_error("path `$rel` is empty, absolute or malformed")
    any(s -> s in ("", ".", ".."), split(rel, '/')) && _bundle_error("path `$rel` escapes or is not normalized")
    return rel
end

function _fileref(table, key)
    ref = _req(table, key, AbstractDict)
    file = _check_relpath(_req(ref, "file", String))
    sha = _req(ref, "sha256", String)
    occursin(r"^[0-9a-f]{64}$", sha) || _bundle_error("`$key` has a malformed sha256")
    _only_keys(ref, ("file", "sha256", "bytes", "header_bytes"), "file reference `$key`")
    bytes = _req_int(ref, "bytes"; min=1)
    header_bytes = _req_int(ref, "header_bytes"; min=16, max=4160)     # mandatory for every NPY payload
    return (; file, sha, bytes, header_bytes)
end

function _resolve(root_real, rel)
    path = root_real
    for segment in split(rel, '/')
        path = joinpath(path, segment)
        islink(path) && _bundle_error("path `$rel` contains a symbolic link")
        ispath(path) || _bundle_error("missing file `$rel`")
    end
    isfile(path) || _bundle_error("`$rel` is not a regular file")
    startswith(realpath(path), root_real * "/") || _bundle_error("path `$rel` resolves outside the bundle root")
    return path
end

"Resolve, then require exact byte count and SHA-256 (streamed)."
function _verified_path(root_real, ref)
    path = _resolve(root_real, ref.file)
    filesize(path) == ref.bytes || _bundle_error("`$(ref.file)`: byte count $(filesize(path)) != manifest $(ref.bytes)")
    bytes2hex(open(sha256, path)) == ref.sha || _bundle_error("`$(ref.file)`: SHA-256 mismatch")
    return path
end

# Bounded minimal NPY (v1/v2) header reader: validates, never a general framework.
function _npy_header(path)
    open(path) do io
        read(io, 6) == b"\x93NUMPY" || _bundle_error("`$(basename(path))` is not an NPY file")
        major = read(io, UInt8); read(io, UInt8)
        len = major == 1 ? Int(ltoh(read(io, UInt16))) : major in (2, 3) ? Int(ltoh(read(io, UInt32))) :
              _bundle_error("unsupported NPY version")
        len ≤ 4096 || _bundle_error("NPY header too large")
        text = String(read(io, len))
        descr = match(r"'descr':\s*'([^']*)'", text)
        fortran = match(r"'fortran_order':\s*(True|False)", text)
        shape = match(r"'shape':\s*\(([^)]*)\)", text)
        (isnothing(descr) || isnothing(fortran) || isnothing(shape)) && _bundle_error("malformed NPY header")
        dims = Tuple(parse(Int, strip(s)) for s in split(shape[1], ',') if !isempty(strip(s)))
        return (; descr = descr[1], fortran = fortran[1] == "True", shape = dims, header_bytes = position(io))
    end
end

"Load a payload after header/size checks: Float64 little-endian, Fortran order, expected shape, finite."
function _load_array(root_real, ref, shape::Tuple)
    all(d -> 1 ≤ d ≤ _MAX_TABLE_ROWS, shape) || _bundle_error("`$(ref.file)`: implausible expected shape $shape")
    path = _verified_path(root_real, ref)
    header = _npy_header(path)
    header.descr == "<f8" || _bundle_error("`$(ref.file)`: dtype `$(header.descr)` is not little-endian Float64")
    header.fortran || _bundle_error("`$(ref.file)`: not Fortran order")
    header.shape == shape || _bundle_error("`$(ref.file)`: shape $(header.shape) != expected $shape")
    ref.header_bytes == header.header_bytes || _bundle_error("`$(ref.file)`: header byte count differs from manifest")
    payload = try
        reduce(Base.checked_mul, shape; init=8)
    catch
        _bundle_error("`$(ref.file)`: payload size overflows")
    end
    filesize(path) == header.header_bytes + payload || _bundle_error("`$(ref.file)`: truncated or padded payload")
    data = try
        npzread(path)
    catch err
        _bundle_error("`$(ref.file)`: unreadable ($(sprint(showerror, err)))")
    end
    (data isa Array{Float64} && size(data) == shape) || _bundle_error("`$(ref.file)`: unexpected array type or shape")
    all(isfinite, data) || _bundle_error("`$(ref.file)`: non-finite values")
    return data
end

# ── Fixed fiducial normalizers (exactly the generator contract) ──────────────
# Generator: As = exp(ln10_10_As)*1e-10 (NOT exp(.)/1e10: different rounding).
_fiducial_phi(k, ln10_10_As, n_s) = (As = exp(ln10_10_As) * 1e-10;
    @. (9 / 25) * 2 * pi^2 * As / k^3 * (k / 0.05)^(n_s - 1))

"""
    FiducialNormalizers

Raw fixed weights at a design's physical sites, from saved fiducial data (never the
sampled cosmology): `qTT = Pmm_fid(k)`, `qT = sqrt(PPhi_fid*Pmm_fid)`,
`qPHI = PPhi_fid`, with `Pmm_fid = exp(natural cubic spline of ln Pmm vs ln k)`.
"""
struct FiducialNormalizers
    qTT::Vector{Float64}
    qT::Vector{Float64}
    qPHI::Vector{Float64}
end

# Final q vectors must be finite and positive AFTER the product/sqrt: positive finite factors
# can still overflow or underflow. Generator formulas are not changed to repair such data.
function _normalizers(pmm::AbstractVector, phi::AbstractVector)
    all(v -> all(x -> isfinite(x) && x > 0, v), (pmm, phi)) ||
        _bundle_error("fiducial normalizers qTT/qT/qPHI are not all finite and positive (overflow, underflow or invalid values)")
    qT = sqrt.(phi .* pmm)
    all(v -> all(x -> isfinite(x) && x > 0, v), (pmm, phi, qT)) ||
        _bundle_error("fiducial normalizers qTT/qT/qPHI are not all finite and positive (overflow, underflow or invalid values)")
    return FiducialNormalizers(Vector{Float64}(pmm), Vector{Float64}(qT), Vector{Float64}(phi))
end

function _fiducial_normalizers(table::Matrix{Float64}, sites::Vector{Float64}, ln10_10_As, n_s)
    spline = AbstractCosmologicalEmulators.CubicSpline(log.(table[:, 2]), log.(table[:, 1]))
    pmm = exp.(spline(log.(sites)))
    return _normalizers(pmm, _fiducial_phi(sites, ln10_10_As, n_s))
end

# ── Types ────────────────────────────────────────────────────────────────────
"Fixed numerical grid of one design: physical sites, authoritative log₁₀ sites, full knot vector."
struct SplineGrid
    sites::Vector{Float64}
    logsites::Vector{Float64}
    knots::Vector{Float64}
end

"Direct (unweighted) spline bank: the eight shared tensors plus the grid. `tensors.signature` is the full v2 signature."
struct DirectSplineBank{A<:AbstractArray{Float64,4}}
    tensors::TensorBank{A}
    grid::SplineGrid
end

"Ratio spline bank: TT (8), T (4), PHI (1) tensors, grid, fixed normalizers and the v2 signature."
struct RatioSplineBank{A<:AbstractArray{Float64,4}}
    TT::NamedTuple{_BANK_GEOMETRIES, NTuple{8, A}}
    T::NamedTuple{(:T_2_00, :T_0_00, :T_2_02, :T_2_20), NTuple{4, A}}
    PHI::A
    grid::SplineGrid
    q::FiducialNormalizers
    signature::BankSignature
    axes::BankAxes
end
RatioSplineBank(TT,T,PHI,grid,q,signature)=RatioSplineBank(TT,T,PHI,grid,q,signature,_production_bank_axes())

bank_signature(b::DirectSplineBank) = b.tensors.signature
bank_signature(b::RatioSplineBank) = b.signature
bank_axes(b::DirectSplineBank)=bank_axes(b.tensors)
bank_axes(b::RatioSplineBank)=b.axes

# ── Signature v2 ─────────────────────────────────────────────────────────────
function _digest_ints!(ctx, v)
    SHA.update!(ctx, reinterpret(UInt8, [htol(UInt64(x)) for x in v]))
end

# Typed, length-prefixed fields in a FIXED order (inputs are allowlisted and type-checked, so
# nothing depends on dictionary order, extra keys or repr()).
function _bundle_signature(representation, count, provenance, axes, grid, fiducial, source_rows, banks)
    ctx = SHA.SHA256_CTX()
    for text in (_BUNDLE_SIGNATURE_FORMAT, BUNDLE_SCHEMA, String(representation))
        _digest_text!(ctx, text)
    end
    _digest_ints!(ctx, (count,))
    for key in _CONVENTION_KEYS[1:end-1]; _digest_text!(ctx, String(key)); _digest_text!(ctx, getfield(BUNDLE_CONVENTION, key)); end
    _digest_floats!(ctx, [BUNDLE_ENDPOINT_EPSILON])
    for key in _PROVENANCE_KEYS
        _digest_text!(ctx, key)
        value = provenance[key]
        value isa Integer ? (_digest_text!(ctx, "int"); _digest_ints!(ctx, (value,))) :
                            (_digest_text!(ctx, "str"); _digest_text!(ctx, value))
    end
    for axis in (axes.ell, axes.chi, axes.R, grid.sites, grid.logsites, grid.knots)
        _digest_floats!(ctx, axis)
    end
    _digest_text!(ctx, fiducial.table_sha256)
    for key in keys(BUNDLE_FIDUCIAL); _digest_text!(ctx, String(key)); _digest_text!(ctx, getfield(BUNDLE_FIDUCIAL, key)); end
    _digest_floats!(ctx, [fiducial.ln10_10_As, fiducial.n_s, fiducial.pivot, fiducial.kmin, fiducial.kmax])
    for row in source_rows
        for text in row; _digest_text!(ctx, text); end
    end
    for (family, geometry, shape, bytes, sha) in banks
        _digest_text!(ctx, family); _digest_text!(ctx, geometry); _digest_ints!(ctx, shape)
        _digest_ints!(ctx, (bytes,)); _digest_text!(ctx, sha)
    end
    words = ltoh.(reinterpret(UInt64, SHA.digest!(ctx)))
    return (words[1], words[2], words[3], words[4])
end

# ── Source map: validated against the existing component table ───────────────
const _SOURCE_NAME = Dict(:ϕTT => "TT", :ϕT => "T", :ϕ => "PHI")
const _RATIO_FAMILY = Dict("TT" => "ratio_TT", "T" => "ratio_T", "PHI" => "ratio_PHI")
_ascii_component(name) = replace(string(name), "ϕ" => "phi")

function _validate_source_map(rows, representation)
    expected = Dict(_ascii_component(name) => (geometry = String(tfield), source = _SOURCE_NAME[src], r1 = r1)
                    for (name, tfield, src, r1) in _PROJECTED_MATTER_COMPONENTS)
    seen = Set{String}(); selected = Vector{NTuple{5,String}}()
    for row in rows
        row isa AbstractDict || _bundle_error("source_map entries must be tables")
        component = _req(row, "component", String)
        component in seen && _bundle_error("duplicate source_map component `$component`")
        push!(seen, component)
        haskey(expected, component) || _bundle_error("unknown source_map component `$component`")
        exp = expected[component]
        source = _req(row, "source", String); orientation = _req(row, "orientation", String)
        source == exp.source || _bundle_error("source_map `$component`: wrong source family `$source`")
        orientation == (exp.r1 ? "R1" : source) || _bundle_error("source_map `$component`: wrong orientation")
        direct = _req(row, "family_direct", AbstractVector); ratio = _req(row, "family_ratio", AbstractVector)
        (length(direct) == 2 && length(ratio) == 2 && all(x -> x isa String, vcat(direct, ratio))) ||
            _bundle_error("source_map `$component`: families must be [family, geometry]")
        direct == ["direct", exp.geometry] || _bundle_error("source_map `$component`: wrong direct bank")
        ratio == [_RATIO_FAMILY[exp.source], exp.geometry] || _bundle_error("source_map `$component`: wrong ratio bank")
        chosen = representation === :direct ? direct : ratio
        push!(selected, (component, source, orientation, chosen[1], chosen[2]))
    end
    seen == Set(keys(expected)) || _bundle_error("source_map must list exactly the 17 projected components")
    return sort!(selected)
end

# ── Loader ───────────────────────────────────────────────────────────────────
"""
    load_tensor_bank(manifest_path; count, representation, axes = production_axes(),
                     expected_k_range = PRODUCTION_K_RANGE)

Cold, fail-closed loader of one count and representation (`:direct` or `:ratio`) from
a portable bundle. The manifest must describe exactly the semantics this runtime
implements (machine IDs and numeric parameters in `BUNDLE_CONVENTION`,
`BUNDLE_FIDUCIAL`, `BUNDLE_ENDPOINT_EPSILON`); anything else is rejected, not hashed.
The k range, axes and basis are checked against TRUSTED expectations (`expected_k_range`,
`axes`, the pinned ACE basis), never against the manifest's own declarations. The whole
manifest is validated structurally; only the selected payloads, meta files and
fiducial table are opened, always with byte count, SHA-256 and finiteness checks.
Returns a [`DirectSplineBank`](@ref) or [`RatioSplineBank`](@ref) carrying the full v2
[`BankSignature`](@ref). `status` and timestamps are outside the identity.
"""
function load_tensor_bank(manifest_path::AbstractString; count::Integer, representation::Symbol,
                          axes::NamedTuple{(:ell, :chi, :R)} = production_axes(),
                          expected_k_range::NTuple{2,Float64} = PRODUCTION_K_RANGE)
    representation in _BUNDLE_REPRESENTATIONS || _bundle_error("representation must be :direct or :ratio")
    isfile(manifest_path) || _bundle_error("manifest not found")
    root = realpath(dirname(abspath(manifest_path)))
    manifest = try
        TOML.parsefile(manifest_path)
    catch err
        _bundle_error("unparseable manifest ($(sprint(showerror, err)))")
    end
    # ---- structural + semantic validation of the whole manifest (no array materialization) ----
    _only_keys(manifest, _TOP_KEYS, "the manifest")
    _req(manifest, "schema", String) == BUNDLE_SCHEMA || _bundle_error("unsupported schema")
    status = _req(manifest, "status", String)
    status in _BUNDLE_STATUSES || _bundle_error("unknown status `$status`")
    convention = _req(manifest, "convention", AbstractDict)
    _only_keys(convention, _CONVENTION_KEYS, "convention")
    for key in keys(BUNDLE_CONVENTION)
        _req(convention, String(key), String) == getfield(BUNDLE_CONVENTION, key) ||
            _bundle_error("unsupported convention.$key `$(convention[String(key)])`; this runtime implements `$(getfield(BUNDLE_CONVENTION, key))`")
    end
    _req(convention, "endpoint_epsilon", Float64) == BUNDLE_ENDPOINT_EPSILON ||
        _bundle_error("unsupported convention.endpoint_epsilon; this runtime's generator uses $BUNDLE_ENDPOINT_EPSILON")
    provenance = _req(manifest, "provenance", AbstractDict)
    _only_keys(provenance, _PROVENANCE_KEYS, "provenance")
    for key in ("generator_sha256", "inputs_manifest_sha256")
        occursin(_HEX64, _req(provenance, key, String)) || _bundle_error("provenance.$key must be 64 lowercase hex characters")
    end
    for key in ("generator_ace_version", "runtime_ace_rev", "basis_compat_evidence")
        isempty(_req(provenance, key, String)) && _bundle_error("provenance.$key is empty")
    end
    _req_int(provenance, "quadrature_N"; min=3)
    range = _req(manifest, "range", AbstractDict)
    _only_keys(range, ("kmin", "kmax"), "range")
    kmin = _req(range, "kmin", Float64); kmax = _req(range, "kmax", Float64)
    (isfinite(kmin) && isfinite(kmax) && 0 < expected_k_range[1] < expected_k_range[2]) || _bundle_error("invalid k range")
    (abs(kmin / expected_k_range[1] - 1) ≤ 4eps() && abs(kmax / expected_k_range[2] - 1) ≤ 4eps()) ||
        _bundle_error("declared k range ($kmin, $kmax) is not the trusted range $expected_k_range")
    axes_table = _req(manifest, "axes", AbstractDict)
    _only_keys(axes_table, ("ell", "chi", "R"), "axes")
    axis_refs = (ell = _fileref(axes_table, "ell"), chi = _fileref(axes_table, "chi"), R = _fileref(axes_table, "R"))
    fid = _req(manifest, "fiducial", AbstractDict)
    _only_keys(fid, _FIDUCIAL_KEYS, "fiducial")
    for key in keys(BUNDLE_FIDUCIAL)
        _req(fid, String(key), String) == getfield(BUNDLE_FIDUCIAL, key) ||
            _bundle_error("unsupported fiducial.$key `$(fid[String(key)])`; this runtime implements `$(getfield(BUNDLE_FIDUCIAL, key))`")
    end
    table_ref = _fileref(fid, "table")
    table_rows = _req_int(fid, "rows"; min=3, max=_MAX_TABLE_ROWS)
    ln10_10_As = _req(fid, "ln10_10_As", Float64); n_s = _req(fid, "n_s", Float64); pivot = _req(fid, "pivot", Float64)
    (isfinite(ln10_10_As) && isfinite(n_s) && pivot == BUNDLE_PIVOT) || _bundle_error("invalid fiducial parameters or pivot (supported pivot: $BUNDLE_PIVOT)")

    entries = _req(manifest, "count", AbstractVector)
    seen_entries = Set{Tuple{Int,String}}(); selected = nothing
    for entry in entries
        entry isa AbstractDict || _bundle_error("count entries must be tables")
        _only_keys(entry, _COUNT_KEYS, "a count entry")
        n = _req_int(entry, "count"; min=4, max=_MAX_COUNT); rep = _req(entry, "representation", String)
        rep in string.(_BUNDLE_REPRESENTATIONS) || _bundle_error("unknown representation `$rep`")
        (n, rep) in seen_entries && _bundle_error("duplicate entry for count $n, $rep")
        push!(seen_entries, (n, rep))
        for key in ("sites_phys", "sites_log10", "knots"); _fileref(entry, key); end
        haskey(entry, "witness") && _fileref(entry, "witness")
        banks = _req(entry, "bank", AbstractVector)
        wanted = rep == "direct" ? _DIRECT_BANKS : _RATIO_BANKS
        listed = Tuple{String,String}[]
        for bank in banks
            bank isa AbstractDict || _bundle_error("bank entries must be tables")
            _only_keys(bank, ("family", "geometry", "file_ref"), "a bank entry")
            family = _req(bank, "family", String); geometry = _req(bank, "geometry", String)
            (family, geometry) in listed && _bundle_error("duplicate bank $family/$geometry in count $n")
            push!(listed, (family, geometry))
            ref = _fileref(bank, "file_ref")
            ref.file == "N$(lpad(n, 3, '0'))/$family/$geometry.npy" ||
                _bundle_error("bank $family/$geometry of count $n has a non-canonical path `$(ref.file)`")
        end
        Set(listed) == Set(wanted) || _bundle_error("count $n ($rep) must list exactly its $(length(wanted)) banks")
        n == count && rep == string(representation) && (selected = entry)
    end
    source_rows = _validate_source_map(_req(manifest, "source_map", AbstractVector), representation)
    isnothing(selected) && _bundle_error("count $count / $representation is not in the bundle (available: $(sort(collect(seen_entries))))")

    # ---- selected data only ----
    ell = vec(_load_array(root, axis_refs.ell, (length(axes.ell),)))
    chi = vec(_load_array(root, axis_refs.chi, (length(axes.chi),)))
    Rax = vec(_load_array(root, axis_refs.R, (length(axes.R),)))
    (ell == axes.ell && chi == axes.chi && Rax == axes.R) || _bundle_error("axes differ from the expected axes")
    (issorted(Rax; lt = <=) && Rax[end] == 1.0) || _bundle_error("R axis must be increasing and end at 1.0")
    _strictly_monotonic(v) = issorted(v; lt = <=) || issorted(v; lt = >=)
    (_strictly_monotonic(ell) && _strictly_monotonic(chi)) || _bundle_error("ell/chi axes are not strictly monotonic")

    n_sites = Int(selected["count"])
    sites = vec(_load_array(root, _fileref(selected, "sites_phys"), (n_sites,)))
    logsites = vec(_load_array(root, _fileref(selected, "sites_log10"), (n_sites,)))
    knots = vec(_load_array(root, _fileref(selected, "knots"), (n_sites + 4,)))
    (issorted(sites; lt = <=) && issorted(logsites; lt = <=)) || _bundle_error("sites must be strictly increasing")
    (abs(sites[1] / expected_k_range[1] - 1) ≤ 4eps() && abs(sites[end] / expected_k_range[2] - 1) ≤ 4eps()) ||
        _bundle_error("site endpoints are not the trusted k range $expected_k_range")
    lo, hi = log10(expected_k_range[1]), log10(expected_k_range[2])
    (abs(logsites[1] - lo) ≤ 4eps(abs(lo)) && abs(logsites[end] - hi) ≤ 4eps(abs(hi))) ||
        _bundle_error("log₁₀ site endpoints are inconsistent with the trusted k range")
    all(abs.(10 .^ logsites ./ sites .- 1) .≤ 4eps()) || _bundle_error("physical and log₁₀ sites disagree beyond roundoff")
    plan = CubicBSplinePlan(logsites, logsites; extrapolation = :throw)
    knot_vector(bspline_basis(plan)) == knots || _bundle_error("knots differ from the not-a-knot basis rebuilt by the runtime ACE")

    table = _load_array(root, table_ref, (table_rows, 2))
    (all(>(0), table) && issorted(@view table[:, 1]; lt = <=)) || _bundle_error("fiducial table must be positive with increasing k")
    (table[1, 1] ≤ sites[1] && table[end, 1] ≥ sites[end]) || _bundle_error("fiducial table does not cover the sites")
    q = _fiducial_normalizers(table, sites, ln10_10_As, n_s)
    if haskey(selected, "witness")
        witness = _load_array(root, _fileref(selected, "witness"), (n_sites, 3))
        recomputed = hcat(q.qTT, q.qT, q.qPHI)
        all(abs.(witness .- recomputed) .≤ 8eps() .* abs.(recomputed)) ||
            _bundle_error("saved normalizer witness disagrees with the recomputation beyond roundoff")
    end

    shape = (length(ell), length(chi), length(Rax), n_sites)
    arrays = Dict{Tuple{String,String},Array{Float64,4}}()
    bank_rows = Tuple{String,String,Tuple{Int,Int,Int,Int},Int,String}[]
    for bank in selected["bank"]
        key = (bank["family"], bank["geometry"])
        ref = _fileref(bank, "file_ref")
        arrays[key] = _load_array(root, ref, shape)
        push!(bank_rows, (key[1], key[2], shape, ref.bytes, ref.sha))
    end
    sort!(bank_rows)
    grid = SplineGrid(sites, logsites, knots)
    fiducial = (; table_sha256 = table_ref.sha, ln10_10_As, n_s, pivot, kmin, kmax)
    digest = _bundle_signature(representation, n_sites, provenance, (; ell, chi, R = Rax), grid,
                               fiducial, source_rows, bank_rows)
    bankaxes=_copy_bank_axes((;ell,chi,R=Rax))
    signature = BankSignature(:bspline_cubic_notaknot, n_sites, digest,
        _geometry_signature(ell,chi,Rax))
    if representation === :direct
        # Raw typed constructor AFTER validation: the inner signature is the full v2 one.
        tensors = TensorBank(arrays[("direct", "T_2_00")], arrays[("direct", "T_0_00")], arrays[("direct", "T_minus2_00")],
                             arrays[("direct", "T_0_02")], arrays[("direct", "T_0_20")], arrays[("direct", "T_2_02")],
                             arrays[("direct", "T_2_20")], arrays[("direct", "T_2_22")], logsites, signature, bankaxes)
        return DirectSplineBank(tensors, grid)
    end
    TT = NamedTuple{_BANK_GEOMETRIES}(Tuple(arrays[("ratio_TT", String(g))] for g in _BANK_GEOMETRIES))
    T = NamedTuple{(:T_2_00, :T_0_00, :T_2_02, :T_2_20)}(Tuple(arrays[("ratio_T", g)] for g in _RATIO_T_GEOMETRIES))
    return RatioSplineBank(TT, T, arrays[("ratio_PHI", "T_2_00")], grid, q, signature, bankaxes)
end

# Phase 4A: cold loader of portable cubic B-spline tensor bundles (synthetic, public-shaped).
using Test, Blast, LinearAlgebra
using AbstractCosmologicalEmulators, FastChebInterp      # same imports as the full suite (name collisions)
isdefined(@__MODULE__, :write_synthetic_bundle) || include("tensor_bundle_fixture.jl")

function with_bundle(f; kwargs...)
    root = mktempdir()
    manifest = write_synthetic_bundle(root; kwargs...)
    return f(root, manifest)
end
# `why` (a regex) pins the reason so a case cannot pass by failing for an unrelated cause.
rejects(root; why=r"", kwargs...) = try; syn_load(root; kwargs...); false; catch e
    e isa ArgumentError && occursin(why, e.msg) || (@warn "unexpected rejection" e; false)
end
function mutated_rejects(mutate; why=r"", kwargs...)
    with_bundle() do root, manifest
        replacement = mutate(root, manifest)
        replacement isa AbstractDict && save_synthetic_manifest(root, replacement)
        rejects(root; why, kwargs...)
    end
end
function flip_value!(path, index)      # change one element in place, same size and header
    bytes = read(path); offset = length(bytes) - 8 * 4 * 5 * 9 * 3 + 8 * (index - 1)
    bytes[offset + 8] ⊻= 0x40; write(path, bytes)   # sign/exponent flip of one Float64
end

@testset "Tensor bundle: fixed fiducial formulas pinned to the generator expression" begin
    k = [1e-3, 2.5e-2, 0.3, 4.0, 10.0]
    As = exp(3.044) * 1e-10
    @test Blast._fiducial_phi(k, 3.044, 0.9649) == @. (9/25)*2*pi^2*As/k^3*(k/.05)^(0.9649-1)
end

@testset "Tensor bundle: supported semantic contract is pinned by literal machine IDs" begin
    @test Blast.BUNDLE_CONVENTION == (coordinate = "log10_k_per_mpc",
        basis = "ace_cubic_bspline_notaknot_degree3_extrapolation_throw",
        coefficient = "bspline_coefficients_collocation_solve", layout = "float64_ell_chi_r_coef_fortran",
        endpoint = "inward_quadrature_nodes_relative_epsilon", direct_source = "unweighted_shared_bank",
        ratio_source = "fixed_q_divide_full_physical_source_before_solve")
    @test Blast.BUNDLE_ENDPOINT_EPSILON == 1.0e-8 && Blast.BUNDLE_PIVOT == 0.05
    @test Blast.BUNDLE_FIDUCIAL.pmm_interpolator == "ln_ln_natural_cubic_spline_ace"
    @test Blast.BUNDLE_FIDUCIAL.q_T == "sqrt_phi_fid_times_pmm_fid" && Blast.BUNDLE_FIDUCIAL.q_PHI == "phi_fid" && Blast.BUNDLE_FIDUCIAL.q_TT == "pmm_fid"
    # The trusted production range is the legacy k_cheb range (constants.jl).
    @test isapprox(10 .^ extrema(Blast.k_cheb) |> collect, collect(Blast.PRODUCTION_K_RANGE); rtol=4eps())
end

@testset "Tensor bundle: loads direct and ratio banks (selected only, index-coded values)" begin
    with_bundle() do root, manifest
        direct = syn_load(root; rep=:direct)
        @test direct isa Blast.DirectSplineBank{Array{Float64,4}}
        @test Blast.bank_signature(direct) === direct.tensors.signature
        @test direct.tensors.signature.basis === :bspline_cubic_notaknot && direct.tensors.signature.n_k == 9
        for g in Blast._BANK_GEOMETRIES
            @test getfield(direct.tensors, g) == syn_tensor("direct", string(g), 9)     # no permutation / family swap
        end
        @test direct.tensors.sites == direct.grid.logsites
        @test isapprox(10 .^ direct.grid.logsites, direct.grid.sites; rtol=1e-14)
        @test length(direct.grid.knots) == 9 + 4

        ratio = syn_load(root; rep=:ratio)
        @test ratio isa Blast.RatioSplineBank{Array{Float64,4}}
        for g in Blast._BANK_GEOMETRIES; @test getproperty(ratio.TT, g) == syn_tensor("ratio_TT", string(g), 9); end
        for g in keys(ratio.T); @test getproperty(ratio.T, g) == syn_tensor("ratio_T", string(g), 9); end
        @test ratio.PHI == syn_tensor("ratio_PHI", "T_2_00", 9)
        pmm, qT, phi = syn_q(syn_table(), ratio.grid.sites)
        @test isapprox(ratio.q.qTT, pmm; rtol=8eps())
        @test isapprox(ratio.q.qT, qT; rtol=8eps())
        @test isapprox(ratio.q.qPHI, phi; rtol=8eps())
        @test ratio.q.qPHI == phi                                                       # exact generator formula
        @test Blast.bank_signature(direct) != Blast.bank_signature(ratio)
        @test Blast.bank_signature(syn_load(root; n=11, rep=:direct)) != Blast.bank_signature(direct)
        @test_throws ArgumentError syn_load(root; n=13)                                  # count absent
        @test_throws ArgumentError Blast.load_tensor_bank(joinpath(root, "bundle.toml"); count=9, representation=:other)
    end
end

@testset "Tensor bundle: only the selected count is opened" begin
    with_bundle() do root, manifest
        for e in manifest["count"], b in e["bank"]
            e["count"] == 11 && rm(joinpath(root, b["file_ref"]["file"]))               # other count's payloads vanish
        end
        for e in manifest["count"]; e["count"] == 11 && rm(joinpath(root, e["knots"]["file"])); end
        @test syn_load(root; n=9, rep=:ratio) isa Blast.RatioSplineBank
        @test rejects(root; n=11, rep=:ratio)
    end
end

@testset "Tensor bundle: signature identity and relocation" begin
    with_bundle() do root, manifest
        sig = Blast.bank_signature(syn_load(root; rep=:ratio))
        @test sig == Blast.bank_signature(syn_load(root; rep=:ratio))
        moved = mktempdir(); cp(root, joinpath(moved, "elsewhere"))                     # relocation as a unit
        @test Blast.bank_signature(syn_load(joinpath(moved, "elsewhere"); rep=:ratio)) == sig
        # Excluded from the identity: status, timestamps and the manifest text itself.
        m = load_synthetic_manifest(root); m["status"] = "certified"; m["created_utc"] = "2030-12-31T00:00:00"
        save_synthetic_manifest(root, m)
        @test Blast.bank_signature(syn_load(root; rep=:ratio)) == sig
        # Legitimate identity changes (all still valid, supported data): they change the signature.
        for change in (m -> m["provenance"]["generator_sha256"] = syn_sha_text("another generator"),
                       m -> m["provenance"]["quadrature_N"] = 262145,
                       m -> m["provenance"]["runtime_ace_rev"] = "another-revision")
            m2 = load_synthetic_manifest(root); change(m2); save_synthetic_manifest(root, m2)
            @test Blast.bank_signature(syn_load(root; rep=:ratio)) != sig
            save_synthetic_manifest(root, m)
        end
        # Fixed fiducial parameters within physical ranges (direct bank: no witness to refresh).
        dsig = Blast.bank_signature(syn_load(root; rep=:direct))
        m2 = load_synthetic_manifest(root); m2["fiducial"]["ln10_10_As"] = 3.1; save_synthetic_manifest(root, m2)
        @test Blast.bank_signature(syn_load(root; rep=:direct)) != dsig
        save_synthetic_manifest(root, m)
        # Dictionary insertion/print order does not matter (canonical typed field order).
        for by in (identity, reverse, s -> hash(s))
            open(io -> TOML.print(io, m; sorted=true, by=by), joinpath(root, "bundle.toml"), "w")
            @test Blast.bank_signature(syn_load(root; rep=:ratio)) == sig
        end
        save_synthetic_manifest(root, m)
        # A same-sized payload with different (self-consistent) content passes integrity but changes identity.
        m3 = load_synthetic_manifest(root)
        ref = syn_bank(m3, 9, "ratio", "ratio_TT", "T_0_20")["file_ref"]
        path = joinpath(root, ref["file"]); bytes = read(path); bytes[end] ⊻= 0x01; write(path, bytes)
        syn_refresh!(root, ref); save_synthetic_manifest(root, m3)
        @test Blast.bank_signature(syn_load(root; rep=:ratio)) != sig
    end
end

@testset "Tensor bundle: rejections (PLAN 6.3 and same-size metadata/payload mismatches)" begin
    @test mutated_rejects((r, m) -> (m["schema"] = "blast.tensor_bank_bundle/2"; m), why=r"schema")
    @test mutated_rejects((r, m) -> (m["status"] = "gold"; m), why=r"status")
    @test mutated_rejects((r, m) -> (delete!(m, "source_map"); m))
    @test mutated_rejects((r, m) -> (delete!(m, "fiducial"); m))
    @test mutated_rejects((r, m) -> (e = syn_entry(m, 9, "direct"); pop!(e["bank"]); m), why=r"exactly its")      # missing family
    @test mutated_rejects((r, m) -> (e = syn_entry(m, 9, "ratio"); filter!(b -> b["family"] != "ratio_T", e["bank"]); m), rep=:ratio)
    @test mutated_rejects((r, m) -> (e = syn_entry(m, 9, "direct"); push!(e["bank"], copy(e["bank"][1])); m), why=r"duplicate bank")   # duplicate ban
    @test mutated_rejects((r, m) -> (push!(m["count"], copy(syn_entry(m, 9, "direct"))); m), why=r"duplicate entry")   # duplicate count entr
    @test mutated_rejects((r, m) -> (push!(m["source_map"], copy(m["source_map"][1])); m), why=r"duplicate source_map")   # duplicate componen
    @test mutated_rejects((r, m) -> (pop!(m["source_map"]); m), why=r"exactly the 17")   # missing componen
    @test mutated_rejects((r, m) -> (m["source_map"][1]["family_direct"] = ["direct", "T_0_00"]; m), why=r"wrong direct bank")   # wrong mappin
    # files
    @test mutated_rejects((r, m) -> rm(joinpath(r, syn_bank(m, 9, "direct", "direct", "T_2_00")["file_ref"]["file"])), why=r"missing file")
    @test mutated_rejects((r, m) -> (p = joinpath(r, syn_bank(m, 9, "direct", "direct", "T_2_00")["file_ref"]["file"]);
                                     write(p, read(p)[1:end-8])), why=r"byte count")   # truncated, manifest unchange
    @test mutated_rejects((r, m) -> (syn_bank(m, 9, "direct", "direct", "T_2_00")["file_ref"]["bytes"] += 8; m))
    @test mutated_rejects((r, m) -> (syn_bank(m, 9, "direct", "direct", "T_2_00")["file_ref"]["header_bytes"] += 64; m))
    @test mutated_rejects((r, m) -> flip_value!(joinpath(r, syn_bank(m, 9, "direct", "direct", "T_2_00")["file_ref"]["file"]), 5), why=r"SHA-256")   # sh
    @test mutated_rejects((r, m) -> (m["axes"]["R"]["sha256"] = "0"^64; m))
    @test mutated_rejects((r, m) -> (m["fiducial"]["table"]["sha256"] = "0"^64; m))
    # paths
    @test mutated_rejects((r, m) -> (syn_bank(m, 9, "direct", "direct", "T_2_00")["file_ref"]["file"] = "/etc/passwd"; m), why=r"absolute")
    @test mutated_rejects((r, m) -> (syn_bank(m, 9, "direct", "direct", "T_2_00")["file_ref"]["file"] = "../x/direct/T_2_00.npy"; m), why=r"escapes")
    @test mutated_rejects((r, m) -> (syn_bank(m, 9, "direct", "direct", "T_2_00")["file_ref"]["file"] = "N009/./direct/T_2_00.npy"; m))
    @test mutated_rejects((r, m) -> (m["axes"]["ell"]["file"] = "meta\\axes_ell.npy"; m))
    @test mutated_rejects((r, m) -> (a = syn_bank(m, 9, "direct", "direct", "T_2_00")["file_ref"];   # swapped same-size families
                                     b = syn_bank(m, 9, "direct", "direct", "T_0_00")["file_ref"];
                                     for k in ("file", "sha256", "bytes", "header_bytes"); a[k], b[k] = b[k], a[k]; end; m), why=r"non-canonical path")
    # symlink component (skipped silently where unsupported)
    with_bundle() do root, m
        try
            real = joinpath(root, "N009", "direct_real"); mv(joinpath(root, "N009", "direct"), real)
            symlink(real, joinpath(root, "N009", "direct"))
            @test rejects(root)
        catch err
            err isa Base.IOError || rethrow()
        end
    end
    # array content, dtype, rank, shape, order
    redo(r, m, rel_shape, writer) = begin
        b = syn_bank(m, 9, "direct", "direct", "T_2_00")["file_ref"]
        writer(joinpath(r, b["file"])); syn_refresh!(r, b); m
    end
    n = 3 * 4 * 5 * 9
    @test mutated_rejects((r, m) -> redo(r, m, n, p -> syn_write_npy(p, syn_tensor("direct", "T_2_00", 9); descr="<f4")))
    @test mutated_rejects((r, m) -> redo(r, m, n, p -> syn_write_npy(p, syn_tensor("direct", "T_2_00", 9)[:, :, :, 1])))
    @test mutated_rejects((r, m) -> redo(r, m, n, p -> syn_write_npy(p, syn_tensor("direct", "T_2_00", 10))))
    @test mutated_rejects((r, m) -> redo(r, m, n, p -> syn_write_npy(p, syn_tensor("direct", "T_2_00", 9); fortran=false)))
    @test mutated_rejects((r, m) -> redo(r, m, n, p -> syn_write_npy(p, permutedims(syn_tensor("direct", "T_2_00", 9), (2, 1, 3, 4)))))
    @test mutated_rejects((r, m) -> redo(r, m, n, p -> (A = syn_tensor("direct", "T_2_00", 9); A[2, 2, 2, 2] = NaN; syn_write_npy(p, A))))
    @test mutated_rejects((r, m) -> redo(r, m, n, p -> (A = syn_tensor("direct", "T_2_00", 9); A[1] = Inf; syn_write_npy(p, A))))
    # axes
    wrong_axes = (; SYN_AXES..., chi = [1.0, 2.5, 4.0, 7.5])
    @test mutated_rejects((r, m) -> nothing; axes=wrong_axes)
    @test mutated_rejects((r, m) -> nothing; axes=(; SYN_AXES..., R = [0.2, 0.5, 0.8, 0.95]))           # 4-R "legacy" shape
    @test mutated_rejects((r, m) -> nothing; axes=Blast.production_axes())                              # synthetic bank is not production
    @test mutated_rejects((r, m) -> (syn_write_npy(joinpath(r, "meta", "axes_R.npy"), [0.2, 0.5, 0.8, 0.95, 0.99]);
                                     syn_refresh!(r, m["axes"]["R"]); m))                                # no R = 1
    # sites, knots, fiducial
    function edit_meta(key, f; rep="direct")
        (r, m) -> begin
            e = syn_entry(m, 9, rep); p = joinpath(r, e[key]["file"]); a = Blast.npzread(p); f(a)
            syn_write_npy(p, a); syn_refresh!(r, e[key]); m
        end
    end
    @test mutated_rejects(edit_meta("sites_phys", a -> (a[4], a[5]) = (a[5], a[4])), why=r"strictly increasing")   # not increasing
    @test mutated_rejects(edit_meta("sites_phys", a -> a[3] = NaN))
    @test mutated_rejects(edit_meta("sites_phys", a -> a[1] = 1.2e-3), why=r"endpoints")                              # wrong endpoint
    @test mutated_rejects(edit_meta("sites_log10", a -> a[4] += 1e-9), why=r"disagree beyond roundoff")             # roundoff-only violated
    @test mutated_rejects(edit_meta("knots", a -> a[6] = nextfloat(a[6])), why=r"knots differ")                     # one-ulp knot change
    @test mutated_rejects(edit_meta("knots", a -> a[1] = prevfloat(a[1])), why=r"knots differ")
    @test mutated_rejects(edit_meta("witness", a -> a[2, 2] *= 1 + 1e-9; rep="ratio"), rep=:ratio, why=r"witness")  # witness disagrees
    function edit_table(f)
        (r, m) -> begin
            p = joinpath(r, "meta", "fiducial_pmm.npy"); a = Blast.npzread(p); f(a)
            syn_write_npy(p, a); syn_refresh!(r, m["fiducial"]["table"]); m
        end
    end
    @test mutated_rejects(edit_table(a -> a[10, 2] = -1.0), rep=:ratio, why=r"positive")
    @test mutated_rejects(edit_table(a -> a[10, 2] = 0.0), rep=:ratio)
    @test mutated_rejects(edit_table(a -> a[10, 2] = NaN), rep=:ratio)
    @test mutated_rejects(edit_table(a -> (a[5, 1], a[6, 1]) = (a[6, 1], a[5, 1])), rep=:ratio)
    @test mutated_rejects((r, m) -> (p = joinpath(r, "meta", "fiducial_pmm.npy"); a = Blast.npzread(p)[1:40, :];  # table too short
                                     syn_write_npy(p, a); syn_refresh!(r, m["fiducial"]["table"]); m["fiducial"]["rows"] = 40; m), rep=:ratio)
    @test mutated_rejects((r, m) -> (m["fiducial"]["pivot"] = 0.06; m), rep=:ratio)
    @test mutated_rejects((r, m) -> (syn_entry(m, 9, "direct")["count"] = 10; m))                      # count/sites mismatch
    # duplicate key in TOML text
    with_bundle() do root, m
        path = joinpath(root, "bundle.toml"); write(path, replace(read(path, String), "status = \"candidate\"" => "status = \"candidate\"\nstatus = \"candidate\""))
        @test rejects(root)
    end
end

@testset "Tensor bundle: unsupported semantics, wrong domain and overflow are rejected (not hashed)" begin
    # Parent reproduction: a manifest claiming another coordinate and a zero formula must NOT load.
    @test mutated_rejects((r, m) -> (m["convention"]["coordinate"] = "ln(k)"; m["fiducial"]["phi_formula"] = "zero everywhere"; m),
                          rep=:ratio, why=r"unsupported convention.coordinate")
    for key in ("coordinate", "basis", "coefficient", "layout", "endpoint", "direct_source", "ratio_source")
        @test mutated_rejects((r, m) -> (m["convention"][key] = "something_else"; m), why=Regex("unsupported convention.$key"))
        @test mutated_rejects((r, m) -> (m["convention"][key] = 1; m), why=Regex("`$key`"))                     # wrong type
    end
    @test mutated_rejects((r, m) -> (m["convention"]["endpoint_epsilon"] = 1.0e-7; m), why=r"endpoint_epsilon")
    @test mutated_rejects((r, m) -> (m["convention"]["endpoint_epsilon"] = "1e-8"; m), why=r"endpoint_epsilon")
    for key in ("pmm_interpolator", "phi_formula", "q_TT", "q_T", "q_PHI")
        @test mutated_rejects((r, m) -> (m["fiducial"][key] = "something_else"; m), why=Regex("unsupported fiducial.$key"))
    end
    @test mutated_rejects((r, m) -> (m["fiducial"]["phi_formula"] = "zero everywhere"; m), rep=:ratio, why=r"unsupported fiducial.phi_formula")
    # Unknown / nested / ambiguous fields are not silently accepted.
    @test mutated_rejects((r, m) -> (m["extra"] = 1; m), why=r"unsupported key `extra`")
    @test mutated_rejects((r, m) -> (m["convention"]["extra"] = "x"; m), why=r"unsupported key")
    @test mutated_rejects((r, m) -> (m["provenance"]["extra"] = Dict("a" => 1); m), why=r"unsupported key")
    @test mutated_rejects((r, m) -> (m["fiducial"]["extra"] = 1; m), why=r"unsupported key")
    @test mutated_rejects((r, m) -> (m["range"]["extra"] = 1.0; m), why=r"unsupported key")
    @test mutated_rejects((r, m) -> (syn_entry(m, 9, "direct")["extra"] = 1; m), why=r"unsupported key")
    @test mutated_rejects((r, m) -> (syn_bank(m, 9, "direct", "direct", "T_2_00")["extra"] = 1; m), why=r"unsupported key")
    @test mutated_rejects((r, m) -> (syn_bank(m, 9, "direct", "direct", "T_2_00")["file_ref"]["extra"] = 1; m), why=r"unsupported key")
    @test mutated_rejects((r, m) -> (m["provenance"]["quadrature_N"] = "131073"; m), why=r"quadrature_N")     # "1" vs 1
    @test mutated_rejects((r, m) -> (m["provenance"]["quadrature_N"] = true; m), why=r"quadrature_N")
    @test mutated_rejects((r, m) -> (m["provenance"]["generator_sha256"] = "g"^64; m), why=r"64 lowercase hex")
    @test mutated_rejects((r, m) -> (m["provenance"]["inputs_manifest_sha256"] = uppercase(syn_sha_text("x")); m), why=r"64 lowercase hex")
    # Self-consistent wrong domain: bounds, sites, knots and fiducial all agree on another range.
    let root = mktempdir()
        write_synthetic_bundle(root; counts=(9,), kmin=2.0 * SYN_KMIN, kmax=SYN_KMAX)
        @test rejects(root; why=r"declared k range")
        @test syn_load(root; expected_k_range=(2.0 * SYN_KMIN, SYN_KMAX)) isa Blast.DirectSplineBank     # only with explicit trust
    end
    let root = mktempdir()
        write_synthetic_bundle(root; counts=(9,))
        @test rejects(root; expected_k_range=Blast.PRODUCTION_K_RANGE, why=r"declared k range")        # synthetic range != production
        m = load_synthetic_manifest(root); m["range"]["kmin"] = 2.0 * SYN_KMIN; save_synthetic_manifest(root, m)
        @test rejects(root; why=r"declared k range")
    end
    # Normalizers: finite positive factors can still overflow/underflow in q_T = sqrt(phi*pmm).
    @test_throws ArgumentError Blast._normalizers([1e200, 1.0], [1e200, 1.0])            # overflow -> Inf
    @test_throws ArgumentError Blast._normalizers([1e-200, 1.0], [1e-200, 1.0])          # underflow -> 0
    @test_throws ArgumentError Blast._normalizers([1.0, 2.0], [1.0, -1.0])
    @test_throws ArgumentError Blast._normalizers([1.0, NaN], [1.0, 1.0])
    @test Blast._normalizers([1.0, 4.0], [1.0, 9.0]).qT == [1.0, 6.0]
    # Sanity bounds before any allocation.
    @test mutated_rejects((r, m) -> (syn_entry(m, 9, "direct")["count"] = 3; m), why=r"out of range")
    @test mutated_rejects((r, m) -> (syn_entry(m, 9, "direct")["count"] = 10^12; m), why=r"out of range")
    @test mutated_rejects((r, m) -> (syn_entry(m, 9, "direct")["count"] = true; m), why=r"non-integer")
    @test mutated_rejects((r, m) -> (m["fiducial"]["rows"] = 2; m), why=r"out of range")
    @test mutated_rejects((r, m) -> (m["fiducial"]["rows"] = 10^15; m), why=r"out of range")
    @test mutated_rejects((r, m) -> (m["axes"]["ell"]["bytes"] = 0; m), why=r"out of range")
    @test mutated_rejects((r, m) -> (delete!(m["axes"]["R"], "header_bytes"); m), why=r"header_bytes")           # mandatory for meta NPY
    @test mutated_rejects((r, m) -> (syn_entry(m, 9, "direct")["knots"]["header_bytes"] = 1; m), why=r"out of range")
    @test mutated_rejects((r, m) -> (syn_entry(m, 9, "direct")["knots"]["header_bytes"] += 64; m), why=r"header byte count")
end

@testset "Tensor bundle: loading leaves the lazy legacy bank untouched (fresh process)" begin
    root = mktempdir(); write_synthetic_bundle(root; counts=(9,))
    code = """
        using Blast
        include($(repr(joinpath(@__DIR__, "tensor_bundle_fixture.jl"))))
        bank = syn_load($(repr(root)); rep=:ratio)
        exit(isassigned(Blast._LEGACY_BANK) ? 1 : 0)
    """
    @test success(`$(Base.julia_cmd()) --project=$(Base.active_project()) -e $code`)
end

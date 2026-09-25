# Blast side of the ARFCAMB validation. Shared by test/test_arfcamb_validation.jl
# and compare_blast_camb.ipynb so both compute exactly the same thing.
#
# Requires `Blast` and `TOML` to be loaded by the caller.

"""
    blast_arf_pieces(cfg, fx; matter=:cb) -> Dict{String, Vector{Float64}}

Compute Blast's C_ℓ for the single-bin ARFCAMB validation. `cfg` is the parsed
`config.toml`, `fx` the loaded fixture (`npzread` of `gaussian_z1.npz`), which
also provides CAMB's linear P(k, z) on Blast's grids, so only kernels and
projections are compared. Keys match the fixture: dd, dd_rsd, dd_mu, gg, kd,
kd_rsd, kd_mu, kg (counts), zz, dz, kz (density-only ARF).

`matter` picks which CAMB P(k) Blast uses: `:cb` (CDM+baryons, the default) or
`:tot` (total matter). CAMB sources galaxy counts with the CDM density but CMB
lensing with the total one, while Blast takes a single P(k); with mν = 0.06 eV
`:tot` biases galaxy autos low by ~1%, `:cb` leaves ~0.1% on κ×g crosses.
"""
function blast_arf_pieces(cfg, fx; matter::Symbol=:cb)
    c = cfg["cosmology"]
    cosmo = Blast.w0waCDMCosmology(ln10Aₛ=log(1e10 * c["As"]), nₛ=c["ns"], h=c["h"],
        ωb=c["ombh2"], ωc=c["omch2"], ωk=0.0, mν=c["mnu"], w0=c["w0"], wa=c["wa"])
    bg = Blast.Background(cosmo)

    z = fx["z"]
    nz = reshape(fx["nz"], 1, :)
    n_bg = length(bg.z)
    bias = fill(cfg["tracer"]["bias"], 1, n_bg)
    zero_bias = zeros(1, n_bg)
    s = fill(cfg["tracer"]["s"], 1, n_bg)
    zbar = fx["zbar"]

    δ   = Blast.NumberCounts(nz=nz, z=z, bias=bias)
    δ0  = Blast.NumberCounts(nz=nz, z=z, bias=zero_bias)  # zero kernel: isolates other terms
    rsd = Blast.RedshiftSpaceDistortions(nz=nz, z=z, growth_rate=bg.f)
    μ   = Blast.MagnificationBias(nz=nz, z=z, s=s)
    arf() = Blast.AngularRedshiftFluctuations(nz=nz, z=z, bias=bias, zbar=zbar)

    probes = Dict(
        "d"      => Blast.GalaxyClustering(δ=δ),
        "d_rsd"  => Blast.GalaxyClustering(δ=δ, RSD=rsd),
        "d_mu"   => Blast.GalaxyClustering(δ=δ, μ=μ),
        "g"      => Blast.GalaxyClustering(δ=δ, RSD=rsd, μ=μ),
        "z"      => Blast.GalaxyClustering(δ=δ0, ARF=arf()),
        "d_z"    => Blast.GalaxyClustering(δ=δ, ARF=arf()),
    )
    foreach(GC -> Blast.evaluate_components!(GC, bg), values(probes))

    CMB = Blast.CMB(κ=Blast.CMBLensing())
    Blast.evaluate_components!(CMB, bg)

    # The full GC (δ+RSD+μ+ARF) allocates every weight tensor any subset needs.
    W, Plans = Blast.SetUp(Blast.GalaxyClustering(δ=δ, RSD=rsd, μ=μ, ARF=arf()), CMB)
    pk, pk_limber = matter === :cb  ? (fx["pk_cb"], fx["pk_limber_cb"]) :
                    matter === :tot ? (fx["pk"], fx["pk_limber_lin"]) :
                    throw(ArgumentError("matter must be :cb or :tot, got :$matter"))
    PS = Blast.prepare_pk_workspace(Plans, pk, pk_limber, pk_limber, bg)
    W = Blast.allocate_compute_w(W, PS)
    Blast.compute_w!(W, PS)

    ℓ = Float64.(fx["ell"])
    auto(key) = Blast.get_Cℓ(ℓ, probes[key], PS, W, bg, Plans)[:, 1, 1]
    cross(key) = Blast.get_Cℓ(ℓ, CMB, probes[key], PS, W, bg, Plans)[:, 1, 1]

    out = Dict{String, Vector{Float64}}(
        "dd" => auto("d"), "dd_rsd" => auto("d_rsd"), "dd_mu" => auto("d_mu"), "gg" => auto("g"),
        "zz" => auto("z"),
        "kd" => cross("d"), "kd_rsd" => cross("d_rsd"), "kd_mu" => cross("d_mu"), "kg" => cross("g"),
        "kz" => cross("z"),
    )
    # GalaxyClustering sums its components into one field: C[δ+ARF] = dd + 2 dz + zz.
    out["dz"] = (auto("d_z") .- out["dd"] .- out["zz"]) ./ 2
    return out
end

"""
    normalized_difference(blast, camb, auto_a, auto_b)

(blast - camb) / sqrt(auto_a * auto_b): a fractional difference that stays
well defined for cross-spectra passing through zero. For an auto-spectrum pass
the same spectrum twice to get the plain relative difference.
"""
normalized_difference(blast, camb, auto_a, auto_b) = (blast .- camb) ./ sqrt.(abs.(auto_a .* auto_b))

# Which CAMB autos normalize each piece's difference (κκ comes from the fixture).
const ARFCAMB_NORMS = Dict(
    "dd" => ("dd", "dd"), "dd_rsd" => ("dd_rsd", "dd_rsd"), "dd_mu" => ("dd_mu", "dd_mu"), "gg" => ("gg", "gg"),
    "zz" => ("zz", "zz"), "dz" => ("dd", "zz"),
    "kd" => ("kk", "dd"), "kd_rsd" => ("kk", "dd_rsd"), "kd_mu" => ("kk", "dd_mu"), "kg" => ("kk", "gg"),
    "kz" => ("kk", "zz"),
)

# Blast side of the fNL validation (issue #38). Computes the galaxy-clustering
# spectra of config.toml's bins at f_NL = 0, +f_NL, -f_NL on CAMB's own linear
# P_cb(k, z) from camb_png.npz, and writes them to the path given as ARGS[1].
#
# Runs against whichever Blast the active project provides, so the same script
# gives the fixed and pre-fix numbers (see README.md):
#
#   julia --project=. validation/png_external/blast_png.jl validation/png_external/blast_fixed.npz
#
# Also writes swift_inputs.npz (background and P_Φ on the shared n(z) table) when
# ARGS[2] == "--swift-inputs", so SwiftCℓ sees Blast's χ(z), H(z) and Φ convention.

using Blast, TOML
import NPZ

cfg = TOML.parsefile(joinpath(@__DIR__, "config.toml"))
fx = NPZ.npzread(joinpath(@__DIR__, "camb_png.npz"))
out_path = ARGS[1]

c = cfg["cosmology"]
cosmo = Blast.w0waCDMCosmology(ln10Aₛ=log(1e10 * c["As"]), nₛ=c["ns"], h=c["h"],
    ωb=c["ombh2"], ωc=c["omch2"], ωk=0.0, mν=c["mnu"], w0=c["w0"], wa=c["wa"])
bg = Blast.Background(cosmo)

z, nz = fx["z"], fx["nz"]
nbins, n_bg = size(nz, 1), length(bg.z)
tr = cfg["tracer"]
bias = repeat([b["bias"] for b in cfg["bins"]], 1, n_bg)
f_NL = tr["f_NL"]

δ = Blast.NumberCounts(nz=nz, z=z, bias=bias)
png(f) = Blast.PrimordialNonGaussianity(nz=nz, z=z, bias=bias, f_NL=f, p=tr["p"])
probes = (c0 = Blast.GalaxyClustering(δ=δ),
          cp = Blast.GalaxyClustering(δ=δ, PNG=png(f_NL)),
          cm = Blast.GalaxyClustering(δ=δ, PNG=png(-f_NL)))
foreach(GC -> Blast.evaluate_components!(GC, bg), probes)

W, Plans = Blast.SetUp(probes.cp)
PS = Blast.prepare_pk_workspace(Plans, fx["pk_cb"], fx["pk_limber_cb"], fx["pk_limber_cb"], bg)
W = Blast.compute_w(W, PS)

ℓ = Float64.(fx["ell"])
spectra = Dict(String(k) => Blast.get_Cℓ(ℓ, GC, PS, W, bg, Plans) for (k, GC) in pairs(probes))
NPZ.npzwrite(out_path, merge(spectra, Dict("ell" => ℓ)))
println("Wrote ", out_path)

if get(ARGS, 2, "") == "--swift-inputs"
    χ_of_z = Blast.akima_interpolation(bg.χ, bg.z, z)
    H_of_z = Blast.akima_interpolation(bg.H, bg.z, z)
    swift_path = joinpath(@__DIR__, "swift_inputs.npz")
    NPZ.npzwrite(swift_path, Dict(
        "z" => z, "chi" => χ_of_z, "H" => H_of_z,     # H in km/s/Mpc
        "k" => fx["k_swift"], "P_phi" => Blast.get_PΦ(fx["k_swift"], bg)))
    println("Wrote ", swift_path)
end

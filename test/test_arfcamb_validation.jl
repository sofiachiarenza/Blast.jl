using Test
using Blast
using NPZ
using TOML

# External validation against ARFCAMB for a single Gaussian bin (z0 = 1, σ = 0.2).
# The fixture (CAMB C_ℓ + CAMB's linear P(k) on Blast's grids) is produced by
# validation/arfcamb/generate_camb_reference.py; see validation/arfcamb/README.md.
include(joinpath(@__DIR__, "..", "validation", "arfcamb", "blast_pieces.jl"))

@testset "ARF vs ARFCAMB (single Gaussian bin)" begin
    cfg = TOML.parsefile(joinpath(@__DIR__, "..", "validation", "arfcamb", "config.toml"))
    fx = npzread(joinpath(@__DIR__, "reference_fixtures", "arfcamb", "gaussian_z1.npz"))
    blast = blast_arf_pieces(cfg, fx)

    # The reference comes from ARFCAMB patched to drop the Newtonian-gauge
    # correction (comoving_density_ev - 3 adotoa) sigma/k from the counts
    # density source, which Blast's NumberCounts kernel does not model (8% on δδ
    # at ℓ = 2 in stock CAMB; kept in the fixture as `*_stock`). The tolerance
    # bounds (Blast - CAMB)/sqrt(C_aa C_bb) over all ℓ; the observed maximum is 0.37%.
    tol = 0.006

    # Density-only ARF: Blast's ARF kernel is (z - zbar) K^δ, so CAMB's
    # ARF-weighted RSD/magnification terms (zz_full etc.) are not comparable.
    for key in ("dd", "dd_rsd", "dd_mu", "gg", "kd", "kd_rsd", "kd_mu", "kg", "zz", "dz", "kz")
        a, b = ARFCAMB_NORMS[key]
        Δ = normalized_difference(blast[key], fx[key], fx[a], fx[b])
        @test maximum(abs, Δ) < tol
    end

    # ARFCAMB's zbar = Σ z W / Σ W on the table grid; for a symmetric Gaussian ≈ z0.
    @test fx["zbar"][1] ≈ cfg["window"]["z0"] atol = 1e-5
end

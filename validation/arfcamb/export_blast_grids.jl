# Dump Blast's fixed P(k) grids so generate_camb_reference.py can evaluate
# CAMB's P(k, z) exactly where Blast expects it.
#
#   julia --project=validation/arfcamb validation/arfcamb/export_blast_grids.jl

using Blast
using NPZ

npzwrite(joinpath(@__DIR__, "grids.npz"),
    z_lin = collect(Blast.z_lin),
    log10k_cheb = collect(Blast.k_cheb),      # log10(k / Mpc^-1)
    z_cheb = collect(Blast.z_cheb),
    log10k_limber = collect(Blast.k_limber))  # log10(k / Mpc^-1)

println("Wrote ", joinpath(@__DIR__, "grids.npz"))

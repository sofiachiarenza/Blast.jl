# GLOBAL CONSTANTS FOR GRIDS AND INTEGRATION
# These values are fixed for the Blast pipeline to ensure consistency across probes.

# ── Physical / unit constants ─────────────────────────────────────────────────
# Conversion factor H0 [km/s/Mpc] = H_0_CONV * h
const H_0_CONV = 100.0

# ── Background redshift grid parameters ───────────────────────────────────────
# Number of fine redshift points used to build the z(χ) inversion spline.
# N=100 is deliberately chosen: akima convergence on bg.H/D/f is ~4e-6 relative
# vs N=1000 (well below cosmic variance), while Background() is ~7.6x faster
# (4.44 ms → 0.58 ms). The gain flows through to the Mooncake tape: T5 MC
# gradient drops proportionally because the r_z-broadcast over fine_z is the
# dominant cost inside Background.  See benchmark/bg_fine_grid_sweep.jl for
# the full accuracy-vs-speed table.
const N_BG_FINE_GRID = 100
# Maximum redshift considered in the background spline.
const Z_MAX_BACKGROUND = 5.0


const full_ℓ_range = reverse(chebpoints(100, 2.0, 2000.0))
const ℓ_nonlimber = full_ℓ_range[full_ℓ_range .< 220]
const ℓ_limber = full_ℓ_range[full_ℓ_range .> 220]

const nχ = 96   # must match T_tildes artifact χ-axis dimension
const χ = Array(LinRange(26.0, 7000.0, nχ))

const _R_nodes = chebpoints(64*2, -1.0, 1.0)
const R = reverse(_R_nodes[_R_nodes .> 0])

const k_cheb = chebpoints(160, Float64(log10(5e-5)), Float64(log10(16)))
const k_limber = chebpoints(256, Float64(log10(1e-4)), Float64(log10(80)))
const z_cheb = chebpoints(49, 0.0, 3.6)
const z_lin = collect(LinRange(0.0, 3.6, 50))

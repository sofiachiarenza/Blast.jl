# GLOBAL CONSTANTS FOR GRIDS AND INTEGRATION
# These values are fixed for the Blast pipeline to ensure consistency across probes.

const full_ℓ_range = reverse(chebpoints(100, 2.0, 2000.0))
const ℓ_nonlimber = full_ℓ_range[full_ℓ_range .< 220]
const ℓ_limber = full_ℓ_range[full_ℓ_range .> 220]

const nχ = 96
const χ = Array(LinRange(26.0, 7000.0, nχ))

const _R_nodes = chebpoints(64*2, -1.0, 1.0)
const R = reverse(_R_nodes[_R_nodes .> 0])

const k_cheb = chebpoints(160, Float64(log10(5e-5)), Float64(log10(16)))
const k_limber = chebpoints(256, Float64(log10(1e-4)), Float64(log10(80)))
const z_cheb = chebpoints(49, 0.0, 3.6)
const z_lin = LinRange(0.0, 3.6, 50)

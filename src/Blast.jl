module Blast

using LoopVectorization
using Tullio
using FastTransforms
using FastChebInterp
using SpecialFunctions
using DataInterpolations
using Interpolations
using StaticArrays
using FFTW
using NPZ
using QuadGK
using Artifacts
using PhysicalConstants
using FastPower

include("cosmo.jl")
include("deprecated.jl")
include("probes.jl")
include("chebcoefs.jl")
include("setup.jl")
include("projected_matter.jl")
include("integrals.jl")
include("cls.jl")
include("utils.jl")
include("limber.jl")

import PhysicalConstants.CODATA2018: c_0

const C_LIGHT = c_0.val * 10^(-3) #speed of light in Km/s


struct T̃{A<:AbstractArray{<:Any,4}}
    T_2_00::A
    T_0_00::A
    T_minus2_00::A
    T_0_02::A
    T_0_20::A
    T_2_02::A
    T_2_20::A
    T_2_22::A
end

function __init__()

    Tdir = artifact"T_tildes"

    global T_tildes = T̃(
        npzread(joinpath(Tdir, "T_2_00.npz")),
        npzread(joinpath(Tdir, "T_0_00.npz")),
        npzread(joinpath(Tdir, "T_minus2_00.npz")),
        npzread(joinpath(Tdir, "T_0_02.npz")),
        npzread(joinpath(Tdir, "T_0_20.npz")),
        npzread(joinpath(Tdir, "T_2_02.npz")),
        npzread(joinpath(Tdir, "T_2_20.npz")),
        npzread(joinpath(Tdir, "T_2_22.npz"))
    )


    global full_ℓ_range = reverse(chebpoints(100, 2, 2000))
    global ℓ_nonlimber = full_ℓ_range[full_ℓ_range .< 220]
    global ℓ_limber = full_ℓ_range[full_ℓ_range .> 220]

    nχ = 96
    nR = 64
    global χ = Array(LinRange(26, 7000, nχ))
    R = chebpoints(nR*2, -1, 1)
    global R = reverse(R[R.>0])
    nR = length(R)

    kmin = 5e-5
    kmax = 16
    n_cheb = 160
    global k_cheb = chebpoints(n_cheb, log10(kmin), log10(kmax))

    global k_limber = chebpoints(256, log10(1e-4), log10(80))
    global z_cheb = chebpoints(49, 0, 3.6)
    global z_lin = LinRange(0,3.6, 50)

end

end # module Blast

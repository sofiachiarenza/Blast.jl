using AbstractCosmologicalEmulators
using Blast
using Reactant
using Test
using NPZ

const ROOT = "/home/marcobonici/Desktop/work/CosmologicalEmulators/blast_play/for_marco_copy"
include(joinpath(ROOT, "performance", "workload.jl"))

function main()
    Reactant.set_default_backend("cpu")
    workload = build_workload(ROOT, 1)
    bg, lp = workload.bg, workload.plans.plan_limber
    loglin, lognon = Blast.reactant_limber_power_products(
        workload.pk_limber_lin, workload.pk_limber_nonlin, bg)
    Kk, Kz = lp.K
    Mk = Blast.reactant_chebyshev_matrix(prepare_chebyshev_plan(
        minimum(lp.nodes[1]), maximum(lp.nodes[1]), Kk;
        size_nd=(Kk + 1,), dim=1))
    Mz = Blast.reactant_chebyshev_matrix(prepare_chebyshev_plan(
        minimum(lp.nodes[2]), maximum(lp.nodes[2]), Kz;
        size_nd=(Kz + 1,), dim=1))
    coeff_ref = (chebyshev_decomposition(lp, loglin),
                 chebyshev_decomposition(lp, lognon))
    coeff_fun = (a, b, mk, mz) -> (
        Blast.reactant_limber_chebyshev_coefficients(a, mk, mz),
        Blast.reactant_limber_chebyshev_coefficients(b, mk, mz))
    coeff_args = map(Reactant.to_rarray, (loglin, lognon, Mk, Mz))
    coeff_comp = Reactant.compile(coeff_fun, coeff_args; sync=true)
    coeff0_r = coeff_comp(coeff_args...)
    coeff0 = map(Array, coeff0_r)
    loglin1 = loglin .+ 1e-5 .* reshape(sin.(1:size(loglin, 1)), :, 1)
    lognon1 = lognon .+ 1e-5 .* reshape(cos.(1:size(lognon, 1)), :, 1)
    coeff_args1 = (Reactant.to_rarray(loglin1), Reactant.to_rarray(lognon1), coeff_args[3:end]...)
    coeff1 = map(Array, coeff_comp(coeff_args1...))
    coeff_ref1 = (chebyshev_decomposition(lp, loglin1),
                  chebyshev_decomposition(lp, lognon1))

    Tz, Tk = Blast.get_limber_coords_polynomials(lp, bg.z), workload.plans.T_k_limber
    grid_ref = (10.0 .^ Blast.limber_eval(coeff_ref[1], Tz, Tk),
                10.0 .^ Blast.limber_eval(coeff_ref[2], Tz, Tk))
    grid_fun = (a, b, tz, tk) -> (
        Blast.reactant_limber_grid_from_coefficients(a, tz, tk),
        Blast.reactant_limber_grid_from_coefficients(b, tz, tk))
    grid_args = (coeff0_r..., Reactant.to_rarray(Tz), Reactant.to_rarray(Tk))
    grid_comp = Reactant.compile(grid_fun, grid_args; sync=true)
    grid0 = map(Array, grid_comp(grid_args...))
    grid_args1 = (Reactant.to_rarray(coeff_ref1[1]), Reactant.to_rarray(coeff_ref1[2]), grid_args[3:end]...)
    grid1 = map(Array, grid_comp(grid_args1...))
    grid_ref1 = (10.0 .^ Blast.limber_eval(coeff_ref1[1], Tz, Tk),
                 10.0 .^ Blast.limber_eval(coeff_ref1[2], Tz, Tk))

    @testset "Reactant full Stages 5-6" begin
        @test coeff0[1] ≈ coeff_ref[1] rtol=1e-11 atol=1e-12
        @test coeff0[2] ≈ coeff_ref[2] rtol=1e-11 atol=1e-12
        @test coeff1[1] ≈ coeff_ref1[1] rtol=1e-11 atol=1e-12
        @test coeff1[2] ≈ coeff_ref1[2] rtol=1e-11 atol=1e-12
        @test maximum(abs, coeff1[1] .- coeff0[1]) > 0
        @test maximum(abs, coeff1[2] .- coeff0[2]) > 0
        @test grid0[1] ≈ grid_ref[1] rtol=1e-10 atol=1e-12
        @test grid0[2] ≈ grid_ref[2] rtol=1e-10 atol=1e-12
        @test grid1[1] ≈ grid_ref1[1] rtol=1e-10 atol=1e-12
        @test grid1[2] ≈ grid_ref1[2] rtol=1e-10 atol=1e-12
        @test maximum(abs, grid1[1] .- grid0[1]) > 0
        @test maximum(abs, grid1[2] .- grid0[2]) > 0
    end
    println("stage5_shapes=", map(size, coeff0))
    println("stage5_maxabs=", map((a,b)->maximum(abs,a.-b), coeff0, coeff_ref))
    println("stage6_shapes=", map(size, grid0))
    println("stage6_maxabs=", map((a,b)->maximum(abs,a.-b), grid0, grid_ref))
end

main()

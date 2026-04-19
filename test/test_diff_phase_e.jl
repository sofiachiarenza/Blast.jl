using Test
using Blast
using Random

Random.seed!(4321)

@testset "Differentiation: Phase E (Chebyshev / Limber helpers)" begin

    # ------------------------------------------------------------------
    # Shared Chebyshev/Limber fixtures
    # ------------------------------------------------------------------
    k_min, k_max = 1e-3, 1.0
    chi_min, chi_max = 0.2, 2.0
    K_k, K_chi = 4, 3

    plan_limber = Blast.prepare_chebyshev_plan((log10(k_min), chi_min), (log10(k_max), chi_max), (K_k, K_chi))

    chi_test = collect(LinRange(chi_min, chi_max, 5))
    ell_test = collect(10.0:10.0:40.0)

    @testset "chebinterp_native w.r.t. coefficients" begin
        c0 = rand(6)
        x_grid = collect(LinRange(0.0, 1.0, 8))

        f(c) = sum(Blast.chebinterp_native(c, x_grid, 0.0, 1.0))
        check_gradient(f, c0; skip_zygote=true)
    end

    @testset "build_coeff w.r.t. vals" begin
        plan_1d = Blast.prepare_chebyshev_plan(0.0, 1.0, 8)
        n_nodes = size(plan_1d.fft_plan, 1)
        vals0 = rand(n_nodes)

        f(vals) = sum(Blast.build_coeff(Blast.cϕ, vals, plan_1d).coefs)
        check_gradient(f, vals0)
    end

    @testset "limber_eval w.r.t. coeff matrix" begin
        T_chi = Blast.get_limber_coords_polynomials(plan_limber, chi_test)
        T_k = Blast.get_limber_k_polynomials(plan_limber, ell_test, chi_test; is_log_k=true)
        c0 = rand(size(T_k, 2), size(T_chi, 2))

        f(c) = sum(Blast.limber_eval(c, T_chi, T_k))
        check_gradient(f, c0)
    end
end

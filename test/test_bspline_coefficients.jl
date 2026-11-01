using Test, Blast, LinearAlgebra
using AbstractCosmologicalEmulators
using ForwardDiff, Mooncake, DifferentiationInterface, ADTypes
# ACE also exports `gradient`; take the AD-interface functions explicitly.
using DifferentiationInterface: prepare_gradient, gradient

function bspline_coefficient_objective(::Type{C}, plan, shape) where {C<:Blast.AbstractCoeff}
    return y -> sum(abs2, Blast.build_coeff(C, reshape(y, shape), plan).coefs)
end

@testset "B-spline source coefficient adapter" begin
    for n in (7, 33, 49, 65, 81, 97, 129)
        sites = -4.3 .+ 5.5 .* collect(range(0.0, 1.0; length=n)).^1.2
        plan = CubicBSplinePlan(sites, sites; extrapolation=:throw)
        A = basis_matrix(bspline_basis(plan), sites)
        values = [exp(0.15x) + 0.03j + 0.07r + 0.01x^3
                  for x in sites, j in 1:2, r in 1:3]
        original = copy(values)
        reference = reshape(A \ reshape(values, n, :), size(values))

        for C in (Blast.cϕTT, Blast.cϕT)
            result = @inferred Blast.build_coeff(C, values, plan)
            @test result isa C
            @test size(result.coefs) == (n, 2, 3)
            @test result.coefs ≈ reference rtol=2e-12 atol=2e-12
            @test A * reshape(result.coefs, n, :) ≈ reshape(values, n, :) rtol=2e-12 atol=2e-12
            @test values == original
            @test_throws DimensionMismatch Blast.build_coeff(C, values[1:end-1, :, :], plan)
        end

        primordial = copy(values[:, 1, 1])
        result = @inferred Blast.build_coeff(Blast.cϕ, primordial, plan)
        @test result.coefs ≈ A \ primordial rtol=2e-12 atol=2e-12
        @test !isapprox(result.coefs, primordial; rtol=1e-8, atol=1e-8)
        @test Blast.build_coeff(Blast.cϕ, primordial, nothing) === nothing
        @test_throws DimensionMismatch Blast.build_coeff(Blast.cϕ, primordial[1:end-1], plan)

        # A cubic polynomial lies in the spline space: the solved coefficients
        # must reproduce it at off-site interior points.
        poly(x) = 0.3 - 0.2x + 0.05x^2 + 0.01x^3
        queries = sites[1] .+ (sites[end] - sites[1]) .* [0.123, 0.377, 0.5, 0.811]
        pc = Blast.build_coeff(Blast.cϕ, poly.(sites), plan).coefs
        @test basis_matrix(bspline_basis(plan), queries) * pc ≈ poly.(queries) rtol=1e-10 atol=1e-10
        # Partition of unity: unit coefficients evaluate to one.
        @test basis_matrix(bspline_basis(plan), queries) * ones(n) ≈ ones(length(queries)) rtol=1e-12

        # Dual-valued sources keep their value and partials through the solve.
        dual = ForwardDiff.Dual.(primordial, 1.0)
        dual_result = Blast.build_coeff(Blast.cϕ, dual, plan)
        @test eltype(dual_result.coefs) <: ForwardDiff.Dual
        @test ForwardDiff.value.(dual_result.coefs) ≈ A \ primordial rtol=2e-12 atol=2e-12
        @test ForwardDiff.partials.(dual_result.coefs, 1) ≈ A \ ones(n) rtol=2e-12 atol=2e-12

        storage = repeat(values; inner=(2, 1, 1))
        strided = @view storage[1:2:end, :, :]
        @test Blast.build_coeff(Blast.cϕTT, strided, plan).coefs ≈ reference rtol=2e-12 atol=2e-12
    end
end

@testset "B-spline source adapter prepared gradients" begin
    sites = [-3.0, -2.7, -1.8, -0.9, 0.0, 0.4, 1.2]
    n = length(sites)
    plan = CubicBSplinePlan(sites, sites; extrapolation=:throw)
    A = basis_matrix(bspline_basis(plan), sites)
    for C in (Blast.cϕ, Blast.cϕTT, Blast.cϕT)
        shape = C === Blast.cϕ ? (n,) : (n, 2, 3)
        x = collect(range(0.2, 1.3; length=prod(shape)))
        objective = bspline_coefficient_objective(C, plan, shape)
        forward = AutoForwardDiff()
        reverse = AutoMooncake()
        forward_prep = prepare_gradient(objective, forward, x)
        reverse_prep = prepare_gradient(objective, reverse, x)
        # Reuse both preparations with changed inputs to expose stale caches.
        for y in (x, 1.03 .* x .+ 0.02)
            c = A \ reshape(y, n, :)
            expected = vec(2 .* (transpose(A) \ c))
            gf = gradient(objective, forward_prep, forward, y)
            gr = gradient(objective, reverse_prep, reverse, y)
            @test gf ≈ expected rtol=2e-10 atol=2e-11
            @test gr ≈ expected rtol=2e-10 atol=2e-11
            @test gr ≈ gf rtol=2e-10 atol=2e-11
        end
    end
end

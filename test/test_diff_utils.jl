# test/test_diff_utils.jl
#
# Unit-level gradient tests for the Akima spline primitives in utils.jl.
# Each testset targets one function and one differentiable input, using the
# four backends defined in test_diff_helpers.jl:
#   - ForwardDiff (forward mode, dual numbers)
#   - Zygote      (reverse mode, ChainRules rrules)
#   - Mooncake    (reverse mode, @from_chainrules patches)
#   - FiniteDifferences (central 5-point, reference ground truth)
#
# Strategy: scalar-ise via sum() then call check_gradient().
# Pre-computed intermediates (m_ref, b_ref, …) are captured as constants so
# each closure has exactly one differentiable argument.

using Test
using Blast

include("test_diff_helpers.jl")

@testset "Differentiation: Akima primitives" begin

    # ------------------------------------------------------------------ #
    # Shared test data                                                     #
    # ------------------------------------------------------------------ #
    x_knots = collect(LinRange(0.0, 5.0, 20))   # knot positions, n = 20
    y_vals  = sin.(x_knots)                       # knot values
    x_query = collect(LinRange(0.5, 4.5, 30))    # query grid,   n_q = 30

    # Pre-compute reference intermediates used as captured constants below.
    # These are plain Float64 arrays — no AD involvement.
    m_ref             = Blast._akima_slopes(y_vals, x_knots)          # length n+3 = 23
    b_ref, c_ref, d_ref = Blast._akima_coefficients(x_knots, m_ref)  # 20, 19, 19

    # ------------------------------------------------------------------ #
    # 1. _akima_slopes(u, t)                                              #
    #    rrule: chainrules.jl line 38                                     #
    #    Mooncake: @from_chainrules line 213                              #
    # ------------------------------------------------------------------ #

    @testset "_akima_slopes w.r.t. u (knot values)" begin
        f(u) = sum(Blast._akima_slopes(u, x_knots))
        check_gradient(f, y_vals)
    end

    @testset "_akima_slopes w.r.t. t (knot positions)" begin
        f(t) = sum(Blast._akima_slopes(y_vals, t))
        check_gradient(f, x_knots)
    end

    # ------------------------------------------------------------------ #
    # 2. _akima_coefficients(t, m)                                        #
    #    rrule: chainrules.jl line 58                                     #
    #    Mooncake: @from_chainrules line 214                              #
    #    Returns a tuple (b, c, d); we sum all three components.          #
    # ------------------------------------------------------------------ #

    @testset "_akima_coefficients w.r.t. m (slopes)" begin
        f(m) = let (b, c, d) = Blast._akima_coefficients(x_knots, m)
            sum(b) + sum(c) + sum(d)
        end
        check_gradient(f, m_ref)
    end

    @testset "_akima_coefficients w.r.t. t (knot positions)" begin
        f(t) = let (b, c, d) = Blast._akima_coefficients(t, m_ref)
            sum(b) + sum(c) + sum(d)
        end
        check_gradient(f, x_knots)
    end

    # ------------------------------------------------------------------ #
    # 3. _akima_eval(u, t, b, c, d, tq)                                  #
    #    rrule: chainrules.jl line 102                                    #
    #    Mooncake: @from_chainrules line 215                              #
    #    Six separate closures — one per differentiable argument.         #
    #    Sizes: u→20, t→20, b→20, c→19, d→19, tq→30                     #
    # ------------------------------------------------------------------ #

    @testset "_akima_eval w.r.t. u (knot values)" begin
        f(u) = sum(Blast._akima_eval(u, x_knots, b_ref, c_ref, d_ref, x_query))
        check_gradient(f, y_vals)
    end

    @testset "_akima_eval w.r.t. t (knot positions)" begin
        f(t) = sum(Blast._akima_eval(y_vals, t, b_ref, c_ref, d_ref, x_query))
        check_gradient(f, x_knots)
    end

    @testset "_akima_eval w.r.t. b" begin
        f(b) = sum(Blast._akima_eval(y_vals, x_knots, b, c_ref, d_ref, x_query))
        check_gradient(f, b_ref)
    end

    @testset "_akima_eval w.r.t. c" begin
        f(c) = sum(Blast._akima_eval(y_vals, x_knots, b_ref, c, d_ref, x_query))
        check_gradient(f, c_ref)
    end

    @testset "_akima_eval w.r.t. d" begin
        f(d) = sum(Blast._akima_eval(y_vals, x_knots, b_ref, c_ref, d, x_query))
        check_gradient(f, d_ref)
    end

    @testset "_akima_eval w.r.t. tq (query points)" begin
        f(tq) = sum(Blast._akima_eval(y_vals, x_knots, b_ref, c_ref, d_ref, tq))
        check_gradient(f, x_query)
    end

    # ------------------------------------------------------------------ #
    # 4. _akima_interpolation(u, t, t_new) — vector variant               #
    #    Chains _akima_slopes → _akima_coefficients → _akima_eval         #
    #    All three individual rrules fire in sequence.                    #
    # ------------------------------------------------------------------ #

    @testset "_akima_interpolation (vector) w.r.t. u" begin
        f(u) = sum(Blast._akima_interpolation(u, x_knots, x_query))
        check_gradient(f, y_vals)
    end

    @testset "_akima_interpolation (vector) w.r.t. t" begin
        f(t) = sum(Blast._akima_interpolation(y_vals, t, x_query))
        check_gradient(f, x_knots)
    end

    @testset "_akima_interpolation (vector) w.r.t. t_new" begin
        f(tq) = sum(Blast._akima_interpolation(y_vals, x_knots, tq))
        check_gradient(f, x_query)
    end

    # ------------------------------------------------------------------ #
    # 5. _akima_interpolation(u, t, t_new) — matrix variant              #
    #    Single fused rrule: chainrules.jl line 129                       #
    #    Mooncake: @from_chainrules line 220                              #
    #    u is (n_knots × n_cols); independent column interpolations.      #
    # ------------------------------------------------------------------ #

    @testset "_akima_interpolation (matrix) w.r.t. U" begin
        Y_mat = hcat(sin.(x_knots), cos.(x_knots), sin.(2 .* x_knots))  # 20×3
        f(U) = sum(Blast._akima_interpolation(U, x_knots, x_query))
        check_gradient(f, Y_mat)
    end

    @testset "_akima_interpolation (matrix) w.r.t. t" begin
        Y_mat = hcat(sin.(x_knots), cos.(x_knots), sin.(2 .* x_knots))
        f(t) = sum(Blast._akima_interpolation(Y_mat, t, x_query))
        check_gradient(f, x_knots)
    end

    @testset "_akima_interpolation (matrix) w.r.t. t_new" begin
        Y_mat = hcat(sin.(x_knots), cos.(x_knots), sin.(2 .* x_knots))
        f(tq) = sum(Blast._akima_interpolation(Y_mat, x_knots, tq))
        check_gradient(f, x_query)
    end

end

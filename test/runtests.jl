"""
Main test runner for Blast.jl

This file orchestrates all tests by including individual test files.
Each test file is self-contained and tests a specific component.

To run all tests:
    julia --project=. -e 'using Pkg; Pkg.test()'

To run specific test files during development:
    julia --project=. test/test_cosmo.jl
"""

using Test
using Blast

# Shared test data and utilities
include("test_fixtures.jl")
global data = load_reference_data()

@testset "Blast.jl Comprehensive Tests" begin
    
    @testset "Utilities" begin
        include("test_utils.jl")
    end

    @testset "Cosmology" begin
        include("test_cosmo.jl")
    end

    @testset "Chebyshev Engine" begin
        include("test_chebcoefs.jl")
    end

    @testset "Probes and Kernels" begin
        include("test_probes.jl")
    end

    @testset "Setup and Workspace" begin
        include("test_setup.jl")
    end

    @testset "C_ℓ Pipeline" begin
        include("test_cls.jl")
    end

    @testset "Limber" begin
        include("test_limber.jl")
    end

    @testset "Automatic Differentiation" begin
        include("test_differentiation.jl")
    end
    
end

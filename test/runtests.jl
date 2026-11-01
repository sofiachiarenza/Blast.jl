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

# Development selectors: `Pkg.test(test_args=["bspline-coefficients"])` and
# `Pkg.test(test_args=["spectral-setup"])` / `["tensor-bundle"]` run only those focused files. The
# default (no arguments) suite below always includes both.
if !isempty(ARGS)
    selectors = Dict("png-transfer-ordering" => ("PNG Transfer-Leg Ordering", "test_png_transfer_ordering.jl"),
                     "arfcamb-validation" => ("ARFCAMB Validation", "test_arfcamb_validation.jl"),
                     "legacy-characterization" => ("Legacy Characterization", "test_legacy_characterization.jl"),
                     "radial-quadrature" => ("Radial Quadrature", "test_radial_quadrature.jl"),
                     "radial-public" => ("Radial Public Path", "test_radial_public.jl"),
                     "radial-ad" => ("Radial AD", "test_radial_ad.jl"),
                     "radial-fixed-ad" => ("Radial Fixed Background AD", "test_radial_fixed_background_ad.jl"),
                     "radial-bundle" => ("Radial Bundle", "test_radial_bundle.jl"),
                     "bspline-coefficients" => ("B-spline Coefficient Adapter", "test_bspline_coefficients.jl"),
                     "spectral-setup" => ("Spectral Setup and Banks", "test_spectral_setup.jl"),
                     "bspline-backend" => ("B-spline Backend", "test_bspline_backend.jl"),
                     "bspline-backend-ad" => ("B-spline Backend AD", "test_bspline_backend_ad.jl"),
                     "ratio-canary" => ("Ratio Canary", "test_ratio_canary.jl"),
                     "bspline-ratio" => ("B-spline Ratio", "test_bspline_ratio.jl"),
                     "bspline-ratio-ad" => ("B-spline Ratio AD", "test_bspline_ratio_ad.jl"),
                     "ratio-public" => ("Ratio Public Correctness", "test_ratio_public_correctness.jl"),
                     "ratio-probe-ad" => ("Ratio Probe AD", "test_ratio_probe_ad.jl"),
                     "probes" => ("Probes and Kernels", "test_probes.jl"),
                     "promotion" => ("PowerSpectrum Promotion", "test_powerspectrum_promotion.jl"),
                     "tensor-bundle" => ("Tensor Bundle Loader", "test_tensor_bundle.jl"),
                     "tensor-bundle-finalizer" => ("Tensor Bundle Finalizer", "test_tensor_bundle_finalizer.jl"))
    (length(ARGS) == 1 && haskey(selectors, only(ARGS))) || error("Unknown test selector(s): $(ARGS)")
    name, file = selectors[only(ARGS)]
    if file == "test_probes.jl"
        include("test_fixtures.jl")
        global data = load_reference_data()
    end
    @testset "$name" begin
        include(file)
    end
    exit(0)
end

# Shared test data and utilities
include("test_fixtures.jl")
include("test_diff_helpers.jl")
global data = load_reference_data()

@testset "Blast.jl Comprehensive Tests" begin
    
    @testset "Utilities" begin
        include("test_utils.jl")
    end

    @testset "Radial Quadrature" begin
        include("test_radial_quadrature.jl")
    end

    @testset "Radial Public Path" begin
        include("test_radial_public.jl")
    end

    @testset "Radial AD" begin
        include("test_radial_ad.jl")
    end

    @testset "Radial Fixed Background AD" begin
        include("test_radial_fixed_background_ad.jl")
    end

    @testset "Cosmology" begin
        include("test_cosmo.jl")
    end

    @testset "Chebyshev Engine" begin
        include("test_chebcoefs.jl")
    end

    @testset "B-spline Coefficient Adapter" begin
        include("test_bspline_coefficients.jl")
    end

    @testset "Spectral Setup and Banks" begin
        include("test_spectral_setup.jl")
    end

    @testset "B-spline Backend" begin
        include("test_bspline_backend.jl")
    end

    @testset "B-spline Backend AD" begin
        include("test_bspline_backend_ad.jl")
    end

    @testset "Ratio Canary" begin
        include("test_ratio_canary.jl")
    end

    @testset "B-spline Ratio" begin
        include("test_bspline_ratio.jl")
    end

    @testset "B-spline Ratio AD" begin
        include("test_bspline_ratio_ad.jl")
    end

    @testset "Ratio Public Correctness" begin
        include("test_ratio_public_correctness.jl")
    end

    @testset "Ratio Probe AD" begin
        include("test_ratio_probe_ad.jl")
    end

    @testset "PowerSpectrum Promotion" begin
        include("test_powerspectrum_promotion.jl")
    end

    @testset "Tensor Bundle Loader" begin
        include("test_tensor_bundle.jl")
    end

    @testset "Radial Bundle" begin
        include("test_radial_bundle.jl")
    end

    @testset "Tensor Bundle Finalizer" begin
        include("test_tensor_bundle_finalizer.jl")
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

    @testset "PNG Transfer-Leg Ordering" begin
        include("test_png_transfer_ordering.jl")
    end

    @testset "Legacy Characterization" begin
        include("test_legacy_characterization.jl")
    end

    @testset "Public Input Guardrails" begin
        include("test_guardrails.jl")
    end

    @testset "3x2pt Type Stability" begin
        include("test_type_stability_3x2.jl")
    end

    @testset "Limber" begin
        include("test_limber.jl")
    end

    @testset "ARFCAMB Validation" begin
        include("test_arfcamb_validation.jl")
    end

@testset "Automatic Differentiation" begin
    include("test_diff_cls.jl")
    include("test_diff_cosmo.jl")
    include("test_diff_setup.jl")
    include("test_diff_phase_e.jl")
    include("test_diff_glue.jl")
end
    
end

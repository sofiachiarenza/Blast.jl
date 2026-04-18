# test/test_diff_helpers.jl
#
# Shared helper for gradient consistency tests.
# Include this file from within a @testset block — the @test calls inside
# check_gradient are attributed to the enclosing testset.

using DifferentiationInterface
using ADTypes
using ForwardDiff
using FiniteDifferences
using Zygote
using Mooncake
using Test

"""
    check_gradient(f, x; kwargs...)

Compute the gradient of the scalar-valued function `f` at input `x` using
up to four backends and verify numerical consistency.

All active AD backends (ForwardDiff, Zygote, Mooncake) are compared pairwise
at `rtol_ad`. Additionally every AD gradient is compared against a central
5-point FiniteDifferences reference at the looser `rtol_fd`.

# Keyword arguments
- `rtol_ad = 1e-9`        — relative tolerance for AD-vs-AD comparisons
- `rtol_fd = 1e-5`        — relative tolerance for AD-vs-FiniteDifferences
- `atol    = 0.0`         — absolute tolerance for all comparisons
- `skip_forward  = false` — exclude ForwardDiff
- `skip_fd       = false` — exclude the FiniteDifferences reference check
- `skip_zygote   = false` — exclude Zygote
- `skip_mooncake = false` — exclude Mooncake

# Returns
A `NamedTuple` mapping each active backend name (`:ForwardDiff`, `:Zygote`,
`:Mooncake`) to the gradient it computed.
"""
function check_gradient(f, x;
                        rtol_ad=1e-9,
                        rtol_fd=1e-5,
                        atol=0.0,
                        skip_forward=false,
                        skip_fd=false,
                        skip_zygote=false,
                        skip_mooncake=false)
    labels   = Symbol[]
    backends = Any[]

    if !skip_forward
        push!(labels, :ForwardDiff)
        push!(backends, AutoForwardDiff())
    end
    if !skip_zygote
        push!(labels, :Zygote)
        push!(backends, AutoZygote())
    end
    if !skip_mooncake
        push!(labels, :Mooncake)
        push!(backends, AutoMooncake())
    end

    grads = [DifferentiationInterface.gradient(f, b, x) for b in backends]

    # ---- Pairwise AD consistency ----------------------------------------
    for i in 1:length(grads), j in (i+1):length(grads)
        @test grads[i] ≈ grads[j]  rtol=rtol_ad  atol=atol
    end

    # ---- AD vs FiniteDifferences reference ------------------------------
    if !skip_fd && !isempty(grads)
        fdm  = central_fdm(5, 1)
        g_fd = DifferentiationInterface.gradient(f, AutoFiniteDifferences(fdm=fdm), x)
        for g in grads
            @test g ≈ g_fd  rtol=rtol_fd  atol=atol
        end
    end

    return NamedTuple{Tuple(labels)}(Tuple(grads))
end

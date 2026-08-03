include(joinpath(@__DIR__, "staged_reactant_gradient.jl"))
using ForwardDiff
using LinearAlgebra
using Random

Random.seed!(8642)
Reactant.set_default_backend("cpu")
nℓ, nχ, nR, nA, nB = 1, 2, 2, 1, 1
w = ntuple(_ -> randn(nℓ, nχ, nR), 17)
kernels = ntuple(_ -> randn(nA, nχ, nR), 6)
integration = (ones(nχ), ones(nR), ones(nχ),)
pref_names = (:δδ, :δRSD, :RSDRSD, :δμ, :μμ, :μRSD, :δfNL, :fNLδ, :fNLRSD,
              :RSDfNL, :μfNL, :fNLμ, :fNLfNL, :γγ, :γIA, :IAγ, :IAIA, :δγ,
              :δIA, :RSDγ, :RSDIA, :μγ, :μIA, :fNLγ, :fNLIA)
pref = NamedTuple{pref_names}(ntuple(_ -> ones(nℓ), length(pref_names)))
C_terms = ntuple(_ -> randn(1, nA, nB), 6)
config = (
    kernels=kernels, integration=integration, C_terms=C_terms,
    finalization=(randn(2), Matrix{Float64}(I, 2, 2), randn(1, 2), ones(1)),
    Δχ=0.5, pref=pref,
)
upstream = ntuple(_ -> ones(1, nA, nB), 3)
grads = staged_endpoint_pullback(w, config, upstream)
loss(wi, i) = begin
    ww = ntuple(j -> j == i ? wi : w[j], 17)
    y = _endpoint_value(ww, config)
    sum(y[1]) + sum(y[2]) + sum(y[3])
end
for i in 1:17
    ref = ForwardDiff.gradient(x -> loss(x, i), w[i])
    println("w", i, "_gradient_error=", maximum(abs, grads[i] .- ref))
end

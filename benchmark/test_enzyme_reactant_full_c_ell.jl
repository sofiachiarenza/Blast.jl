using Blast
using Enzyme
using ForwardDiff
using LinearAlgebra
using Reactant
using Random

function make_prefactor()
    names = (
        :δδ, :δRSD, :RSDRSD, :δμ, :μμ, :μRSD, :δfNL, :fNLδ,
        :fNLRSD, :RSDfNL, :μfNL, :fNLμ, :fNLfNL, :γγ, :γIA,
        :IAγ, :IAIA, :δγ, :δIA, :RSDγ, :RSDIA, :μγ, :μIA,
        :fNLγ, :fNLIA,
    )
    return NamedTuple{names}(ntuple(_ -> randn(1), length(names)))
end

function main()
    Reactant.set_default_backend("cpu")
    Random.seed!(24680)
    nb = 2
    nχ, nR, nℓ, nout = 4, 3, 1, 2
    w0 = ntuple(_ -> randn(nℓ, nχ, nR), 17)
    kernels = ntuple(_ -> randn(nb, nχ, nR), 6)
    integration = (rand(nχ), rand(nR), rand(nχ), 0.6)
    pref = make_prefactor()
    final_arrays = (
        randn(nℓ, nb, nb), randn(2, nb, nb),
        randn(nℓ, nb, nb), randn(2, nb, nb),
        randn(nℓ, nb, nb), randn(2, nb, nb),
    )
    ell2 = rand(nℓ + 2)
    transform = Matrix{Float64}(I, nℓ + 2, nℓ + 2)
    T_eval = randn(nout, nℓ + 2)
    invell2 = rand(nout)
    Δχ = 0.6

    args_host = (
        w0..., kernels..., integration[1], integration[2], integration[3],
        final_arrays..., ell2, transform, T_eval, invell2, Δχ, pref,
    )
    loss(args...) = begin
        r = Blast.reactant_full_3x2pt(args...)
        sum(r[1]) + sum(r[2]) + sum(r[3])
    end
    reference = ForwardDiff.gradient(x -> loss(x, args_host[2:end]...), args_host[1])

    wR = map(Reactant.to_rarray, w0)
    otherR = map(Reactant.to_rarray, args_host[18:end-4])
    # Keep the static prefactor and scalar Δχ in the closure; all array inputs
    # that carry the W/C_l graph remain explicit Reactant arguments.
    grad_dynamic = let pref=pref, Δχ=Δχ, integration=integration, final_arrays=final_arrays,
                       ell2=ell2, transform=transform, T_eval=T_eval, invell2=invell2
        function (w1, rest...)
            loss_dynamic(w1, rest...) = begin
                r = Blast.reactant_full_3x2pt(
                    w1, rest[1:16]...,
                    rest[17:22]...,
                    rest[23:25]...,
                    rest[26:31]...,
                    rest[32], rest[33], rest[34], Δχ, pref,
                )
                sum(r[1]) + sum(r[2]) + sum(r[3])
            end
            Enzyme.gradient(Reverse, loss_dynamic, w1, Const(rest...))[1]
        end
    end
    # The closure form above is intentionally explicit about the active W
    # argument; build the compiled function with all arrays in positional order.
    all_arrays = (
        wR..., map(Reactant.to_rarray, kernels)...,
        Reactant.to_rarray(integration[1]), Reactant.to_rarray(integration[2]), Reactant.to_rarray(integration[3]),
        map(Reactant.to_rarray, final_arrays)...,
        Reactant.to_rarray(ell2), Reactant.to_rarray(transform), Reactant.to_rarray(T_eval), Reactant.to_rarray(invell2),
    )
    function grad_fun(w1, w2, w3, w4, w5, w6, w7, w8, w9, w10, w11, w12, w13, w14, w15, w16, w17,
                     k1, k2, k3, k4, k5, k6, wi, wr, ch,
                     c1, c2, c3, c4, c5, c6, ell2_, M_, Te_, inv2_)
        loss_dynamic(w1_, args...) = begin
            r = Blast.reactant_full_3x2pt(
                w1_, args[1:16]...,
                args[17:22]..., args[23:25]..., args[26:31]...,
                args[32], args[33], args[34], args[35], Δχ, pref,
            )
            sum(r[1]) + sum(r[2]) + sum(r[3])
        end
        return Enzyme.gradient(
            Reverse, loss_dynamic, w1,
            Const(w2), Const(w3), Const(w4), Const(w5), Const(w6), Const(w7), Const(w8),
            Const(w9), Const(w10), Const(w11), Const(w12), Const(w13), Const(w14), Const(w15), Const(w16), Const(w17),
            Const(k1), Const(k2), Const(k3), Const(k4), Const(k5), Const(k6), Const(wi), Const(wr), Const(ch),
            Const(c1), Const(c2), Const(c3), Const(c4), Const(c5), Const(c6), Const(ell2_), Const(M_), Const(Te_), Const(inv2_),
        )[1]
    end
    compiled = Reactant.@compile sync=true grad_fun(all_arrays...)
    result = compiled(all_arrays...)
    Reactant.synchronize(result)
    println("Reactant=", pkgversion(Reactant))
    println("gradient_shape=", size(reference))
    println("gradient_error=", maximum(abs, Array(result) .- reference))
    println("gradient_norm=", maximum(abs, reference))
end

main()

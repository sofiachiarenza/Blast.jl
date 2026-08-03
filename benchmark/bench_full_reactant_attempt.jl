using BenchmarkTools
using Blast
using NPZ
using Reactant

const ROOT = "/home/marcobonici/Desktop/work/CosmologicalEmulators/blast_play/for_marco_copy"
include(joinpath(ROOT, "performance", "workload.jl"))

function baseline(workload)
    return full_power_to_cls(workload; pass_plans=true)
end

function main()
    Reactant.set_default_backend("cpu")
    Blast.FFTW.set_num_threads(8)
    workload = build_workload(ROOT, 8)

    println("Reactant=", pkgversion(Reactant))
    println("Blast=", pathof(Blast))
    println("grid=(nχ=", length(Blast.χ), ", nR=", length(Blast.R), ")")

    println("building ordinary full reference")
    reference = baseline(workload)
    println("ordinary sizes=", (size(reference.cl_gg), size(reference.cl_gs), size(reference.cl_ss)))
    println("ordinary sums=", (sum(reference.cl_gg), sum(reference.cl_gs), sum(reference.cl_ss)))

    baseline_call = () -> baseline(workload)
    baseline_trial = @benchmark $baseline_call() samples=3 evals=1
    println("ordinary_full_min_ms=", minimum(baseline_trial).time / 1e6)
    println("ordinary_full_alloc_bytes=", minimum(baseline_trial).memory)

    # Keep all configuration and precomputed kernel data outside the traced
    # argument list. The three P(k,z) arrays are the dynamic inputs.
    full_dynamic = let plans=workload.plans,
                       bg=workload.bg,
                       W_template=workload.W_template,
                       galaxy=workload.galaxy,
                       lensing=workload.lensing,
                       ell=workload.ell
        (pk, pk_lin, pk_nonlin) -> begin
            spectrum = Blast.prepare_pk_workspace(plans, pk, pk_lin, pk_nonlin, bg)
            W = Blast.compute_w(W_template, spectrum)
            cl_gg = Blast.get_Cℓ(ell, galaxy, spectrum, W, bg, plans)
            cl_ss = Blast.get_Cℓ(ell, lensing, spectrum, W, bg, plans)
            cl_gs = Blast.get_Cℓ(ell, galaxy, lensing, spectrum, W, bg, plans)
            return (; cl_gg, cl_gs, cl_ss)
        end
    end

    pk_r = Reactant.to_rarray(workload.pk)
    pk_lin_r = Reactant.to_rarray(workload.pk_limber_lin)
    pk_nonlin_r = Reactant.to_rarray(workload.pk_limber_nonlin)

    println("compiling full Reactant wrapper")
    try
        compiled = Reactant.compile(
            full_dynamic,
            (pk_r, pk_lin_r, pk_nonlin_r);
            sync=true,
        )
        result = compiled(pk_r, pk_lin_r, pk_nonlin_r)
        Reactant.synchronize(result)
        println("reactant_full_compiled=true")
        println("reactant_sizes=", (size(Array(result.cl_gg)), size(Array(result.cl_gs)), size(Array(result.cl_ss))))
        println("reactant_sums=", (sum(Array(result.cl_gg)), sum(Array(result.cl_gs)), sum(Array(result.cl_ss))))
        compiled_call = () -> compiled(pk_r, pk_lin_r, pk_nonlin_r)
        compiled_trial = @benchmark $compiled_call() samples=3 evals=1
        println("reactant_full_min_ms=", minimum(compiled_trial).time / 1e6)
        println("reactant_full_alloc_bytes=", minimum(compiled_trial).memory)
    catch err
        println("reactant_full_compiled=false")
        showerror(stdout, err, catch_backtrace())
        println()
    end
end

main()

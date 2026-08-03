using AbstractCosmologicalEmulators
using BenchmarkTools
using Blast
using Reactant

const ROOT = "/home/marcobonici/Desktop/work/CosmologicalEmulators/blast_play/for_marco_copy"
include(joinpath(ROOT, "performance", "workload.jl"))

function max_tuple_error(actual, expected)
    return maximum(maximum(abs, Array(a) .- b) for (a, b) in zip(actual, expected))
end

function main()
    Reactant.set_default_backend("cpu")
    Blast.FFTW.set_num_threads(8)
    workload = build_workload(ROOT, 8)
    ordinary_spectrum = prepare_spectrum(workload)
    ordinary_W = Blast.compute_w(workload.W_template, ordinary_spectrum)

    plan_1d = prepare_chebyshev_plan(
        minimum(Blast.k_cheb), maximum(Blast.k_cheb), length(Blast.k_cheb) - 1;
        size_nd=(length(Blast.k_cheb),), dim=1,
    )
    transform = Blast.reactant_chebyshev_matrix(plan_1d)
    transform_r = Reactant.to_rarray(transform)

    prepare_fun = let bg=workload.bg
        (pk, transform_) -> Blast.reactant_prepare_nonlimber_spectrum(pk, bg, transform_)
    end

    pk0 = workload.pk
    pk1 = 0.91 .* workload.pk .+ 0.13
    pk0r = Reactant.to_rarray(pk0)
    pk1r = Reactant.to_rarray(pk1)
    prepare_compiled = Reactant.compile(prepare_fun, (pk0r, transform_r); sync=true)
    flat0 = prepare_compiled(pk0r, transform_r)
    flat1 = prepare_compiled(pk1r, transform_r)
    flat0_arrays = map(Array, flat0)
    flat1_arrays = map(Array, flat1)

    ordinary1 = Blast.prepare_pk_workspace(
        workload.plans, pk1, workload.pk_limber_lin,
        workload.pk_limber_nonlin, workload.bg,
    )
    ordinary0 = (
        ordinary_spectrum.cϕTT.coefs,
        ordinary_spectrum.cϕT.coefs,
        ordinary_spectrum.cϕ.coefs,
    )
    ordinary1_arrays = (
        ordinary1.cϕTT.coefs,
        ordinary1.cϕT.coefs,
        ordinary1.cϕ.coefs,
    )

    println("Reactant=", pkgversion(Reactant))
    println("grid=(nχ=", length(Blast.χ), ", nR=", length(Blast.R), ")")
    println("prepare_parity_error=", max_tuple_error(flat0_arrays, ordinary0))
    println("prepare_new_input_error=", max_tuple_error(flat1_arrays, ordinary1_arrays))
    println("prepare_compiled_delta=", maximum(maximum(abs, a .- b) for (a, b) in zip(flat1_arrays, flat0_arrays)))
    println("prepare_reference_delta=", maximum(maximum(abs, a .- b) for (a, b) in zip(ordinary1_arrays, ordinary0)))

    # Use dynamic T_lijk arguments. They are not captured constants, so the W
    # compilation test can also change a precomputed kernel after tracing.
    T_original = (
        Blast.T_tildes.T_2_00,
        Blast.T_tildes.T_minus2_00,
        Blast.T_tildes.T_0_00,
        Blast.T_tildes.T_0_02,
        Blast.T_tildes.T_0_20,
        Blast.T_tildes.T_2_02,
        Blast.T_tildes.T_2_20,
        Blast.T_tildes.T_2_22,
    )
    T_lfirst = map(T -> permutedims(T, (4, 1, 2, 3)), T_original)
    T_r = map(Reactant.to_rarray, T_lfirst)
    c_r = (
        Reactant.to_rarray(flat0_arrays[1]),
        Reactant.to_rarray(flat0_arrays[2]),
        Reactant.to_rarray(flat0_arrays[3] .* ones(1, size(T_original[1], 2))),
    )
    c1_r = (
        Reactant.to_rarray(flat1_arrays[1]),
        Reactant.to_rarray(flat1_arrays[2]),
        Reactant.to_rarray(flat1_arrays[3] .* ones(1, size(T_original[1], 2))),
    )

    w_compiled = Reactant.compile(
        Blast.reactant_compute_w,
        (c_r..., T_r...);
        sync=true,
    )
    w0 = map(Array, w_compiled(c_r..., T_r...))
    w1 = map(Array, w_compiled(c1_r..., T_r...))

    ordinary_fields = (
        ordinary_W.w_2_00_ϕTT.w,
        ordinary_W.w_minus2_00_ϕTT.w,
        ordinary_W.w_0_00_ϕTT.w,
        ordinary_W.w_0_02_ϕTT.w,
        ordinary_W.w_0_20_ϕTT.w,
        ordinary_W.w_2_02_ϕTT.w,
        ordinary_W.w_2_20_ϕTT.w,
        ordinary_W.w_2_22_ϕTT.w,
        ordinary_W.w_2_00_ϕT.w,
        ordinary_W.w_2_00_ϕT_R1.w,
        ordinary_W.w_0_00_ϕT.w,
        ordinary_W.w_0_00_ϕT_R1.w,
        ordinary_W.w_2_02_ϕT.w,
        ordinary_W.w_2_02_ϕT_R1.w,
        ordinary_W.w_2_20_ϕT.w,
        ordinary_W.w_2_20_ϕT_R1.w,
        ordinary_W.w_2_00_ϕ.w,
    )
    ordinary_W1 = Blast.compute_w(
        workload.W_template,
        ordinary1,
    )
    ordinary_fields1 = (
        ordinary_W1.w_2_00_ϕTT.w,
        ordinary_W1.w_minus2_00_ϕTT.w,
        ordinary_W1.w_0_00_ϕTT.w,
        ordinary_W1.w_0_02_ϕTT.w,
        ordinary_W1.w_0_20_ϕTT.w,
        ordinary_W1.w_2_02_ϕTT.w,
        ordinary_W1.w_2_20_ϕTT.w,
        ordinary_W1.w_2_22_ϕTT.w,
        ordinary_W1.w_2_00_ϕT.w,
        ordinary_W1.w_2_00_ϕT_R1.w,
        ordinary_W1.w_0_00_ϕT.w,
        ordinary_W1.w_0_00_ϕT_R1.w,
        ordinary_W1.w_2_02_ϕT.w,
        ordinary_W1.w_2_02_ϕT_R1.w,
        ordinary_W1.w_2_20_ϕT.w,
        ordinary_W1.w_2_20_ϕT_R1.w,
        ordinary_W1.w_2_00_ϕ.w,
    )
    println("compute_w_parity_error=", max_tuple_error(w0, ordinary_fields))
    println("compute_w_new_input_error=", max_tuple_error(w1, ordinary_fields1))
    println("compute_w_compiled_delta=", maximum(maximum(abs, a .- b) for (a, b) in zip(w1, w0)))
    println("compute_w_reference_delta=", maximum(maximum(abs, a .- b) for (a, b) in zip(ordinary_fields1, ordinary_fields)))

    prep_call = () -> prepare_compiled(pk0r, transform_r)
    w_call = () -> w_compiled(c_r..., T_r...)
    println("prepare_reactant_ms=", minimum(@benchmark $prep_call() samples=3 evals=1).time / 1e6)
    println("compute_w_reactant_ms=", minimum(@benchmark $w_call() samples=3 evals=1).time / 1e6)
end

main()

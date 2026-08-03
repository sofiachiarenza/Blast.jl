using AbstractCosmologicalEmulators
using Blast
using Enzyme
using ForwardDiff
using Reactant
using Random

const ROOT = "/home/marcobonici/Desktop/work/CosmologicalEmulators/blast_play/for_marco_copy"
include(joinpath(ROOT, "performance", "workload.jl"))

function test_chebyshev_and_grid()
    Random.seed!(1234)
    nk, nz, nl, nout = 5, 4, 7, 6
    values = randn(nk, nz)
    Mk = randn(3, nk)
    Mz = randn(4, nz)
    Tz = randn(nz, 4)
    Tk = randn(nout, 3, nz)
    loss_coeff(v) = sum(Blast.reactant_limber_chebyshev_coefficients(v, Mk, Mz))
    ref_coeff = ForwardDiff.gradient(loss_coeff, values)
    loss_eval(c) = sum(Blast.reactant_limber_grid_from_coefficients(c, Tz, Tk))
    ref_eval = ForwardDiff.gradient(loss_eval, randn(3, 4))

    grad_coeff(v, Mk_, Mz_) = Enzyme.gradient(
        Reverse, (x, a, b) -> sum(Blast.reactant_limber_chebyshev_coefficients(x, a, b)),
        v, Const(Mk_), Const(Mz_),
    )[1]
    grad_eval(c, Tz_, Tk_) = Enzyme.gradient(
        Reverse, (x, a, b) -> sum(Blast.reactant_limber_grid_from_coefficients(x, a, b)),
        c, Const(Tz_), Const(Tk_),
    )[1]
    vals_r, Mk_r, Mz_r = map(Reactant.to_rarray, (values, Mk, Mz))
    coeff_compiled = Reactant.@compile sync=true grad_coeff(vals_r, Mk_r, Mz_r)
    coeff_result = coeff_compiled(vals_r, Mk_r, Mz_r)
    c = randn(3, 4)
    Tz_r, Tk_r = map(Reactant.to_rarray, (Tz, Tk))
    c_r = Reactant.to_rarray(c)
    eval_ref = ForwardDiff.gradient(loss_eval, c)
    eval_compiled = Reactant.@compile sync=true grad_eval(c_r, Tz_r, Tk_r)
    eval_result = eval_compiled(c_r, Tz_r, Tk_r)
    Reactant.synchronize(coeff_result)
    Reactant.synchronize(eval_result)
    println("coeff_gradient_error=", maximum(abs, Array(coeff_result) .- ref_coeff))
    println("grid_gradient_error=", maximum(abs, Array(eval_result) .- eval_ref))
end

function test_power_products()
    workload = build_workload(ROOT, 8)
    bg = workload.bg
    pk = workload.pk_limber_lin
    pknl = workload.pk_limber_nonlin
    loss(p) = sum(Blast.reactant_limber_power_products(p, pknl, bg)[1])
    reference = ForwardDiff.gradient(loss, pk)
    loss_dynamic(p, pknl_) = sum(Blast.reactant_limber_power_products(p, pknl_, bg)[1])
    grad_fun(p, pknl_) = Enzyme.gradient(Reverse, loss_dynamic, p, Const(pknl_))[1]
    p_r, pknl_r = Reactant.to_rarray(pk), Reactant.to_rarray(pknl)
    compiled = Reactant.@compile sync=true grad_fun(p_r, pknl_r)
    result = compiled(p_r, pknl_r)
    Reactant.synchronize(result)
    println("power_product_gradient_error=", maximum(abs, Array(result) .- reference))
    println("power_product_gradient_norm=", maximum(abs, reference))
end

Reactant.set_default_backend("cpu")
test_chebyshev_and_grid()
test_power_products()

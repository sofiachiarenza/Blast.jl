using Blast
using NPZ

function build_realistic_workload(root::AbstractString, fftw_threads::Integer)
    Blast.FFTW.set_num_threads(fftw_threads)
    data_directory = joinpath(root, "data")
    nz_data = npzread(joinpath(data_directory, "lsst_nz.npz"))
    nz_g = nz_data["nz_clustering"]
    z_g = nz_data["z_clustering"]
    nz_s = nz_data["nz_lensing"]
    z_s = nz_data["z_lensing"]
    z_g == z_s || error("clustering and lensing n(z) grids differ")

    ell = npzread(joinpath(data_directory, "ell_n5k.npz"))
    pk = npzread(joinpath(data_directory, "pk_nonlimber_161.npz"))["pk"]
    pk_limber_lin = npzread(joinpath(data_directory, "pk_limber_lin.npz"))["pk"]
    pk_limber_nonlin = npzread(joinpath(data_directory, "pk_limber_nonlin.npz"))["pk"]

    cosmo = characterization_cosmology()
    bg = Blast.Background(cosmo)
    bias_values = [1.376695, 1.451179, 1.528404, 1.607983, 1.689579,
                   1.772899, 1.857700, 1.943754, 2.030887, 2.118943]
    slope_values = [0.1648, 0.2496, 0.2708, 0.3300, 0.3880,
                    0.2960, 0.3580, 0.3960, 0.4320, 0.5680]
    bias = bias_values .* ones(size(nz_g, 1), length(bg.z))
    slope = slope_values .* ones(size(nz_g, 1), length(bg.z))

    # Opt-in (BLAST_REALISTIC_WITH_PNG=1): add PNG so all 17 projected-matter
    # components are active. Same constructor pattern as the characterization
    # workload. Default (unset) workload is unchanged.
    with_png = get(ENV, "BLAST_REALISTIC_WITH_PNG", "0") == "1"
    png = with_png ? Blast.PrimordialNonGaussianity(nz=nz_g, z=z_g, bias=bias, f_NL=7.0, p=1.0) : nothing
    galaxy = Blast.GalaxyClustering(
        δ=Blast.NumberCounts(nz=nz_g, z=z_g, bias=bias),
        RSD=Blast.RedshiftSpaceDistortions(nz=nz_g, z=z_g, growth_rate=bg.f),
        μ=Blast.MagnificationBias(nz=nz_g, z=z_g, s=slope),
        PNG=png,
    )
    lensing = Blast.WeakLensing(
        γ=Blast.CosmicShear(nz=nz_s, z=z_s),
        IA=Blast.IntrinsicAlignment(nz=nz_s, z=z_s, A=1.72),
    )
    Blast.evaluate_components!(galaxy, bg)
    Blast.evaluate_components!(lensing, bg)
    weights_template, plans = Blast.SetUp(galaxy, lensing)
    active_fields = [name for name in fieldnames(typeof(weights_template))
                     if !isnothing(getfield(weights_template, name))]
    with_png && length(active_fields) != 17 &&
        error("PNG workload must activate all 17 projected components; got $(length(active_fields))")
    return (; ell, pk, pk_limber_lin, pk_limber_nonlin, bg, galaxy, lensing,
            weights_template, plans, with_png, active_fields)
end

function prepare_realistic_spectrum(workload, x; bg=workload.bg)
    return Blast.prepare_pk_workspace(
        workload.plans,
        x[1] .* workload.pk,
        x[2] .* workload.pk_limber_lin,
        x[3] .* workload.pk_limber_nonlin,
        bg,
    )
end

function project_realistic(workload, spectrum, weights)
    cl_gg = Blast.get_Cℓ(workload.ell, workload.galaxy, spectrum, weights,
                         workload.bg, workload.plans)
    cl_gl = Blast.get_Cℓ(workload.ell, workload.galaxy, workload.lensing,
                         spectrum, weights, workload.bg, workload.plans)
    cl_ll = Blast.get_Cℓ(workload.ell, workload.lensing, spectrum, weights,
                         workload.bg, workload.plans)
    return (; cl_gg, cl_gl, cl_ll)
end

"""
Validation-only: `bg` with every stored vector lifted to the ForwardDiff Dual
type of `like` with ZERO partials (`zero(like) .+ v`), built through
Background's raw constructor. `cosmo` and `z_of_χ` are the identical objects, no
ODE is re-solved, and every physical array is bit-exact. `get_As`/`get_ns`
convert the fixed Float64 cosmology to the Dual element type, so the primordial
coefficients (PNG's cϕ) promote to Dual with zero amplitude derivative. Used
only by the independent ForwardDiff reference; the Float64/Mooncake objective
never sees it.
"""
function lift_background(bg, like)
    lift(v) = zero(like) .+ v
    T = typeof(zero(like) + one(eltype(bg.z)))
    return Blast.Background{T, typeof(bg.z_of_χ), typeof(bg.cosmo)}(
        bg.cosmo, lift(bg.z), lift(bg.χ), lift(bg.H), lift(bg.D), lift(bg.D_norm), lift(bg.f), bg.z_of_χ)
end

objective(result, normalization) =
    normalization[1] * sum(abs2, result.cl_gg) +
    normalization[2] * sum(abs2, result.cl_gl) +
    normalization[3] * sum(abs2, result.cl_ll)

function make_realistic_objective(workload, x0)
    spectrum = prepare_realistic_spectrum(workload, x0)
    weights = Blast.compute_w(workload.weights_template, spectrum)
    result = project_realistic(workload, spectrum, weights)
    normalization = (inv(sum(abs2, result.cl_gg)),
                     inv(sum(abs2, result.cl_gl)),
                     inv(sum(abs2, result.cl_ll)))
    functional = let workload=workload, normalization=normalization
        x -> begin
            spectrum_x = prepare_realistic_spectrum(workload, x)
            weights_x = Blast.compute_w(workload.weights_template, spectrum_x)
            return objective(project_realistic(workload, spectrum_x, weights_x), normalization)
        end
    end
    # Independent ForwardDiff reference of the SAME 3-amplitude objective: the
    # stored background is lifted exactly to zero-partial Duals so primordial
    # coefficients (PNG's cϕ) promote to Dual like the Pmm-derived ones. The
    # fixed probe kernels were built on the Float64 background; projection uses it.
    functional_fd = let workload=workload, normalization=normalization
        x -> begin
            bg_x = lift_background(workload.bg, first(x))
            spectrum_x = prepare_realistic_spectrum(workload, x; bg=bg_x)
            weights_x = Blast.compute_w(workload.weights_template, spectrum_x)
            return objective(project_realistic(workload, spectrum_x, weights_x), normalization)
        end
    end
    workspace = Blast.allocate_compute_w(workload.weights_template, spectrum)
    reusable = let workload=workload, workspace=workspace, normalization=normalization
        x -> begin
            spectrum_x = prepare_realistic_spectrum(workload, x)
            Blast.compute_w!(workspace, spectrum_x)
            return objective(project_realistic(workload, spectrum_x, workspace), normalization)
        end
    end
    return (; functional, functional_fd, reusable, workspace, normalization)
end

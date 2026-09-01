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

    galaxy = Blast.GalaxyClustering(
        δ=Blast.NumberCounts(nz=nz_g, z=z_g, bias=bias),
        RSD=Blast.RedshiftSpaceDistortions(nz=nz_g, z=z_g, growth_rate=bg.f),
        μ=Blast.MagnificationBias(nz=nz_g, z=z_g, s=slope),
    )
    lensing = Blast.WeakLensing(
        γ=Blast.CosmicShear(nz=nz_s, z=z_s),
        IA=Blast.IntrinsicAlignment(nz=nz_s, z=z_s, A=1.72),
    )
    Blast.evaluate_components!(galaxy, bg)
    Blast.evaluate_components!(lensing, bg)
    weights_template, plans = Blast.SetUp(galaxy, lensing)
    return (; ell, pk, pk_limber_lin, pk_limber_nonlin, bg, galaxy, lensing,
            weights_template, plans)
end

function prepare_realistic_spectrum(workload, x)
    return Blast.prepare_pk_workspace(
        workload.plans,
        x[1] .* workload.pk,
        x[2] .* workload.pk_limber_lin,
        x[3] .* workload.pk_limber_nonlin,
        workload.bg,
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

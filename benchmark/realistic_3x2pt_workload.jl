using Blast
using NPZ

function build_realistic_no_png_workload(root::AbstractString, fftw_threads::Integer)
    Blast.FFTW.set_num_threads(fftw_threads)

    nz_data = npzread(joinpath(root, "data", "lsst_nz.npz"))
    nz_g = nz_data["nz_clustering"]
    z_g = nz_data["z_clustering"]
    nz_s = nz_data["nz_lensing"]
    z_s = nz_data["z_lensing"]
    z_g == z_s || error("clustering and lensing n(z) grids differ")

    ell = npzread(joinpath(root, "data", "ell_n5k.npz"))
    pk = npzread(joinpath(root, "data", "pk_nonlimber_161.npz"))["pk"]
    pk_limber_lin = npzread(joinpath(root, "data", "pk_limber_lin.npz"))["pk"]
    pk_limber_nonlin = npzread(joinpath(root, "data", "pk_limber_nonlin.npz"))["pk"]

    cosmo = Blast.w0waCDM()
    bg = Blast.Background(cosmo)

    bias_values = [
        1.376695, 1.451179, 1.528404, 1.607983, 1.689579,
        1.772899, 1.857700, 1.943754, 2.030887, 2.118943,
    ]
    bias = ones(size(nz_g, 1), length(bg.z)) .* bias_values
    slope_values = [
        0.1648, 0.2496, 0.2708, 0.3300, 0.3880,
        0.2960, 0.3580, 0.3960, 0.4320, 0.5680,
    ]
    slope = slope_values .* ones(size(nz_g, 1), length(bg.z))

    galaxy = Blast.GalaxyClustering(
        δ=Blast.NumberCounts(nz=nz_g, z=z_g, bias=bias),
        RSD=Blast.RedshiftSpaceDistortions(nz=nz_g, z=z_g, growth_rate=bg.f),
        μ=Blast.MagnificationBias(nz=nz_g, z=z_g, s=slope),
    )
    Blast.evaluate_components!(galaxy, bg)

    Ωm = Blast.get_Ωm(bg.cosmo)
    A_IA = @. -1.72 * 0.0134 * Ωm / bg.D
    A_IA_matrix = ones(size(nz_s, 1), length(bg.z)) .* A_IA'
    lensing = Blast.WeakLensing(
        γ=Blast.CosmicShear(nz=nz_s, z=z_s),
        IA=Blast.IntrinsicAlignment(nz=nz_s, z=z_s, A_IA=A_IA_matrix),
    )
    Blast.evaluate_components!(lensing, bg)

    W_template, plans = Blast.SetUp(galaxy, lensing)
    return (;
        root, ell, pk, pk_limber_lin, pk_limber_nonlin,
        bg, galaxy, lensing, W_template, plans,
    )
end

function prepare_scaled_spectrum(workload, x)
    return Blast.prepare_pk_workspace(
        workload.plans,
        x[1] .* workload.pk,
        x[2] .* workload.pk_limber_lin,
        x[3] .* workload.pk_limber_nonlin,
        workload.bg,
    )
end

function project_all(workload, spectrum, W)
    cl_gg = Blast.get_Cℓ(
        workload.ell, workload.galaxy, spectrum, W,
        workload.bg, workload.plans,
    )
    cl_ss = Blast.get_Cℓ(
        workload.ell, workload.lensing, spectrum, W,
        workload.bg, workload.plans,
    )
    cl_gs = Blast.get_Cℓ(
        workload.ell, workload.galaxy, workload.lensing, spectrum, W,
        workload.bg, workload.plans,
    )
    return (; cl_gg, cl_gs, cl_ss)
end

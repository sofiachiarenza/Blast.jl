using Blast

"""Construct the exact ACE cosmology represented by Blast's legacy defaults.

The physical-density fields are written explicitly so this workload survives
the removal of Blast's `w0waCDM(H0, Ωm, Ωb, ...)` convenience constructor.
This is a characterization cosmology, not a statement that the legacy
big-Ω constructor had the desired total-matter semantics.
"""
function characterization_cosmology()
    h = 0.6727
    Ωb_legacy = 0.0492
    Ωc_legacy = 0.3156 - Ωb_legacy
    As = 2.12107e-9
    return Blast.w0waCDMCosmology(
        ln10Aₛ=log(1.0e10 * As),
        nₛ=0.9645,
        h=h,
        ωb=Ωb_legacy * h^2,
        ωc=Ωc_legacy * h^2,
        ωk=0.0,
        mν=0.06,
        w0=-1.0,
        wa=0.0,
    )
end

function characterization_source_distributions()
    z = collect(LinRange(0.0, 3.6, 101))
    profile(z, center, width) = z^2 * exp(-0.5 * ((z - center) / width)^2)
    nz_g = reduce(vcat, (permutedims(profile.(z, center, width)) for
                         (center, width) in ((0.55, 0.18), (1.05, 0.24))))
    nz_s = reduce(vcat, (permutedims(profile.(z, center, width)) for
                         (center, width) in ((0.75, 0.22), (1.35, 0.30))))
    return (; z, nz_g, nz_s)
end

function characterization_pmm()
    pmm(k, z) = 2.0e4 * (k / 0.2)^0.96 /
                (1 + (k / 0.2)^2)^2 / (1 + z)^2
    nonlinear_boost(k, z) = 1 + 0.45 * (k / (0.35 + k))^2 / (1 + z)

    k_nonlimber = 10.0 .^ Blast.k_cheb
    pk_nonlimber = [pmm(k, z) for z in Blast.z_lin, k in k_nonlimber]

    k_limber = 10.0 .^ Blast.k_limber
    pk_limber_linear = [pmm(k, z) for z in Blast.z_cheb, k in k_limber]
    pk_limber_nonlinear = [pmm(k, z) * nonlinear_boost(k, z)
                           for z in Blast.z_cheb, k in k_limber]
    return (; pk_nonlimber, pk_limber_linear, pk_limber_nonlinear)
end

function build_characterization_probes(bg)
    sources = characterization_source_distributions()
    n_bg = length(bg.z)

    bias = Matrix{eltype(bg.z)}(undef, 2, n_bg)
    slope = similar(bias)
    @views begin
        bias[1, :] .= 1.20 .+ 0.10 .* bg.z
        bias[2, :] .= 1.55 .+ 0.08 .* bg.z
        slope[1, :] .= 0.24 .+ 0.015 .* bg.z
        slope[2, :] .= 0.34 .+ 0.010 .* bg.z
    end

    galaxy = Blast.GalaxyClustering(
        δ=Blast.NumberCounts(nz=sources.nz_g, z=sources.z, bias=bias),
        RSD=Blast.RedshiftSpaceDistortions(
            nz=sources.nz_g,
            z=sources.z,
            growth_rate=bg.f,
        ),
        μ=Blast.MagnificationBias(nz=sources.nz_g, z=sources.z, s=slope),
        PNG=Blast.PrimordialNonGaussianity(
            nz=sources.nz_g,
            z=sources.z,
            bias=bias,
            f_NL=7.0,
            p=1.0,
        ),
    )
    lensing = Blast.WeakLensing(
        γ=Blast.CosmicShear(nz=sources.nz_s, z=sources.z),
        IA=Blast.IntrinsicAlignment(nz=sources.nz_s, z=sources.z, A=1.72),
    )
    cmb = Blast.CMB(
        κ=Blast.CMBLensing(),
        ISW=Blast.IntegratedSachsWolfe(growth_rate=bg.f),
    )

    galaxy = Blast.prepare_probe(galaxy, bg)
    lensing = Blast.prepare_probe(lensing, bg)
    cmb = Blast.prepare_probe(cmb, bg)
    return (; galaxy, lensing, cmb)
end

function project_characterization_spectra(state, spectrum, weights)
    (; ell, bg, probes, plans) = state
    (; galaxy, lensing, cmb) = probes
    cl_gg = Blast.get_Cℓ(ell, galaxy, spectrum, weights, bg, plans)
    cl_gl = Blast.get_Cℓ(ell, galaxy, lensing, spectrum, weights, bg, plans)
    cl_gc = Blast.get_Cℓ(ell, galaxy, cmb, spectrum, weights, bg, plans)
    cl_ll = Blast.get_Cℓ(ell, lensing, spectrum, weights, bg, plans)
    cl_lc = Blast.get_Cℓ(ell, lensing, cmb, spectrum, weights, bg, plans)
    return (; cl_gg, cl_gl, cl_gc, cl_ll, cl_lc)
end

function build_characterization_state()
    cosmo = characterization_cosmology()
    bg = Blast.Background(cosmo)
    probes = build_characterization_probes(bg)
    pmm = characterization_pmm()
    weights_template, plans = Blast.SetUp(probes.galaxy, probes.lensing, probes.cmb)
    ell = copy(Blast.full_ℓ_range)
    return (; cosmo, bg, probes, pmm, weights_template, plans, ell)
end

function prepare_characterization_spectrum(state; scales=(1.0, 1.0, 1.0))
    return Blast.prepare_pk_workspace(
        state.plans,
        scales[1] .* state.pmm.pk_nonlimber,
        scales[2] .* state.pmm.pk_limber_linear,
        scales[3] .* state.pmm.pk_limber_nonlinear,
        state.bg,
    )
end

function run_characterization(state=build_characterization_state(); scales=(1.0, 1.0, 1.0))
    spectrum = prepare_characterization_spectrum(state; scales)
    weights = Blast.compute_w(state.weights_template, spectrum)
    spectra = project_characterization_spectra(state, spectrum, weights)
    return (; state..., spectrum, weights, spectra)
end

function _three_indices(n)
    return unique((1, cld(n, 2), n))
end

function _coefficient_sample(coefficients)
    ndims(coefficients) == 1 && return copy(coefficients)
    indices = ntuple(d -> d == 1 ? axes(coefficients, d) : _three_indices(size(coefficients, d)),
                     ndims(coefficients))
    return copy(coefficients[indices...])
end

function _weight_sample(component)
    w = component.w
    indices = (axes(w, 1), _three_indices(size(w, 2)), _three_indices(size(w, 3)))
    return copy(w[indices...])
end

"""Return the arrays frozen by the legacy characterization fixture.

Large coefficient and projected-matter tensors are represented by deterministic
slices. Final spectra and all background/component kernels are stored in full.
"""
function characterization_arrays(result)
    (; bg, probes, spectrum, weights, spectra) = result
    arrays = Dict{String, Any}(
        "background_chi" => bg.χ,
        "background_z" => bg.z,
        "background_H" => bg.H,
        "background_D" => bg.D,
        "background_f" => bg.f,
        "galaxy_delta_nz_norm" => probes.galaxy.δ.nz_norm,
        "galaxy_delta_kernel" => probes.galaxy.δ.Kernel,
        "galaxy_rsd_kernel" => probes.galaxy.RSD.Kernel,
        "galaxy_magnification_kernel" => probes.galaxy.μ.Kernel,
        "galaxy_png_kernel" => probes.galaxy.PNG.Kernel,
        "lensing_shear_nz_norm" => probes.lensing.γ.nz_norm,
        "lensing_shear_kernel" => probes.lensing.γ.Kernel,
        "lensing_ia_kernel" => probes.lensing.IA.Kernel,
        "cmb_lensing_kernel" => probes.cmb.κ.Kernel,
        "cmb_isw_kernel" => probes.cmb.ISW.Kernel,
        "power_delta_limber" => spectrum.Pδ_limber,
        "power_delta_limber_correction" => spectrum.ΔP_limber,
        "coeff_phi_tt_sample" => _coefficient_sample(spectrum.cϕTT.coefs),
        "coeff_phi_t_sample" => _coefficient_sample(spectrum.cϕT.coefs),
        "coeff_phi" => _coefficient_sample(spectrum.cϕ.coefs),
        "ell" => result.ell,
        "cl_gg" => spectra.cl_gg,
        "cl_gl" => spectra.cl_gl,
        "cl_gc" => spectra.cl_gc,
        "cl_ll" => spectra.cl_ll,
        "cl_lc" => spectra.cl_lc,
    )
    for field in fieldnames(typeof(weights))
        component = getfield(weights, field)
        isnothing(component) && continue
        arrays["weight_$(field)_sample"] = _weight_sample(component)
    end
    return arrays
end

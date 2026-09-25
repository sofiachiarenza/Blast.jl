"""
Generate the ARFCAMB reference fixture for the Blast ARF validation.

Runs ARFCAMB (https://github.com/chmATiac/ARFCAMB, branch version_March11th_2025)
for a single Gaussian bin, with two builds:

  * patched (`patches/arfcamb_no_gauge_term.patch`): the Newtonian-gauge
    correction (comoving_density_ev - 3 adotoa) sigma/k is removed from the
    counts density source, leaving b * synchronous-gauge delta_c, which is what
    Blast's NumberCounts kernel models. This is the strict reference.
  * stock: unmodified ARFCAMB, kept to show the size of that relativistic term.

and writes test/reference_fixtures/arfcamb/gaussian_z1.npz with the linear
matter P(k, z) on Blast's grids (so Blast consumes CAMB's own power spectrum and
only projections/kernels are tested) and the C_ell pieces

    dd, dd_rsd, dd_mu, gg            counts auto, SourceTerms toggled
    kd, kd_rsd, kd_mu, kg            CMB kappa x counts, same toggles
    zz, dz, kz                       ARF (density only) auto / x counts / x kappa
    zz_full, dz_full, kz_full        ARF with density+RSD+lensing (informative:
                                     Blast's ARF kernel is density-only)
    kk                               CMB kappa auto (normalizes kappa crosses)
    <piece>_stock                    same pieces from the stock build

Usage (from the repo root, inside the `arfcamb` conda env; see README.md):

    julia --project=validation/arfcamb validation/arfcamb/export_blast_grids.jl
    python validation/arfcamb/generate_camb_reference.py \\
        --patched ~/Desktop/PhD/ARFCAMB_nogauge --stock ~/Desktop/PhD/ARFCAMB

Each build runs in its own subprocess, since two compiled CAMB libraries cannot
share one Python process.
"""
import argparse
import os
import subprocess
import sys
import tempfile
import time
import tomllib

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "..", "test", "reference_fixtures", "arfcamb", "gaussian_z1.npz")
PATCH_MARKER = "[BLAST-VALIDATION] Newtonian-gauge correction"

CL_KEYS = ("dd", "dd_rsd", "dd_mu", "gg", "kd", "kd_rsd", "kd_mu", "kg",
           "zz", "dz", "kz", "zz_full", "dz_full", "kz_full", "kk")

camb = None  # imported from the requested build by import_camb()


def import_camb(build_dir):
    global camb
    sys.path.insert(0, os.path.abspath(os.path.expanduser(build_dir)))
    import camb as _camb
    camb = _camb
    return camb


def is_patched(build_dir):
    with open(os.path.join(os.path.expanduser(build_dir), "fortran", "equations.f90")) as f:
        return PATCH_MARKER in f.read()


def load_config():
    with open(os.path.join(HERE, "config.toml"), "rb") as f:
        return tomllib.load(f)


def window_table(cfg):
    w = cfg["window"]
    z = np.linspace(w["zmin"], w["zmax"], w["nz"])
    W = np.exp(-0.5 * ((z - w["z0"]) / w["sigma"]) ** 2) / (np.sqrt(2 * np.pi) * w["sigma"])
    # Same discrete estimator ARFCAMB uses internally (SourceWindows.f90, SetTable).
    zbar = np.sum(z * W) / np.sum(W)
    return z, W, zbar


def base_params(cfg, lmax):
    c = cfg["cosmology"]
    pars = camb.set_params(
        H0=100 * c["h"], ombh2=c["ombh2"], omch2=c["omch2"], mnu=c["mnu"],
        num_massive_neutrinos=1, omk=0.0, As=c["As"], ns=c["ns"],
        w=c["w0"], wa=c["wa"], lmax=lmax + 200, lens_potential_accuracy=1,
    )
    pars.NonLinear = camb.model.NonLinear_none
    pars.Want_CMB = True
    pars.Want_Cls = True
    pars.Want_cl_2D_array = True  # needed for window cross-spectra
    pars.DoLensing = False  # we want unlensed phi and source windows, not lensed CMB
    # Converged: AccuracyBoost 2 -> 3, lAccuracyBoost 2 -> 5, Transfer.high_precision
    # and accurate_massive_neutrinos all change C_ell by < 0.02% at every ell.
    pars.Accuracy.AccuracyBoost = 2
    pars.Accuracy.lSampleBoost = 50  # every ell computed exactly, no interpolation
    pars.Accuracy.lAccuracyBoost = 2
    return pars


def set_windows(pars, cfg, z, W):
    b, s = cfg["tracer"]["bias"], cfg["tracer"]["s"]
    counts = camb.sources.SplinedSourceWindow(source_type="counts", z=z, W=W, bias=b, dlog10Ndm=s)
    arf = camb.sources.SplinedSourceWindow(source_type="arf", z=z, W=W, bias=b, dlog10Ndm=s)
    pars.SourceWindows = [counts, arf]  # W1 = counts, W2 = ARF


def set_terms(pars, density, redshift, lensing):
    st = pars.SourceTerms
    st.limber_windows = False
    st.limber_phi_lmin = 10000  # keep phi x window non-Limber over our whole ell range
    st.counts_density = density
    st.counts_redshift = redshift
    st.counts_lensing = lensing
    st.counts_velocity = False
    st.counts_radial = False
    st.counts_timedelay = False
    st.counts_ISW = False
    st.counts_potential = False
    st.counts_evolve = False


def run_cls(cfg, z, W, lmax, density, redshift, lensing):
    pars = base_params(cfg, lmax)
    set_windows(pars, cfg, z, W)
    set_terms(pars, density, redshift, lensing)
    t0 = time.time()
    res = camb.get_results(pars)
    cls = res.get_source_cls_dict(raw_cl=True, lmax=lmax)
    print(f"  CAMB run (density={density}, redshift={redshift}, lensing={lensing}) "
          f"took {time.time() - t0:.1f}s", flush=True)
    ell = np.arange(lmax + 1)
    phi_to_kappa = ell * (ell + 1) / 2.0  # kappa_lm = l(l+1)/2 phi_lm
    return {
        "11": cls["W1xW1"], "12": cls["W1xW2"], "22": cls["W2xW2"],
        "k1": phi_to_kappa * cls["PxW1"], "k2": phi_to_kappa * cls["PxW2"],
        "kk": phi_to_kappa**2 * cls["PxP"],
    }


def all_cls(cfg, z, W):
    lmin, lmax = cfg["ell"]["lmin"], cfg["ell"]["lmax"]
    d = run_cls(cfg, z, W, lmax, True, False, False)
    dr = run_cls(cfg, z, W, lmax, True, True, False)
    dm = run_cls(cfg, z, W, lmax, True, False, True)
    full = run_cls(cfg, z, W, lmax, True, True, True)
    out = dict(
        dd=d["11"], dd_rsd=dr["11"], dd_mu=dm["11"], gg=full["11"],
        kd=d["k1"], kd_rsd=dr["k1"], kd_mu=dm["k1"], kg=full["k1"],
        zz=d["22"], dz=d["12"], kz=d["k2"],
        zz_full=full["22"], dz_full=full["12"], kz_full=full["k2"],
        kk=d["kk"],
    )
    return {k: np.asarray(v)[lmin:lmax + 1] for k, v in out.items()}


def eval_on_grid(interp, z, k):
    """interp.P on (z, k) in whatever order Blast stores them (Chebyshev grids are descending)."""
    iz, ik = np.argsort(z), np.argsort(k)
    pk_sorted = interp.P(z[iz], k[ik])
    pk = np.empty_like(pk_sorted)
    pk[np.ix_(iz, ik)] = pk_sorted
    return pk


def matter_pk(cfg, grids):
    """Linear P(k, z) in Mpc^3 on Blast's (z, k) grids, k in Mpc^-1."""
    pars = base_params(cfg, 100)
    pars.WantTransfer = True  # (this ARFCAMB version rejects Want_CMB=False without windows)
    k_nl = 10 ** grids["log10k_cheb"]
    k_l = 10 ** grids["log10k_limber"]
    kmax = 1.2 * max(k_nl.max(), k_l.max())
    pars.set_matter_power(redshifts=np.linspace(0.0, 3.7, 38)[::-1], kmax=kmax, nonlinear=False)
    res = camb.get_results(pars)
    out = {}
    for var, tag in (("delta_tot", ""), ("delta_nonu", "_cb")):
        interp = res.get_matter_power_interpolator(
            nonlinear=False, hubble_units=False, k_hunit=False, extrap_kmax=10 * kmax,
            var1=var, var2=var)
        out["pk" + tag] = eval_on_grid(interp, grids["z_lin"], k_nl)
        out["pk_limber" + tag] = eval_on_grid(interp, grids["z_cheb"], k_l)
    return out


def worker(build_dir, out_path, with_pk):
    """Run every CAMB computation with one build and save the raw results."""
    import_camb(build_dir)
    print(f"camb {camb.__version__} from {camb.__file__}", flush=True)
    cfg = load_config()
    z, W, _ = window_table(cfg)
    res = all_cls(cfg, z, W)
    if with_pk:
        print("  P(k, z) on Blast grids ...", flush=True)
        res.update(matter_pk(cfg, np.load(os.path.join(HERE, "grids.npz"))))
    np.savez(out_path, **res)


def run_worker(build_dir, out_path, with_pk):
    cmd = [sys.executable, "-u", __file__, "--worker", build_dir, out_path] + (["--with-pk"] if with_pk else [])
    subprocess.run(cmd, check=True)
    return dict(np.load(out_path))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--patched", help="ARFCAMB build with patches/arfcamb_no_gauge_term.patch applied")
    ap.add_argument("--stock", help="unmodified ARFCAMB build")
    ap.add_argument("--worker", nargs=2, metavar=("BUILD", "OUT"), help=argparse.SUPPRESS)
    ap.add_argument("--with-pk", action="store_true", help=argparse.SUPPRESS)
    args = ap.parse_args()

    if args.worker:
        worker(*args.worker, args.with_pk)
        return
    if not (args.patched and args.stock):
        ap.error("both --patched and --stock are required")
    if not is_patched(args.patched):
        sys.exit(f"{args.patched} does not carry patches/arfcamb_no_gauge_term.patch")
    if is_patched(args.stock):
        sys.exit(f"{args.stock} is patched; --stock must be unmodified ARFCAMB")

    cfg = load_config()
    z, W, zbar = window_table(cfg)
    print(f"zbar = {zbar:.10f}")

    with tempfile.TemporaryDirectory() as tmp:
        print("Patched build (strict reference) ...")
        patched = run_worker(args.patched, os.path.join(tmp, "patched.npz"), with_pk=True)
        print("Stock build (informative) ...")
        stock = run_worker(args.stock, os.path.join(tmp, "stock.npz"), with_pk=False)

    lmin, lmax = cfg["ell"]["lmin"], cfg["ell"]["lmax"]
    fixture = dict(
        ell=np.arange(lmin, lmax + 1), z=z, nz=W, zbar=np.array([zbar]),
        bias=np.array([cfg["tracer"]["bias"]]), s=np.array([cfg["tracer"]["s"]]),
        pk=patched["pk"], pk_limber_lin=patched["pk_limber"], pk_limber_nonlin=patched["pk_limber"],
        pk_cb=patched["pk_cb"], pk_limber_cb=patched["pk_limber_cb"],
    )
    for key in CL_KEYS:
        fixture[key] = patched[key]
        fixture[key + "_stock"] = stock[key]

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    np.savez(OUT, **fixture)
    print(f"Wrote {os.path.normpath(OUT)}")


if __name__ == "__main__":
    main()

"""
ARFCAMB side of the Blast fNL validation (issue #38).

Runs the no-gauge ARFCAMB build (validation/arfcamb/patches/arfcamb_no_gauge_term.patch)
with counts density only, so each window sources b * delta_c plus ARFCAMB's
built-in local-f_NL scale-dependent bias, (b - 1) 3 delta_c f_NL Omega_m H0^2 /
(c^2 k^2 T D) written through CAMB's own phi(k, a). Every bin pair is computed
at f_NL = 0, +f_NL and -f_NL. Writes camb_png.npz with

    ell, z, nz                 multipoles, shared n(z) table (nbins, nz)
    cl_0, cl_p, cl_m           raw C_ell, shape (nbins, nbins, nell)
    pk_cb, pk_limber_cb        linear P_cb(k, z) on Blast's grids (Blast input)
    z_swift, k_swift, pk_swift the same P_cb on SwiftCl's (z, k) grid

Usage, from the repo root inside the `arfcamb` conda env:

    python validation/png_external/generate_camb_png.py --build ~/Desktop/PhD/ARFCAMB_nogauge
"""
import argparse
import os
import sys
import time
import tomllib

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "arfcamb"))
import generate_camb_reference as arf  # noqa: E402  (shared CAMB settings and helpers)

OUT = os.path.join(HERE, "camb_png.npz")
K_SWIFT = np.logspace(-5, 3, 700)  # SwiftCl's k grid in the paper notebooks


def load_config():
    with open(os.path.join(HERE, "config.toml"), "rb") as f:
        return tomllib.load(f)


def window_table(cfg):
    w = cfg["window"]
    z = np.linspace(w["zmin"], w["zmax"], w["nz"])
    nz = np.array([np.exp(-0.5 * ((z - b["z0"]) / b["sigma"]) ** 2) / (np.sqrt(2 * np.pi) * b["sigma"])
                   for b in cfg["bins"]])
    return z, nz


def run_cls(cfg, z, nz, f_nl):
    lmin, lmax = cfg["ell"]["lmin"], cfg["ell"]["lmax"]
    pars = arf.base_params(cfg, lmax)
    pars.f_NL = f_nl
    pars.SourceWindows = [
        arf.camb.sources.SplinedSourceWindow(source_type="counts", z=z, W=W, bias=b["bias"], dlog10Ndm=0.0)
        for W, b in zip(nz, cfg["bins"])
    ]
    arf.set_terms(pars, density=True, redshift=False, lensing=False)
    t0 = time.time()
    cls = arf.camb.get_results(pars).get_source_cls_dict(raw_cl=True, lmax=lmax)
    print(f"  f_NL = {f_nl:+g}: {time.time() - t0:.1f}s", flush=True)
    n = len(nz)
    return np.array([[cls[f"W{i + 1}xW{j + 1}"][lmin:lmax + 1] for j in range(n)] for i in range(n)])


def matter_pk(cfg, grids, z_swift):
    """Linear P_cb(k, z) in Mpc^3, k in Mpc^-1, on Blast's grids and SwiftCl's grid."""
    pars = arf.base_params(cfg, 100)
    pars.WantTransfer = True
    k_nl = 10 ** grids["log10k_cheb"]
    k_l = 10 ** grids["log10k_limber"]
    kmax = 1.2 * max(k_nl.max(), k_l.max())
    pars.set_matter_power(redshifts=np.linspace(0.0, 3.7, 38)[::-1], kmax=kmax, nonlinear=False)
    interp = arf.camb.get_results(pars).get_matter_power_interpolator(
        nonlinear=False, hubble_units=False, k_hunit=False, extrap_kmax=10 * K_SWIFT.max(),
        var1="delta_nonu", var2="delta_nonu")
    return dict(
        pk_cb=arf.eval_on_grid(interp, grids["z_lin"], k_nl),
        pk_limber_cb=arf.eval_on_grid(interp, grids["z_cheb"], k_l),
        pk_swift=interp.P(z_swift, K_SWIFT),
    )


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--build", required=True, help="ARFCAMB build carrying the no-gauge patch")
    args = ap.parse_args()
    if not arf.is_patched(args.build):
        sys.exit(f"{args.build} does not carry validation/arfcamb/patches/arfcamb_no_gauge_term.patch")

    arf.import_camb(args.build)
    print(f"camb {arf.camb.__version__} from {arf.camb.__file__}", flush=True)
    cfg = load_config()
    z, nz = window_table(cfg)
    f = cfg["tracer"]["f_NL"]

    out = dict(ell=np.arange(cfg["ell"]["lmin"], cfg["ell"]["lmax"] + 1), z=z, nz=nz,
               z_swift=z, k_swift=K_SWIFT)
    out["cl_0"] = run_cls(cfg, z, nz, 0.0)
    out["cl_p"] = run_cls(cfg, z, nz, f)
    out["cl_m"] = run_cls(cfg, z, nz, -f)
    print("  P(k, z) ...", flush=True)
    out.update(matter_pk(cfg, np.load(os.path.join(HERE, "..", "arfcamb", "grids.npz")), z))
    np.savez(OUT, **out)
    print(f"Wrote {os.path.relpath(OUT)}")


if __name__ == "__main__":
    main()

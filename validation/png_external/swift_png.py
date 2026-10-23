"""
SwiftCl side of the Blast fNL validation (issue #38). Needs swiftcl >= 0.1.3,
which restores the growth factor missing from 0.1.2's PNG kernel.

SwiftCl integrates over (chi_1, chi_2) with P(k; z1, z2) = P(k, z_eff) D(z1, k)
D(z2, k), each kernel carrying its own D, so the transfer function sits on the
non-PNG leg without any folding. Inputs come from the other two scripts so all
three codes share one P(k):

    camb_png.npz       CAMB linear P_cb(k, z) on the shared n(z) table
    swift_inputs.npz   Blast's chi(z), H(z) and P_Phi(k) (Phi = 3/5 zeta)

T(z, k) = sqrt(P(k, z) / P_Phi(k)) is delta_m / Phi, Blast's get_Tm. SwiftCl's
w_fnl carries an explicit 2/3, so b_fNL = (3/2) b_Phi = 3 delta_c (b - p)
reproduces Blast's b_Phi f_NL Phi field (Dalal et al. 2008 normalization).

Usage, from the repo root inside the `swift_env` conda env:

    python validation/png_external/swift_png.py
"""
import os
import tomllib

import jax

jax.config.update("jax_platform_name", "cpu")
jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp  # noqa: E402
import numpy as np  # noqa: E402
import swiftcl  # noqa: E402
from swiftcl.cl import ClComp  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, "swift_cache")  # SwiftCl stores its FFTLog coefficients under <path>/SwiftCl/g_r
DELTA_C = 1.686


def main():
    with open(os.path.join(HERE, "config.toml"), "rb") as f:
        cfg = tomllib.load(f)
    fx = np.load(os.path.join(HERE, "camb_png.npz"))
    sw = np.load(os.path.join(HERE, "swift_inputs.npz"))
    print(f"swiftcl {swiftcl.__version__ if hasattr(swiftcl, '__version__') else ''} from {swiftcl.__file__}")

    z, nz = sw["z"], fx["nz"]
    k, pk = sw["k"], fx["pk_swift"]  # pk: (nz, nk)
    z_eff = 1.0
    P_eff = np.array([np.interp(z_eff, z, pk[:, i]) for i in range(len(k))])
    D = np.sqrt(pk / P_eff)                      # D(z, k) relative to z_eff
    T = np.sqrt(pk / sw["P_phi"][None, :])       # delta_m / Phi on the same table

    c = cfg["cosmology"]
    h = c["h"]
    H0 = 100 * h
    Om = (c["ombh2"] + c["omch2"] + c["mnu"] / 93.14) / h**2  # unused by the g and f_NL kernels
    tr = cfg["tracer"]
    f = tr["f_NL"]
    bias = [b["bias"] for b in cfg["bins"]]
    b_fnl = [1.5 * 2 * DELTA_C * (b - tr["p"]) for b in bias]

    ell = jnp.array(fx["ell"])
    nbins = nz.shape[0]
    arr = {key: jnp.array(val) for key, val in
           (("chi", sw["chi"]), ("H", sw["H"]), ("D", D), ("T", T), ("P", P_eff))}
    out = {"ell": np.asarray(fx["ell"])}
    for name, f_nl, contr in (("c0", 0.0, "g"), ("cp", f, "g,f_NL"), ("cm", -f, "g,f_NL")):
        cl = np.zeros((len(ell), nbins, nbins))
        for i in range(nbins):
            for j in range(i, nbins):
                n1 = jnp.array(np.column_stack([z, nz[i]]))
                n2 = jnp.array(np.column_stack([z, nz[j]]))
                comp = ClComp(z1=[z[0], z[-1]], z2=[z[0], z[-1]], contr1=contr, contr2=contr,
                              l=ell, k=jnp.array(k), jit=False, path=CACHE, N_interp=400)
                kw = dict(n1=n1, n2=n2, chis1=arr["chi"], chis2=arr["chi"], D1=arr["D"], D2=arr["D"],
                          P=arr["P"], H1=arr["H"], H2=arr["H"], H0=H0, O_m=Om, b_g1=bias[i], b_g2=bias[j])
                if contr != "g":
                    kw.update(b_fNL1=b_fnl[i], b_fNL2=b_fnl[j], f_NL1=f_nl, f_NL2=f_nl, T1=arr["T"], T2=arr["T"])
                cl[:, i, j] = cl[:, j, i] = np.asarray(comp.C_l(**kw))
                print(f"  {name} bins ({i + 1},{j + 1}) done", flush=True)
        out[name] = cl
    path = os.path.join(HERE, "swift_png.npz")
    np.savez(path, **out)
    print(f"Wrote {os.path.relpath(path)}")


if __name__ == "__main__":
    os.makedirs(os.path.join(CACHE, "SwiftCl", "g_r"), exist_ok=True)
    main()

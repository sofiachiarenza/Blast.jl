"""
Compare Blast's fNL pieces before and after the issue #38 fix with ARFCAMB_nogauge
and SwiftCl. Each code's spectra at f_NL = 0, +f, -f give

    D = C(0)                              density only (shared-input sanity check)
    A = (C(+f) - C(-f)) / (2 f)           delta x PNG + PNG x delta, per unit f_NL
    B = (C(+f) + C(-f) - 2 C(0)) / (2 f^2) PNG x PNG, per unit f_NL^2

The ordering fix acts on A only. Prints a table per bin pair and ell range and
writes png_comparison.png.

    python validation/png_external/compare_png.py
"""
import os
import tomllib

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RANGES = ((2, 10), (11, 50), (51, 200))


def pieces(path, f, transpose=False):
    d = np.load(path)
    get = (lambda k: np.moveaxis(d[k], -1, 0)) if transpose else (lambda k: d[k])  # -> (nell, i, j)
    c0, cp, cm = get("c0" if "c0" in d else "cl_0"), get("cp" if "cp" in d else "cl_p"), get("cm" if "cm" in d else "cl_m")
    return dict(D=c0, A=(cp - cm) / (2 * f), B=(cp + cm - 2 * c0) / (2 * f**2))


def main():
    with open(os.path.join(HERE, "config.toml"), "rb") as fh:
        cfg = tomllib.load(fh)
    f = cfg["tracer"]["f_NL"]
    ell = np.load(os.path.join(HERE, "camb_png.npz"))["ell"]

    codes = {"camb": pieces(os.path.join(HERE, "camb_png.npz"), f, transpose=True)}
    for name, file in (("swift", "swift_png.npz"), ("old", "blast_old.npz"), ("fixed", "blast_fixed.npz")):
        path = os.path.join(HERE, file)
        if os.path.exists(path):
            codes[name] = pieces(path, f)
    refs = [r for r in ("camb", "swift") if r in codes]
    nb = codes["camb"]["D"].shape[1]
    bins = [f"z{b['z0']:g}" for b in cfg["bins"]]

    def rel(x, ref, i, j, lo, hi):
        m = (ell >= lo) & (ell <= hi)
        # normalized by the reference autos, so crosses near zero stay well defined
        norm = np.sqrt(np.abs(ref[m, i, i] * ref[m, j, j]))
        return np.median(np.abs(x[m, i, j] - ref[m, i, j]) / norm)

    for piece, label in (("D", "density only"), ("A", "delta x PNG (fix acts here)"), ("B", "PNG x PNG")):
        print(f"\n=== {piece}: {label} — median |X - ref| / sqrt(ref_ii ref_jj) ===")
        cols = [(c, r) for r in refs for c in ("old", "fixed", "swift") if c in codes and c != r]
        print(f"{'pair':>12} {'ell':>8} " + " ".join(f"{c + ' vs ' + r:>16}" for c, r in cols))
        for i in range(nb):
            for j in range(i, nb):
                for lo, hi in RANGES:
                    vals = [rel(codes[c][piece], codes[r][piece], i, j, lo, hi) for c, r in cols]
                    print(f"{bins[i] + 'x' + bins[j]:>12} {f'{lo}-{hi}':>8} " + " ".join(f"{v:16.2e}" for v in vals))

    for c in ("old", "fixed"):
        if c in codes:
            A = codes[c]["A"]
            asym = np.max(np.abs(A - np.swapaxes(A, 1, 2))) / np.max(np.abs(A))
            print(f"\nBlast {c}: max |A_ij - A_ji| / max |A| = {asym:.2e}")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return
    pairs = [(i, j) for i in range(nb) for j in range(i, nb)]
    fig, axes = plt.subplots(2, len(pairs), figsize=(3.2 * len(pairs), 6), sharex=True, squeeze=False)
    styles = {"old": dict(ls="--", color="C3"), "fixed": dict(ls="-", color="C0"), "swift": dict(ls=":", color="C2")}
    for col, (i, j) in enumerate(pairs):
        for row, piece in enumerate(("A", "B")):
            ax = axes[row, col]
            ref = codes["camb"][piece]
            norm = np.sqrt(np.abs(ref[:, i, i] * ref[:, j, j]))
            for c in ("old", "fixed", "swift"):
                if c in codes:
                    ax.plot(ell, (codes[c][piece][:, i, j] - ref[:, i, j]) / norm, label=c, **styles[c])
            ax.axhline(0, color="k", lw=0.5)
            ax.set_xscale("log")
            ax.set_title(f"{piece}  {bins[i]}x{bins[j]}", fontsize=9)
            if col == 0:
                ax.set_ylabel(r"$(X - \mathrm{CAMB}) / \sqrt{\mathrm{CAMB}_{ii}\mathrm{CAMB}_{jj}}$", fontsize=8)
            if row == 1:
                ax.set_xlabel(r"$\ell$")
    axes[0, 0].legend(fontsize=8)
    fig.suptitle("Blast fNL pieces vs ARFCAMB_nogauge (A: δ×PNG, B: PNG×PNG)", fontsize=10)
    fig.tight_layout()
    out = os.path.join(HERE, "png_comparison.png")
    fig.savefig(out, dpi=130)
    print(f"\nWrote {os.path.relpath(out)}")


if __name__ == "__main__":
    main()

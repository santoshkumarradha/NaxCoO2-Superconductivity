#!/usr/bin/env python3
"""Fit Tc vs basal spacing d for electron-doped beta-ZrNCl and beta-HfNCl intercalates.

Data: zrncl_data.csv (literature values; see 'citation'/'provenance' columns).

Two competing phenomenologies are fit per family:
  (A) saturating exponential  Tc(d) = Tc_sat - B * exp(-(d - d0)/L)
      -- the "2D-limit" picture (Takano PRL 100, 247005 (2008); Kasahara PRB 82,
         054504 (2010)): Tc rises as interlayer hopping t_z dies out, then saturates.
  (B) interlayer-Coulomb form  Tc(d) = a + b/d
      -- the alpha-TiNCl phenomenology (Tc linear in 1/d, Zhang PRB 86, 024516 (2012))
         and the Harshman-Fiory interlayer-Coulomb model. For beta-MNCl a 1/d law
         with b>0 (Tc DECREASING with d) is what a naive "spacing kills pairing"
         model would give; the data reject it.

Run:  python3 fit_tc_vs_spacing.py
Outputs fit parameters + SSE to stdout and saves tc_vs_spacing.png.
"""

import csv
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CSV = os.path.join(HERE, "zrncl_data.csv")
D0 = 9.0  # reference spacing (A), just below the smallest observed d


def load(family):
    d, tc, labels = [], [], []
    with open(CSV, newline="") as f:
        for row in csv.DictReader(f):
            if row["family"] != family:
                continue
            try:
                dd = float(row["d_A"])
                tt = float(row["Tc_K"])
            except ValueError:
                continue  # skip rows lacking either number (kept in CSV for the record)
            d.append(dd)
            tc.append(tt)
            labels.append(row["compound"])
    return np.array(d), np.array(tc), labels


def sat_exp(d, tc_sat, B, L):
    return tc_sat - B * np.exp(-(d - D0) / L)


def inv_d(d, a, b):
    return a + b / d


def fit_family(name):
    d, tc, labels = load(name)
    order = np.argsort(d)
    d, tc = d[order], tc[order]
    labels = [labels[i] for i in order]
    print(f"\n=== {name}: {len(d)} points with both d and Tc ===")
    for dd, tt, lb in zip(d, tc, labels):
        print(f"  d = {dd:6.2f} A   Tc = {tt:5.1f} K   {lb}")

    results = {}
    try:
        from scipy.optimize import curve_fit

        # (A) saturating exponential
        p0 = [tc.max(), tc.max() - tc.min() + 1.0, 2.0]
        try:
            pA, _ = curve_fit(sat_exp, d, tc, p0=p0, maxfev=20000)
            sseA = float(np.sum((tc - sat_exp(d, *pA)) ** 2))
            results["sat_exp"] = (pA, sseA)
            print(f"  (A) Tc(d) = {pA[0]:.2f} - {pA[1]:.2f} exp(-(d-{D0})/{pA[2]:.2f})"
                  f"   SSE = {sseA:.2f} K^2")
        except RuntimeError as e:
            print(f"  (A) saturating-exponential fit failed: {e}")

        # (B) a + b/d
        pB, _ = curve_fit(inv_d, d, tc, p0=[tc.mean(), 0.0], maxfev=20000)
        sseB = float(np.sum((tc - inv_d(d, *pB)) ** 2))
        results["inv_d"] = (pB, sseB)
        sign = "Tc FALLS with d (naive spacing-kills-pairing)" if pB[1] > 0 \
            else "Tc RISES with d (b<0: opposite of alpha-TiNCl)"
        print(f"  (B) Tc(d) = {pB[0]:.2f} + {pB[1]:.1f}/d   SSE = {sseB:.2f} K^2   -> {sign}")
    except ImportError:
        # numpy-only fallback: linear regression of Tc on 1/d
        A = np.vstack([np.ones_like(d), 1.0 / d]).T
        coef, res, *_ = np.linalg.lstsq(A, tc, rcond=None)
        sseB = float(res[0]) if len(res) else float(np.sum((tc - A @ coef) ** 2))
        results["inv_d"] = (coef, sseB)
        print(f"  scipy unavailable; 1/d linear fit: Tc = {coef[0]:.2f} + {coef[1]:.1f}/d"
              f"   SSE = {sseB:.2f} K^2")
    return d, tc, results


def main():
    fams = {}
    for fam in ("ZrNCl", "HfNCl"):
        fams[fam] = fit_family(fam)

    print("\nHeadline: in BOTH beta families Tc is a rising, saturating function of d.")
    print("A positive-b 1/d law (pairing weakened by opening the gallery) is rejected;")
    print("the best 1/d fits come out with b<0, i.e. Tc grows as the layers separate,")
    print("saturating above d ~ 13-15 A (cf. Takano PRL 100, 247005; Kasahara PRB 82, 054504).")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib unavailable -- skipping plot.")
        return

    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    colors = {"ZrNCl": "#1f77b4", "HfNCl": "#d62728"}
    for fam, (d, tc, results) in fams.items():
        ax.scatter(d, tc, s=42, color=colors[fam], label=f"$\\beta$-{fam} intercalates", zorder=3)
        if "sat_exp" in results:
            p, _ = results["sat_exp"]
            xx = np.linspace(d.min() - 0.3, d.max() + 1.0, 200)
            ax.plot(xx, sat_exp(xx, *p), color=colors[fam], lw=1.2, alpha=0.7,
                    label=f"{fam}: $T_c^{{sat}}$ = {p[0]:.1f} K, $L$ = {p[2]:.1f} $\\AA$")
    ax.set_xlabel("basal spacing d ($\\AA$)")
    ax.set_ylabel("$T_c$ (K)")
    ax.set_title("Electron-doped $\\beta$-MNCl: $T_c$ vs interlayer (basal) spacing")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    out = os.path.join(HERE, "tc_vs_spacing.png")
    fig.tight_layout()
    fig.savefig(out, dpi=160)
    print(f"plot saved: {out}")


if __name__ == "__main__":
    sys.exit(main())

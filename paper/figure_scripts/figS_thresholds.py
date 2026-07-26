#!/usr/bin/env python
"""
Figure S (threshold ordering) - figS_thresholds.
The Landau stiffness alpha(c) of the alkali off-centring mode for the three
1x1 ions Li, Na, K.  Where alpha crosses zero the single well turns into a
double well; that critical spacing c* is the design threshold, and it scales
with the ion radius:  c*_Li <= 5.5 A  <  c*_Na = 6.40 A  <  c*_K = 6.75 A.
Opening the K threshold is the paper's testable prediction.

Data (all read-only, reusing the existing parses):
 * Li 1x1, Na 1x1 alpha(c) -> theory/results/well_fits_v4.csv
 * K  1x1 alpha(c)         -> theory/results_extra_analysis/summary_extra.json
                             (set E, parsed by analyze_extra.py from
                              runpod/results_extra)
 * zero crossings c*_Li/Na/K -> the same set-E summary.

House style: pastel fill for the alpha<0 double-well domain, physics labelled
inside; flat fills only, no gradients.  Vertical drop-lines mark each c*.
"""
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.interpolate import PchipInterpolator

import _style as S

S.use_house_style()

REPO = Path(__file__).resolve().parents[2]
fits = S.load_well_fits()
E = json.loads((REPO / "theory" / "results_extra_analysis"
                / "summary_extra.json").read_text())["E"]

C_K = "#6a4c93"   # sober violet, the new (K) ion

# alpha(c) series: (element, cell, colour, marker, label, c_data, alpha_data)
def wf(el):
    w = fits[(fits.element == el) & (fits.cell == "1x1")].sort_values("c")
    return w["c"].values, w["alpha"].values

cLi, aLi = wf("Li")
cNa, aNa = wf("Na")
cK = np.array([r["c"] for r in E["rows"]])
aK = np.array([r["alpha"] for r in E["rows"]])

SERIES = [
    ("Li", S.C_AQUA, "s", r"Li $1{\times}1$", cLi, aLi, E["alpha_Li55"]),
    ("Na", S.C_BLUE, "o", r"Na $1{\times}1$", cNa, aNa, E["cstar_Na"]),
    ("K",  C_K,      "D", r"K $1{\times}1$",  cK,  aK,  E["cstar_K"]),
]
CSTAR = {"Li": 5.5, "Na": E["cstar_Na"], "K": E["cstar_K"]}  # c*_Li<=5.5 bound

XLO, XHI = 5.4, 10.1
YLO, YHI = -1.15, 3.5

fig, ax = plt.subplots(figsize=(3.5, 3.45))

# double-well domain (alpha < 0): amber pastel, labelled inside
ax.axhspan(YLO, 0.0, color=S.FILL_AMBER, lw=0, zorder=0)
ax.text(8.5, -1.0, "double well  (soft mode)", ha="center", va="center",
        color=S.INK_AMBER, fontsize=8.0, zorder=1.5)
ax.text(7.6, 3.05, "single well  (stiff)", ha="left", va="center",
        color=S.C_SEC, fontsize=8.0, zorder=1.5)
ax.axhline(0, color=S.C_INK, lw=0.8, zorder=2)

cg = np.linspace(5.5, 9.9, 400)
for el, col, mk, lab, cc, aa, cstar in SERIES:
    p = PchipInterpolator(cc, aa)
    ax.plot(cg, p(cg), color=col, lw=1.4, zorder=3)
    ax.plot(cc, aa, mk, ms=6, color=col, mfc="white", mew=1.5, zorder=4,
            label=lab)
    # vertical drop-line from the alpha=0 axis to the bottom at c*
    cs = CSTAR[el]
    ax.plot([cs, cs], [YLO, 0.0], color=col, lw=1.1, ls=(0, (2, 2)),
            zorder=3.5)
    ax.plot(cs, 0.0, mk, ms=5.5, color=col, mfc=col, mew=0, zorder=5)

# c* labels, colour-coded, ordered Li<Na<K
ax.text(5.58, -1.06, r"$c^*_{\mathrm{Li}}\!\leq\!5.5$", color=S.C_AQUA,
        fontsize=7.2, ha="left", va="bottom", zorder=6)
ax.text(CSTAR["Na"] - 0.05, 0.12, r"$c^*_{\mathrm{Na}}$" "\n" r"$6.40$",
        color=S.C_BLUE, fontsize=7.2, ha="right", va="bottom", zorder=6)
ax.text(CSTAR["K"] + 0.07, 0.12, r"$c^*_{\mathrm{K}}$" "\n" r"$6.75$",
        color=C_K, fontsize=7.2, ha="left", va="bottom", zorder=6)

ax.set_xlim(XLO, XHI)
ax.set_ylim(YLO, YHI)
ax.set_xlabel(r"CoO$_2$$-$CoO$_2$ spacing  $c$  (Å)")
ax.set_ylabel(r"Landau stiffness  $\alpha$  (eV/Å$^2$)")
ax.legend(loc="upper right", handletextpad=0.4, labelspacing=0.3,
          borderaxespad=0.7, fontsize=7.6)
S.thin_spines(ax)

fig.subplots_adjust(left=0.155, right=0.965, top=0.965, bottom=0.145)
S.save(fig, "figS_thresholds")
print("figS_thresholds written")

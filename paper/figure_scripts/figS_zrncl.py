#!/usr/bin/env python
"""
Figure S (null-hypothesis contrast) - figS_zrncl.
Two side-by-side panels, both Tc (K) vs an interlayer/gallery spacing (Å):

 (a) THE NITRIDE (beta-M NCl, M = Zr, Hf): literature Tc vs basal spacing d,
     from reanalysis/zrncl_data.csv.  A layered host + rigidly-ionized donor
     with NO soft donor ion: Tc rises modestly then SATURATES and never falls,
     right out to the 19-22 A propylene-carbonate / THF-bilayer swellings.
     1/d guide curves (ZRNCL.md fits) trace the monotonic-saturating trend.

 (b) THE HYDRATE (our result): the computed Allen-Dynes Tc(c) dome (reused
     verbatim from theory/results/coupling_tc_vs_c.csv, exactly as fig2.py
     panel (c) builds it) with its gamma-sensitivity band and the 4.5 K
     experimental star at 9.9 A.  A narrow, NON-monotonic pairing window.

Takeaway (caption): a gallery WITHOUT a soft donor ion gives a monotonic,
saturating Tc(spacing); the hydrate's nonmonotonic window is the anomaly that
demands the soft ion.  The two panels share the same y-concept and styling so
the shape contrast (broad blue "SC everywhere" vs narrow amber window) lands.

House style: pastel region fills name a physical domain, physics labelled INSIDE
in matching ink; flat fills only, no gradients.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

import _style as S

S.use_house_style()

REPO = Path(__file__).resolve().parents[2]

# sober family colours (not tied to the main Na/Li/K series meaning)
C_ZR = "#5b6b7a"   # slate  : ZrNCl family
C_HF = "#b07a3c"   # ochre  : HfNCl family

# ----------------------------------------------------------------------
fig, (axA, axB) = plt.subplots(1, 2, figsize=(6.7, 3.15))

# ======================================================================
# (a) beta-MNCl : Tc vs basal spacing d
# ======================================================================
z = pd.read_csv(REPO / "reanalysis" / "zrncl_data.csv").dropna(
    subset=["d_A", "Tc_K"])

DLO, DHI = 8.6, 23.2
# broad "carriers present / SC at every spacing" domain (blue = the style's
# gas-on paramagnetic-metal / carriers-present colour)
axA.axvspan(DLO, DHI, color=S.FILL_BLUE, lw=0, zorder=0)
axA.text(0.5, 0.085, "superconducting at every spacing (no window)",
         transform=axA.transAxes, ha="center", va="center",
         color=S.INK_BLUE, fontsize=7.6, zorder=1.5)

# 1/d guide curves (fits quoted in ZRNCL.md): monotonic, saturating
dd = np.linspace(9.2, DHI, 200)
axA.plot(dd, 15.9 - 24.4 / dd, color=C_ZR, lw=1.1, ls=(0, (5, 3)), zorder=2)
dd2 = np.linspace(9.3, 19.2, 200)
axA.plot(dd2, 31.1 - 90.2 / dd2, color=C_HF, lw=1.1, ls=(0, (5, 3)), zorder=2)
# direct curve labels (no legend)
axA.text(20.4, 27.0, r"$\beta$-HfNCl", color=C_HF, fontsize=8, ha="left",
         va="center", zorder=5)
axA.text(20.4, 16.3, r"$\beta$-ZrNCl", color=C_ZR, fontsize=8, ha="left",
         va="center", zorder=5)

for fam, col in [("ZrNCl", C_ZR), ("HfNCl", C_HF)]:
    s = z[z.family == fam]
    bare = s[s.cointercalant == "none"]
    solv = s[s.cointercalant != "none"]
    # bare alkali (narrow gallery) = filled ; solvent-cointercalated
    # (swollen gallery) = open -> the spacing knob, turned hard, does nothing
    axA.plot(bare.d_A, bare.Tc_K, "o", ms=6, color=col, mfc=col, mew=0,
             zorder=4)
    axA.plot(solv.d_A, solv.Tc_K, "o", ms=6, color=col, mfc="white",
             mew=1.4, zorder=4)

axA.set_xlim(DLO, DHI)
axA.set_ylim(0, 29)
axA.set_xlabel(r"basal / gallery spacing  $d$  (Å)")
axA.set_ylabel(r"$T_c$  (K)")
# minimal fill-key: filled = bare ion, open = solvent-swollen gallery
axA.plot(10.05, 8.8, "o", ms=6, color=S.C_SEC, mfc=S.C_SEC, mew=0, zorder=4)
axA.text(10.5, 8.8, "bare ion", va="center", fontsize=6.9, color=S.C_SEC)
axA.plot(10.05, 6.8, "o", ms=6, color=S.C_SEC, mfc="white", mew=1.4, zorder=4)
axA.text(10.5, 6.8, "solvated (swollen)", va="center", fontsize=6.9,
         color=S.C_SEC)
S.thin_spines(axA)
S.panel_label(axA, "a")

# ======================================================================
# (b) the hydrate : our Allen-Dynes Tc(c) dome  (verbatim data path)
# ======================================================================
cp = pd.read_csv(REPO / "theory" / "results" / "coupling_tc_vs_c.csv"
                 ).sort_values("c_A")
c = cp["c_A"].values
tc0 = cp["Tc_K"].values
g1 = cp["Tc_K_g1"].values
g2 = cp["Tc_K_g2"].values
lo = np.minimum.reduce([tc0, g1, g2])
hi = np.maximum.reduce([tc0, g1, g2])

SC_LO, LAM_C = 6.35, 6.44          # SC window (amber) / lambda=2 wall
CLO, CHI = 5.5, 10.1

axB.axvspan(CLO, SC_LO, color=S.FILL_NONE, lw=0, zorder=0)
axB.axvspan(SC_LO, LAM_C, color=S.FILL_AMBER, lw=0, zorder=0)
axB.axvspan(LAM_C, CHI, color=S.FILL_RED, lw=0, zorder=0)

axB.fill_between(c, lo, hi, color=S.C_BLUE, alpha=0.22, lw=0, zorder=2)
axB.plot(c, tc0, color=S.C_BLUE, lw=1.8, zorder=4,
         label=r"$T_c$ (Allen$-$Dynes)")
# experimental SC anchor
axB.plot(9.9, 4.5, marker="*", ms=15, color=S.C_INK, mec="white", mew=0.6,
         zorder=6)
axB.text(9.9, 5.9, "4.5 K\n(9.9 Å)", ha="center", va="bottom", fontsize=7.0,
         color=S.C_INK, zorder=6)

# region labels INSIDE fills
axB.text(0.5 * (CLO + SC_LO), 0.5, "no 2DEG", transform=axB.get_xaxis_transform(),
         ha="center", va="center", color=S.INK_NONE, fontsize=7.8, zorder=3)
axB.text(6.395, 0.52, "SC window", rotation=90, transform=axB.get_xaxis_transform(),
         ha="center", va="center", color=S.INK_AMBER, fontsize=7.2, zorder=5)
axB.text(0.5 * (LAM_C + CHI), 0.30, r"polaronic ($\lambda\!>\!2$)",
         transform=axB.get_xaxis_transform(), ha="center", va="center",
         color=S.INK_RED, fontsize=7.8, zorder=3)

axB.set_xlim(CLO, CHI)
axB.set_ylim(0, 16)
axB.set_xlabel(r"CoO$_2$$-$CoO$_2$ spacing  $c$  (Å)")
axB.set_ylabel(r"$T_c$  (K)")
axB.legend(loc="upper right", handletextpad=0.5, borderaxespad=0.6,
           fontsize=7.4)
S.thin_spines(axB)
S.panel_label(axB, "b")

fig.subplots_adjust(left=0.085, right=0.985, top=0.95, bottom=0.145,
                    wspace=0.28)
S.save(fig, "figS_zrncl")
print("figS_zrncl written")

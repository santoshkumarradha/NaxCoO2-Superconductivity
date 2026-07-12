#!/usr/bin/env python
"""
Figure 2 - SINGLE-COLUMN variant (for judging).

Same three stacked broken-axis panels as fig2.py [alpha(c), N(0), Tc dome],
but re-laid for one journal column (~3.37 in wide) so the figure occupies far
less page area than the full-width figure*.  Fonts are kept at their absolute
point sizes, so against the smaller canvas the labels read proportionally
larger.  Output: figures/fig2_col.{pdf,png}.  Nothing in main.tex is touched.
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy.interpolate import PchipInterpolator

import _style as S

S.use_house_style()
# nudge a few sizes for the compact canvas
import matplotlib as mpl
mpl.rcParams.update({
    "font.size": 8.5,
    "axes.labelsize": 8.5,
    "xtick.labelsize": 7.5,
    "ytick.labelsize": 7.5,
})

res = S.load_results()
fits = S.load_well_fits()
cp = S.load_coupling().sort_values("c_A")

SC_LO, LAM_C = 6.35, 6.44
XLO, XHI = 5.4, 10.1
L0, L1 = 5.4, 7.1
R0, R1 = 8.2, 10.1
CG = np.linspace(5.5, 9.9, 400)
C_ANN = S.C_INK

# ----- single-column canvas --------------------------------------------
fig = plt.figure(figsize=(3.37, 5.05))
gs = GridSpec(3, 2, width_ratios=[2.0, 1.0], height_ratios=[1, 1, 1],
              hspace=0.14, wspace=0.045,
              left=0.17, right=0.985, top=0.99, bottom=0.075)

axAL = fig.add_subplot(gs[0, 0]); axAR = fig.add_subplot(gs[0, 1], sharey=axAL)
axBL = fig.add_subplot(gs[1, 0], sharex=axAL); axBR = fig.add_subplot(gs[1, 1], sharex=axAR, sharey=axBL)
axCL = fig.add_subplot(gs[2, 0], sharex=axAL); axCR = fig.add_subplot(gs[2, 1], sharex=axAR, sharey=axCL)
LEFT = [axAL, axBL, axCL]
RIGHT = [axAR, axBR, axCR]
ROWS = [(axAL, axAR), (axBL, axBR), (axCL, axCR)]

axAL.set_xlim(L0, L1); axAR.set_xlim(R0, R1)
for axl, axr in ROWS:
    for ax in (axl, axr):
        ax.axvspan(XLO, SC_LO, color=S.FILL_NONE, lw=0, zorder=0)
        ax.axvspan(SC_LO, LAM_C, color=S.FILL_AMBER, lw=0, zorder=0)
        ax.axvspan(LAM_C, XHI, color=S.FILL_RED, lw=0, zorder=0)
        for c, _lab, _col, _sc in S.ANCHORS:
            ax.axvline(c, color=S.C_MUT, ls=(0, (5, 3)), lw=0.8, alpha=0.85, zorder=1)


def break_marks(axl, axr):
    d = 0.7
    kw = dict(marker=[(-1, -d), (1, d)], markersize=6, linestyle="none",
              color=S.C_SEC, mec=S.C_SEC, mew=0.9, clip_on=False, zorder=6)
    axl.plot([1, 1], [0, 1], transform=axl.transAxes, **kw)
    axr.plot([0, 0], [0, 1], transform=axr.transAxes, **kw)


# ---------------------------------------------------------------- (a) alpha
for (el, cell), st in S.SERIES.items():
    w = fits[(fits.element == el) & (fits.cell == cell)].sort_values("c")
    p = PchipInterpolator(w["c"].values, w["alpha"].values)
    for ax in (axAL, axAR):
        ax.plot(CG, p(CG), color=st["color"], lw=1.3, zorder=3)
    axAL.plot(w["c"], w["alpha"], st["marker"], ms=4.5, color=st["color"],
              mfc="white", mew=1.3, label=st["label"], zorder=4)
    axAR.plot(w["c"], w["alpha"], st["marker"], ms=4.5, color=st["color"],
              mfc="white", mew=1.3, zorder=4)
for ax in (axAL, axAR):
    ax.axhline(0, color=S.C_INK, lw=0.8, zorder=2)
axAL.set_ylim(-1.2, 3.4)
axAL.set_ylabel(r"$\alpha$  (eV/Å$^2$)")
axAL.legend(loc="upper right", ncol=1, handletextpad=0.3, labelspacing=0.25,
            borderaxespad=0.4, handlelength=1.3, fontsize=6.4)
S.panel_label(axAL, "a")

# ------------------------------------------------- (b) 2DEG turn-on: N(0)
lab = {("Na", "1x1"): (9.85, 1.78, "right", "bottom"),
       ("Li", "1x1"): (9.85, 0.92, "right", "top"),
       ("Na", "s3"):  (8.30, 4.05, "left", "center")}
for (el, cell), st in S.SERIES.items():
    r = res[(res.element == el) & (res.cell == cell)].sort_values("c")
    for ax in (axBL, axBR):
        ax.plot(r["c"], r["N0"], st["marker"] + "-", ms=4.5, lw=1.4,
                color=st["color"], mfc=st["color"], mew=0, zorder=4)
    lx, ly, lha, lva = lab[(el, cell)]
    axBR.text(lx, ly, st["label"], color=C_ANN, fontsize=6.6, ha=lha, va=lva)
for ax in (axBL, axBR):
    ax.axhline(0.10, color=S.C_MUT, ls=":", lw=1.0, zorder=2)
axBR.text(8.30, 0.36, r"$N(0){=}0.1$ gate", ha="left", va="bottom",
          fontsize=6.6, color=S.C_SEC)
axBL.set_ylim(-0.12, 4.3)
axBL.set_ylabel(r"DFT $N(0)$" "\n" r"(eV$^{-1}$ cell$^{-1}$)", fontsize=7.8)
S.panel_label(axBL, "b")

# ---------------------------------------------------------------- (c) Tc dome
c = cp["c_A"].values
tc0 = cp["Tc_K"].values
g1 = cp["Tc_K_g1"].values
g2 = cp["Tc_K_g2"].values
lo = np.minimum.reduce([tc0, g1, g2])
hi = np.maximum.reduce([tc0, g1, g2])
for k, ax in enumerate((axCL, axCR)):
    ax.fill_between(c, lo, hi, color=S.C_BLUE, alpha=0.22, lw=0, zorder=2,
                    label=("ansatz range" if k == 0 else None))
    ax.plot(c, tc0, color=S.C_BLUE, lw=1.6, zorder=4,
            label=(r"$T_c$ (A$-$D)" if k == 0 else None))
axCL.axvline(LAM_C, color=S.EDGE_STONER, lw=1.0, zorder=1.5)
axCR.plot(9.9, 4.5, marker="*", ms=13, color=C_ANN, mec="white", mew=0.6, zorder=6)
axCL.text(0.5 * (XLO + SC_LO), 0.5, "no 2DEG",
          transform=axCL.get_xaxis_transform(), ha="center", va="center",
          color=S.INK_NONE, fontsize=7.2, zorder=3)
axCL.text(6.395, 19.5, "SC window", rotation=90, ha="center", va="center",
          color=S.INK_AMBER, fontsize=6.8, zorder=5)
axCL.text(6.80, 0.90, r"polaronic", transform=axCL.get_xaxis_transform(),
          ha="left", va="top", color=S.INK_RED, fontsize=7.2, zorder=3)
axCL.set_ylim(0, 26)
axCL.set_ylabel(r"$T_c$  (K)")
axCL.legend(loc="upper left", handletextpad=0.4, borderaxespad=0.4,
            labelspacing=0.25, handlelength=1.3, fontsize=6.4)
S.panel_label(axCL, "c")

# ----- spines / ticks / breaks -----------------------------------------
for ax in LEFT:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(top=False, right=False)
for ax in RIGHT:
    ax.spines["top"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.tick_params(top=False, left=False, labelleft=False)

axAL.set_xticks([5.5, 6.0, 6.5, 7.0])
axAR.set_xticks([8.5, 9.0, 9.5, 10.0])
for ax in (axAL, axBL, axAR, axBR):
    ax.tick_params(labelbottom=False)

for axl, axr in ROWS:
    break_marks(axl, axr)

fig.supxlabel(r"CoO$_2$$-$CoO$_2$ spacing  $c$  (Å)", fontsize=8, y=0.008)

S.save(fig, "fig2_col")
print("fig2_col written")

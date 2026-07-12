#!/usr/bin/env python
"""
Figure 2 - The money plot (paper centrepiece).
Three stacked panels sharing the CoO2-CoO2 spacing axis c:
 (a) Landau alpha(c) for the three series; the alpha=0 crossing / SC window shaded.
 (b) Carrier turn-on: DFT N(0) (left) and Bader areal density n_2D (right).
 (c) Tc(c) dome (Allen-Dynes) with the electron-phonon gamma-sensitivity band
     (11-23 K), and the lambda>2 "polaronic - SC killed" region hatched.
Experimental anchors: 5.5 & 6.9 A (grey, not SC), 9.9 A (red star, SC 4.5 K).
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import PchipInterpolator

import _style as S

S.use_house_style()

res = S.load_results()
fits = S.load_well_fits()
cp = S.load_coupling().sort_values("c_A")

SC_LO, SC_HI = 6.35, 6.45        # SC window (task)
LAM_C = 6.44                     # first c with lambda > 2 (polaronic)
XLO, XHI = 5.4, 10.0             # panel x-range
CG = np.linspace(5.5, 9.9, 400)

# House style (spin_phase_diagram.pdf): pastel region fills name domains and
# the physics is labelled INSIDE each fill in matching darker ink; data curves
# and annotation text ride on top.  amber = pairing window, red = polaronic
# (SC killed), neutral = mechanism off.
C_ANN = S.C_INK                  # annotation text: black
C_ANN2 = S.C_SEC                 # secondary annotation text: ~40% grey

fig, axes = plt.subplots(3, 1, figsize=(6.9, 7.4), sharex=True)
axA, axB, axC = axes


def shade_window(ax):
    # the pairing window as an amber strip, echoing panel (c)
    ax.axvspan(SC_LO, SC_HI, color=S.FILL_AMBER, lw=0, zorder=0)


def anchors(ax, star=False):
    for c, _lab, col, is_sc in S.ANCHORS:
        ax.axvline(c, color=S.C_MUT, ls=(0, (5, 3)), lw=0.8,
                   alpha=0.8, zorder=1)


# ---------------------------------------------------------------- (a) alpha
for (el, cell), st in S.SERIES.items():
    w = fits[(fits.element == el) & (fits.cell == cell)].sort_values("c")
    p = PchipInterpolator(w["c"].values, w["alpha"].values)
    axA.plot(CG, p(CG), color=st["color"], lw=1.4, zorder=3)
    axA.plot(w["c"], w["alpha"], st["marker"], ms=6, color=st["color"],
             mfc="white", mew=1.5, label=st["label"], zorder=4)
axA.axhline(0, color=S.C_INK, lw=0.8, zorder=2)
shade_window(axA)
anchors(axA)
axA.annotate(r"$c^\ast_{\mathrm{Na}}\!\approx\!6.40$ Å", xy=(6.40, 0),
             xytext=(7.15, 1.35), fontsize=8, color=C_ANN,
             arrowprops=dict(arrowstyle="->", color=C_ANN2, lw=0.7))
axA.set_ylim(-1.2, 3.4)
axA.set_ylabel(r"$\alpha$  (eV/Å$^2$)")
axA.legend(loc="upper right", ncol=1, handletextpad=0.4, fontsize=8)
S.thin_spines(axA)
S.panel_label(axA, "a")
axA.text(SC_LO - 0.02, 3.05, "double well\n" r"($\alpha<0$)", ha="right",
         va="top", fontsize=8, color=S.C_SEC)

# ------------------------------------------------- (b) 2DEG turn-on: N(0)
# Single axis (no dual scale): DFT N(0) is the direct carrier-turn-on measure.
# Bader n_2D tracks the same transition and is tabulated in the README; it is
# omitted here (it falls with dilution and would blur the turn-on message).
lab = {("Na", "1x1"): (9.82, 1.78, "right"), ("Li", "1x1"): (9.82, 0.92, "right"),
       ("Na", "s3"): (8.55, 3.66, "right")}
for (el, cell), st in S.SERIES.items():
    r = res[(res.element == el) & (res.cell == cell)].sort_values("c")
    axB.plot(r["c"], r["N0"], st["marker"] + "-", ms=6, lw=1.6,
             color=st["color"], mfc=st["color"], mew=0, zorder=4)
    lx, ly, lha = lab[(el, cell)]
    axB.text(lx, ly, st["label"], color=C_ANN, fontsize=8,
             ha=lha, va="center")
axB.axhline(0.10, color=S.C_MUT, ls=":", lw=1.0, zorder=2)
axB.text(7.35, 0.28, r"$N(0){=}0.1$ gate", ha="left", va="bottom",
         fontsize=7.6, color=C_ANN2)
shade_window(axB)
anchors(axB)
axB.set_ylim(-0.12, 4.3)
axB.set_ylabel(r"DFT $N(0)$" "\n" r"(states eV$^{-1}$ cell$^{-1}$ spin$^{-1}$)")
S.thin_spines(axB)
S.panel_label(axB, "b")

# ---------------------------------------------------------------- (c) Tc dome
c = cp["c_A"].values
tc0 = cp["Tc_K"].values
g1 = cp["Tc_K_g1"].values
g2 = cp["Tc_K_g2"].values
lo = np.minimum.reduce([tc0, g1, g2])
hi = np.maximum.reduce([tc0, g1, g2])
def draw_dome(ax, lw_line=1.8, band_lab=None, line_lab=None):
    ax.fill_between(c, lo, hi, color=S.C_BLUE, alpha=0.22, lw=0, zorder=2,
                    label=band_lab)
    ax.plot(c, tc0, color=S.C_BLUE, lw=lw_line, zorder=4, label=line_lab)


# three labelled domains (the reference-figure region-fill style):
#   neutral  no 2DEG (carriers off, no pairing)  [XLO, SC_LO]
#   amber    the SC window (the narrow dome)      [SC_LO, LAM_C]
#   red      polaronic, lambda>2, SC killed       [LAM_C, XHI]
S.vspan_region(axC, XLO, SC_LO, S.FILL_NONE, "no 2DEG", ytext=0.5,
               ha="center", fontsize=8.5)
axC.axvspan(SC_LO, LAM_C, color=S.FILL_AMBER, lw=0, zorder=0)
S.vspan_region(axC, LAM_C, XHI, S.FILL_RED, r"polaronic ($\lambda>2$)",
               ytext=0.93, ha="center", fontsize=8.5)
axC.axvline(LAM_C, color=S.EDGE_STONER, lw=1.0, zorder=1)

draw_dome(axC, band_lab=r"$\gamma$ band (11$-$23 K)",
          line_lab=r"$T_c$ (Allen$-$Dynes)")
# point to the hairline amber window (too narrow to label in place)
axC.annotate("SC window", xy=(SC_LO, 6.0), xytext=(6.62, 19.0),
             fontsize=8.5, color=S.INK_AMBER, ha="left", va="center",
             arrowprops=dict(arrowstyle="->", color=S.INK_AMBER, lw=0.9))
anchors(axC, star=True)
# experimental SC anchor: distinct star marker at 9.9 A, Tc = 4.5 K (black)
axC.plot(9.9, 4.5, marker="*", ms=16, color=C_ANN, mec="white", mew=0.6,
         zorder=6)
axC.annotate("expt. " r"$T_c=4.5$ K", xy=(9.82, 4.7),
             xytext=(9.15, 10.5), fontsize=8, color=C_ANN, ha="center",
             arrowprops=dict(arrowstyle="->", color=C_ANN2, lw=0.9))
axC.set_ylim(0, 26)
axC.set_ylabel(r"$T_c$  (K)")
axC.set_xlabel(r"CoO$_2$$-$CoO$_2$ spacing  $c$  (Å)")
axC.legend(loc="upper left", handletextpad=0.5, borderaxespad=0.4, fontsize=8)
S.thin_spines(axC)
S.panel_label(axC, "c")

# --- zoom inset: the narrow SC dome is otherwise a hairline ---
# placed over the (data-free) polaronic dead zone, right of the real dome
# spike. Clean-inset rules: fully opaque white background, a complete square
# frame on all four sides, drawn ABOVE the panel shading, and NO internal
# shading (the parent's grey bands would only clutter the zoom); the lambda=2
# wall is a single thin grey line, as in the parent.
axz = axC.inset_axes([0.40, 0.30, 0.33, 0.54])
axz.set_facecolor("white")
axz.set_zorder(6)
axz.patch.set_alpha(1.0)
draw_dome(axz, lw_line=1.6)
axz.axvline(LAM_C, color=S.C_MUT, lw=0.7)
axz.set_xlim(6.28, 6.52)
axz.set_ylim(0, 25)
axz.text(0.5, 0.965, r"dome (zoom), peak $T_c\!\approx\!14$ K",
         transform=axz.transAxes, ha="center", va="top", fontsize=7.5,
         color=C_ANN)
axz.tick_params(labelsize=7, length=2)
axz.set_xticks([6.3, 6.4, 6.5])
axz.set_yticks([0, 10, 20])
for s in axz.spines.values():           # full square border
    s.set_visible(True)
    s.set_color(S.C_SEC)
    s.set_linewidth(0.8)
axC.indicate_inset_zoom(axz, edgecolor=S.C_MUT, lw=0.7, alpha=0.8)

for ax in axes:
    ax.set_xlim(5.4, 10.0)

fig.subplots_adjust(hspace=0.10, left=0.095, right=0.905, top=0.985,
                    bottom=0.115)
S.save(fig, "fig2")
print("fig2 written")

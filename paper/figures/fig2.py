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
CG = np.linspace(5.5, 9.9, 400)

fig, axes = plt.subplots(3, 1, figsize=(6.9, 7.4), sharex=True)
axA, axB, axC = axes


def shade_window(ax):
    ax.axvspan(SC_LO, SC_HI, color=S.C_YELL, alpha=0.14, lw=0, zorder=0)


def anchors(ax, star=False):
    for c, _lab, col, is_sc in S.ANCHORS:
        ax.axvline(c, color=col, ls=(0, (5, 3)), lw=1.0,
                   alpha=0.85 if is_sc else 0.55, zorder=1)


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
             xytext=(7.15, 1.35), fontsize=8, color=S.C_BLUE,
             arrowprops=dict(arrowstyle="->", color=S.C_BLUE, lw=0.7))
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
    axB.text(lx, ly, st["label"], color=st["color"], fontsize=8,
             ha=lha, va="center", fontweight="bold")
axB.axhline(0.10, color=S.C_RED, ls=":", lw=1.0, zorder=2)
axB.text(7.35, 0.30, "2DEG gate  " r"$N(0){=}0.1$", ha="left", va="bottom",
         fontsize=7.6, color=S.C_RED)
shade_window(axB)
anchors(axB)
axB.set_ylim(-0.12, 4.3)
axB.set_ylabel(r"DFT $N(0)$" "\n" r"(states eV$^{-1}$ cell$^{-1}$ spin$^{-1}$)")
S.thin_spines(axB)
S.panel_label(axB, "b")
axB.annotate("2DEG switches on", xy=(6.55, 0.9), xytext=(6.75, 2.3),
             fontsize=8, color=S.C_INK, ha="left",
             arrowprops=dict(arrowstyle="->", color=S.C_INK, lw=0.8))

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


draw_dome(axC, band_lab=r"$\gamma$ band (11$-$23 K)",
          line_lab=r"$T_c$ (Allen$-$Dynes)")
# polaronic region: lambda > 2 -> SC killed
axC.axvspan(LAM_C, 9.9, facecolor="none", edgecolor=S.C_CRIT, lw=0.0,
            hatch="////", alpha=0.5, zorder=1)
axC.text(9.6, 21.5, "polaronic self-trapping\n" r"($\lambda>2$) $-$ SC killed",
         ha="right", va="center", fontsize=8, color=S.C_CRIT)
shade_window(axC)
anchors(axC, star=True)
# experimental SC anchor: star at 9.9 A, Tc = 4.5 K
axC.plot(9.9, 4.5, marker="*", ms=16, color=S.C_RED, mec="white", mew=0.6,
         zorder=6)
axC.annotate("Takada 2003\n" r"expt. $T_c=4.5$ K", xy=(9.82, 4.7),
             xytext=(9.15, 10.5), fontsize=8, color=S.C_RED, ha="center",
             arrowprops=dict(arrowstyle="->", color=S.C_RED, lw=0.9))
axC.set_ylim(0, 26)
axC.set_ylabel(r"$T_c$  (K)")
axC.set_xlabel(r"CoO$_2$$-$CoO$_2$ spacing  $c$  (Å)")
axC.legend(loc="upper left", handletextpad=0.5, borderaxespad=0.4, fontsize=8)
S.thin_spines(axC)
S.panel_label(axC, "c")

# --- zoom inset: the narrow SC dome is otherwise a hairline ---
# placed over the (data-free) polaronic dead zone, right of the real dome
# spike, with an opaque white background so the hatch never shows through.
axz = axC.inset_axes([0.40, 0.32, 0.33, 0.52])
axz.set_facecolor("white")
axz.set_zorder(5)
axz.patch.set_alpha(1.0)
draw_dome(axz, lw_line=1.6)
axz.axvspan(LAM_C, 6.6, facecolor="none", edgecolor=S.C_CRIT, lw=0.0,
            hatch="////", alpha=0.5)
axz.axvspan(SC_LO, SC_HI, color=S.C_YELL, alpha=0.14, lw=0)
axz.set_xlim(6.28, 6.52)
axz.set_ylim(0, 25)
axz.set_title(r"dome (zoom), peak $T_c\!\approx\!14$ K", fontsize=7.5, pad=3)
axz.tick_params(labelsize=7, length=2)
axz.set_xticks([6.3, 6.4, 6.5])
axz.set_yticks([0, 10, 20])
for s in ("top", "right"):
    axz.spines[s].set_visible(False)
axC.indicate_inset_zoom(axz, edgecolor=S.C_MUT, lw=0.7, alpha=0.8)

axC.text(0.5, -0.30,
         r"$9.9$ Å reconciliation requires water softening of the Na well "
         r"(prediction, untested)",
         transform=axC.transAxes, ha="center", va="top", fontsize=8,
         color=S.C_SEC, style="italic")

for ax in axes:
    ax.set_xlim(5.4, 10.0)

fig.subplots_adjust(hspace=0.10, left=0.095, right=0.905, top=0.985,
                    bottom=0.115)
S.save(fig, "fig2")
print("fig2 written")

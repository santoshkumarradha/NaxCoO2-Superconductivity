#!/usr/bin/env python
"""
Figure 5 - Mode softening.
Effective zero-point frequency omega_eff(c) (log scale) for the three series
(theory/results/quantum_wells_v4.csv, physical alkali masses): the alkali mode
collapses from ~30 meV (stiff, single well) to ~10 mueV (soft double well) as
the gallery opens, entering a quantum-paraelectric window. Inset: the quantized
levels E0/E1/E2 in the c = 6.9 A Na double well (E0,E1 a near-degenerate
tunnelling doublet).
"""
import numpy as np
import matplotlib.pyplot as plt

import _style as S

S.use_house_style()

qw = S.load_quantum()
fits = S.load_well_fits()

fig, ax = plt.subplots(figsize=(3.4, 3.5))

# quantum-paraelectric (soft-mode) window: omega < 1 meV.  Amber pastel (the
# reference-figure "favourable window" colour), labelled inside in matching ink.
ax.axhspan(5e-3, 1.0, color=S.FILL_AMBER, lw=0, zorder=0)
ax.axhline(1.0, color=S.EDGE_SSTAR, lw=0.9, ls=(0, (4, 3)), zorder=1)
ax.text(5.42, 0.0075, "soft-mode window", ha="left", va="bottom",
        fontsize=8, color=S.INK_AMBER)

for (el, cell), st in S.SERIES.items():
    q = qw[(qw.element == el) & (qw.cell == cell) & (qw.mass == el)]
    q = q.sort_values("c")
    ax.semilogy(q["c"], q["hw_eff_meV"], st["marker"] + "-", ms=5.5, lw=1.4,
                color=st["color"], mfc=st["color"], mew=0, zorder=4,
                label=st["label"])

ax.set_xlim(5.3, 10.0)
ax.set_ylim(5e-3, 60)
ax.set_xlabel(r"CoO$_2$$-$CoO$_2$ spacing  $c$  (Å)")
ax.set_ylabel(r"$\hbar\omega_{\mathrm{eff}}=\hbar^2/2M\langle\delta^2"
              r"\rangle_0$  (meV)")
ax.legend(loc="upper right", handletextpad=0.4, fontsize=7.6)
S.thin_spines(ax)

# ---------------- inset: levels in the 6.9 A Na double well -------------
# centre-right, dropped below the upper-right legend. Clean-inset rules: a
# padded opaque white card behind the WHOLE inset (axes + tick labels + axis
# labels + title, so the window shading never shows through anywhere), plus a
# complete square frame on the inset itself.
r = fits[(fits.element == "Na") & (fits.cell == "1x1") & (fits.c == 6.9)].iloc[0]
lv = qw[(qw.element == "Na") & (qw.cell == "1x1") & (qw.c == 6.9)
        & (qw.mass == "Na")].iloc[0]
E0, E1, E2 = lv["E0_meV"], lv["E1_meV"], lv["E2_meV"]

import matplotlib.patches as mpatches
# The card sits in the data-free pocket below the legend and above the
# low-lying (c >= 7.5 A) curve segments; verified against the log-scale
# curve positions so NO data point or curve segment is occluded.  Enlarged
# so the inset tick/axis labels read at >= 7 pt.
card = mpatches.FancyBboxPatch(
    (0.455, 0.285), 0.535, 0.475, transform=ax.transAxes,
    boxstyle="square,pad=0", facecolor="white", edgecolor="none",
    zorder=5)
ax.add_patch(card)

axi = ax.inset_axes([0.555, 0.37, 0.42, 0.35])
axi.set_facecolor("white")
axi.patch.set_alpha(1.0)
axi.set_zorder(6)
dd = np.linspace(-0.95, 0.95, 300)
V = S.vpoly(dd, r["a_q"], r["b_q"], r["g_q"]) * 1e3
axi.plot(dd, V, color=S.C_INK, lw=1.2, zorder=3)


def level(E, col, lab, lab_dy=3.0):
    # draw the level across its classically allowed width; label at the LEFT
    # end just above the line, inside the well mouth where the curve is clear.
    within = np.where(V <= E)[0]
    if len(within):
        x0, x1 = dd[within[0]], dd[within[-1]]
    else:
        x0, x1 = -0.2, 0.2
    axi.plot([x0, x1], [E, E], color=col, lw=1.6, zorder=4)
    axi.text(0.0, E + lab_dy, lab, ha="center", va="bottom",
             fontsize=7.5, color=col, zorder=5)


level(E2, S.C_YELL, r"$E_2$", lab_dy=2.5)
level(E0, S.C_BLUE, r"$E_0,E_1$", lab_dy=2.5)   # tunnelling doublet
axi.set_ylim(V.min() - 8, 42)
axi.set_xlim(-1.0, 1.0)
axi.set_xticks([-0.5, 0.0, 0.5])
axi.tick_params(labelsize=7, length=2.4)
axi.set_xlabel(r"$\delta$ (Å)", fontsize=7.5, labelpad=1.5)
axi.set_ylabel("meV", fontsize=7.5, labelpad=1.5)
axi.set_title(r"$c=6.9$ Å Na double well", fontsize=7.5, pad=3)
for s in axi.spines.values():           # full square border
    s.set_visible(True)
    s.set_color(S.C_SEC)
    s.set_linewidth(0.8)

fig.subplots_adjust(left=0.155, right=0.97, top=0.97, bottom=0.135)
S.save(fig, "fig5")
print("fig5 written")

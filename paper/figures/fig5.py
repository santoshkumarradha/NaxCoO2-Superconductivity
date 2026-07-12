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

# quantum-paraelectric (soft-mode) window: omega < 1 meV  (neutral shading)
ax.axhspan(1e-3, 1.0, color="#e2e1db", alpha=0.9, lw=0, zorder=0)
ax.axhline(1.0, color=S.C_MUT, lw=0.7, ls=(0, (4, 3)), alpha=0.8, zorder=1)
ax.text(5.42, 0.0068, "quantum-paraelectric\n(soft-mode) window", ha="left",
        va="bottom", fontsize=7.6, color=S.C_SEC)

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
# centre-right, dropped below the upper-right legend; opaque white background
# masks the window shading beneath it.
r = fits[(fits.element == "Na") & (fits.cell == "1x1") & (fits.c == 6.9)].iloc[0]
lv = qw[(qw.element == "Na") & (qw.cell == "1x1") & (qw.c == 6.9)
        & (qw.mass == "Na")].iloc[0]
E0, E1, E2 = lv["E0_meV"], lv["E1_meV"], lv["E2_meV"]

axi = ax.inset_axes([0.45, 0.37, 0.50, 0.37])
axi.set_facecolor("white")
axi.patch.set_alpha(1.0)
dd = np.linspace(-0.95, 0.95, 300)
V = S.vpoly(dd, r["a_q"], r["b_q"], r["g_q"]) * 1e3
axi.plot(dd, V, color=S.C_INK, lw=1.1, zorder=3)


def level(E, col, lab, lw=1.3):
    # draw the level across its classically allowed width
    within = np.where(V <= E)[0]
    if len(within):
        x0, x1 = dd[within[0]], dd[within[-1]]
    else:
        x0, x1 = -0.2, 0.2
    axi.plot([x0, x1], [E, E], color=col, lw=lw, zorder=4)
    axi.text(0.98, E, lab, transform=axi.get_yaxis_transform(), ha="right",
             va="bottom", fontsize=7, color=S.C_INK)


level(E2, S.C_YELL, r"$E_2$")
level(E0, S.C_BLUE, r"$E_0,E_1$")   # near-degenerate tunnelling doublet
axi.set_ylim(V.min() - 8, 40)
axi.set_xlim(-1.0, 1.0)
axi.tick_params(labelsize=6.5, length=2)
axi.set_xlabel(r"$\delta$ (Å)", fontsize=7, labelpad=1)
axi.set_ylabel("meV", fontsize=7, labelpad=1)
axi.set_title(r"$c=6.9$ Å Na double well", fontsize=7, pad=3)
for s in ("top", "right"):
    axi.spines[s].set_visible(False)

fig.subplots_adjust(left=0.155, right=0.97, top=0.97, bottom=0.135)
S.save(fig, "fig5")
print("fig5 written")

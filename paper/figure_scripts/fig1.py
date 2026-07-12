#!/usr/bin/env python
"""
Figure 1 - The transition.
(a) E(delta) for Na 1x1 at c = 5.5/6.9/8.4/9.9 A; (b) same for Na sqrt3.
Points = DFT totals (runpod jobs/*/pw.out); lines = quartic/sextic fits
(theory/results/well_fits_v4.csv). Well minima marked; the single->double-well
crossover annotated; open markers + down-arrows flag the unterminated (c>=8.4)
adsorption-regime scans whose depths are only lower bounds.
"""
import numpy as np
import matplotlib.pyplot as plt

import _style as S

S.use_house_style()

CS = [5.5, 6.9, 8.4, 9.9]
PANELS = [("Na", "1x1", S.RAMP_BLUE, (-620, 360), "a"),
          ("Na", "s3", S.RAMP_YELL, (-460, 320), "b")]

fits = S.load_well_fits()
scans = S.load_scans()

fig, axes = plt.subplots(2, 1, figsize=(3.4, 4.9), sharex=True)

for ax, (el, cell, ramp, ylim, letter) in zip(axes, PANELS):
    for c, col in zip(CS, ramp):
        key = (el, cell, c)
        if key not in scans:
            continue
        d, E = scans[key]                       # meV, referenced to delta=0
        r = fits[(fits.element == el) & (fits.cell == cell)
                 & (fits.c == c)].iloc[0]
        unterm = not r["min_bracketed"] and c >= 8.4   # scan still falling
        # --- fit curve over the sampled |delta| range (symmetrised) ---
        dmax = d.max()
        xx = np.linspace(-dmax, dmax, 400)
        yy = S.vpoly(xx, r["a_q"], r["b_q"], r["g_q"]) * 1e3
        ax.plot(xx, yy, color=col, lw=1.3,
                ls="--" if r["refit_sextic"] else "-", zorder=2)
        # --- DFT points, symmetrised about delta = 0 ---
        ds = np.concatenate([-d[::-1], d])
        Es = np.concatenate([E[::-1], E])
        mfc = "white" if unterm else col
        ax.plot(ds, Es, "o", ms=3.6, color=col, mfc=mfc, mew=0.9,
                zorder=4, label=f"{c:.1f}")
        # --- well minimum (only if bracketed by the scan) ---
        if r["min_bracketed"] and r["depth_eV"] > 1e-3:
            d0 = r["d0_A"]
            y0 = S.vpoly(d0, r["a_q"], r["b_q"], r["g_q"]) * 1e3
            for sgn in (-1, 1):
                ax.plot(sgn * d0, y0, marker="|", ms=7, color=col,
                        mew=1.3, zorder=5)
        # --- unterminated: down-arrow at the last sampled point ---
        if unterm:
            for sgn in (-1, 1):
                ax.annotate("", xy=(sgn * dmax, E[-1] - 0.16 * (ylim[1] - ylim[0])),
                            xytext=(sgn * dmax, E[-1]),
                            arrowprops=dict(arrowstyle="-|>", color=col, lw=1.0),
                            zorder=5)
    ax.axhline(0, color=S.C_MUT, lw=0.6, ls=(0, (4, 3)), zorder=1)
    ax.set_ylim(*ylim)
    ax.set_xlim(-1.15, 1.15)
    S.thin_spines(ax)
    ax.set_ylabel(r"$E(\delta)-E(0)$  (meV)")
    S.panel_label(ax, letter)
    lab = S.SERIES[(el, cell)]["label"]
    ax.text(0.5, 0.94, lab, transform=ax.transAxes, ha="center", va="top",
            fontsize=8, color=S.C_INK)
    leg = ax.legend(title=r"$c$ (Å)", loc="lower center", ncol=4,
                    handletextpad=0.3, columnspacing=0.9, borderpad=0.3,
                    handlelength=0.9)
    leg.get_title().set_fontsize(7)

axes[-1].set_xlabel(r"off-centre displacement  $\delta$  (Å)")

# single -> double-well crossover annotation (alpha changes sign 5.5 -> 6.9 A)
axes[0].annotate("single well\n" r"($\alpha>0$)", xy=(0.0, 250),
                 xytext=(0.62, 250), fontsize=6.6, color=S.C_SEC,
                 ha="center", va="center",
                 arrowprops=dict(arrowstyle="->", color=S.C_SEC, lw=0.7))
axes[0].annotate("double well\n" r"($\alpha<0$)", xy=(-0.69, -110),
                 xytext=(-0.62, -330), fontsize=6.6, color=S.C_SEC,
                 ha="center", va="center",
                 arrowprops=dict(arrowstyle="->", color=S.C_SEC, lw=0.7))

fig.subplots_adjust(hspace=0.08, left=0.17, right=0.97, top=0.97, bottom=0.09)
S.save(fig, "fig1")
print("fig1 written")

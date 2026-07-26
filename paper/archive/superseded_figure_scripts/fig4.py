#!/usr/bin/env python
"""
Figure 4 - Spin turn-on.
Cell-integrated |m|(c) for Li 1x1 and Na 1x1 (LSDA spin scans,
spin_analysis/magnetization_v2.csv): the local moment switches on exactly at
the alpha<0 well transition. Markers = mean over the delta scan, whiskers =
min/max over delta. The alkali gallery moment is antiparallel to Co (noted).
The Na sqrt3 cell is over-polarised in LSDA (m ~ const ~ 2.3 muB, unphysical)
and is shown only as a shaded caveat band.
"""
import numpy as np
import matplotlib.pyplot as plt

import _style as S

S.use_house_style()

mag = S.load_magnetization()


def series_stats(el, cell):
    d = mag[(mag.alkali == el) & (mag.cell == cell)]
    cs, mean, lo, hi = [], [], [], []
    for c, g in d.groupby("c_A"):
        v = g["abs_mag_muB"].values
        cs.append(c); mean.append(v.mean()); lo.append(v.min()); hi.append(v.max())
    o = np.argsort(cs)
    return (np.array(cs)[o], np.array(mean)[o], np.array(lo)[o], np.array(hi)[o])


fig, ax = plt.subplots(figsize=(3.4, 3.05))

C_TR = 6.40             # well transition / moment turn-on (c*_Na)

# --- two pastel domains along c (reference-figure region-fill style):
#     neutral  no local moment (paramagnetic)   [5.3, C_TR]
#     red      magnetic: gallery moment on       [C_TR, 10.0]
S.vspan_region(ax, 5.3, C_TR, S.FILL_NONE)
S.vspan_region(ax, C_TR, 10.0, S.FILL_RED, "magnetic", ytext=0.60,
               ha="center", fontsize=8.5)
ax.axvline(C_TR, color=S.EDGE_STONER, lw=1.0, zorder=1)
ax.text(5.95, 1.4, "no moment", ha="center", va="center", rotation=90,
        fontsize=8, color=S.INK_NONE)
ax.text(C_TR - 0.06, 0.10, "well transition", ha="right", va="bottom",
        rotation=90, fontsize=6.8, color=S.INK_RED)

# --- s3 LSDA over-polarisation caveat band (top), kept neutral ---
s3 = mag[(mag.alkali == "Na") & (mag.cell == "s3")]["abs_mag_muB"]
ax.axhspan(s3.min(), s3.max(), color=S.FILL_NONE, lw=0, zorder=0.5)
ax.text(9.85, s3.max() - 0.04, r"$\sqrt{3}$ (excluded)", ha="right", va="top",
        fontsize=7.0, color=S.C_SEC)

for el in ("Na", "Li"):
    st = S.SERIES[(el, "1x1")]
    cs, mean, lo, hi = series_stats(el, "1x1")
    ax.fill_between(cs, lo, hi, color=st["color"], alpha=0.16, lw=0, zorder=2)
    ax.plot(cs, mean, st["marker"] + "-", ms=5.5, lw=1.4, color=st["color"],
            mfc=st["color"], mew=0, zorder=4, label=st["label"])

ax.set_xlim(5.3, 10.0)
ax.set_ylim(-0.05, 2.82)
ax.set_xlabel(r"CoO$_2$$-$CoO$_2$ spacing  $c$  (Å)")
ax.set_ylabel(r"cell moment  $\langle|m|\rangle$  ($\mu_B$)")
leg = ax.legend(loc="upper left", bbox_to_anchor=(0.16, 0.80),
                handletextpad=0.4, fontsize=7.6, frameon=True,
                facecolor="white", edgecolor=S.C_GRID, framealpha=0.92)
leg.get_frame().set_linewidth(0.6)
S.thin_spines(ax)

fig.subplots_adjust(left=0.145, right=0.97, top=0.97, bottom=0.145)
S.save(fig, "fig4")
print("fig4 written")

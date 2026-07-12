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

C_SHADE = "#e2e1db"     # neutral light-grey fill (sobriety: no decorative colour)

# --- s3 LSDA over-polarisation caveat band (top) ---
s3 = mag[(mag.alkali == "Na") & (mag.cell == "s3")]["abs_mag_muB"]
ax.axhspan(s3.min(), s3.max(), color=C_SHADE, alpha=0.9, lw=0, zorder=0)
ax.text(9.85, s3.max(), r"Na $\sqrt{3}$: LSDA over-polarised," "\n"
        r"$|m|\!\approx\!$const (excluded)", ha="right", va="top",
        fontsize=7.0, color=S.C_SEC)

# --- well-transition band (consistent with Fig. 2) ---
ax.axvspan(6.35, 6.45, color=C_SHADE, alpha=0.9, lw=0, zorder=0)
ax.text(6.31, 1.55, "well transition", ha="right", va="center", rotation=90,
        fontsize=7.2, color=S.C_SEC)

for el in ("Na", "Li"):
    st = S.SERIES[(el, "1x1")]
    cs, mean, lo, hi = series_stats(el, "1x1")
    ax.fill_between(cs, lo, hi, color=st["color"], alpha=0.16, lw=0, zorder=2)
    ax.plot(cs, mean, st["marker"] + "-", ms=5.5, lw=1.4, color=st["color"],
            mfc=st["color"], mew=0, zorder=4, label=st["label"])

ax.annotate("moment switches on\n" r"at $\alpha<0$", xy=(6.99, 0.34),
            xytext=(7.95, 0.06), fontsize=7.2, color=S.C_SEC, ha="center",
            arrowprops=dict(arrowstyle="->", color=S.C_SEC, lw=0.7))
ax.text(9.85, 1.66, r"gallery moment"
        "\n" r"antiparallel to Co"
        "\n" r"($m_{\mathrm{alkali}}\!<\!0$)", fontsize=7.2, color=S.C_MUT,
        va="center", ha="right", style="italic", linespacing=1.4)

ax.set_xlim(5.3, 10.0)
ax.set_ylim(-0.05, 2.82)
ax.set_xlabel(r"CoO$_2$$-$CoO$_2$ spacing  $c$  (Å)")
ax.set_ylabel(r"cell moment  $\langle|m|\rangle$  ($\mu_B$)")
ax.legend(loc="lower left", bbox_to_anchor=(0.035, 0.62), handletextpad=0.4,
          fontsize=7.6)
S.thin_spines(ax)

fig.subplots_adjust(left=0.145, right=0.97, top=0.97, bottom=0.145)
S.save(fig, "fig4")
print("fig4 written")

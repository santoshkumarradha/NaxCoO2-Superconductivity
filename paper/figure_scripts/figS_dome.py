#!/usr/bin/env python
"""
Figure S (experimental dome + magnetic neighbour) - figS_dome.
The real superconducting dome of Na_xCoO2.yH2O plotted against its physically
correct abscissa -- the titrated Co valence -- with the magnetic phase that
borders it.

Data: reanalysis/experimental_dome.csv.  Following EXPERIMENTAL_DOME.md, the
only rows carrying a *redox-titrated* Co valence (the defensible charge axis;
naive 4-x is superseded, see the .md) are Milne et al., PRL 93, 247007 (2004);
those 8 points are plotted.  Marker fill distinguishes the deuteration of the
sample (D2O / H2O / mixed) -- an incidental label, not a variable of interest.

Two pastel domains (house style, physics labelled inside, flat fills):
 * amber  = the SC pairing window (optimal Co valence 3.24-3.35, Milne);
 * red    = the magnetic-ordered phase that borders the dome on the
            over-oxidised (high-valence) side, where Tc collapses above ~3.36
            and NaxCoO2.yH2O gives way to magnetic order (Sakurai 2015; the
            boundary is schematic/abstract-sourced, see caption).

Takeaway (caption): magnetism borders the dome -- with real, titrated data.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

import _style as S

S.use_house_style()

REPO = Path(__file__).resolve().parents[2]
e = pd.read_csv(REPO / "reanalysis" / "experimental_dome.csv")
m = e[(e.source == "Milne2004") & e.Co_valence_titrated.notna()].copy()

VLO, VHI = 3.18, 3.50
V_SC0, V_MAG = 3.205, 3.36        # SC-window / magnetic-order boundary

fig, ax = plt.subplots(figsize=(3.5, 3.35))

# ---- domains (labelled inside) --------------------------------------
S.vspan_region(ax, VLO, V_SC0, S.FILL_NONE, "underdoped", ytext=0.09,
               fontsize=7.4, rotation=90)
S.vspan_region(ax, V_SC0, V_MAG, S.FILL_AMBER, "SC pairing window",
               ytext=0.09, fontsize=7.8)
ax.axvspan(V_MAG, VHI, color=S.FILL_RED, lw=0, zorder=0)
ax.text(3.487, 0.5, "magnetic order (Sakurai 2015)", rotation=90,
        transform=ax.get_xaxis_transform(), ha="center", va="center",
        color=S.INK_RED, fontsize=7.6, zorder=1.5)

# ---- dome guide: quadratic through the titrated points --------------
vv = m.Co_valence_titrated.values
tt = m.Tc_K.values
cf = np.polyfit(vv, tt, 2)
xg = np.linspace(3.205, 3.47, 200)
ax.plot(xg, np.polyval(cf, xg), color=S.C_SEC, lw=1.1, ls=(0, (5, 3)),
        zorder=3)

# ---- data, marker fill by deuteration -------------------------------
styles = {"D2O": ("o", S.C_BLUE, "D$_2$O"),
          "H2O": ("s", "white", "H$_2$O"),
          "(H+D)2O": ("D", "white", "(H+D)$_2$O")}
seen = set()
for _, r in m.iterrows():
    mk, fc, lab = styles[r.hydration]
    ax.plot(r.Co_valence_titrated, r.Tc_K, mk, ms=6.5, color=S.C_BLUE,
            mfc=(S.C_BLUE if fc != "white" else "white"), mew=1.4, zorder=5,
            label=(lab if lab not in seen else None))
    seen.add(lab)

ax.set_xlim(VLO, VHI)
ax.set_ylim(0, 5.6)
ax.set_xlabel(r"titrated Co valence  $v_{\mathrm{Co}}$")
ax.set_ylabel(r"$T_c$  (K)")
ax.legend(loc="upper left", title="Milne 2004", handletextpad=0.4,
          labelspacing=0.3, borderaxespad=0.7, fontsize=7.0,
          ncol=3, columnspacing=0.8, handlelength=1.0)
leg = ax.get_legend()
leg.get_title().set_fontsize(7.2)
S.thin_spines(ax)

fig.subplots_adjust(left=0.155, right=0.965, top=0.965, bottom=0.145)
S.save(fig, "figS_dome")
print("figS_dome written")

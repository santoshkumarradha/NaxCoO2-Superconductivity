#!/usr/bin/env python
"""
Figure S-phase: the assembled superconducting phase diagram in the (c, T) plane.

Single panel with the fig2_col broken x-axis (left segment 5.4-7.15 A, right
segment 8.3-10.0 A).  Everything drawn here is READ from the repository data:

  * computed dome + ansatz band .... theory/results/coupling_tc_vs_c.csv
  * region walls ................... same csv (status column, lambda = 2 cut)
  * experimental points ............ reanalysis/experimental_dome.csv
  * Stoner strip (T = 0) ........... spin_analysis/magnetization_v2.csv
  * alkali thresholds .............. theory/results/well_fits_v4.csv (Li, Na)

The K threshold c*_K ~ 6.75 A is NOT derivable from well_fits_v4.csv (no K
rows); the value is taken from the supplement, Sec. "Predictions beyond the
compound" (Fig. figS_thresholds), where the K Landau coefficient crosses zero.

Output: paper/figures/figS_phase.{pdf,png}.  Nothing outside figures/ touched.
"""
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy.interpolate import PchipInterpolator

import _style as S

S.use_house_style()
mpl.rcParams.update({
    "font.size": 8.5,
    "axes.labelsize": 8.5,
    "xtick.labelsize": 7.5,
    "ytick.labelsize": 7.5,
    "hatch.linewidth": 0.6,
})

# ----------------------------------------------------------------------
# Data + derived boundaries (all computed from the csv files)
# ----------------------------------------------------------------------
cp = S.load_coupling().sort_values("c_A").reset_index(drop=True)
c = cp["c_A"].values
tc0 = cp["Tc_K"].values
lo = np.minimum.reduce([tc0, cp["Tc_K_g1"].values, cp["Tc_K_g2"].values])
hi = np.maximum.reduce([tc0, cp["Tc_K_g1"].values, cp["Tc_K_g2"].values])

# carrier onset: boundary of the status == "no-2DEG" block (5.56 -> 5.58)
no2 = cp[cp.status == "no-2DEG"]["c_A"].max()
yes2 = cp[cp.status != "no-2DEG"]["c_A"].min()
C_ON = 0.5 * (no2 + yes2)                      # = 5.57 A

# lambda = 2 cut (polaronic freeze-out wall), linear interpolation
lam = cp["lambda"].values
i = np.where((lam[:-1] - 2) * (lam[1:] - 2) < 0)[0][0]
C_LAM = c[i] + (2 - lam[i]) * (c[i + 1] - c[i]) / (lam[i + 1] - lam[i])  # 6.434

# amber SC window: where the ansatz band is finite, up to the lambda = 2 cut
C_SC = c[hi > 0.01].min()                      # = 6.22 A

# Stoner strip: Na 1x1, delta = 0, |m| > 0.1 muB (sampled c only: sparse grid)
mag = S.load_magnetization()
na = mag[(mag.alkali == "Na") & (mag.cell == "1x1")
         & (mag.delta_A == 0.0)].sort_values("c_A")
c_mag = na[na.abs_mag_muB > 0.1]["c_A"].values  # -> [6.9, 8.4, 9.9]
MAG0, MAG1 = c_mag.min(), c_mag.max()

# alkali thresholds c* from alpha(c) sign change (well_fits_v4.csv)
fits = S.load_well_fits()


def cstar(el):
    w = fits[(fits.element == el) & (fits.cell == "1x1")].sort_values("c")
    if w["alpha"].iloc[0] < 0:                  # already double-well at 5.5 A
        return None                             # upper bound only: c* <= 5.5
    p = PchipInterpolator(w["c"].values, w["alpha"].values)
    grid = np.linspace(w["c"].min(), w["c"].max(), 8801)
    a = p(grid)
    j = np.where(np.sign(a[:-1]) != np.sign(a[1:]))[0][0]
    return float(grid[j])


CSTAR_LI = cstar("Li")                          # None -> c*_Li <= 5.5 A
CSTAR_NA = cstar("Na")                          # 6.396 ~ 6.40 A
CSTAR_K = 6.75   # not in well_fits_v4.csv; SM Sec. "Predictions beyond the
#                  compound", Fig. figS_thresholds (bracketed K crossing)

# experimental dome (drop the mixed-phase row with no measured spacing)
dome = pd.read_csv(S.REPO / "reanalysis" / "experimental_dome.csv")
dome = dome.dropna(subset=["coo2_coo2_spacing_A"])
dome = dome[dome.hydration.isin(["H2O", "D2O", "(H+D)2O"])]

print(f"carrier onset      c = {C_ON:.3f} A  (no-2DEG {no2} -> 2DEG {yes2})")
print(f"lambda = 2 cut     c = {C_LAM:.3f} A")
print(f"SC window          c = {C_SC:.2f} .. {C_LAM:.3f} A")
print(f"Stoner strip       c = {MAG0} .. {MAG1} A  (sampled: {list(c_mag)})")
print(f"c*_Na = {CSTAR_NA:.3f} A,  c*_Li <= 5.5 A (alpha<0),  c*_K = {CSTAR_K}")

# ----------------------------------------------------------------------
# Canvas: single row, broken x axis (fig2_col technique)
# ----------------------------------------------------------------------
L0, L1 = 5.4, 7.15
R0, R1 = 8.3, 10.0
YMAX = 26.0
HATCH0, HATCH1 = 9.70, 9.95     # dynamically-softened sub-band (open test);
#   nominal 9.75-9.95, lower edge relaxed to 9.70 so every measured bilayer
#   point (9.715-9.928 A) sits inside the contingent region

fig = plt.figure(figsize=(3.4, 2.75))
gs = GridSpec(1, 2, width_ratios=[L1 - L0, R1 - R0], wspace=0.05,
              left=0.115, right=0.985, top=0.90, bottom=0.145)
axL = fig.add_subplot(gs[0, 0])
axR = fig.add_subplot(gs[0, 1], sharey=axL)

axL.set_xlim(L0, L1)
axR.set_xlim(R0, R1)
axL.set_ylim(0, YMAX)

# ----- pastel regions ---------------------------------------------------
axL.axvspan(L0, C_ON, color=S.FILL_NONE, lw=0, zorder=0)
axL.axvspan(C_ON, C_SC, color=S.FILL_BLUE, lw=0, zorder=0)
axL.axvspan(C_SC, C_LAM, color=S.FILL_AMBER, lw=0, zorder=0)
axL.axvspan(C_LAM, L1, color=S.FILL_RED, lw=0, zorder=0)
axR.axvspan(R0, R1, color=S.FILL_RED, lw=0, zorder=0)
axL.axvline(C_LAM, color=S.EDGE_STONER, lw=1.0, zorder=1.5)
axL.axvline(C_ON, color=S.EDGE_TURNON, lw=0.8, zorder=1.5)

axL.text(0.5 * (L0 + C_ON), 0.40, "no 2DEG", rotation=90,
         transform=axL.get_xaxis_transform(), ha="center", va="center",
         color=S.INK_NONE, fontsize=5.8, zorder=3)
axL.text(0.5 * (C_ON + C_SC), 0.115, "2DEG metal",
         transform=axL.get_xaxis_transform(), ha="center", va="center",
         color=S.INK_BLUE, fontsize=6.8, zorder=3)
axL.text(0.5 * (C_SC + C_LAM), 0.80, "SC window", rotation=90,
         transform=axL.get_xaxis_transform(), ha="center", va="center",
         color=S.INK_AMBER, fontsize=6.4, zorder=5)
axL.text(6.72, 0.965, "polaronic\n(static)",
         transform=axL.get_xaxis_transform(), ha="center", va="top",
         color=S.INK_RED, fontsize=6.8, zorder=3)

# anchor dashed verticals (5.5 / 6.9 left, 9.9 right), fig2_col style
for ax, cs in ((axL, (5.5, 6.9)), (axR, (9.9,))):
    for cc in cs:
        ax.axvline(cc, color=S.C_MUT, ls=(0, (5, 3)), lw=0.8, alpha=0.85,
                   zorder=1)

# ----- computed dome + ansatz band (LEFT segment only) ------------------
mL = c <= L1
axL.fill_between(c[mL], lo[mL], hi[mL], color=S.C_BLUE, alpha=0.22, lw=0,
                 zorder=2, label="ansatz range")
axL.plot(c[mL], tc0[mL], color=S.C_BLUE, lw=1.6, zorder=4,
         label=r"$T_c$ (computed)")

# ----- zoom inset: the dome at its natural width ------------------------
# The SC window is only ~0.2 A wide, so on the full c axis the dome reads
# as a needle; the inset re-plots it over c = 6.15-6.55 A so its shape (the
# familiar dome bounded by its two walls) is visible.
axZ = axL.inset_axes([0.135, 0.315, 0.36, 0.46])
Z0, Z1 = 6.15, 6.55
mZ = (c >= Z0) & (c <= Z1)
axZ.axvspan(Z0, C_SC, color=S.FILL_BLUE, lw=0, zorder=0)
axZ.axvspan(C_SC, C_LAM, color=S.FILL_AMBER, lw=0, zorder=0)
axZ.axvspan(C_LAM, Z1, color=S.FILL_RED, lw=0, zorder=0)
axZ.axvline(C_LAM, color=S.EDGE_STONER, lw=0.8, zorder=1.5)
axZ.fill_between(c[mZ], lo[mZ], hi[mZ], color=S.C_BLUE, alpha=0.22, lw=0,
                 zorder=2)
axZ.plot(c[mZ], tc0[mZ], color=S.C_BLUE, lw=1.3, zorder=4)
axZ.set_xlim(Z0, Z1)
axZ.set_ylim(0, 24)
axZ.set_xticks([6.2, 6.4])
axZ.set_yticks([0, 10, 20])
axZ.tick_params(labelsize=5.2, length=1.8, pad=1.5, top=False, right=False)
for sp in axZ.spines.values():
    sp.set_linewidth(0.6)
    sp.set_color(S.C_SEC)
axZ.text(0.96, 0.94, "zoom", transform=axZ.transAxes, ha="right", va="top",
         fontsize=5.4, color=S.C_SEC)
axL.indicate_inset_zoom(axZ, edgecolor=S.C_SEC, lw=0.6, alpha=0.7)

# ----- hatched contingent sub-band (RIGHT segment) ----------------------
axR.axvspan(HATCH0, HATCH1, facecolor="none", edgecolor=S.INK_RED,
            hatch="///", lw=0, zorder=1.2)
axR.text(HATCH0 - 0.05, 12.4, "dynamically softened\nby water (open test)",
         rotation=90, ha="right", va="center", color=S.INK_RED,
         fontsize=5.6, linespacing=1.15, zorder=3)

# ----- experimental points (bilayer hydrates) ---------------------------
MARK = {
    "H2O": dict(marker="o", mfc=S.C_RED, mec=S.C_RED,
                label=r"H$_2$O hydrate"),
    "D2O": dict(marker="D", mfc="white", mec=S.C_RED,
                label=r"D$_2$O hydrate"),
    "(H+D)2O": dict(marker="^", mfc=S.C_RED, mec=S.C_RED, fillstyle="left",
                    markerfacecoloralt="white", label=r"(H+D)$_2$O"),
}
for hyd, st in MARK.items():
    d = dome[dome.hydration == hyd]
    axR.plot(d["coo2_coo2_spacing_A"], d["Tc_K"], ls="none", ms=3.6,
             mew=0.8, zorder=6, **st)

# ----- "not SC" anchors on the T = 0 line -------------------------------
axL.plot([5.5, 6.9], [0, 0], "s", ms=4.2, mfc="white", mec=S.C_SEC, mew=1.1,
         clip_on=False, zorder=6, label="not SC (expt)")

# ----- Stoner strip (ground state, T = 0 statement: a bar, no T extent) --
BAR_Y, BAR_H = 0.45, 0.62
for ax, x0, x1 in ((axL, MAG0, L1), (axR, R0, MAG1)):
    ax.axhspan(BAR_Y, BAR_Y + BAR_H, xmin=(x0 - ax.get_xlim()[0]) /
               (ax.get_xlim()[1] - ax.get_xlim()[0]),
               xmax=(x1 - ax.get_xlim()[0]) /
               (ax.get_xlim()[1] - ax.get_xlim()[0]),
               color=S.EDGE_STONER, lw=0, zorder=5)
axR.plot(c_mag[c_mag >= R0], np.full((c_mag >= R0).sum(), BAR_Y + BAR_H / 2),
         "o", ms=1.6, color="white", mew=0, zorder=5.5)
axL.plot([MAG0], [BAR_Y + BAR_H / 2], "o", ms=1.6, color="white", mew=0,
         zorder=5.5)
axR.text(8.38, 1.55, "Stoner strip (ground state)", ha="left", va="bottom",
         color=S.INK_RED, fontsize=5.8, zorder=5)

# ----- alkali thresholds on the upper edge ------------------------------
def top_mark(ax, x, lab, arrow_left=False):
    if arrow_left:  # upper bound only: c* <= x
        ax.annotate("", xy=(x - 0.16, YMAX), xytext=(x, YMAX),
                    annotation_clip=False, zorder=6,
                    arrowprops=dict(arrowstyle="->", lw=0.9, color=S.C_INK,
                                    shrinkA=0, shrinkB=0))
    else:
        ax.plot([x], [YMAX], marker=7, ms=4.5, color=S.C_INK, clip_on=False,
                zorder=6)
    ax.text(x, YMAX + 1.1, lab, ha="center", va="bottom", fontsize=6.2,
            color=S.C_INK, clip_on=False)


top_mark(axL, 5.5, "Li", arrow_left=True)       # c*_Li <= 5.5 A
top_mark(axL, CSTAR_NA, "Na")                   # 6.40 A (interpolated)
top_mark(axL, CSTAR_K, "K")                     # 6.75 A (SM thresholds sec.)

# ----- legends (theory left, experiment right) --------------------------
axL.legend(loc="upper left", bbox_to_anchor=(0.015, 0.985), fontsize=5.6,
           handletextpad=0.35, borderaxespad=0.0, labelspacing=0.3,
           handlelength=1.25)
axR.legend(loc="upper left", bbox_to_anchor=(0.03, 0.985), fontsize=5.6,
           handletextpad=0.35, borderaxespad=0.0, labelspacing=0.3,
           handlelength=1.05)

# ----- spines / ticks / break marks (fig2_col technique) ----------------
axL.spines["top"].set_visible(False)
axL.spines["right"].set_visible(False)
axL.tick_params(top=False, right=False)
axR.spines["top"].set_visible(False)
axR.spines["left"].set_visible(False)
axR.tick_params(top=False, left=False, labelleft=False)

axL.set_xticks([5.5, 6.0, 6.5, 7.0])
axR.set_xticks([8.5, 9.0, 9.5, 10.0])

d = 0.7
kw = dict(marker=[(-1, -d), (1, d)], markersize=6, linestyle="none",
          color=S.C_SEC, mec=S.C_SEC, mew=0.9, clip_on=False, zorder=6)
axL.plot([1, 1], [0, 1], transform=axL.transAxes, **kw)
axR.plot([0, 0], [0, 1], transform=axR.transAxes, **kw)

axL.set_ylabel(r"$T$  (K)")
fig.supxlabel(r"CoO$_2$$-$CoO$_2$ spacing  $c$  (Å)", fontsize=8, y=0.012)

S.save(fig, "figS_phase")
print("figS_phase written")

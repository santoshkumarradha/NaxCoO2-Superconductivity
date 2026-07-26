#!/usr/bin/env python
"""
Figure 10 - The gallery-2DEG spin phase diagram (Stoner model on the SSH4
gallery band), regenerated as a paper figure in the house style and anchored to
the v4 DFT magnetisation data.

Physics + model: spin_analysis/stoner_model.py (imported, not duplicated) --
the Radha-Lambrecht SSH4 surface-charge model with a Stoner criterion on the
gallery band.  The (x, c) plane is tiled into four pastel domains:
  neutral  no 2DEG (gallery band empty)
  red      polarized 2DEG (magnetic order; kills singlet SC; equal-spin f-wave)
  amber    paramagnon-boosted paramagnet (Stoner S > S*)
  blue     plain paramagnetic 2DEG (phonon singlet + paramagnon channel)
The physics is labelled INSIDE each fill in matching darker ink; the two side
panels are cuts at x = 0.35 sharing the c axis.

DFT anchors (this figure's addition over the standalone model):
  * Na 1x1 (x=1) local moment ONSET from magnetization_v2.csv -- |m| = 0 at
    c = 5.5 A, turns on by 6.9 A: it switches on exactly where the model's
    carrier turn-on line sits at full coverage, i.e. the polarized strip hugs
    the turn-on line (drawn as a red onset bar at x = 1).
  * the sqrt3 x sqrt3 (x = 1/3) cell is LSDA over-polarized and is flagged as an
    excluded caveat (the wells themselves shift < 6% under the spin treatment).
The gallery moment is antiparallel to Co (stated in the polarized-region label).

Aesthetic: _style.py house rcParams (serif, thin spines) + the region-fill
vocabulary; identical fills/inks to fig2/fig4/fig5.
"""
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "spin_analysis"))
import _style as S
import stoner_model as SM

S.use_house_style()

FILLS = [S.FILL_NONE, S.FILL_RED, S.FILL_AMBER, S.FILL_BLUE]   # region codes 0..3
P = SM.P

# ---------------------------------------------------------------- model grid
x = np.linspace(0.12, 1.03, 460)
c = np.linspace(5.5, 10.6, 460)
X, C, q, F, Sf, region, w_alk = SM.classify(x, c)

two_delta = 2.0 * P["delta_Li"] * P["delta_Na_over_Li"]
gamma_co = 6.0 * P["t_Co_xy"]
lv_on = two_delta - gamma_co
lv_st = P["I_s"] * P["nu"]
lv_ss = lv_st * P["S_star"] / (P["S_star"] - 1.0)

# ---------------------------------------------------------------- DFT anchors
mag = S.load_magnetization()
na = mag[(mag.alkali == "Na") & (mag.cell == "1x1")].groupby("c_A")["abs_mag_muB"].mean()
# onset bracket: last c with |m|~0 and first c with |m|>0
c_zero = float(na.index[np.where(na.values < 0.02)[0].max()])   # 5.5
c_on_dft = float(na.index[np.where(na.values > 0.05)[0].min()])  # 6.9

# ---------------------------------------------------------------- figure
fig = plt.figure(figsize=(7.5, 3.4))
gs = fig.add_gridspec(2, 3, width_ratios=[1.68, 1.0, 1.05],
                      hspace=0.15, wspace=0.50)
ax = fig.add_subplot(gs[:, 0])          # (a) phase diagram
axq = fig.add_subplot(gs[0, 1])         # (b) q(c) cut
axs = fig.add_subplot(gs[1, 1], sharex=axq)   # (b) S(c) cut
axm = fig.add_subplot(gs[:, 2])         # (c) DFT moment turn-on

# --- (a) phase diagram --------------------------------------------------
cmap = ListedColormap(FILLS)
ax.pcolormesh(X, C, region, cmap=cmap, vmin=-0.5, vmax=3.5, shading="auto",
              rasterized=True)

for lv, col, ls in ((lv_on, S.EDGE_TURNON, "-"), (lv_st, S.EDGE_STONER, "-"),
                    (lv_ss, S.EDGE_SSTAR, "--")):
    ax.contour(X, C, w_alk, levels=[lv], colors=col, linewidths=1.2,
               linestyles=ls)

# Stoner-enhancement contours inside the plain-PM (blue) region
S_plot = np.where(region == 3, Sf, np.nan)
cs = ax.contour(X, C, S_plot, levels=[2.0], colors=["#5598e7"], linewidths=0.8)
ax.clabel(cs, fmt=lambda v: f"S={v:g}", fontsize=7.0, colors=S.C_SEC)

# minimal region tags (phase name only; full reading is in the caption)
ax.text(0.235, 6.5, "no 2DEG", color=S.INK_NONE, ha="center", fontsize=8)
ax.text(0.50, 6.05, "polarized\n2DEG", color=S.INK_RED, ha="center",
        fontsize=7.6, linespacing=1.15)
ax.text(0.66, 7.35, "$S>S^\\ast$", color=S.INK_AMBER, ha="center", fontsize=7.6)
ax.text(0.66, 9.0, "paramagnetic\n2DEG", color=S.INK_BLUE, ha="center",
        fontsize=7.6, linespacing=1.15)

# experimental anchors at x = 0.35 (short tags; detail in caption)
ax.axvline(0.35, color=S.C_SEC, lw=0.7, ls=":")
ax.plot([0.35], [6.9], marker="x", ms=7, mew=1.8, color=S.C_INK, zorder=5)
ax.annotate("6.9 Å", (0.35, 6.9), xytext=(0.155, 6.75), fontsize=7,
            color=S.C_INK)
ax.plot([0.35], [9.9], marker="o", ms=6, mfc="#104281", mec="white", zorder=5)
ax.annotate("9.9 Å", (0.35, 9.9), xytext=(0.155, 9.9), fontsize=7,
            color=S.C_INK)

# DFT magnetisation anchor: Na 1x1 moment onset bar in the empty right space
xd = 1.0
ax.annotate("", xy=(xd, c_on_dft), xytext=(xd, c_zero + 0.18),
            arrowprops=dict(arrowstyle="-|>", color=S.INK_RED, lw=2.0))
ax.plot([xd], [c_zero + 0.18], marker="o", ms=4.2, mfc="white", mec=S.INK_RED,
        mew=1.2, zorder=6)
ax.plot([xd], [c_on_dft], marker="o", ms=4.8, color=S.INK_RED, mec="white",
        mew=0.6, zorder=6)
ax.text(0.90, 7.25, "DFT $|m|$\nonset", fontsize=7.0, color=S.INK_RED,
        ha="center", va="bottom", linespacing=1.15)

ax.set_xlabel("Na content $x$")
ax.set_ylabel(r"CoO$_2$ plane spacing $c$ (Å)  [gallery opening]")
ax.set_xlim(x[0], x[-1])
ax.set_ylim(c[0], c[-1])
S.panel_title(ax, "(a) Gallery-2DEG spin phase diagram")

ax2 = ax.twinx()                     # same c, mapped to the model half-width
ax2.set_ylim(ax.get_ylim())
ticks = np.arange(6, 11)
ax2.set_yticks(ticks)
ax2.set_yticklabels([f"{g:.1f}" for g in SM.gamma0_of_c(ticks)])
ax2.set_ylabel(r"$\Gamma_0(c)$: full-coverage gallery half-width (eV)",
               fontsize=7.2, color=S.C_SEC)
ax2.tick_params(colors=S.C_SEC, labelsize=7)

# boundary-stroke legend
handles = [Line2D([], [], color=S.EDGE_TURNON, lw=1.2,
                  label=r"carrier turn-on  $\Gamma=2\delta$"),
           Line2D([], [], color=S.EDGE_STONER, lw=1.2,
                  label=r"Stoner  $I_s N(0)=1$"),
           Line2D([], [], color=S.EDGE_SSTAR, lw=1.2, ls="--",
                  label=f"$S=S^\\ast={P['S_star']:g}$")]
ax.legend(handles=handles, loc="upper right", fontsize=7.0, frameon=True,
          facecolor="white", edgecolor=S.C_GRID, framealpha=0.9,
          handlelength=1.4, borderpad=0.4).get_frame().set_linewidth(0.6)

# --- (b) cuts at x = 0.35 -----------------------------------------------
f35 = float(SM.f_dilution(0.35, P["a_Na"], P["lam_Na"]))
cc = np.linspace(5.5, 10.6, 500)
w35 = SM.gamma0_of_c(cc) * f35
q35 = SM.q_carrier(w35, gamma_co, two_delta)
with np.errstate(divide="ignore"):
    F35 = P["I_s"] * P["nu"] / w35
    S35 = np.where(F35 < 1, 1.0 / (1.0 - F35), np.nan)

b = SM.boundaries_cut(f35)
for a_ in (axq, axs):
    a_.axvspan(5.5, b["on"], color=S.FILL_NONE)
    a_.axvspan(b["on"], b["st"], color=S.FILL_RED)
    a_.axvspan(b["st"], b["ss"], color=S.FILL_AMBER)
    a_.axvspan(b["ss"], 10.6, color=S.FILL_BLUE)
    a_.set_xlim(5.5, 10.6)
    a_.grid(color=S.C_GRID, lw=0.5)
    a_.set_axisbelow(True)
    S.thin_spines(a_)
    for c_ref, mk in ((6.9, "x"), (9.9, "o")):
        yv = np.interp(c_ref, cc, q35 if a_ is axq else
                       np.nan_to_num(S35, nan=0.0))
        a_.plot([c_ref], [yv], marker=mk, color=S.C_INK, ms=5,
                mfc="none" if mk == "o" else S.C_INK, mew=1.4, zorder=5)

axq.plot(cc, q35, color=S.C_BLUE, lw=1.8)
axq.set_ylabel("q per Na", fontsize=8)
axq.tick_params(labelbottom=False)
S.panel_title(axq, "(b) cut at $x=0.35$", fontsize=8.5, pad=3)

axs.plot(cc, S35, color=S.C_BLUE, lw=1.8)
axs.set_ylim(0, 8)
axs.set_ylabel("Stoner $S$", fontsize=8)
axs.set_xlabel("c (Å)", fontsize=8)

# --- (c) DFT cell-moment turn-on (folded in from the standalone fig4) -----
# Cell-integrated |m|(c) for Na and Li 1x1 (LSDA spin scans): the local moment
# switches on exactly at the alpha<0 well transition -- the DFT counterpart of
# (a)'s red polarized strip.  Same region-fill vocabulary; c-axis matched to (b).
C_TR = 6.40                              # well transition / moment turn-on
S.vspan_region(axm, 5.5, C_TR, S.FILL_NONE)
S.vspan_region(axm, C_TR, 10.6, S.FILL_RED, "magnetic", ytext=0.62,
               ha="center", fontsize=8.0)
axm.axvline(C_TR, color=S.EDGE_STONER, lw=1.0, zorder=1)
axm.text(5.98, 1.35, "no moment", ha="center", va="center", rotation=90,
         fontsize=7.6, color=S.INK_NONE)
axm.text(C_TR - 0.07, 0.12, "well transition", ha="right", va="bottom",
         rotation=90, fontsize=7.0, color=S.INK_RED)

# sqrt3 LSDA over-polarisation caveat band (kept neutral; excluded)
s3 = mag[(mag.alkali == "Na") & (mag.cell == "s3")]["abs_mag_muB"]
axm.axhspan(s3.min(), s3.max(), color=S.FILL_NONE, lw=0, zorder=0.5)
axm.text(10.45, s3.max() - 0.05, r"$\sqrt{3}$ (excluded)", ha="right", va="top",
         fontsize=7.2, color=S.C_SEC)

for el in ("Na", "Li"):
    st = S.SERIES[(el, "1x1")]
    d = mag[(mag.alkali == el) & (mag.cell == "1x1")]
    cs_, mean_, lo_, hi_ = [], [], [], []
    for cv, g in d.groupby("c_A"):
        v = g["abs_mag_muB"].values
        cs_.append(cv); mean_.append(v.mean()); lo_.append(v.min()); hi_.append(v.max())
    o = np.argsort(cs_)
    cs_ = np.array(cs_)[o]; mean_ = np.array(mean_)[o]
    lo_ = np.array(lo_)[o]; hi_ = np.array(hi_)[o]
    axm.fill_between(cs_, lo_, hi_, color=st["color"], alpha=0.16, lw=0, zorder=2)
    axm.plot(cs_, mean_, st["marker"] + "-", ms=5.0, lw=1.4, color=st["color"],
             mfc=st["color"], mew=0, zorder=4, label=st["label"])

axm.set_xlim(5.5, 10.6)                  # c-axis matched to (b)
axm.set_ylim(-0.05, 2.82)
axm.set_xlabel("c (Å)", fontsize=8)
axm.set_ylabel(r"cell moment  $\langle|m|\rangle$  ($\mu_B$)", fontsize=8)
leg = axm.legend(loc="upper left", bbox_to_anchor=(0.03, 0.78),
                 handletextpad=0.4, fontsize=7.2, frameon=True,
                 facecolor="white", edgecolor=S.C_GRID, framealpha=0.92)
leg.get_frame().set_linewidth(0.6)
S.thin_spines(axm)
S.panel_title(axm, "(c) DFT $|m|$ turn-on", fontsize=8.5, pad=3)

fig.subplots_adjust(left=0.068, right=0.965, top=0.90, bottom=0.13)
S.save(fig, "fig10")
print("fig10 written; DFT Na1x1 onset bracket c =", c_zero, "->", c_on_dft)

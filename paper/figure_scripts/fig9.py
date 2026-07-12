#!/usr/bin/env python
"""
Figure 9 (supplementary) - The gallery band at the physical composition x=1/3.

Na sqrt3 x sqrt3 cell, c = 5.5 vs 9.9 A, same Kohn-Sham spectral construction as
Fig. 7 (grayscale total intensity, HWHM 40 meV; blue Na fatband from projwfc.x).
This confirms the x=1 result at the real sodium coverage: no Na weight at E_F for
the dry spacing, a strongly Na-character band through E_F at 9.9 A.

HONEST FOLDING NOTE.  These bands are those of the SUPERCELL, whose BZ is 1/3 the
area of the 1x1 BZ and rotated 30 deg; every band is a fold of the primitive
zone.  We do NOT have the plane-wave coefficients needed for a weight-resolved
unfolding, so we do not draw an unfolded 1x1 dispersion.  Instead we mark the
known folding relations from the pw.in header: the 1x1 K point folds onto the
supercell Gamma and the 1x1 M point folds onto the supercell M.  The reader must
read every band here as a supercell (folded) band, not a primitive one.

Aesthetic-sobriety rule: blue marks Na character (data) only; intensity is
neutral grayscale, E_F / ticks / folding marks are neutral ink.
"""
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap, PowerNorm

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "theory" / "spectral"))
import _style as S
import spectral as SP

S.use_house_style()

GAMMA = 0.040
EWIN = (-3.0, 3.0)
JOBS = [("bands_Na_s3_c5.5", r"$c=5.5$ Å"),
        ("bands_Na_s3_c9.9", r"$c=9.9$ Å")]
GRAY = LinearSegmentedColormap.from_list("tot", ["#ffffff", S.C_INK])

data = {job: SP.load_spinsummed(job) for job, _ in JOBS}
kfine = None
omega = np.linspace(*EWIN, 700)
imgs, vmax = {}, 0.0
for job, _ in JOBS:
    d = data[job]
    if kfine is None:
        kfine = np.linspace(d["kdist"][0], d["kdist"][-1], 500)
    ones = np.ones_like(d["E"])
    imgs[job] = SP.spectral_image(d["kdist"], d["E"], ones, kfine, omega, GAMMA)
    vmax = max(vmax, np.percentile(imgs[job], 99.5))

fig, axes = plt.subplots(1, 2, figsize=(5.2, 3.35), sharey=True)
extent = [kfine[0], kfine[-1], omega[0], omega[-1]]
# supercell path labels; Gamma also carries the folded 1x1 K
LAB = [r"$\Gamma$", "M", "K", r"$\Gamma$"]

for ax, (job, title) in zip(axes, JOBS):
    d = data[job]
    ax.imshow(imgs[job], origin="lower", extent=extent, aspect="auto",
              cmap=GRAY, norm=PowerNorm(0.62, vmin=0, vmax=vmax),
              interpolation="bilinear", zorder=1)
    segs, lws = [], []
    for b in range(d["E"].shape[0]):
        eb, wb = d["E"][b], d["Na"][b]
        if eb.min() > EWIN[1] or eb.max() < EWIN[0]:
            continue
        for i in range(len(d["kdist"]) - 1):
            w = 0.5 * (wb[i] + wb[i + 1])
            if w < 0.04:
                continue
            segs.append([(d["kdist"][i], eb[i]), (d["kdist"][i + 1], eb[i + 1])])
            lws.append(0.3 + 6.0 * w)
    if segs:
        ax.add_collection(LineCollection(segs, linewidths=lws, colors=S.C_BLUE,
                                         alpha=0.9, capstyle="round", zorder=3))
    ax.axhline(0.0, color=S.C_SEC, lw=0.8, ls=(0, (5, 3)), zorder=2)
    ax.set_xticks(d["ticks"])
    ax.set_xticklabels(LAB)
    for xt in d["ticks"][1:-1]:
        ax.axvline(xt, color=S.C_GRID, lw=0.6, zorder=1.5)
    ax.set_xlim(kfine[0], kfine[-1])
    ax.set_ylim(*EWIN)
    ax.set_title(title, pad=4)
    ax.tick_params(top=False)

axes[0].set_ylabel(r"$E - E_F$  (eV)")
axes[1].text(0.97, 0.0, r"$E_F$", transform=axes[1].get_yaxis_transform(),
             ha="right", va="bottom", fontsize=7.5, color=S.C_SEC)
for j, ax in enumerate(axes):
    S.panel_label(ax, "ab"[j], x=0.0, y=1.03)

# folding annotations under the Gamma / M ticks (supercell -> which 1x1 point)
for ax in axes:
    tk = data[JOBS[0][0]]["ticks"]
    ax.annotate(r"$+\,1{\times}1\,K$", xy=(tk[0], EWIN[0]),
                xytext=(tk[0] + 0.03, EWIN[0] + 0.18), fontsize=6.0,
                color=S.C_SEC, ha="left", va="bottom")
    ax.annotate(r"$+\,1{\times}1\,M$", xy=(tk[1], EWIN[0]),
                xytext=(tk[1], EWIN[0] + 0.18), fontsize=6.0,
                color=S.C_SEC, ha="center", va="bottom")

# Na-weight key (white bbox; blue swatch = data), placed in the dry panel
kx = axes[0]
kx.add_patch(plt.Rectangle((0.03, 0.03, ), 0.40, 0.30, transform=kx.transAxes,
             facecolor="white", edgecolor=S.C_GRID, lw=0.5, zorder=5))
for frac, yy in ((1.0, 0.255), (0.5, 0.165), (0.15, 0.08)):
    kx.plot([0.07, 0.19], [yy, yy], transform=kx.transAxes, color=S.C_BLUE,
            lw=0.3 + 6.0 * frac, solid_capstyle="round", zorder=6)
    kx.text(0.22, yy, f"{frac:.2f}", transform=kx.transAxes, fontsize=6.3,
            va="center", color=S.C_INK, zorder=6)
kx.text(0.05, 0.305, "Na weight", transform=kx.transAxes, fontsize=6.6,
        va="center", color=S.C_INK, zorder=6)

fig.text(0.5, 0.005, r"supercell $\Gamma$ carries the folded $1{\times}1\,K$; "
         r"all bands are $\sqrt{3}{\times}\sqrt{3}$ (folded)", ha="center",
         fontsize=6.2, color=S.C_SEC)
fig.subplots_adjust(left=0.10, right=0.985, top=0.90, bottom=0.13, wspace=0.07)
S.save(fig, "fig9")
print("fig9 written; s3 ef/mag:",
      {j: (round(data[j]["ef"], 2), round(data[j]["mag"], 2)) for j, _ in JOBS})

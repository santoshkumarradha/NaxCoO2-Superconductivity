#!/usr/bin/env python
"""
Figure 7 - The gallery band appears (the ARPES-comparable reviewer figure).

Three panels sharing the E - E_F axis: Na 1x1 at c = 5.5, 6.9, 9.9 A.  Each
panel is a Kohn-Sham spectral map along Gamma-M-K-Gamma: every band (both
spins, summed) is Lorentzian-broadened (HWHM Gamma = 40 meV) onto an (k, omega)
grid and shown as GRAYSCALE total intensity (the ARPES-like image).  On top,
the Na-orbital projected weight from projwfc.x is drawn as a BLUE fatband
(linewidth proportional to Na character).  The eye should read it instantly:
at 5.5 A the gallery is empty and no band reaches E_F; by 9.9 A a strongly
Na-character band has descended through E_F -> the interlayer 2DEG.

Data: runpod/results_bands/bands_Na_1x1_c{5.5,6.9,9.9} (pw.bands.*.dat.gnu +
pw.proj.projwfc_{up,down}; E_F from pw.out).  No self-energy: honest DFT
spectral weight.  Aesthetic-sobriety rule: colour (blue) marks Na character
only; total intensity is neutral grayscale, E_F/ticks/labels are neutral ink.
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

GAMMA = 0.040                     # Lorentzian HWHM, eV (40 meV)
EWIN = (-3.0, 3.0)                # energy window, E - E_F (eV)
JOBS = [("bands_Na_1x1_c5.5", r"$c=5.5$ Å"),
        ("bands_Na_1x1_c6.9", r"$c=6.9$ Å"),
        ("bands_Na_1x1_c9.9", r"$c=9.9$ Å")]

# grayscale-total ramp: white -> ink (high spectral weight = dark)
GRAY = LinearSegmentedColormap.from_list("tot", ["#ffffff", S.C_INK])

# ------------------------------------------------------------------ load
data = {job: SP.load_spinsummed(job) for job, _ in JOBS}

kfine = None
omega = np.linspace(*EWIN, 700)
imgs = {}
vmax = 0.0
for job, _ in JOBS:
    d = data[job]
    if kfine is None:
        kfine = np.linspace(d["kdist"][0], d["kdist"][-1], 500)
    ones = np.ones_like(d["E"])
    tot = SP.spectral_image(d["kdist"], d["E"], ones, kfine, omega, GAMMA)
    imgs[job] = tot
    vmax = max(vmax, np.percentile(tot, 99.5))

# ------------------------------------------------------------------ plot
fig, axes = plt.subplots(1, 3, figsize=(7.1, 3.25), sharey=True)
extent = [kfine[0], kfine[-1], omega[0], omega[-1]]

for ax, (job, title) in zip(axes, JOBS):
    d = data[job]
    # grayscale total spectral weight (ARPES-like intensity)
    ax.imshow(imgs[job], origin="lower", extent=extent, aspect="auto",
              cmap=GRAY, norm=PowerNorm(0.62, vmin=0, vmax=vmax),
              interpolation="bilinear", zorder=1)

    # blue Na fatband overlay: linewidth prop. to Na projected weight
    segs, lws = [], []
    for b in range(d["E"].shape[0]):
        eb = d["E"][b]
        wb = d["Na"][b]
        if eb.min() > EWIN[1] or eb.max() < EWIN[0]:
            continue
        for i in range(len(d["kdist"]) - 1):
            w = 0.5 * (wb[i] + wb[i + 1])
            if w < 0.04:              # skip negligible Na character
                continue
            segs.append([(d["kdist"][i], eb[i]),
                         (d["kdist"][i + 1], eb[i + 1])])
            lws.append(0.3 + 6.0 * w)
    if segs:
        lc = LineCollection(segs, linewidths=lws, colors=S.C_BLUE,
                            alpha=0.9, capstyle="round", zorder=3)
        ax.add_collection(lc)

    # E_F line
    ax.axhline(0.0, color=S.C_SEC, lw=0.8, ls=(0, (5, 3)), zorder=2)

    # high-symmetry ticks
    ax.set_xticks(d["ticks"])
    ax.set_xticklabels(SP.KPATH_LABELS)
    for xt in d["ticks"][1:-1]:
        ax.axvline(xt, color=S.C_GRID, lw=0.6, zorder=1.5)
    ax.set_xlim(kfine[0], kfine[-1])
    ax.set_ylim(*EWIN)
    ax.set_title(title, pad=4)
    ax.tick_params(top=False)

axes[0].set_ylabel(r"$E - E_F$  (eV)")
axes[2].text(0.97, 0.0, r"$E_F$", transform=axes[2].get_yaxis_transform(),
             ha="right", va="bottom", fontsize=7.5, color=S.C_SEC)
for j, ax in enumerate(axes):
    S.panel_label(ax, "abc"[j], x=0.0, y=1.03)

# a compact Na-weight linewidth key (neutral text, blue swatch = data).
# White bbox keeps it legible over the dense 9.9 A bands (neutral, not coloured).
key_ax = axes[2]
key_ax.add_patch(plt.Rectangle((0.58, 0.665), 0.40, 0.325,
                 transform=key_ax.transAxes, facecolor="white",
                 edgecolor=S.C_GRID, lw=0.5, zorder=5))
for frac, yy in ((1.0, 0.895), (0.5, 0.805), (0.15, 0.72)):
    key_ax.plot([0.62, 0.74], [yy, yy], transform=key_ax.transAxes,
                color=S.C_BLUE, lw=0.3 + 6.0 * frac, solid_capstyle="round",
                zorder=6)
    key_ax.text(0.77, yy, f"{frac:.2f}", transform=key_ax.transAxes,
                fontsize=6.3, va="center", color=S.C_INK, zorder=6)
key_ax.text(0.60, 0.955, "Na weight", transform=key_ax.transAxes,
            fontsize=6.6, va="center", color=S.C_INK, zorder=6)

fig.subplots_adjust(left=0.075, right=0.99, top=0.90, bottom=0.10, wspace=0.06)
S.save(fig, "fig7")
print("fig7 written; E_F/mag:",
      {j: (round(data[j]["ef"], 2), round(data[j]["mag"], 2)) for j, _ in JOBS})

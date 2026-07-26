#!/usr/bin/env python
"""
Figure S (valence-thickness channel) - figS_zscan.
How the CoO2-sheet thickness -- parametrised by the oxygen height z_O above the
Co plane, a structural proxy for the Co valence -- tunes the alkali well.

Data: runpod/results_zscan/*/pw.out (set K): vacuum sqrt3xsqrt3 Na_{1/3}CoO2
E(delta) scans at fixed c = 9.9 A for z_O = 0.90 and 1.02 A, re-parsed here from
the raw total energies exactly as _style.load_scans does.  The z_O = 0.96 A
baseline is the existing set-G vacuum scan (theory/results_extra_analysis/
summary_extra.json).  Every scan is the vacuum (unterminated) Na potential, so
alpha < 0 throughout -- the honest quantity these SCF jobs support is the SLOPE
d(alpha)/d(z_O), not an absolute stiffness.

 (a) the three vacuum E(delta) curves, coloured by z_O: a taller CoO2 sheet
     DEEPENS the well.
 (b) alpha (quadratic coefficient of each E(delta) fit) vs z_O, with the linear
     fit; the slope d(alpha)/d(z_O) is the microscopic basis of the XAS/NQR
     valence-tuning prediction.

House style: pastel fill for the alpha<0 (softening) domain in (b), flat only.
"""
import re
import json
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

import _style as S

S.use_house_style()

REPO = Path(__file__).resolve().parents[2]
ZDIR = REPO / "runpod" / "results_zscan"
RY_EV = S.RY_EV
_E_RE = re.compile(r"^!\s+total energy\s+=\s+(-?[\d.]+)\s+Ry", re.M)

DELTAS = [0.00, 0.15, 0.30, 0.50, 0.75]


def scan(zo):
    """Vacuum E(delta) (meV, rel. delta=0) for a set-K z_O from raw pw.out."""
    E = []
    for d in DELTAS:
        txt = (ZDIR / f"zscan_Na_s3_c9.9_zO{zo:.2f}_d{d:.2f}" / "pw.out").read_text()
        assert "JOB DONE" in txt
        E.append(float(_E_RE.findall(txt)[-1]) * RY_EV)
    E = (np.array(E) - E[0]) * 1e3
    return np.array(DELTAS), E


# set-G vacuum baseline (z_O = 0.96) from the existing parse
G = json.loads((REPO / "theory" / "results_extra_analysis"
                / "summary_extra.json").read_text())["G"]
curves = {0.90: scan(0.90),
          0.96: (np.array(G["d"]), np.array(G["E_vac_meV"])),
          1.02: scan(1.02)}


def alpha_of(dd, EE):
    """V = a d^2 + b d^4 (eV, d in A); return a (= alpha)."""
    a, _b = np.linalg.lstsq(np.column_stack([dd**2, dd**4]), EE / 1e3,
                            rcond=None)[0]
    return a


zos = np.array(sorted(curves))
alphas = np.array([alpha_of(*curves[z]) for z in zos])
slope, intc = np.polyfit(zos, alphas, 1)

# mono-hue ramp (thin sheet -> thick sheet): light -> dark blue
RAMP = {0.90: S.RAMP_BLUE[0], 0.96: S.RAMP_BLUE[1], 1.02: S.RAMP_BLUE[3]}

fig, (axA, axB) = plt.subplots(1, 2, figsize=(6.7, 3.15))

# ================================================================= (a)
for z in zos:
    dd, EE = curves[z]
    xx = np.linspace(0, 0.75, 200)
    a, b = np.linalg.lstsq(np.column_stack([dd**2, dd**4]), EE / 1e3,
                           rcond=None)[0]
    axA.plot(xx, (a * xx**2 + b * xx**4) * 1e3, color=RAMP[z], lw=1.4, zorder=3)
    axA.plot(dd, EE, "o", ms=5, color=RAMP[z], mfc="white", mew=1.4, zorder=4,
             label=fr"$z_{{\mathrm{{O}}}}={z:.2f}$ Å")
axA.set_xlim(0, 0.78)
axA.set_xlabel(r"Na off-centre displacement  $\delta$  (Å)")
axA.set_ylabel(r"$E(\delta)-E(0)$  (meV)")
axA.legend(loc="lower left", handletextpad=0.4, labelspacing=0.3,
           borderaxespad=0.7, fontsize=7.2, title="sheet thickness")
axA.get_legend().get_title().set_fontsize(7.2)
S.thin_spines(axA)
S.panel_label(axA, "a")

# ================================================================= (b)
YLO, YHI = -0.40, -0.29
axB.axhspan(YLO, 0.0, color=S.FILL_AMBER, lw=0, zorder=0)  # softening domain
zg = np.linspace(0.87, 1.05, 50)
axB.plot(zg, slope * zg + intc, color=S.C_SEC, lw=1.2, zorder=2)
for z, a in zip(zos, alphas):
    axB.plot(z, a, "o", ms=7, color=RAMP[z], mfc=RAMP[z], mew=0, zorder=4)
# slope quoted as a minimal curve label
axB.text(0.995, -0.312, fr"$d\alpha/dz_{{\mathrm{{O}}}}={slope:.2f}$"
         "\n" r"eV Å$^{-2}$ / Å", ha="left", va="top", color=S.C_SEC,
         fontsize=7.6, zorder=5)
axB.text(0.885, -0.312, "thicker sheet\n$\\rightarrow$ softer well",
         ha="left", va="top", color=S.INK_AMBER, fontsize=7.4, zorder=5)
axB.set_xlim(0.87, 1.05)
axB.set_ylim(YLO, YHI)
axB.set_xlabel(r"oxygen height  $z_{\mathrm{O}}$  (Å)   [Co-valence proxy]")
axB.set_ylabel(r"Landau stiffness  $\alpha$  (eV/Å$^2$)")
S.thin_spines(axB)
S.panel_label(axB, "b")

fig.subplots_adjust(left=0.085, right=0.985, top=0.95, bottom=0.145,
                    wspace=0.28)
S.save(fig, "figS_zscan")
print(f"figS_zscan written  (alpha: {dict(zip(zos, alphas.round(4)))}, "
      f"slope={slope:.3f})")

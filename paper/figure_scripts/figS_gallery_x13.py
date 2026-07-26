#!/usr/bin/env python
"""
Figure S (gallery-band onset at x = 1/3) - figS_gallery_x13.

Referee-response figure: at the ordered sqrt3 x sqrt3 Na_{1/3}CoO2 cell, where
exactly does the interlayer two-dimensional electron gas (2DEG) form, how much
charge does the sodium hold, and how does this connect to the double well?

All quantities are re-parsed from the existing raw outputs; no new calculations.

 (a) Na-projected gallery band along Gamma-M-K-Gamma (spin up) at c = 5.5
     (band 56) and c = 9.9 A (band 50), referenced to each Fermi level.  The
     band drops ~5 eV as the gallery opens and its bottom crosses E_F only at
     9.9 A; in-plane widths 2.92 eV (5.5 A) and 1.52 eV (9.9 A) quantify the
     k-parallel broadening that the surface-model picture requires.
 (b) the same band at Gamma vs c for delta = 0 and 0.75 A: the Gamma pocket
     opens between 8.4 and 9.9 A (occupation 0.92 at 9.9 A, delta = 0).  There
     is no finite-delta threshold at 9.9 A: the pocket exists at delta = 0 and
     is PARTIALLY EMPTIED as Na moves off centre (occupancy 0.92 -> 0.74), the
     negative feedback that bounds the double-well minimum.
 (c) Bader charge donated by Na (9 - q_Na, semicore-inclusive PAW valence) vs
     delta for the four spacings: ~0.85-0.94 e transferred to the CoO2 sheet
     while the gallery is closed, dropping to ~0.65 e once the pocket is
     occupied -- the missing ~0.25 e is the 2DEG itself, residing in the Na
     (gallery) basin.
 (d) E(delta) wells at the four spacings (meV, rel. delta = 0): the bifurcation
     deepens hand in hand with the pocket opening.

Data: runpod/results_v4/jobs/Na_s3_c*_d*/pw.out + ACF.dat (set B), and
runpod/results_bands/bands_Na_s3_c{5.5,9.9}/ (set I1 fatbands + projwfc).
Assembled table: runpod/results_v4/gallery_band_vs_delta_s3.json.
"""
import json
import re
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

import _style as S

S.use_house_style()

REPO = Path(__file__).resolve().parents[2]
V4 = REPO / "runpod" / "results_v4"
BANDS = REPO / "runpod" / "results_bands"
RY_EV = S.RY_EV

CS = ["5.5", "6.9", "8.4", "9.9"]
C_COL = {"5.5": S.C_SEC, "6.9": S.C_AQUA, "8.4": S.C_BLUE, "9.9": S.C_RED}
K_RE = re.compile(r"k\s*=\s*(-?[\d.]+)\s+(-?[\d.]+)\s+(-?[\d.]+).*bands \(ev\)")
NUM_RE = re.compile(r"-?\d+\.\d+")


def gamma_block(path):
    """Last spin-up Gamma eigenvalue+occupation block and E_F from pw.out."""
    lines = Path(path).read_text().splitlines()
    ef = [float(x) for x in
          re.findall(r"Fermi energy is\s*(-?[\d.]+)", "\n".join(lines))]
    ef = ef[-1]
    spin, res = "up", {}
    i = 0
    while i < len(lines):
        ln = lines[i]
        if "SPIN UP" in ln:
            spin = "up"
        elif "SPIN DOWN" in ln:
            spin = "down"
        m = K_RE.search(ln)
        if m and all(abs(float(g)) < 1e-4 for g in m.groups()):
            eigs, occ, j = [], [], i + 1
            while j < len(lines):
                lj = lines[j]
                if K_RE.search(lj) or "SPIN" in lj:
                    break
                if "occupation numbers" in lj:
                    k2 = j + 1
                    while k2 < len(lines):
                        lk = lines[k2].strip()
                        if lk == "" and occ:
                            break
                        occ += NUM_RE.findall(lk)
                        k2 += 1
                    j = k2
                    break
                eigs += NUM_RE.findall(lj)
                j += 1
            if eigs and occ:
                res[spin] = ([float(x) for x in eigs],
                             [float(x) for x in occ])
            i = j
        i += 1
    return res["up"], ef


def fatband(cdir, band_index, na_states=range(55, 60)):
    """(k, E, Na weight per k) for one band from a set-I1 bands directory."""
    blocks, cur = [], []
    for ln in (cdir / "pw.bands.up.dat.gnu").read_text().splitlines():
        if ln.strip() == "":
            if cur:
                blocks.append(np.array(cur))
                cur = []
        else:
            cur.append([float(x) for x in ln.split()[:2]])
    if cur:
        blocks.append(np.array(cur))
    blk = blocks[band_index - 1]
    hdr = re.compile(r"^\s*(\d+)\s+(\d+)\s+([A-Za-z]+)\s+(\S+)\s+")
    kbr = re.compile(r"^\s*(\d+)\s+(\d+)\s+([-\d.]+)\s*$")
    w, cur = {}, None
    for ln in (cdir / "pw.proj.projwfc_up").read_text().splitlines():
        m = hdr.match(ln)
        if m and not kbr.match(ln):
            st = int(m.group(1))
            cur = st if st in set(na_states) else None
            continue
        m = kbr.match(ln)
        if m and cur is not None:
            k, b, wt = int(m.group(1)), int(m.group(2)), float(m.group(3))
            if b == band_index:
                w[k] = w.get(k, 0.0) + wt
    wt = np.array([w.get(k, 0.0) for k in range(1, len(blk) + 1)])
    ef = float(re.findall(r"Fermi energy is\s*(-?[\d.]+)",
                          (cdir / "pw.out").read_text())[-1])
    return blk[:, 0], blk[:, 1] - ef, wt


tab = json.loads((V4 / "gallery_band_vs_delta_s3.json").read_text())

fig, axs = plt.subplots(1, 4, figsize=(13.4, 3.3))
a, b, cax, d = axs

# (a) Na-projected dispersion at 5.5 vs 9.9 ----------------------------------
curves = {}
for cval, band in (("5.5", 56), ("9.9", 50)):
    k, e, w = fatband(BANDS / f"bands_Na_s3_c{cval}", band)
    x = (k - k.min()) / (k.max() - k.min())
    curves[cval] = (x, e, w)
    a.plot(x, e, color=C_COL[cval], lw=1.3)
    a.scatter(x, e, s=26 * w + 2, color=C_COL[cval], alpha=0.45,
              edgecolors="none")
a.axhline(0, color=S.C_MUT, lw=0.8, ls=(0, (4, 3)))
x55, e55, _ = curves["5.5"]
x99, e99, _ = curves["9.9"]
a.text(0.635, 9.15, r"5.5 $\AA$", color=S.C_SEC, fontsize=8.2,
       ha="center")
a.text(0.335, 2.35, r"9.9 $\AA$", color=S.C_RED, fontsize=8.2,
       ha="center")
a.set_xticks([0, 1 / 3, 2 / 3, 1])
a.set_xticklabels([r"$\Gamma$", "M", "K", r"$\Gamma$"])
a.set_ylabel(r"$\varepsilon_{\rm gal}(k)-E_F$ (eV)")
a.set_ylim(-1.2, 10.4)
S.thin_spines(a, hide=())

# (b) Gamma energy and occupation vs c ---------------------------------------
cc = np.array([float(x) for x in CS])
for delta, mk, lab in ((0.00, "o", r"$\delta=0$"), (0.75, "s", r"$\delta=0.75$ $\AA$")):
    ee = [tab[cv][int(delta / 0.25)]["gal_e"] for cv in CS]
    b.plot(cc, ee, marker=mk, ms=4.5, color=S.C_BLUE if delta == 0 else S.C_YELL,
           lw=1.2, label=lab)
b.axhline(0, color=S.C_MUT, lw=0.8, ls=(0, (4, 3)))
b.axvspan(8.4, 9.9, color=S.C_RED, alpha=0.07)
b.text(9.13, 1.32, "2DEG", ha="center", fontsize=8.2, color=S.C_RED)
b.set_xlabel(r"gallery spacing $c$ ($\AA$)")
b.set_ylabel(r"$\varepsilon_{\rm gal}(\Gamma)-E_F$ (eV)")
b.legend(frameon=False, fontsize=7.6, loc="lower right")
b.set_xticks(cc)
S.thin_spines(b, hide=())

# (c) Bader donated charge vs delta ------------------------------------------
endlab = {"5.5": -0.020, "6.9": 0.032, "8.4": -0.036, "9.9": 0.0}
for cv in CS:
    ds = [r["delta"] for r in tab[cv]]
    dn = [r["donated"] for r in tab[cv]]
    cax.plot(ds, dn, marker="o", ms=4, lw=1.2, color=C_COL[cv])
    cax.text(0.775, dn[-1] + endlab[cv], rf"{cv} $\AA$", color=C_COL[cv],
             fontsize=7.8, va="center", clip_on=False)
cax.set_xlabel(r"Na displacement $\delta$ ($\AA$)")
cax.set_ylabel(r"charge donated to CoO$_2$ ($e$)")
cax.set_xlim(-0.03, 0.78)
cax.set_ylim(0.55, 1.0)
S.thin_spines(cax, hide=())

# (d) E(delta) wells ----------------------------------------------------------
for cv in CS[1:]:
    ds = np.array([r["delta"] for r in tab[cv]])
    E = np.array([r["E_ry"] for r in tab[cv]]) * RY_EV
    d.plot(ds, (E - E[0]) * 1e3, marker="o", ms=4, lw=1.2, color=C_COL[cv],
           label=rf"{cv} $\AA$")
d.axhline(0, color=S.C_MUT, lw=0.8, ls=(0, (4, 3)))
d.text(0.04, 0.945, "5.5 $\AA$: $+2.56$ eV at $\delta=0.75$ $\AA$",
       transform=d.transAxes, fontsize=7.6, color=S.C_SEC, va="top")
d.set_xlabel(r"Na displacement $\delta$ ($\AA$)")
d.set_ylabel(r"$E(\delta)-E(0)$ (meV)")
d.set_ylim(-260, 120)
d.legend(frameon=False, fontsize=7.4, loc="upper right", handlelength=1.3)
S.thin_spines(d, hide=())

for ax, lt in zip(axs, "abcd"):
    S.panel_label(ax, lt)

fig.tight_layout(w_pad=1.6)
S.save(fig, "figS_gallery_x13")
print("wrote figures/figS_gallery_x13.{pdf,png}")

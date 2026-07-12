#!/usr/bin/env python
"""
Figure 6 - The rigid-cage water demonstration (Sec. IV).

(a) E(delta) for the off-centre Na displacement at c = 9.9 A in the
    sqrt3 x sqrt3 (x=1/3) cell, WITH a rigid 4-H2O solvation cage (hydrate)
    vs the SAME cell with the water deleted (vacuum). Raw QE total energies
    (runpod/results_extra/jobs_extra/Na_s3{hyd,vac}_c9.9_d*). The vacuum curve
    runs away monotonically (Na chemisorbs onto a CoO2 sheet: unstable, no
    bound well); the hydrate curve is a BOUND double well with minima near
    +/-0.30 A and a ~0.20 eV central barrier -> a real, confined c-axis mode.

(b) The hydrate unit cell rendered FROM THE ACTUAL pw.in COORDINATES
    (Na_s3hyd_c9.9_d0.30/pw.in), side view along c, no generative imagery.

House style from _style.py. Aesthetic-sobriety rule: colour only on the two
data curves in (a) and on atomic species in (b); every decorative element
(guides, shading, annotations) is neutral grey/black.
"""
import re
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch

import _style as S

S.use_house_style()

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
JOBS = REPO / "runpod" / "results_extra" / "jobs_extra"
BOHR = 0.52917721067
RY_EV = S.RY_EV
E_RE = re.compile(r"^!\s+total energy\s+=\s+(-?[\d.]+)\s+Ry", re.M)


def total_E(job):
    txt = (JOBS / job / "pw.out").read_text()
    return float(E_RE.findall(txt)[-1]) * RY_EV


# ---------------------------------------------------------------- data (a)
dd = np.array([0.0, 0.15, 0.30, 0.50, 0.75])
Eh = np.array([total_E(f"Na_s3hyd_c9.9_d{d:.2f}") for d in dd])
Ev = np.array([total_E(f"Na_s3vac_c9.9_d{d:.2f}") for d in dd])
Eh = (Eh - Eh[0]) * 1e3       # meV rel. centre
Ev = (Ev - Ev[0]) * 1e3

# ---------------------------------------------------------------- data (b)
def parse_pwin(job):
    """Return {species: Nx3 cartesian (A)} and (a, c) from an ibrav=4 pw.in."""
    txt = (JOBS / job / "pw.in").read_text().splitlines()
    a = c = None
    for ln in txt:
        m = re.search(r"celldm\(1\)\s*=\s*([\d.]+)", ln)
        if m:
            a = float(m.group(1)) * BOHR
        m = re.search(r"celldm\(3\)\s*=\s*([\d.]+)", ln)
        if m:
            c3 = float(m.group(1))
    c = a * c3
    a1 = np.array([a, 0, 0])
    a2 = np.array([-a / 2, a * np.sqrt(3) / 2, 0])
    a3 = np.array([0, 0, c])
    M = np.vstack([a1, a2, a3])
    atoms = {}
    read = False
    for ln in txt:
        if ln.strip().startswith("ATOMIC_POSITIONS"):
            read = True
            continue
        if read:
            p = ln.split()
            if len(p) < 4 or p[0] in ("K_POINTS",):
                if p and p[0] == "K_POINTS":
                    break
                continue
            sp = p[0]
            frac = np.array([float(p[1]), float(p[2]), float(p[3])])
            xyz = frac @ M
            atoms.setdefault(sp, []).append(xyz)
    return {k: np.array(v) for k, v in atoms.items()}, a, c


atoms, aL, cL = parse_pwin("Na_s3hyd_c9.9_d0.30")

# species -> (facecolor, edgecolor, radius pt, label). Muted CPK; atoms=data.
SPEC = {
    "Co": ("#3b4a6b", "#243049", 10.5, "Co"),
    "O":  ("#d98a7a", "#9c5344", 7.0, "O"),
    "Na": ("#c9a63a", "#8f7212", 11.5, "Na"),
    "H":  ("#d9d7cf", "#8a8880", 4.2, "H"),
}

# ======================================================================
fig, (axA, axB) = plt.subplots(
    1, 2, figsize=(7.0, 3.3), gridspec_kw={"width_ratios": [1.28, 1.0]})

# ---------------------------------------------------------------- (a)
dsym = np.concatenate([-dd[:0:-1], dd])       # drop the duplicated 0
Ehs = np.concatenate([Eh[:0:-1], Eh])
Evs = np.concatenate([Ev[:0:-1], Ev])

# smooth guides (monotone PCHIP through the symmetric points), neutral weight
from scipy.interpolate import PchipInterpolator
xf = np.linspace(-0.75, 0.75, 300)
axA.plot(xf, PchipInterpolator(dsym, Ehs)(xf), "-", color=S.C_BLUE, lw=1.4,
         zorder=3)
axA.plot(xf, PchipInterpolator(dsym, Evs)(xf), "-", color=S.C_SEC, lw=1.4,
         zorder=3)
axA.plot(dsym, Ehs, "o", ms=5.0, color=S.C_BLUE, mec="white", mew=0.6,
         zorder=4, label="hydrate (rigid 4-H$_2$O cage)")
axA.plot(dsym, Evs, "s", ms=4.6, color=S.C_SEC, mec="white", mew=0.6,
         zorder=4, label="vacuum (water deleted)")

axA.axhline(0, color=S.C_MUT, lw=0.6, zorder=1)
# mark the bound minimum and the runaway, both in neutral ink
axA.annotate("bound well\n$\\delta_0\\approx0.30$ Å", (0.30, Eh[2]),
             xytext=(0.30, -95), ha="center", fontsize=7.0, color=S.C_INK,
             arrowprops=dict(arrowstyle="-", color=S.C_MUT, lw=0.7))
axA.annotate("runs away\n(Na adsorbs)", (0.72, Ev[-1] + 4),
             xytext=(0.40, -232), ha="center", fontsize=7.0, color=S.C_INK,
             arrowprops=dict(arrowstyle="->", color=S.C_MUT, lw=0.7))
axA.set_xlabel(r"Na off-centre displacement  $\delta$  (Å)")
axA.set_ylabel(r"$E(\delta)-E(0)$  (meV)")
axA.set_xlim(-0.80, 0.80)
axA.set_ylim(-250, 55)
axA.legend(loc="upper center", fontsize=7.0, handletextpad=0.4,
           borderaxespad=0.3)
S.thin_spines(axA)
S.panel_label(axA, "a", x=-0.16)

# ---------------------------------------------------------------- (b)
# side view: horizontal = x (A), vertical = z (A, the c-axis stacking).
# draw two in-plane periodic images (shift by +a1) so the layers read as sheets
a1 = np.array([aL, 0, 0])
def scat(P, sp):
    fc, ec, r, _ = SPEC[sp]
    for shift in (0.0, 1.0):
        Q = P + shift * a1
        axB.scatter(Q[:, 0], Q[:, 2], s=r**2, facecolor=fc, edgecolor=ec,
                    linewidth=0.5, zorder={"H": 2, "O": 3, "Co": 4,
                                           "Na": 5}[sp])

# water O-H bonds (within the primary + shifted image)
Ow = atoms["O"][atoms["O"][:, 2] > 0.20 * cL]   # water O (above the lower slab)
Hs = atoms.get("H", np.empty((0, 3)))
for shift in (0.0, 1.0):
    for o in Ow:
        for h in Hs:
            for hs in (0.0, 1.0):
                oo, hh = o + shift * a1, h + hs * a1
                if np.linalg.norm(oo - hh) < 1.1:
                    axB.plot([oo[0], hh[0]], [oo[2], hh[2]], "-",
                             color="#b9b7ae", lw=0.8, zorder=1)
for sp in ("Co", "O", "H", "Na"):
    if sp in atoms:
        scat(atoms[sp], sp)

# CoO2-sheet guide lines (z of the two Co planes = 0 and c) and the mid-plane
for zc, lab in ((0.0, None), (cL, None)):
    axB.axhline(zc, color=S.C_GRID, lw=0.7, zorder=0)
axB.axhline(cL / 2, color=S.C_MUT, lw=0.6, ls=(0, (4, 3)), zorder=0)
axB.text(axB_x := (-0.5 * aL), cL / 2, " Na plane\n ($c/2$)", fontsize=6.3,
         color=S.C_SEC, va="center", ha="left")
axB.text(-0.5 * aL, 0.02 * cL, r" CoO$_2$", fontsize=6.3, color=S.C_SEC,
         va="bottom")

# tidy legend of species (small dots), neutral text
handles = [plt.Line2D([], [], marker="o", ls="", mfc=SPEC[s][0],
                      mec=SPEC[s][1], mew=0.5,
                      ms=np.sqrt(SPEC[s][2] ** 2) / 2.6, label=SPEC[s][3])
           for s in ("Co", "O", "Na", "H")]
axB.legend(handles=handles, loc="upper right", fontsize=6.6, ncol=4,
           handletextpad=0.15, columnspacing=0.7, borderaxespad=0.2,
           bbox_to_anchor=(1.02, 1.06))

axB.set_xlabel(r"$x$ (Å)")
axB.set_ylabel(r"$z$ along $c$ (Å)")
axB.set_xlim(-0.7 * aL, 1.7 * aL)
axB.set_ylim(-0.9, cL + 0.9)
axB.set_aspect("equal", adjustable="box")
S.thin_spines(axB)
S.panel_label(axB, "b", x=-0.12)

fig.subplots_adjust(left=0.085, right=0.985, top=0.95, bottom=0.145,
                    wspace=0.32)
S.save(fig, "fig6")
print("fig6 written")

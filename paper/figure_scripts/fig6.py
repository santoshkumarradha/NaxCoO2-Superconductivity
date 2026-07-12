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


MOB = REPO / "runpod" / "results_mobile"
BND = REPO / "runpod" / "results_bands"


def totE(base, job):
    return float(E_RE.findall((base / job / "pw.out").read_text())[-1]) * RY_EV


def total_E(job):
    return totE(JOBS, job)


# ---------------------------------------------------------------- data (a)
# Three treatments of the SAME sqrt3 hydrate cell at c=9.9 A, each E(delta)
# referenced to its OWN delta=0 (the well SHAPE, not the ~1.9 eV relaxation
# offset between a frozen and a relaxed cage):
#   vacuum    (water deleted)                 -> monotone runaway, open markers
#   rigid     (gas-phase cage frozen)         -> bounded well, min -199 meV at 0.30
#   adiabatic (mobile: the 4 H2O BFGS-relaxed -> bounded AND DEEPER, min -380 meV;
#              at each pinned Na, 100 steps)      the shell FOLLOWS Na, not flattens
dd = np.array([0.0, 0.15, 0.30, 0.50, 0.75])
Eh = np.array([total_E(f"Na_s3hyd_c9.9_d{d:.2f}") for d in dd])
Ev = np.array([total_E(f"Na_s3vac_c9.9_d{d:.2f}") for d in dd])
Eh = (Eh - Eh[0]) * 1e3       # meV rel. centre (rigid cage)
Ev = (Ev - Ev[0]) * 1e3       # meV rel. centre (vacuum)

dm = np.array([0.0, 0.15, 0.30, 0.50, 0.75, 1.00])
Em_abs = np.array([totE(MOB, f"mobile_Na_s3hyd_c9.9_d{d:.2f}") for d in dm])
Em = (Em_abs - Em_abs[0]) * 1e3   # meV rel. centre (adiabatic, mobile water)

# 60-step-capped adiabatic relaxation at delta=0.30: a DIFFERENT, incompletely
# reorganized H-bond configuration that stops 394 meV above the 100-step
# minimum -> an explicit sample of the multi-minimum spread of the 4-water
# network at fixed Na (referenced to the same relaxed delta=0 baseline).
E60_30 = (totE(BND, "relax_Na_s3hyd_c9.9_d0.30") - Em_abs[0]) * 1e3

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
dmsym = np.concatenate([-dm[:0:-1], dm])
Ems = np.concatenate([Em[:0:-1], Em])

from scipy.interpolate import PchipInterpolator
xf = np.linspace(-0.75, 0.75, 400)
rig_i = PchipInterpolator(dsym, Ehs)(xf)
adi_i = PchipInterpolator(dmsym, Ems)(xf)

# configurational-range band: between the frozen-cage well (upper) and the deep
# adiabatic well (lower) -- the range of E(delta) surfaces the moving cage can
# present to Na.  Pastel amber = the softening/prediction domain (house style).
axA.fill_between(xf, rig_i, adi_i, color=S.FILL_AMBER, lw=0, zorder=1)
axA.text(-0.52, -250, "cage\nconfigurational\nrange $\\gtrsim 0.4$ eV",
         ha="center", va="center", fontsize=6.6, color=S.INK_AMBER, zorder=2)
axA.axhline(0, color=S.C_MUT, lw=0.6, zorder=1)

# smooth guides (monotone PCHIP through the mirror-symmetric points)
xf2 = np.linspace(-1.0, 1.0, 500)
axA.plot(xf2, PchipInterpolator(dmsym, Ems)(xf2), "-", color=S.C_RED, lw=1.5,
         zorder=3)
axA.plot(xf, rig_i, "-", color=S.C_BLUE, lw=1.4, zorder=3)
axA.plot(xf, PchipInterpolator(dsym, Evs)(xf), color=S.C_SEC, lw=1.1,
         ls=(0, (4, 2)), zorder=3)

axA.plot(dmsym, Ems, "D", ms=4.3, color=S.C_RED, mec="white", mew=0.5,
         zorder=5, label="adiabatic (mobile 4-H$_2$O, 100 steps)")
axA.plot(dsym, Ehs, "o", ms=4.8, color=S.C_BLUE, mec="white", mew=0.5,
         zorder=5, label="rigid cage (frozen 4-H$_2$O)")
axA.plot(dsym, Evs, "s", ms=4.2, mfc="white", mec=S.C_SEC, mew=1.0,
         zorder=5, label="vacuum (water deleted)")

# 60-step outlier: an explicit multi-minimum sample, and the vertical spread it
# opens above the adiabatic minimum (evidence for the >=0.4 eV cage landscape)
axA.annotate("", xy=(0.30, Em[2] + 6), xytext=(0.30, E60_30 - 6),
             arrowprops=dict(arrowstyle="<->", color=S.C_INK, lw=0.8),
             zorder=4)
axA.text(0.35, 0.5 * (Em[2] + E60_30), "$394$ meV", rotation=90, va="center",
         ha="left", fontsize=6.4, color=S.C_INK, zorder=6)
axA.plot([0.30, -0.30], [E60_30, E60_30], "*", ms=9.0, mfc="white",
         mec=S.C_INK, mew=0.9, zorder=6)
axA.annotate("step-capped\nlocal minimum", (0.30, E60_30), xytext=(0.66, -55),
             ha="center", fontsize=6.4, color=S.C_INK,
             arrowprops=dict(arrowstyle="-", color=S.C_MUT, lw=0.6), zorder=6)
axA.annotate("runs away", (0.72, Ev[-1] + 4), xytext=(0.55, -250),
             ha="center", fontsize=6.6, color=S.C_SEC,
             arrowprops=dict(arrowstyle="->", color=S.C_MUT, lw=0.6))

axA.set_xlabel(r"Na off-centre displacement  $\delta$  (Å)")
axA.set_ylabel(r"$E(\delta)-E(0)$  (meV)")
axA.set_xlim(-1.05, 1.05)
axA.set_ylim(-410, 60)
axA.legend(loc="lower center", fontsize=6.6, handletextpad=0.4,
           borderaxespad=0.4, labelspacing=0.3)
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

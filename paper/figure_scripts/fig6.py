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
import matplotlib.colors as mcolors
from matplotlib.patches import Circle
from matplotlib.lines import Line2D

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

# Species -> flat CPK-ish fill + ionic/covalent-proportioned radius (A) + key
# marker size.  Na large, O medium, Co small, H smallest (per house palette).
SPEC = {
    "Na": dict(color="#d9b23a", r=0.62, ms=9.0, label="Na"),
    "O":  dict(color="#d76b5a", r=0.44, ms=7.5, label="O"),
    "Co": dict(color="#3b5b8c", r=0.38, ms=7.0, label="Co"),
    "H":  dict(color="#d8d6cc", r=0.24, ms=5.5, label="H"),
}
PAPER = "#fcfcfb"


def _blend(color, other, t):
    """Flat de-saturation toward the paper colour (no gradient): mix a fraction
    t of `other` into `color`.  Used to push far-depth atoms slightly back."""
    c = np.array(mcolors.to_rgb(color))
    o = np.array(mcolors.to_rgb(other))
    return tuple((1 - t) * c + t * o)


# ======================================================================
fig, (axA, axB) = plt.subplots(
    1, 2, figsize=(7.0, 2.95), gridspec_kw={"width_ratios": [1.95, 1.0]})

# ================================================================= (a)
dsym = np.concatenate([-dd[:0:-1], dd])       # drop the duplicated 0
Ehs = np.concatenate([Eh[:0:-1], Eh])
Evs = np.concatenate([Ev[:0:-1], Ev])
dmsym = np.concatenate([-dm[:0:-1], dm])
Ems = np.concatenate([Em[:0:-1], Em])

from scipy.interpolate import PchipInterpolator
xf = np.linspace(-0.75, 0.75, 400)
xf2 = np.linspace(-1.0, 1.0, 500)
rig_i = PchipInterpolator(dsym, Ehs)(xf)
adi_i = PchipInterpolator(dmsym, Ems)(xf)
adi_full = PchipInterpolator(dmsym, Ems)(xf2)
vac_i = PchipInterpolator(dsym, Evs)(xf)

# amber configurational-range band = spread of E(delta) surfaces the moving cage
# can present between the frozen-cage well (upper) and the deep adiabatic well.
axA.fill_between(xf, rig_i, adi_i, color=S.FILL_AMBER, lw=0, zorder=1)
axA.text(-0.46, -235, "cage range", ha="center", va="center", fontsize=7.4,
         color=S.INK_AMBER, zorder=2)
axA.axhline(0, color=S.C_MUT, lw=0.6, zorder=1)

# smooth guides (monotone PCHIP through the mirror-symmetric points)
axA.plot(xf2, adi_full, "-", color=S.C_RED, lw=1.6, zorder=3)
axA.plot(xf, rig_i, "-", color=S.C_BLUE, lw=1.5, zorder=3)
axA.plot(xf, vac_i, ls=(0, (4, 2)), color=S.C_SEC, lw=1.1, zorder=3)

# data markers (legend labels short; caption carries the reading)
axA.plot(dsym, Evs, "s", ms=4.0, mfc="white", mec=S.C_SEC, mew=1.0,
         zorder=5, label="vacuum")
axA.plot(dsym, Ehs, "o", ms=4.6, color=S.C_BLUE, mec="white", mew=0.6,
         zorder=5, label="rigid cage")
axA.plot(dmsym, Ems, "D", ms=4.2, color=S.C_RED, mec="white", mew=0.6,
         zorder=5, label="adiabatic water")

# the 60-step-capped local minimum: an explicit high-lying multi-minimum sample
# (its meaning is spelled out in the caption -- no in-plot annotation here)
axA.plot([0.30, -0.30], [E60_30, E60_30], "*", ms=8.5, mfc="white",
         mec=S.C_INK, mew=0.8, zorder=6)

axA.set_xlabel(r"Na off-centre displacement  $\delta$  (Å)")
axA.set_ylabel(r"$E(\delta)-E(0)$  (meV)")
axA.set_xlim(-1.05, 1.05)
axA.set_ylim(-405, 80)
axA.legend(loc="upper center", ncol=3, fontsize=7.3, columnspacing=1.4,
           handletextpad=0.4, borderaxespad=0.35)
S.thin_spines(axA)
S.panel_label(axA, "a", x=-0.13)

# ================================================================= (b)
# Flat vector crystal drawing (side view, projection along the in-plane y):
# horizontal = x, vertical = z (the c-axis stack).  The full periodic crystal
# is generated from the pw.in basis with hexagonal lattice translations
# n1*a1 + n2*a2 + n3*a3, then clipped to a window ~2.5 cells wide and one
# in-plane depth period deep, so the CoO2 sheets appear as complete O-Co-O
# sandwiches (both apical O rows), the gallery as the Na + 4-H2O layer, and
# every drawn bond is a real 3D near-neighbour distance.  Frame OFF; the
# projected unit cell is a thin dashed rectangle; single c-arrow scale cue.
a1 = np.array([aL, 0.0, 0.0])
a2 = np.array([-aL / 2, aL * np.sqrt(3) / 2, 0.0])
a3 = np.array([0.0, 0.0, cL])

XW = (-1.3, 11.3)                # drawing window in x (about 2.5 cells)
ZW = (-1.35, cL + 1.35)          # z window: sheets complete with both O rows
YW = (-0.10, aL * np.sqrt(3) / 2 - 0.10)   # one depth period (no duplicates)

kept = []                        # (pos3, species)
for sp in ("Co", "O", "Na", "H"):
    for P in atoms.get(sp, []):
        for n1 in range(-2, 4):
            for n2 in range(-2, 3):
                for n3 in (-1, 0, 1):
                    q = P + n1 * a1 + n2 * a2 + n3 * a3
                    if (XW[0] <= q[0] <= XW[1] and ZW[0] <= q[2] <= ZW[1]
                            and YW[0] <= q[1] < YW[1]):
                        kept.append((q, sp))
pos = np.array([q for q, _ in kept])
spc = [sp for _, sp in kept]
ymin, ymax = pos[:, 1].min(), pos[:, 1].max()

is_water_O = np.array([sp == "O" and 0.2 * cL < p[2] < 0.8 * cL
                       for p, sp in kept])

# --- bonds (real 3D distances), drawn under the atoms -------------------
BOND_INK = "#8f8d84"
n = len(kept)
for i in range(n):
    for j in range(i + 1, n):
        pi, pj = pos[i], pos[j]
        d3 = np.linalg.norm(pi - pj)
        pair = {spc[i], spc[j]}
        dmean = 0.5 * (pi[1] + pj[1] - 2 * ymin) / (ymax - ymin + 1e-9)
        if pair == {"Co", "O"} and d3 < 2.15:            # CoO6 sheet bond
            col = _blend(BOND_INK, PAPER, 0.45 * dmean)
            axB.plot([pi[0], pj[0]], [pi[2], pj[2]], "-", color=col,
                     lw=1.5, solid_capstyle="round", zorder=2)
        elif pair == {"O", "H"} and d3 < 1.15:           # water O-H capsule
            col = _blend(BOND_INK, PAPER, 0.45 * dmean)
            axB.plot([pi[0], pj[0]], [pi[2], pj[2]], "-", color=col,
                     lw=2.4, solid_capstyle="round", zorder=2.2)
        elif pair == {"Na", "O"} and d3 < 3.15 and \
                (is_water_O[i] or is_water_O[j]):        # Na...O contact
            axB.plot([pi[0], pj[0]], [pi[2], pj[2]], ls=(0, (2, 2)),
                     color=S.C_MUT, lw=0.8, zorder=1.6)

# --- atoms, back (far y) -> front (near), flat de-saturation for depth ---
order = sorted(range(n), key=lambda i: -pos[i][1])
for k, i in enumerate(order):
    x, y, z = pos[i]
    d = (y - ymin) / (ymax - ymin + 1e-9)       # 0 near, 1 far
    fc = _blend(SPEC[spc[i]]["color"], PAPER, 0.20 * d)
    axB.add_patch(Circle((x, z), SPEC[spc[i]]["r"], facecolor=fc,
                         edgecolor="white", linewidth=0.9, zorder=5 + 0.01 * k))

# --- projected unit cell (a1 x a3 face at y=0), thin dashed ink ----------
axB.plot([0, aL, aL, 0, 0], [0, 0, cL, cL, 0], ls=(0, (4, 2.5)),
         color=S.C_INK, lw=0.9, zorder=9)
# label sits in the empty band between the sheet apical O and the water H
axB.text(0.5 * aL, 1.98, "unit cell", ha="center", va="center", fontsize=7.0,
         color=S.C_INK, zorder=9)

# --- layer labels (right edge, muted ink) --------------------------------
axB.text(XW[1] + 0.25, 0.0, r"CoO$_2$", fontsize=7.2, color=S.C_SEC,
         va="center", ha="left")
axB.text(XW[1] + 0.25, cL, r"CoO$_2$", fontsize=7.2, color=S.C_SEC,
         va="center", ha="left")
axB.text(XW[1] + 0.25, cL / 2, "Na +\n4 H$_2$O", fontsize=7.2, color=S.C_SEC,
         va="center", ha="left", linespacing=1.2)

# --- single scale cue: c-spacing double arrow at the left ----------------
xarr = -2.35
axB.annotate("", xy=(xarr, cL), xytext=(xarr, 0.0),
             arrowprops=dict(arrowstyle="<->", color=S.C_SEC, lw=0.9,
                             shrinkA=0, shrinkB=0))
axB.text(xarr - 0.30, cL / 2, r"$c\simeq 9.9$ Å", rotation=90, va="center",
         ha="right", fontsize=7.2, color=S.C_SEC)

axB.set_xlim(-3.6, XW[1] + 2.6)
axB.set_ylim(ZW[0] - 0.15, ZW[1] + 0.15)
axB.set_aspect("equal", adjustable="box")
axB.axis("off")

# species key: labelled circles in one tidy row below the panel
handles = [Line2D([], [], marker="o", ls="", mfc=SPEC[s]["color"],
                  mec="white", mew=0.8, ms=SPEC[s]["ms"], label=SPEC[s]["label"])
           for s in ("Na", "O", "Co", "H")]
axB.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 0.015),
           ncol=4, frameon=False, handletextpad=0.25, columnspacing=1.3,
           fontsize=8, borderaxespad=0.0)
S.panel_label(axB, "b", x=0.02, y=0.98)

fig.subplots_adjust(left=0.088, right=0.99, top=0.94, bottom=0.155, wspace=0.04)
S.save(fig, "fig6")
print("fig6 written")

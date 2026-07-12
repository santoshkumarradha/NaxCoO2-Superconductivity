#!/usr/bin/env python
"""
Figure 8 - Fermi-surface map of the opened gallery (Na 1x1, c = 9.9 A).

A(k, E_F) in the k_x-k_y plane (k_z = 0), built from the dense uniform-grid
NSCF eigenvalues (runpod/results_v4/jobs/nscf_Na_1x1_c9.9/pw.out, a 24x24x8
Gamma-centred mesh).  The irreducible k_z=0 wedge is expanded to the full
Brillouin zone with the 12 C6v operations (A(k) is a scalar invariant, so this
is exact, not an interpolation guess); each band is then interpolated onto a
fine grid and Lorentzian-broadened at E_F (HWHM 50 meV) and summed over bands
and both spins.  The hexagonal 1x1 BZ is drawn with Gamma, M, K marked.

What ARPES sees in real Na_xCoO2: a single large a1g hole barrel around Gamma;
the six e_g' hole pockets predicted by LDA along Gamma-K are NEVER observed.
The caption reports what THIS calculation shows (checked against the map below).

Aesthetic-sobriety rule: colour is the spectral intensity (data) on a single-hue
sequential map; BZ frame, high-symmetry marks and annotations are neutral ink.
"""
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, PowerNorm
from scipy.interpolate import griddata

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "theory" / "spectral"))
import _style as S
import spectral as SP

S.use_house_style()

JOB = "nscf_Na_1x1_c9.9"
ETA = 0.050                       # Lorentzian HWHM at E_F, eV (50 meV)
NG = 400                          # display grid resolution

# ------------------------------------------------------------------ data
kxy, Eb, ef = SP.load_nscf_kz0(JOB)
K, E = SP.hex_symmetrize(kxy, Eb)           # full-zone scattered set

# hexagonal BZ geometry (2pi/alat cart): K corners at radius R_K, angle 0,60..
# K = (B1+B2)/3 in cartesian -> |K| sets the hexagon corner radius.
R_K = float(np.hypot(*((SP.B1 + SP.B2) / 3.0)))
ang = np.deg2rad(np.arange(0, 360, 60))
hexK = np.c_[R_K * np.cos(ang), R_K * np.sin(ang)]
R_M = 0.5 * np.hypot(*SP.B2)                 # M at edge midpoint

# fine display grid, masked to the hexagon
gx = np.linspace(-R_K * 1.02, R_K * 1.02, NG)
GX, GY = np.meshgrid(gx, gx)

# interpolate each band's energy onto the grid, accumulate Lorentzian at E_F
A = np.zeros_like(GX)
for n in range(E.shape[1]):
    En = griddata(K, E[:, n], (GX, GY), method="cubic")
    d = GX * 0 + (ef - En)
    A += (ETA / np.pi) / (d * d + ETA * ETA)
A = np.nan_to_num(A)


def in_hexagon(px, py, verts):
    """point-in-polygon for the convex BZ hexagon."""
    inside = np.ones(px.shape, bool)
    n = len(verts)
    for i in range(n):
        x1, y1 = verts[i]
        x2, y2 = verts[(i + 1) % n]
        cross = (x2 - x1) * (py - y1) - (y2 - y1) * (px - x1)
        inside &= cross >= -1e-9
    return inside


mask = in_hexagon(GX, GY, hexK)
A = np.where(mask, A, np.nan)

# ------------------------------------------------------------------ plot
BLUE = LinearSegmentedColormap.from_list(
    "fs", ["#ffffff", "#cfe0f5", S.C_BLUE, "#0d366b"])

fig, ax = plt.subplots(figsize=(4.0, 3.7))
im = ax.imshow(A, origin="lower",
               extent=[gx[0], gx[-1], gx[0], gx[-1]], aspect="equal",
               cmap=BLUE, norm=PowerNorm(0.6, vmin=0, vmax=np.nanpercentile(A, 99.5)),
               interpolation="bilinear", zorder=1)

# BZ hexagon outline
hx = np.vstack([hexK, hexK[:1]])
ax.plot(hx[:, 0], hx[:, 1], "-", color=S.C_INK, lw=1.0, zorder=3)

# experimental a1g barrel from ARPES on the anhydrous parent Na0.73CoO2
# (Geck 2007): k_F = 0.65-0.70 1/Ang.  a=2.888 Ang -> 2pi/a = 2.175 1/Ang, so
# k_F = 0.30-0.32 in these (2pi/a) units.  Dashed = experimental REFERENCE
# (neutral), distinct from our filled Na-2DEG barrel (blue = our data).
A_ANG = 5.45752905 * 0.52917721        # in-plane a (Ang)
kf_exp = np.array([0.65, 0.70]) / (2 * np.pi / A_ANG)   # -> 2pi/a units
th = np.linspace(0, 2 * np.pi, 200)
for r in kf_exp:
    ax.plot(r * np.cos(th), r * np.sin(th), ls=(0, (4, 3)), color=S.C_SEC,
            lw=0.9, zorder=3)
rr = kf_exp.mean()
ax.annotate(r"ARPES $a_{1g}\,k_F$" "\n" r"(Na$_{0.73}$, dry)",
            xy=(rr * np.cos(np.deg2rad(215)), rr * np.sin(np.deg2rad(215))),
            xytext=(-0.55, -0.40), fontsize=6.0, color=S.C_SEC,
            ha="center", va="center",
            arrowprops=dict(arrowstyle="-", color=S.C_SEC, lw=0.7))

# high-symmetry points
ax.plot(0, 0, "o", ms=3.2, color=S.C_INK, zorder=4)
ax.text(0.02, 0.02, r"$\Gamma$", fontsize=9, color=S.C_INK, zorder=4)
ax.plot(*hexK[0], "o", ms=3.0, color=S.C_INK, zorder=4)
ax.text(hexK[0, 0] + 0.015, hexK[0, 1] + 0.02, "K", fontsize=9,
        color=S.C_INK, ha="left", va="bottom", zorder=4)
Mpt = np.array([R_M * np.cos(np.deg2rad(30)), R_M * np.sin(np.deg2rad(30))])
ax.plot(*Mpt, "o", ms=3.0, color=S.C_INK, zorder=4)
ax.text(Mpt[0] + 0.02, Mpt[1], "M", fontsize=9, color=S.C_INK,
        ha="left", va="center", zorder=4)

# colourbar (intensity = data)
cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
cb.set_label(r"$A(\mathbf{k},E_F)$  (arb. units)", fontsize=7.6)
cb.ax.tick_params(labelsize=6.5)
cb.outline.set_linewidth(0.6)

ax.set_xlabel(r"$k_x$  ($2\pi/a$)")
ax.set_ylabel(r"$k_y$  ($2\pi/a$)")
lim = R_K * 1.08
ax.set_xlim(-lim, lim)
ax.set_ylim(-lim, lim)
ax.set_xticks([-0.5, 0.0, 0.5])
ax.set_yticks([-0.5, 0.0, 0.5])

# neutral annotations: name the single sheet, flag the empty K corners
ax.annotate("electron barrel\n(Na 2DEG)", xy=(0.0, 0.36), xytext=(-0.30, 0.50),
            fontsize=6.8, color=S.C_INK, ha="center", va="center",
            arrowprops=dict(arrowstyle="-", color=S.C_SEC, lw=0.7))
ax.annotate(r"no $e_g'$ pocket", xy=hexK[0] * 0.93, xytext=(0.30, -0.52),
            fontsize=6.8, color=S.C_INK, ha="center", va="center",
            arrowprops=dict(arrowstyle="->", color=S.C_SEC, lw=0.7))

fig.subplots_adjust(left=0.14, right=0.99, top=0.97, bottom=0.12)
S.save(fig, "fig8")

# ------------------------------------------------------------------ report
# quantify sheets: intensity around Gamma (barrel) vs near K (e_g' check)
def ring_intensity(r0, r1):
    R = np.hypot(GX, GY)
    sel = mask & (R >= r0) & (R < r1)
    return np.nanmean(np.where(sel, A, np.nan))

print(f"fig8 written. E_F={ef:.3f} eV")
print(f"  mean A: Gamma barrel(0-0.35)={ring_intensity(0,0.35):.1f}  "
      f"mid(0.35-0.55)={ring_intensity(0.35,0.55):.1f}  "
      f"near-K(0.55-0.67)={ring_intensity(0.55,0.67):.1f}")

#!/usr/bin/env python
"""
Figure 3 - what a magnetic field does to the ordered phase.

(a) Clausius-Clapeyron construction.  The upper stability edge of the ordered
    phase moves as dTc/dB = Delta M / Delta S, so a field series measures the
    latent magnetization per unit latent entropy.  Lines are drawn for a
    reference latent entropy and three values of Delta M; the ordered phase,
    the one that carries the LARGER magnetization, is stabilized by the field.
(b) The collective-moment headcount, on the corrected criterion.  What the
    Clausius-Clapeyron shift needs is not that a Zeeman energy beat kB T ln 2,
    but that the unit actually DELIVER a latent magnetization per cobalt,

        m(N) = rho tanh(N mu_B B / kB T),   rho = 0.48 mu_B/Co,

    at least M_FACTOR = 2.1 times the latent entropy sigma it is asked to
    carry.  Each sigma is therefore a horizontal requirement line and its
    crossing with m(N) is the minimum headcount: N ~ 10 for sigma = 0.05 kB and
    N ~ 21 for sigma = 0.10 kB at 150 K.  Because m saturates at rho, no
    headcount whatever can serve sigma above rho / 2.1 = 0.229 kB, which is the
    ceiling drawn across the panel.  A single spin delivers 0.01 mu_B/Co, two
    orders below the weakest requirement, so it is excluded outright.
(c) Discrimination between the candidate field mechanisms.  Every model
    amplitude is normalized at 5 T, the one field at which anything is known,
    so the 1-9 T series separates them with no free parameters.  The Brillouin
    shape is drawn at the paper's own headcount, N = 20, and at N = 100: at
    N = 20 the response is nearly linear out to 9 T, so measurable curvature in
    a field series would itself indicate clusters larger than the minimum
    headcount.  The shaded band is the projected precision of the proposed
    series, so that the panel advertises only the discrimination it can
    actually deliver; see PROJECTED PRECISION below.

Model curves only.  No experimental point is plotted anywhere in this figure;
the numbers on the axes are model parameters.

----------------------------------------------------------------------------
PROJECTED PRECISION (panel c)
----------------------------------------------------------------------------
BAND_REL is the half-width of the band on the normalized amplitude ratio, set
to 4.5 percent: the 2-sigma band that just separates the Brillouin N = 20 shape
from the linear reference at 9 T, which corresponds to +-3.2 percent per-point
precision on the excess amplitude.  It is a projection, not a measured error.
The band is drawn as y_ref (1 +- BAND_REL) about the N = 20 Brillouin curve,
the paper's own headcount and therefore the shape a field series would be
testing against.

The 5 T normalization does NOT enter the band.  All four predictions are
divided by the SAME measured 5 T amplitude, so its uncertainty is common mode:
it slides the four curves together and cannot change which shape is selected.
What discriminates is the per-field precision alone.
"""
from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

import _style as S

S.use_house_style()

PNG_DIR = S.FIGDIR            # inspection PNG beside the vector PDF

G_FACTOR = 2.0

# Projected per-field relative precision on the excess amplitude of the
# proposed 1-9 T series (panel c).  See PROJECTED PRECISION in the docstring;
# this single number sets the band width.
BAND_REL = 0.045
BAND_FIELDS = [1.0, 2.0, 3.0, 5.0, 7.0, 9.0]   # the proposed field series

fig, (axa, axb, axc) = plt.subplots(1, 3, figsize=(7.0, 2.35))

# ======================================================================
# (a) Clausius-Clapeyron:  Tc(B) = Tc(0) + (Delta M / Delta S) B
# ======================================================================
B = np.linspace(0.0, 9.0, 200)
DS_KB = 0.2                       # latent entropy per Co, in units of kB
DM_LIST = [0.2, 0.4, 0.5]         # latent magnetization per Co, in mu_B
REF = 1                           # the line that carries the shading

# dTc/dB [K/T] = (muB/kB) (Delta M / muB) / (Delta S / kB), with
# muB/kB = 0.6717 K/T.  The illustrative Delta S = 0.2 kB per Co respects the
# bound Delta S <= 0.24 kB per Co that follows from Delta M <= 0.5 muB per Co at
# x = 1/2, and the heavy contour, Delta M = 0.4 muB per Co, then carries
# dTc/dB = 1.34 K/T, the measured shift.  The top contour sits at the fully
# polarized ceiling of 0.5 muB per Co.
slopes = [S.MUB_OVER_KB * dm / DS_KB for dm in DM_LIST]   # K / T
lines = [S.TC0 + s * B for s in slopes]

YLO, YHI = 236.0, 260.0
axa.fill_between(B, YLO, lines[REF], color=S.FILL_RED, lw=0, zorder=0)
axa.fill_between(B, lines[REF], YHI, color=S.FILL_BLUE, lw=0, zorder=0)
axa.fill_betweenx([YLO, YHI], 9.0, 10.6, color="white", lw=0, zorder=1)

for i, (dm, ln) in enumerate(zip(DM_LIST, lines)):
    heavy = i == REF
    axa.plot(B, ln, color=S.C_RED if heavy else S.C_SEC,
             lw=1.7 if heavy else 1.0,
             ls="-" if heavy else (0, (4, 2)), zorder=4)
    axa.text(9.25, ln[-1], f"${dm:g}$", color=S.C_RED if heavy else S.C_SEC,
             fontsize=7.0, ha="left", va="center")
axa.text(9.25, 259.5, r"$\Delta M$" "\n" r"($\mu_{\mathrm{B}}$/Co)",
         color=S.C_INK, fontsize=7.0, ha="left", va="top", linespacing=1.25)

axa.text(0.45, 255.6, "parent P", color=S.INK_BLUE, fontsize=7.5,
         ha="left", va="center", zorder=5)
axa.text(0.45, 238.3, "ordered O\n(larger $M$)", color=S.INK_RED, fontsize=7.5,
         ha="left", va="center", linespacing=1.3, zorder=5)
axa.text(8.75, 236.9,
         r"$\dfrac{dT_{\mathrm{c}}}{dB}=\dfrac{\Delta M}{\Delta S}$,"
         "\n" r"$\Delta S=0.2\,k_{\mathrm{B}}$/Co",
         color=S.C_INK, fontsize=7.0, ha="right", va="bottom",
         linespacing=1.6, zorder=5)

axa.set_xlim(0.0, 10.6)
axa.set_ylim(YLO, YHI)
axa.set_xticks([0, 2, 4, 6, 8])
axa.set_yticks([240, 245, 250, 255])
axa.set_xlabel(r"magnetic field  $B$  (T)")
axa.set_ylabel(r"collapse temperature  $T_{\mathrm{c}}(B)$  (K)")
S.panel_label(axa, "a")

# ======================================================================
# (b) headcount: magnetization DELIVERED per Co against what is required
# ======================================================================
# A field-coupled unit of N aligned moments polarizes as tanh(N muB B / kB T),
# so the latent magnetization it can deliver per cobalt is
#
#     m(N) = rho tanh(N muB B / kB T),   rho = 0.48 muB/Co the saturation value.
#
# The Clausius-Clapeyron shift needs m >= M_FACTOR sigma for a latent entropy
# sigma, so each sigma sets a horizontal requirement line and the crossing sets
# the minimum headcount.  m saturates at rho as N -> infinity, which caps the
# largest latent entropy any headcount can serve at rho / M_FACTOR.
B_HEAD = 5.0
RHO_MUB = 0.48          # muB per Co, saturation magnetization of the unit
M_FACTOR = 2.1          # required m in units of the latent entropy sigma
T_LO, T_HI = 150.0, 183.0
SIGMAS = [0.05, 0.10]   # kB per Co
NMAX = 60.0

N = np.linspace(0.0, NMAX, 600)


def m_delivered(n, T):
    """Latent magnetization per Co delivered by a unit of n aligned moments."""
    return RHO_MUB * np.tanh(np.asarray(n, float) * S.MU_B_EV * B_HEAD
                             / (S.KB_EV * T))


def n_required(m, T):
    """Headcount at which m_delivered first reaches m (inf if unreachable)."""
    if m >= RHO_MUB:
        return np.inf
    return (np.arctanh(m / RHO_MUB) * S.KB_EV * T / (S.MU_B_EV * B_HEAD))


M_REQ = [M_FACTOR * s for s in SIGMAS]                 # 0.105, 0.21 muB/Co
N_CROSS = [n_required(m, T_LO) for m in M_REQ]         # ~10 and ~21 at 150 K
SIGMA_CEIL = RHO_MUB / M_FACTOR                        # 0.229 kB per Co

# excluded: below the weakest requirement nothing works, the single spin least
axb.axvspan(0.0, N_CROSS[0], color=S.FILL_NONE, lw=0, zorder=0)
# the window of requirements spanned by sigma = 0.05 - 0.10 kB
axb.axhspan(M_REQ[0], M_REQ[1], color=S.FILL_AMBER, lw=0, zorder=1)
for m in M_REQ:
    axb.plot([0.0, NMAX], [m, m], color=S.INK_AMBER, lw=1.0, ls=(0, (4, 2)),
             zorder=3)
# the N -> infinity ceiling
axb.plot([0.0, NMAX], [RHO_MUB, RHO_MUB], color=S.C_MUT, lw=0.9,
         ls=(0, (1.4, 1.4)), zorder=3)

for T, col, ls, lw in [(T_LO, S.C_RED, "-", 1.7),
                       (T_HI, S.INK_RED, (0, (3.4, 1.5)), 1.4)]:
    axb.plot(N, m_delivered(N, T), color=col, lw=lw, ls=ls, zorder=4)

for n_c, m in zip(N_CROSS, M_REQ):
    axb.plot([n_c], [m], "o", ms=4.2, mfc="white", mec=S.C_RED, mew=1.1,
             zorder=6)
axb.plot([1.0], [m_delivered(1.0, T_LO)], "o", ms=2.6, mfc=S.C_SEC,
         mec="none", zorder=6)

axb.text(12.5, 0.450, rf"excluded: $N\lesssim{N_CROSS[0]:.0f}$",
         color=S.C_SEC, fontsize=7.0, ha="left", va="top", zorder=6)
axb.annotate(rf"$N=1$: {m_delivered(1.0, T_LO):.3f}$\,\mu_{{\mathrm{{B}}}}$/Co",
             xy=(1.0, m_delivered(1.0, T_LO)), xytext=(16.0, 0.040),
             color=S.C_SEC, fontsize=7.0, ha="left", va="bottom", zorder=7,
             arrowprops=dict(arrowstyle="-", lw=0.6, color=S.C_MUT,
                             shrinkA=2, shrinkB=2))
axb.text(N_CROSS[0] - 0.6, 0.115, rf"$N\simeq{N_CROSS[0]:.0f}$",
         color=S.C_INK, fontsize=7.5, ha="right", va="bottom", zorder=7)
axb.text(N_CROSS[1] - 1.2, 0.235, rf"$N\simeq{N_CROSS[1]:.0f}$",
         color=S.C_INK, fontsize=7.5, ha="right", va="bottom", zorder=7)

axb.text(NMAX - 1.0, M_REQ[1] + 0.005, rf"$\sigma={SIGMAS[1]:.2f}\,k_{{\mathrm{{B}}}}$",
         color=S.INK_AMBER, fontsize=7.0, ha="right", va="bottom", zorder=6)
axb.text(NMAX - 1.0, M_REQ[0] + 0.005, rf"$\sigma={SIGMAS[0]:.2f}\,k_{{\mathrm{{B}}}}$",
         color=S.INK_AMBER, fontsize=7.0, ha="right", va="bottom", zorder=6)
axb.text(NMAX - 1.0, M_REQ[1] - 0.012, rf"required $m={M_FACTOR:g}\,\sigma$",
         color=S.INK_AMBER, fontsize=7.0, ha="right", va="top", zorder=6)
axb.text(NMAX - 1.0, RHO_MUB + 0.007,
         rf"$m\to\rho={RHO_MUB:g}$ ($\sigma={SIGMA_CEIL:.2f}\,k_{{\mathrm{{B}}}}$)",
         color=S.C_SEC, fontsize=7.0, ha="right", va="bottom", zorder=6)

axb.annotate(rf"$T={T_LO:g}$ K", xy=(52.0, m_delivered(52.0, T_LO)),
             xytext=(46.0, 0.437), color=S.C_RED, fontsize=7.0,
             ha="center", va="bottom", zorder=6,
             arrowprops=dict(arrowstyle="-", lw=0.6, color=S.C_RED,
                             alpha=0.8, shrinkA=2, shrinkB=2))
axb.annotate(rf"$T={T_HI:g}$ K", xy=(48.0, m_delivered(48.0, T_HI)),
             xytext=(51.0, 0.288), color=S.INK_RED, fontsize=7.0,
             ha="center", va="top", zorder=6,
             arrowprops=dict(arrowstyle="-", lw=0.6, color=S.INK_RED,
                             alpha=0.8, shrinkA=2, shrinkB=2))

axb.set_xlim(0.0, NMAX)
axb.set_ylim(0.0, 0.52)
axb.set_xticks([0, 20, 40, 60])
axb.set_yticks([0.0, 0.1, 0.2, 0.3, 0.4, 0.5])
axb.set_xlabel(r"aligned moments in the field-coupled unit  $N$")
axb.set_ylabel(r"delivered  $m$  ($\mu_{\mathrm{B}}$/Co)")
S.panel_label(axb, "b")

# ======================================================================
# (c) shapes of Delta R(B), all normalized at the 5 T anchor
# ======================================================================
BB = np.linspace(0.0, 9.0, 400)
B_ANCH = 5.0
T_AMP = 180.0
# Two headcounts are drawn.  N = 20 is the paper's own bound (panel b): at that
# size the 0-9 T response is nearly linear, so any measurable curvature in a
# field series already implies a cluster larger than the minimum headcount.
# N = 100 shows what a genuinely saturating collective moment would look like.
N_SMALL, N_LARGE = 20.0, 100.0


def brillouin(b, n):
    """Saturating response of a collective moment of n aligned spins at T_AMP."""
    x = n * G_FACTOR * S.MU_B_EV * b / (2.0 * S.KB_EV * T_AMP)
    return np.tanh(x)


def percolation(b, b0=4.0, w=1.0):
    """Connectivity threshold: nothing, then everything."""
    return 0.5 * (1.0 + np.tanh((b - b0) / w))


MODELS = [
    (lambda b: brillouin(b, N_SMALL), S.C_RED, "-",
     rf"Brillouin ($N={N_SMALL:.0f}$)", 1.30),
    (lambda b: brillouin(b, N_LARGE), S.INK_RED, (0, (3.4, 1.5)),
     rf"Brillouin ($N={N_LARGE:.0f}$)", 1.30),
    (lambda b: b**2, S.C_BLUE, (0, (5, 2)), r"$\propto B^{2}$", 3.05),
    (lambda b: percolation(b), S.C_AQUA, (0, (1.4, 1.4)),
     "percolation", 0.72),
]

axc.axvline(B_ANCH, color=S.C_MUT, lw=0.8, ls=(0, (2, 2)), zorder=2)
axc.plot([B_ANCH], [1.0], "o", ms=4.0, mfc="white", mec=S.C_SEC, mew=1.0,
         zorder=6)
axc.text(B_ANCH + 0.22, 0.10, "anchored at 5 T", color=S.C_SEC, fontsize=7.0,
         ha="left", va="bottom", zorder=6)

# Projected precision of the proposed series, drawn about the N = 20 Brillouin
# shape (the paper's own headcount) and lying behind every model curve.
y_ref = brillouin(BB, N_SMALL) / brillouin(B_ANCH, N_SMALL)
axc.fill_between(BB, (1.0 - BAND_REL) * y_ref, (1.0 + BAND_REL) * y_ref,
                 color=S.FILL_RED, lw=0, zorder=1)
for b in BAND_FIELDS:
    yb = float(np.interp(b, BB, y_ref))
    axc.plot([b, b], [(1.0 - BAND_REL) * yb, (1.0 + BAND_REL) * yb],
             color=S.INK_RED, lw=0.7, alpha=0.55, zorder=3)

for fn, col, ls, lab, _ in MODELS:
    y = np.asarray(fn(BB), float) / float(fn(B_ANCH))
    axc.plot(BB, y, color=col, lw=1.7, ls=ls, zorder=4, label=lab)

axc.set_xlim(0.0, 9.3)
axc.set_ylim(0.0, 3.45)
axc.set_xticks([0, 2, 4, 6, 8])
axc.set_yticks([0, 1, 2, 3])
axc.set_xlabel(r"magnetic field  $B$  (T)")
axc.set_ylabel(r"$\Delta R(B)\,/\,\Delta R(5\,\mathrm{T})$")
handles, labels = axc.get_legend_handles_labels()
handles.append(Patch(facecolor=S.FILL_RED, edgecolor="none"))
labels.append(rf"projected $\pm{BAND_REL * 100:g}\%$ (1–9 T)")
legc = axc.legend(handles, labels, loc="upper left",
                  bbox_to_anchor=(0.02, 0.885),
                  title="field-response predictions", fontsize=7.0,
                  handlelength=1.7, labelspacing=0.3, borderaxespad=0.0)
legc.get_title().set_fontsize(7.0)
legc.get_title().set_color(S.C_SEC)
legc._legend_box.align = "left"
S.panel_label(axc, "c")

fig.subplots_adjust(left=0.072, right=0.988, top=0.975, bottom=0.175,
                    wspace=0.34)
out = S.save(fig, "fig3_field", png_dir=PNG_DIR)
print(f"wrote {out}")
print(f"  dTc/dB = {slopes[REF]:.2f} K/T for DM = {DM_LIST[REF]} muB, "
      f"DS = {DS_KB} kB")
for s, m, n_c in zip(SIGMAS, M_REQ, N_CROSS):
    print(f"  sigma = {s:.2f} kB -> m = {m:.3f} muB/Co: N = {n_c:.1f} at "
          f"{T_LO:g} K, {n_required(m, T_HI):.1f} at {T_HI:g} K")
print(f"  single spin delivers {m_delivered(1.0, T_LO):.4f} muB/Co "
      f"(tanh = {m_delivered(1.0, T_LO) / RHO_MUB:.4f});  ceiling rho = "
      f"{RHO_MUB:g} muB/Co = sigma {SIGMA_CEIL:.3f} kB")
for n in (N_SMALL, N_LARGE):
    r9 = brillouin(9.0, n) / brillouin(B_ANCH, n)
    print(f"  N = {n:.0f}: DR(9 T)/DR(5 T) = {r9:.2f}  (linear would be 1.80)")

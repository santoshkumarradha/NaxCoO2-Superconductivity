#!/usr/bin/env python
"""
Figure 1 - the mechanism.

(a) Landau free energy f(phi) = a(T) phi^2/2 - b phi^3/3 + c phi^4/4 at three
    temperatures.  The cubic invariant makes the transition first order: above
    Tc the ordered minimum exists but lies above the disordered one, at Tc the
    two are degenerate and the barrier separating them is maximal relative to
    the (vanishing) driving force -- the 2D nucleation barrier
    W* = pi sigma^2/|Delta f| diverges there -- and below Tc the ordered
    minimum deepens while W* collapses.
(b) The protocol in the (T, phi) plane: the cool-down rides the phi = 0 branch
    because nucleation is bypassed, the warm-up unfreezes at T_on and the
    ordered phase collapses thermodynamically at Tc.

Everything drawn here is model output.  No experimental trace is used; the
published scalars T_on ~ 146 K and Tc ~ 240 K enter only as model parameters.
"""
from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch

import _style as S

S.use_house_style()

PNG = ("/private/tmp/claude-501/-Users-santoshkumar/"
       "57bc6bbc-f3db-4267-a42f-5531905aedd5/scratchpad")

fig, (axa, axb) = plt.subplots(2, 1, figsize=(3.375, 4.45))

# ======================================================================
# (a) Landau free energy
# ======================================================================
TEMPS = [(252.0, S.C_BLUE, "-",         r"$252$ K  $(T>T_{\mathrm{c}})$"),
         (240.0, S.C_SEC,  (0, (5, 2)), r"$240$ K  $(T=T_{\mathrm{c}})$"),
         (210.0, S.C_RED,  "-",         r"$210$ K  $(T<T_{\mathrm{c}})$")]

phi = np.linspace(-0.03, 0.92, 900)
axa.axhline(0.0, color=S.C_GRID, lw=0.7, zorder=0)

handles = []
for T, col, ls, lab in TEMPS:
    axa.plot(phi, S.landau_f(phi, T), color=col, lw=1.6, ls=ls, zorder=3)
    pb, po = S.landau_extrema(T)
    axa.plot(po, S.landau_f(po, T), "o", ms=3.6, mfc=col, mec="white",
             mew=0.7, zorder=5)
    handles.append(Line2D([], [], color=col, ls=ls, lw=1.6, label=lab))

# the disordered minimum, common to all three curves
axa.plot(0.0, 0.0, "o", ms=3.6, mfc=S.C_INK, mec="white", mew=0.7, zorder=6)

axa.set_xlim(-0.03, 0.92)
axa.set_ylim(-1.30, 1.05)
axa.set_xticks([0.0, 0.2, 0.4, 0.6, 0.8])
axa.set_yticks([-1.0, -0.5, 0.0, 0.5, 1.0])
axa.set_xlabel(r"order parameter  $\phi$")
axa.set_ylabel(r"$f(\phi)-f(0)$   (meV per Co)")

leg = axa.legend(handles=handles, loc="lower left",
                 bbox_to_anchor=(0.005, -0.012), fontsize=7.0,
                 handlelength=1.6, labelspacing=0.32, borderaxespad=0.35)
leg.set_zorder(9)

# --- phase identity of the two minima -------------------------------------
axa.annotate("P  disordered", xy=(0.004, -0.02), xytext=(0.075, -0.34),
             fontsize=7.5, color=S.C_INK, ha="left", va="center",
             arrowprops=dict(arrowstyle="-", lw=0.7, color=S.C_MUT,
                             shrinkA=2, shrinkB=3))
axa.text(0.755, -1.06, "O  ordered", fontsize=7.5, color=S.C_RED,
         ha="center", va="top")

# --- the nucleation barrier on the degenerate (T = Tc) curve --------------
pb_c = S.landau_extrema(240.0)[0]
fb_c = S.landau_f(pb_c, 240.0)
axa.add_patch(FancyArrowPatch((pb_c, 0.0), (pb_c, fb_c),
                              arrowstyle="<|-|>", mutation_scale=6,
                              lw=0.9, color=S.C_INK, shrinkA=0, shrinkB=0,
                              zorder=8))
axa.annotate(r"$W^{*}\propto 1/|\Delta f|$" "\n"
             r"diverges at $T_{\mathrm{c}}$",
             xy=(pb_c, 0.55 * fb_c), xytext=(0.115, 0.755), fontsize=7.0,
             color=S.C_INK, ha="left", va="center", linespacing=1.3,
             arrowprops=dict(arrowstyle="-", lw=0.7, color=S.C_MUT,
                             shrinkA=3, shrinkB=3))

S.panel_label(axa, "a")

# ======================================================================
# (b) protocol in the (T, phi) plane
# ======================================================================
TW = np.arange(80.0, 300.0 + 1e-9, 0.02)
TCOOL = np.arange(300.0, 80.0 - 1e-9, -0.02)
phi_w = S.warming_phi(TW, 1.0)
phi_c = S.cooling_phi(TCOOL, 1.0)
Ton = S.onset_temperature(TW, phi_w)

# the master equation is integrated on a 0.02 K grid; curves are drawn on
# every 10th point, which is still ~60 points across the 12 K rise and keeps
# the vector file small.
D = slice(None, None, 10)
axb.plot(TW[D], S.phi_eq(TW[D]), color=S.C_MUT, lw=1.0, ls=(0, (1.6, 1.6)),
         zorder=2)
axb.text(93.0, 0.855, r"$\phi_{\mathrm{eq}}(T)$", color=S.C_SEC,
         fontsize=7.5, ha="left", va="center")

axb.fill_between(TW[D], 0.0, phi_w[D], color=S.FILL_RED, lw=0, zorder=1)
axb.text(0.5 * (Ton + S.TC0) + 5.0, 0.45, "ordered\nphase O", color=S.INK_RED,
         fontsize=7.5, ha="center", va="center", linespacing=1.3, zorder=4)

# the two branches coincide wherever phi = 0, so the cool-down is drawn as a
# wider grey underlay and the warm-up as a narrower red line on top: the halo
# makes the shared stretches legible without displacing either curve.
axb.plot(TCOOL[D], phi_c[D], color=S.C_SEC, lw=2.8, zorder=5,
         solid_capstyle="butt")
axb.plot(TW[D], phi_w[D], color=S.C_RED, lw=1.5, zorder=6)

# direction arrows: cooling to the left, warming to the right, never sharing
# a stretch of the phi = 0 line where the two branches coincide.
for x, dx, y, col in [(276.0, -20.0, 0.0, S.C_SEC),
                      (122.0, -20.0, 0.0, S.C_SEC)]:
    axb.add_patch(FancyArrowPatch((x, y), (x + dx, y), arrowstyle="-|>",
                                  mutation_scale=7, lw=2.0, color=col,
                                  shrinkA=0, shrinkB=0, zorder=7))
for x, dx, y in [(198.0, 20.0, 1.0), (283.0, 14.0, 0.0)]:
    axb.add_patch(FancyArrowPatch((x, y), (x + dx, y), arrowstyle="-|>",
                                  mutation_scale=7, lw=1.6, color=S.C_RED,
                                  shrinkA=0, shrinkB=0, zorder=7))

axb.text(86.0, 0.42, "cooling:\nnucleation\nbypassed", color=S.C_SEC,
         fontsize=7.5, ha="left", va="center", linespacing=1.3, zorder=8)
axb.text(196.0, 1.055, "warming", color=S.C_RED, fontsize=7.5,
         ha="center", va="bottom", zorder=8)

for x, lab, sub, col, ha, xoff in [
        (Ton, r"$T_{\mathrm{on}}$", "kinetic", S.C_RED, "right", -4.0),
        (S.TC0, r"$T_{\mathrm{c}}$", "thermodynamic", S.C_SEC, "left", 4.0)]:
    axb.axvline(x, color=col, lw=0.8, ls=(0, (2.5, 2.0)), zorder=3)
    axb.text(x + xoff, 1.36, f"{lab}\n{sub}", color=col, fontsize=7.0,
             ha=ha, va="top", linespacing=1.35, zorder=8)

axb.set_xlim(80.0, 300.0)
axb.set_ylim(-0.06, 1.40)
axb.set_yticks([0.0, 0.5, 1.0])
axb.set_xticks([100, 150, 200, 250, 300])
axb.set_xlabel(r"temperature  $T$  (K)")
axb.set_ylabel(r"transformed fraction  $\phi$")
S.panel_label(axb, "b")

fig.subplots_adjust(left=0.155, right=0.975, top=0.985, bottom=0.086,
                    hspace=0.38)
out = S.save(fig, "fig1_mechanism", png_dir=PNG)
print(f"wrote {out}   (T_on = {Ton:.1f} K)")

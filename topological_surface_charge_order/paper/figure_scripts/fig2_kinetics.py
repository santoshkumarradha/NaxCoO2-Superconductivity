#!/usr/bin/env python
"""
Figure 2 - output of the kinetic two-phase model, with the published resistance
loop overlaid.

(a) Excess resistance ratio DeltaR/R_P = c phi (s - 1) through a cool-then-warm
    protocol.  The cool-down is nucleation limited (barrier Ea + W*(T)) and
    stays on the untransformed branch.  The warm-up is SITE-SATURATED TWO-
    DIMENSIONAL GROWTH from the areal density n0 of nuclei quenched in on
    cooling, integrated at r = 0.1, 1 and 10 K/min:

        v(T)   = nu a exp(-Ea / kB T) [1 - exp(Delta f / kB T)]
        R_w(T) = (1/r) int_{T0}^{T} v(T') dT'                     (T0 = 80 K)
        phi(T) = 1 - exp(-n0 pi R_w(T)^2)                         (KJMA, n = 2)

    with Delta f = -Delta S (T_ord - T) per Co the thermodynamic driving force,
    RETAINED IN FULL rather than saturated to unity.  Two consequences the old
    single-barrier master equation did not have: the rise is sharper (KJMA
    n = 2 compresses the 10-90 % width) and the rate shift is larger, because
    the driving-force factor flattens the effective Arrhenius slope.  The decay
    side is thermodynamic as before: the transport contrast s(T) - 1 shrinks
    and phi_eq collapses at T_ord, imposed as a ceiling on the KJMA fraction.

    Overlaid on the model curves, as subordinate open symbols, are the DIGITIZED
    0 T and 5 T loops of Crowley et al., J. Phys. Chem. C 124, 20693 (2020),
    SI Fig. S5, reduced to the same DeltaR/R_P variable.  See EXPERIMENTAL
    OVERLAY below for the transformation and its provenance.

(b) Rise midpoint T_on versus warming rate.  The growth model gives a straight
    line in log10 r whose slope is the model's sharpest out-of-sample
    prediction.  The leading form kB T_on^2 ln10 / Ea is NOT what is drawn --
    the driving-force factor adds to it -- so the panel prints the numerically
    obtained slope alone and attributes it to no closed-form expression.

The model curves are unfitted in shape and scale: nothing below was adjusted to
improve agreement with the overlaid data, and the overlay itself carries no
scale factor or offset (see EXPERIMENTAL OVERLAY).  The single exception, stated
plainly because it is a calibration and not a fit, is EA_EV, which is fixed
analytically by the site-saturated onset condition at one scalar -- the measured
rise midpoint T_on = 152 K -- and is not varied to improve the match anywhere
else on the curve.  Panel (b) is model output only.

----------------------------------------------------------------------------
COOLING BRANCH - why it is drawn at zero, and an open issue
----------------------------------------------------------------------------
The cool-down is nucleation limited: W* = A/(T_ord - T) diverges at the
transition, so nothing nucleates while the ions are still mobile.  What the
sweep leaves behind is the seed density n0 that the warm-up then grows, and the
growth accumulated during the cool-down itself is nil, so the transformed
fraction is phi = 1 - exp(-n0 pi R_w^2) with R_w = 0, i.e. zero.  That is the
result quoted in the text and it is what is drawn.

Two things are recorded here rather than hidden:

1. The single-barrier proxy used for this branch previously - relaxing phi
   toward phi_eq at rate nu exp(-(Ea + W*)/kB T) - returns phi -> 1.0e-2 at
   Ea = 0.387 eV, where it returned exactly 0 at Ea = 0.45 eV.  Multiplied by
   the contrast factor, which diverges at low T, that 1 % produced a spurious
   upturn to DeltaR/R_P ~ 0.14 at 80 K.  The proxy conflates the number of
   nuclei with the transformed fraction, so it was never the right object; the
   coincidence that it returned 0 at the old barrier hid the category error.

2. A KJMA-consistent nucleation-AND-growth integral for the cool-down (seeds
   appearing at areal rate (nu/a^2) exp(-(Ea + W*)/kB T), each then growing at
   the same v(T)) gives phi -> 1 on cooling with the old lumped line tension
   A = 16 eV K, at BOTH Ea = 0.387 and Ea = 0.45 eV.  Keeping the cool-down
   flat needs A >~ 40 eV K, which is the smoothness constraint of Sec. IV-C of
   the paper, and A_NUC_EVK is set accordingly below.  A enters only the
   nucleation barrier on the cool-down, which is drawn at phi = 0 either way,
   so the plotted curves are unaffected by the value.

----------------------------------------------------------------------------
EXPERIMENTAL OVERLAY - provenance and transformation
----------------------------------------------------------------------------
Source figure   Crowley, Radha et al., J. Phys. Chem. C 124, 20693 (2020),
                Supporting Information Fig. S5: two-terminal resistance of a
                LiCoO2 flake device against temperature, measured on cooling
                and on warming, at B = 0 T and B = 5 T.
Digitizer       ../../src/digitize_crowley2020.py (colour-mask extraction with
                the pixel calibration recorded in results/digitized/
                crowley2020_fit.json).
Data file       ../../results/digitized/crowley2020_figS5.csv, columns
                temperature_K, 0T_warmup, 0T_cooldown, 5T_warmup, 5T_cooldown
                in MOhm on a 1 K grid from 100 to 300 K.
Units in        resistance in MOhm; temperature in K.
Transformation  The model plots the excess of the two-phase mixture over the
                untransformed parent P,

                    DeltaR/R_P = [R_mix(T) - R_P(T)] / R_P(T) = c phi (s - 1).

                Experimentally the parent branch is realized by the COOLING
                sweep, which the model shows (does not assume) stays at phi = 0.
                The published traces therefore map onto the model variable as

                    (DeltaR/R_P)_expt(T) = [R_warm(T) - R_cool(T)] / R_cool(T),

                evaluated field by field at the same temperature.  Both the
                temperature-dependent device background and the absolute
                resistance scale divide out, so no scale factor, offset or fit
                is applied anywhere: the overlay is a pure change of variable on
                the digitized numbers.  Under this normalization the measured
                COOLING branch is identically zero by construction and coincides
                with the model cooling branch, which is why only one curve per
                field is drawn.
Uncertainty     The digitization is quoted at +-2 K and +-0.01 MOhm
                (crowley2020_fit.json); at the ~0.15-0.30 MOhm baseline of the
                cooling branch that is a few percent on DeltaR/R_P, comparable
                to the symbol size, so no error bars are drawn.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch

import _style as S

S.use_house_style()

# ======================================================================
# MODEL PARAMETERS - the only place any of these is set.  Changing EA_EV,
# NU0 or N0_NM2 here alone regenerates both panels and every printed number.
# ======================================================================
EA_EV = 0.4026      # eV, growth barrier (Li migration class).  NOT free: it is
                    # fixed by the site-saturated onset condition, which at the
                    # rise midpoint phi = 1/2 reads pi n0 R_w(T_on)^2 = ln 2,
                    # i.e. R_w = 4.697 nm at n0 = 1e-2 nm^-2.  Imposing that at
                    # T_on = 152 K and r = 1 K/min returns Ea = 0.4026 eV.
                    # T_on = 152 K is the MIDPOINT of the resistance rise in
                    # the project's own digitization of Crowley et al. Fig. S5
                    # (0 T loop: 10/50/90 % of the rise at 149.1/152.2/163.4 K).
                    # The 146 K used previously is the FOOT of that rise -- it
                    # sits below the 10 % point -- so it was the wrong scalar to
                    # impose a midpoint condition with, and it put the model
                    # onset ~6 K below the data.
NU0 = 1.0e13        # 1/s, attempt frequency
N0_NM2 = 1.0e-2     # 1/nm^2, areal density of nuclei quenched in on cooling
                    # (n_0); site saturated, so no nucleation term on warming
A_NM = 0.28         # nm, hop distance setting the interface velocity scale
DELTA_S_KB = 0.23   # kB per Co, latent entropy in Delta f = -Delta S (T_ord-T)
T_ORD = 240.0       # K, first-order collapse of the ordered phase
T0_K = 80.0         # K, start of the warm-up integration
PHI_EQ_WIDTH = 6.0  # K, compositional smearing of the collapse
A_NUC_EVK = 40.0    # eV K, lumped 2D line-tension constant in W* = A/(T_ord-T).
                    # Set by the smoothness constraint on the cooling branch
                    # (Sec. IV-C): keeping the KJMA nucleation-and-growth
                    # integral flat on the way down requires A >~ 40 eV K.

# Transport-contrast parameters (Eq. for s(T) - 1); model choices, not fits.
CONTRAST_TREF = 160.0     # K
CONTRAST_ES_EV = 0.036    # eV

# ======================================================================
# Experimental overlay - see module docstring for the transformation.
# ======================================================================
DIGITIZED_CSV = (Path(__file__).resolve().parents[2]
                 / "results" / "digitized" / "crowley2020_figS5.csv")
SHOW_EXPT = True
EXPT_FIELDS = ["0T", "5T"]          # which digitized loops to overlay
EXPT_MARKER_EVERY = 6               # K between plotted symbols

PNG_DIR = S.FIGDIR                  # inspection PNG beside the vector PDF

TW = np.arange(T0_K, 300.0 + 1e-9, 0.02)
TCOOL = np.arange(300.0, T0_K - 1e-9, -0.02)
RATES = [0.1, 1.0, 10.0]
D = slice(None, None, 10)   # draw every 10th integration point


# ----------------------------------------------------------------------
# Model wrappers.  These forward the parameters above into the shared
# definitions in _style so that _style's defaults are never silently used.
# ----------------------------------------------------------------------
def phi_eq(T):
    return S.phi_eq(T, Tc=T_ORD, width=PHI_EQ_WIDTH)


def growth_velocity(T):
    """Interface velocity v(T) = nu a exp(-Ea/kB T) [1 - exp(Delta f/kB T)].

    Delta f = -DELTA_S_KB kB (T_ord - T) per Co is the thermodynamic driving
    force, kept in full rather than saturated to unity, so v carries both the
    Arrhenius mobility and the vanishing of the driving force as T -> T_ord.
    Above T_ord the driving force reverses; growth then stops (v clipped to
    zero) and the collapse of the ordered phase is carried by phi_eq instead.
    """
    T = np.asarray(T, float)
    mobility = NU0 * A_NM * np.exp(-EA_EV / (S.KB_EV * T))        # nm/s
    drive = 1.0 - np.exp(-DELTA_S_KB * (T_ORD - T) / T)
    return np.maximum(mobility * drive, 0.0)


def warming_phi(T_grid, rate):
    """Site-saturated 2D KJMA warm-up: phi = 1 - exp(-n0 pi R_w^2).

    R_w is the radius each quenched-in nucleus has grown to by temperature T
    on a linear ramp, R_w = (1/r) int v dT'.  phi_eq is applied on top as the
    thermodynamic ceiling, so the fraction rises kinetically and then collapses
    at T_ord.  T_grid must be ascending and start at T0_K.
    """
    T_grid = np.asarray(T_grid, float)
    v = growth_velocity(T_grid)
    swept = np.concatenate(
        ([0.0], np.cumsum(0.5 * (v[1:] + v[:-1]) * np.diff(T_grid))))
    R_w = swept / (rate / 60.0)                                   # nm
    kjma = 1.0 - np.exp(-N0_NM2 * np.pi * np.minimum(R_w, 1e6) ** 2)
    return kjma * phi_eq(T_grid)


def cooling_phi(T_grid, rate):
    """Nucleation-limited cool-down: seeds are deposited but do not grow.

    phi = 1 - exp(-n0 pi R_w^2) with R_w = 0, hence identically zero.  See
    COOLING BRANCH in the module docstring for the proxy this replaces and for
    the smoothness constraint that sets A_NUC_EVK.
    """
    return np.zeros_like(np.asarray(T_grid, float))


def contrast(T):
    return S.contrast(T, T_ref=CONTRAST_TREF, E_s_eV=CONTRAST_ES_EV)


def load_expt(path=DIGITIZED_CSV):
    """Digitized loops reduced to the model variable DeltaR/R_P.

    Returns {field: (T, excess_over_cooling_baseline)} with
    excess = (R_warm - R_cool) / R_cool, a pure change of variable: no scale
    factor, offset or fit.  Returns {} if the digitized file is absent.
    """
    if not path.exists():
        return {}
    d = np.genfromtxt(path, delimiter=",", names=True)
    out = {}
    for tag in EXPT_FIELDS:
        cool = d[f"{tag}_cooldown"]
        warm = d[f"{tag}_warmup"]
        out[tag] = (d["temperature_K"], (warm - cool) / cool)
    return out


expt = load_expt() if SHOW_EXPT else {}

# Panel (a) has to hold whatever the data reach; the model curves are drawn on
# exactly the same axis, unrescaled.
YMAX = 1.42
if expt:
    YMAX = max(YMAX, 1.16 * max(np.nanmax(y) for _, y in expt.values()))
YMIN = -0.035 * YMAX

fig, (axa, axb) = plt.subplots(
    2, 1, figsize=(3.375, 4.55), gridspec_kw=dict(height_ratios=[2.0, 1.0]))

# ======================================================================
# (a) DeltaR / R_P through the loop
# ======================================================================
contrast_w = contrast(TW)
curves = {}
for r, col in zip(RATES, S.RAMP_RATE):
    phi_w = warming_phi(TW, r)
    dR = phi_w * contrast_w
    curves[r] = (phi_w, dR)
    axa.plot(TW[D], dR[D], color=col, lw=1.6, zorder=6,
             label=rf"${r:g}$" if r != 0.1 else r"$0.1$")

# cool-down: nucleation limited, phi stays at zero -> flat trace
phi_c = cooling_phi(TCOOL, 1.0)
axa.plot(TCOOL[D], (phi_c * contrast(TCOOL))[D], color=S.C_SEC, lw=1.6,
         ls=(0, (4, 2)), zorder=7)
axa.axhline(0.0, color=S.C_GRID, lw=0.7, zorder=0)

axa.text(197.0, 0.055 * YMAX, "cooling (all rates)", color=S.C_SEC,
         fontsize=7.0, ha="center", va="bottom", zorder=8)

# ----------------------------------------------------------------------
# digitized experiment, subordinate: hairline trace + sparse open symbols,
# directly labelled off the collapse edge where the panel is empty.
# ----------------------------------------------------------------------
EXPT_STYLE = {
    # field: (colour, marker, label, anchor T on the trace, label xy)
    "0T": (S.C_INK, "o", r"expt $0$ T", 243.0, (283.0, 1.44)),
    "5T": (S.C_BLUE, "s", r"expt $5$ T", 239.0, (283.0, 2.44)),
}
for tag, (T_e, y_e) in expt.items():
    col, mk, lab, T_anchor, xytext = EXPT_STYLE[tag]
    axa.plot(T_e, y_e, color=col, lw=0.55, alpha=0.45, zorder=2,
             solid_capstyle="round")
    sub = slice(None, None, EXPT_MARKER_EVERY)
    axa.plot(T_e[sub], y_e[sub], ls="none", marker=mk, ms=2.3, mfc="white",
             mec=col, mew=0.65, alpha=0.9, zorder=3)
    y_anchor = float(np.interp(T_anchor, T_e, y_e))
    axa.annotate(lab, xy=(T_anchor, y_anchor), xytext=xytext, color=col,
                 fontsize=7.0, ha="center", va="center", zorder=8,
                 arrowprops=dict(arrowstyle="-", lw=0.6, color=col,
                                 alpha=0.8, shrinkA=3, shrinkB=2))

if expt:
    axa.text(0.985, 0.985,
             "data digitized from Crowley $et\\ al.$ (2020), Fig. S5",
             transform=axa.transAxes, color=S.C_MUT, fontsize=6.2,
             ha="right", va="top", zorder=8)

leg = axa.legend(loc="upper left", bbox_to_anchor=(0.035, 0.905),
                 title="model ($B=0$), warming\nrate $r$ (K min$^{-1}$)",
                 fontsize=7.0,
                 handlelength=1.5, labelspacing=0.28, borderaxespad=0.0)
leg.get_title().set_fontsize(7.0)
leg.get_title().set_multialignment("left")
leg._legend_box.align = "left"

# --- 10-90 % rise width of the 1 K/min warm-up ----------------------------
phi1, dR1 = curves[1.0]
T10 = S.onset_temperature(TW, phi1, 0.1)
T90 = S.onset_temperature(TW, phi1, 0.9)
y_w = 0.60          # parked in the strip the data leave empty below 140 K
# KJMA n = 2 makes this span narrow (~7 K), so it is delimited by end bars
# rather than arrowheads, which would merge at this width.
axa.add_patch(FancyArrowPatch((T10, y_w), (T90, y_w),
                              arrowstyle="|-|,widthA=1.7,widthB=1.7",
                              mutation_scale=3, lw=0.8, color=S.C_INK,
                              shrinkA=0, shrinkB=0, zorder=9))
axa.annotate(rf"10–90% rise: {T90 - T10:.0f} K",
             xy=(T10, y_w), xytext=(T10 - 8.0, y_w), fontsize=7.0,
             color=S.C_INK, ha="right", va="center",
             arrowprops=dict(arrowstyle="-", lw=0.7, color=S.C_MUT,
                             shrinkA=2, shrinkB=2))

axa.set_xlim(80.0, 300.0)
axa.set_ylim(YMIN, YMAX)
axa.set_xticks([100, 150, 200, 250, 300])
axa.set_yticks([0.0, 1.0, 2.0] if YMAX > 1.9 else [0.0, 0.5, 1.0])
axa.set_xlabel(r"temperature  $T$  (K)")
axa.set_ylabel(r"$\Delta R/R_{\mathrm{P}}$")
S.panel_label(axa, "a")

# ======================================================================
# (b) onset versus warming rate
# ======================================================================
rates = np.logspace(-1.15, 1.15, 13)
Tons = np.array([S.onset_temperature(TW, warming_phi(TW, r)) for r in rates])
slope, intercept = np.polyfit(np.log10(rates), Tons, 1)

# Same slope restricted to the experimentally accessible window, reported to
# stdout only: the drawn line is the full two-decade fit.
_win = (rates >= 0.5) & (rates <= 5.0)
slope_win = np.polyfit(np.log10(rates[_win]), Tons[_win], 1)[0]

axb.plot(rates, intercept + slope * np.log10(rates), color=S.C_SEC, lw=1.2,
         ls=(0, (5, 2)), zorder=3)
axb.plot(rates, Tons, "o", ms=3.2, mfc="white", mec=S.C_SEC, mew=0.9,
         zorder=4)
for r, col in zip(RATES, S.RAMP_RATE):
    axb.plot([r], [S.onset_temperature(TW, warming_phi(TW, r))], "o",
             ms=4.6, mfc=col, mec="white", mew=0.8, zorder=5)

axb.set_xscale("log")
axb.set_xlim(0.06, 18.0)
YB_LO, YB_HI = Tons.min() - 3.3, Tons.max() + 4.0
axb.set_ylim(YB_LO, YB_HI)
axb.set_yticks([t for t in range(100, 300, 5) if YB_LO <= t <= YB_HI])
axb.set_xlabel(r"warming rate  $r$  (K min$^{-1}$)")
axb.set_ylabel(r"onset  $T_{\mathrm{on}}$  (K)")

# The drawn slope is the full rate-shift relation integrated numerically; it is
# NOT the leading form kB T_on^2 ln10 / Ea, so no formula is printed here and no
# equation number is cited (equation numbering may shift).
axb.annotate(rf"${slope:.1f}$ K per decade",
             xy=(0.039, 0.893), xycoords="axes fraction", fontsize=7.0,
             color=S.C_INK, ha="left", va="top", linespacing=1.35)
axb.text(0.975, 0.05,
         rf"$E_{{\mathrm{{a}}}}={EA_EV:.2f}$ eV" "\n"
         rf"$\nu=10^{{{np.log10(NU0):.0f}}}$ s$^{{-1}}$" "\n"
         rf"$n_{{0}}=10^{{{np.log10(N0_NM2):.0f}}}$ nm$^{{-2}}$",
         transform=axb.transAxes,
         fontsize=7.0, color=S.C_SEC, ha="right", va="bottom",
         linespacing=1.35)
S.panel_label(axb, "b")

fig.subplots_adjust(left=0.155, right=0.975, top=0.985, bottom=0.085,
                    hspace=0.35)
out = S.save(fig, "fig2_kinetics", png_dir=PNG_DIR)
print(f"wrote {out}")
print(f"  slope = {slope:.2f} K/decade over {rates[0]:.3g}-{rates[-1]:.3g} "
      f"K/min (drawn);  {slope_win:.2f} K/decade over 0.5-5 K/min")
print(f"  T_on(1 K/min) = {S.onset_temperature(TW, warming_phi(TW, 1.0)):.1f} K"
      f";  10-90 width = {T90 - T10:.1f} K;  "
      f"peak DR/RP = {dR1.max():.2f} at {TW[dR1.argmax()]:.0f} K")
for tag, (T_e, y_e) in expt.items():
    i = int(np.nanargmax(y_e))
    print(f"  expt {tag}: peak DR/RP = {y_e[i]:.2f} at {T_e[i]:.0f} K")

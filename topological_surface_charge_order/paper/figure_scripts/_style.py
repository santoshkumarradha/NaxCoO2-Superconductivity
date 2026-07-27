"""
House style for the kinetic two-phase transport paper (PRB single column).

Inherits the conventions of the group's previous paper
(../../../paper/figure_scripts/_style.py): STIX serif text with stix
mathtext, ticks in on all four sides, no grid, thin coloured spines, and the
validated colourblind-safe categorical palette.  Only two things are tightened
here relative to that file: axes.linewidth 0.7 -> 0.8 and the base font drops
to 8 pt so that every glyph in the figure set sits in the 7-9 pt band.

Everything drawn by the scripts that import this module is MODEL OUTPUT, with
one declared exception: fig2_kinetics.py overlays the digitized Crowley et al.
(2020) Fig. S5 resistance loops on panel (a), read from
results/digitized/crowley2020_figS5.csv and reduced to DeltaR/R_P by a change
of variable that is documented in that script's docstring.  Nothing anywhere in
this directory is FITTED to an experimental trace: no model parameter is
adjusted to improve agreement.  Published scalar values (onset near 146 K,
collapse near 240 K) enter only as model parameters and axis ranges.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

# ----------------------------------------------------------------------
# Paths.  Scripts live in paper/figure_scripts/, rendered vector PDFs go to
# paper/figures/ (what the .tex consumes).  Inspection PNGs never land there.
# ----------------------------------------------------------------------
HERE = Path(__file__).resolve().parent
FIGDIR = HERE.parent / "figures"

# ----------------------------------------------------------------------
# Physical constants (CODATA), in eV-based units so every energy in the
# figures is directly readable in meV.
# ----------------------------------------------------------------------
KB_EV = 8.617333262e-5          # eV / K
MU_B_EV = 5.7883818060e-5       # eV / T
MUB_OVER_KB = MU_B_EV / KB_EV   # 0.6717 K/T -- the Clausius-Clapeyron factor:
                                # a latent magnetization of 1 muB per Co against
                                # a latent entropy of 1 kB per Co tilts a
                                # first-order line by 0.6717 K/T.  (The inverse
                                # ratio KB_EV/MU_B_EV = 1.4887 is in T/K and must
                                # never be used as a K/T conversion.)

# ----------------------------------------------------------------------
# Palette: carried over verbatim from the previous paper's validated
# categorical set (worst adjacent CVD deltaE = 21.6).  Series identity is
# never carried by colour alone: every curve is either directly labelled or
# given a distinct dash pattern.
# ----------------------------------------------------------------------
C_BLUE = "#2a78d6"   # field / thermodynamic constructions
C_AQUA = "#1baf7a"   # third categorical series
C_YELL = "#eda100"   # predictions, slowest warming rate
C_RED = "#e34948"    # the ordered, resistive, magnetic phase; warm-up branch
C_INK = "#0b0b0b"
C_SEC = "#52514e"
C_MUT = "#898781"
C_GRID = "#e1e0d9"

# Colour grammar used by every figure here: the WARM-UP branch is drawn in
# warm ink (amber -> red -> maroon, an ordinal ramp for the sweep rate) and the
# COOL-DOWN branch in neutral grey.  Series identity is additionally carried by
# dash pattern and by direct labels, never by hue alone.
RAMP_RATE = ["#eda100", "#e34948", "#7c1f1e"]

# Region fills + matching darker label inks (previous paper's vocabulary).
FILL_NONE = "#ebebe8"    # neutral: untransformed / empty
FILL_RED = "#f5c4c4"     # the ordered, resistive phase O
FILL_AMBER = "#f8e2b4"   # prediction / untested window
FILL_BLUE = "#cde2fb"    # parent metal P

INK_NONE = C_SEC
INK_RED = "#7c1f1e"
INK_AMBER = "#6d4c00"
INK_BLUE = "#184f95"

_FILL_INK = {FILL_NONE: INK_NONE, FILL_RED: INK_RED,
             FILL_AMBER: INK_AMBER, FILL_BLUE: INK_BLUE}


def use_house_style() -> None:
    """STIX serif, ticks in on all sides, hairline spines, embedded fonts."""
    mpl.rcParams.update({
        "font.family": "serif",
        "font.serif": ["STIXGeneral", "STIX Two Text", "Times New Roman",
                       "Times", "DejaVu Serif"],
        "mathtext.fontset": "stix",
        "font.size": 8.0,
        "axes.titlesize": 8.5,
        "axes.labelsize": 8.0,
        "xtick.labelsize": 7.5,
        "ytick.labelsize": 7.5,
        "legend.fontsize": 7.0,
        "axes.linewidth": 0.8,
        "axes.grid": False,
        "axes.axisbelow": True,
        "lines.linewidth": 1.4,
        "lines.solid_capstyle": "round",
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "xtick.minor.width": 0.55,
        "ytick.minor.width": 0.55,
        "xtick.major.size": 3.0,
        "ytick.major.size": 3.0,
        "xtick.minor.size": 1.7,
        "ytick.minor.size": 1.7,
        "xtick.top": True,
        "ytick.right": True,
        "legend.frameon": False,
        "legend.handlelength": 1.5,
        "legend.handletextpad": 0.5,
        "legend.labelspacing": 0.3,
        "legend.borderpad": 0.2,
        "axes.edgecolor": C_SEC,
        "text.color": C_INK,
        "axes.labelcolor": C_INK,
        "xtick.color": C_SEC,
        "ytick.color": C_SEC,
        "figure.dpi": 150,
        "savefig.dpi": 400,
        "pdf.fonttype": 42,     # editable vector text
        "ps.fonttype": 42,
    })


def panel_label(ax, letter, x=0.012, y=0.975, **kw):
    """Bold (a)/(b)/(c) in the panel's own top-left corner, same place always."""
    ax.text(x, y, f"({letter})", transform=ax.transAxes, fontsize=9,
            fontweight="bold", va="top", ha="left", color=C_INK,
            zorder=20, **kw)


def outer_panel_label(ax, letter, x=-0.20, y=1.04, **kw):
    """Panel letter parked outside the axes, for panels whose corner is busy."""
    ax.text(x, y, f"({letter})", transform=ax.transAxes, fontsize=9,
            fontweight="bold", va="bottom", ha="left", color=C_INK,
            zorder=20, **kw)


def save(fig, stem, png_dir=None):
    """Vector PDF into paper/figures/; optional inspection PNG elsewhere."""
    FIGDIR.mkdir(parents=True, exist_ok=True)
    out = FIGDIR / f"{stem}.pdf"
    fig.savefig(out)
    if png_dir is not None:
        png_dir = Path(png_dir)
        png_dir.mkdir(parents=True, exist_ok=True)
        fig.savefig(png_dir / f"{stem}.png", dpi=400)
    plt.close(fig)
    return out


# ======================================================================
# Shared model definitions.  Every figure in this directory draws on these,
# so fig1's Landau functional and fig2's kinetics cannot drift apart.
# ======================================================================

# Published scalars used only as model parameters / axis anchors.
T_ON_REF = 146.0     # K, warm-up unfreezing midpoint at r = 1 K/min
TC0 = 240.0          # K, zero-field first-order collapse of the ordered phase
EA_EV = 0.45         # eV, effective barrier (Li migration class)
NU0 = 1.0e13         # 1/s, attempt frequency


def landau_f(phi, T, T0=60.0, Tc=TC0, scale_meV=100.0):
    """Landau free energy f(phi) = a(T) phi^2/2 - b phi^3/3 + c phi^4/4.

    Reduced coefficients b = c = 1; a(T) = a_deg (T - T0)/(Tc - T0) with
    a_deg = 2 b^2 / (9 c), the exact degeneracy value, so that f(phi_O) = f(0)
    at T = Tc by construction: the cubic invariant makes the transition first
    order and the two minima are degenerate on the transition line.  Returned
    in meV per Co after multiplication by `scale_meV`.
    """
    a = (2.0 / 9.0) * (np.asarray(T, float) - T0) / (Tc - T0)
    phi = np.asarray(phi, float)
    return scale_meV * (a * phi**2 / 2.0 - phi**3 / 3.0 + phi**4 / 4.0)


def landau_a(T, T0=60.0, Tc=TC0):
    return (2.0 / 9.0) * (np.asarray(T, float) - T0) / (Tc - T0)


def landau_extrema(T, T0=60.0, Tc=TC0):
    """(phi_barrier, phi_ordered) of f, or (nan, nan) if only phi = 0 exists.

    Stationary points of f solve a - b phi + c phi^2 = 0.
    """
    a = float(landau_a(T, T0, Tc))
    disc = 1.0 - 4.0 * a
    if disc <= 0.0:
        return np.nan, np.nan
    r = np.sqrt(disc)
    return 0.5 * (1.0 - r), 0.5 * (1.0 + r)


def k_arrhenius(T, Ea=EA_EV, nu=NU0):
    """Growth-limited rate k(T) = nu exp(-Ea / kB T)."""
    return nu * np.exp(-Ea / (KB_EV * np.asarray(T, float)))


def nucleation_barrier(T, A_eVK=16.0, Tc=TC0):
    """2D KJMA nucleation barrier W*(T) = pi sigma^2 / |Delta f| = A/(Tc - T).

    Diverges at Tc (no undercooling, no driving force) and falls with
    undercooling.  A is the single lumped line-tension parameter.
    """
    T = np.asarray(T, float)
    return np.where(T < Tc - 1e-9, A_eVK / np.maximum(Tc - T, 1e-9), np.inf)


def phi_eq(T, Tc=TC0, width=6.0):
    """Equilibrium ordered fraction: 1 below Tc, 0 above, first-order collapse
    smeared over `width` by the compositional spread in local Tc."""
    return 1.0 / (1.0 + np.exp((np.asarray(T, float) - Tc) / width))


def sweep(T_grid, rate_K_per_min, phi0, rate_fn, eq_fn=phi_eq):
    """Integrate dphi/dt = k(T)[phi_eq(T) - phi] along a linear ramp.

    T_grid is the ordered list of temperatures actually visited (ascending for
    a warm-up, descending for a cool-down); `rate_K_per_min` is the magnitude
    of |dT/dt|.  Because the master equation is linear in phi, each step is
    solved exactly for frozen coefficients,

        phi_{i+1} = phi_eq + (phi_i - phi_eq) exp(-k dt),

    which is unconditionally stable: no stiff solver is needed even where
    k ~ 1e9 / s.  Returns phi on T_grid.
    """
    T_grid = np.asarray(T_grid, float)
    r_s = rate_K_per_min / 60.0                       # K / s
    dt = np.abs(np.diff(T_grid)) / r_s                # s per step
    Tmid = 0.5 * (T_grid[:-1] + T_grid[1:])
    k = np.asarray(rate_fn(Tmid), float)
    pe = np.asarray(eq_fn(Tmid), float)
    phi = np.empty_like(T_grid)
    phi[0] = phi0
    decay = np.exp(-np.clip(k * dt, 0.0, 700.0))
    for i in range(T_grid.size - 1):
        phi[i + 1] = pe[i] + (phi[i] - pe[i]) * decay[i]
    return phi


def warming_phi(T_grid, rate_K_per_min, Ea=EA_EV, nu=NU0, phi0=0.0):
    """Warm-up branch: growth of quenched-in nuclei, barrier Ea only.

    The quenched state carries subcritical/critical nuclei frozen in at
    maximal undercooling, so the warm-up is growth-limited (W* << Ea there)
    and the effective single-barrier master equation applies.
    """
    return sweep(T_grid, rate_K_per_min, phi0,
                 lambda T: k_arrhenius(T, Ea, nu))


def cooling_phi(T_grid, rate_K_per_min, Ea=EA_EV, nu=NU0, A_eVK=16.0,
                phi0=0.0):
    """Cool-down branch: nucleation-limited, barrier Ea + W*(T).

    Near Tc the driving force vanishes and W* diverges, so nothing nucleates
    while the ions are still mobile; by the time undercooling has cut W* down,
    exp(-Ea/kB T) has frozen the transport of Li.  phi therefore stays at 0
    over the whole cool-down: this is a calculation, not an assumption.
    """
    def k_nuc(T):
        W = nucleation_barrier(T, A_eVK)
        return nu * np.exp(-np.minimum((Ea + W) / (KB_EV * T), 700.0))
    return sweep(T_grid, rate_K_per_min, phi0, k_nuc)


def onset_temperature(T_grid, phi, level=0.5):
    """Rise midpoint: first crossing of phi = `level` on a warm-up branch."""
    i = int(np.argmax(phi >= level))
    if i == 0:
        return np.nan
    t0, t1, p0, p1 = T_grid[i - 1], T_grid[i], phi[i - 1], phi[i]
    return t0 + (level - p0) * (t1 - t0) / (p1 - p0)


def contrast(T, T_ref=160.0, E_s_eV=0.036, amp=1.0):
    """Geometry-weighted transport contrast c [s(T) - 1] of the ordered phase.

    A gapped phase measured against an increasingly conductive parent loses
    contrast on warming; parametrized by a single effective scale E_s and
    normalized to `amp` at T_ref so that DeltaR/R_P is of order unity where
    phi first saturates.

    E_s_eV = 0.036: illustrative shallow decay slope (model choice).
    """
    T = np.asarray(T, float)
    return amp * np.exp(E_s_eV / (KB_EV * T) - E_s_eV / (KB_EV * T_ref))

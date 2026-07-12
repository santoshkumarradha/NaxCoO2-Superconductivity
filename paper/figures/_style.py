"""
Shared house style + data loaders for the Na_xCoO2.yH2O publication figures.

Every figure (fig1..fig5) imports from here so that the series -> color mapping,
rcParams, and data provenance are identical across the whole figure set.

Palette: dataviz reference palette (validated, colorblind-safe).
  worst adjacent CVD deltaE for the 3-series categorical set = 21.6 (>= 12 target).
  aqua/yellow fall below the 3:1 contrast floor -> the "relief rule" is honoured
  everywhere by (i) a distinct marker per series and (ii) a legend / direct labels,
  so series identity is never carried by colour alone.

Nothing in theory/, runpod/, spin_analysis/ is modified: we only READ from them.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

# ----------------------------------------------------------------------
# Paths (paper/figures/_style.py -> repo root is parents[2])
# ----------------------------------------------------------------------
HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
V4 = REPO / "runpod" / "results_v4"
THEORY_RES = REPO / "theory" / "results"
SPIN = REPO / "spin_analysis"

RY_EV = 13.605693122994  # CODATA Ry -> eV

# ----------------------------------------------------------------------
# Palette (dataviz reference instance, light surface #fcfcfb)
# ----------------------------------------------------------------------
C_BLUE = "#2a78d6"   # Na 1x1  (x=1) : the mechanism series
C_AQUA = "#1baf7a"   # Li 1x1  (x=1)
C_YELL = "#eda100"   # Na sqrt3 x sqrt3 (x=1/3)
C_RED = "#e34948"    # RESERVED: superconducting experimental anchor (9.9 A)
C_INK = "#0b0b0b"    # primary ink
C_SEC = "#52514e"    # secondary ink / neutral (grey) anchor
C_MUT = "#898781"    # muted axis/label
C_GRID = "#e1e0d9"   # hairline
C_GOOD = "#0ca30c"   # status: explained
C_WARN = "#eda100"   # status: prediction / untested  (amber)
C_CRIT = "#d03b3b"   # status: not explained

# Fixed series identity -> (colour, marker, label). Used in EVERY figure.
SERIES = {
    ("Na", "1x1"): dict(color=C_BLUE, marker="o",
                        label=r"Na $1{\times}1$  ($x{=}1$)"),
    ("Li", "1x1"): dict(color=C_AQUA, marker="s",
                        label=r"Li $1{\times}1$  ($x{=}1$)"),
    ("Na", "s3"):  dict(color=C_YELL, marker="^",
                        label=r"Na $\sqrt{3}{\times}\sqrt{3}$  ($x{=}1/3$)"),
}

# Ordinal ramps for the four c-values within a single fig-1 panel (magnitude = c).
# One hue per series so a panel's overall colour still reads as that series.
RAMP_BLUE = ["#9ec5f4", "#3987e5", "#1c5cab", "#0d366b"]   # Na 1x1 panel
RAMP_YELL = ["#f6d489", "#eda100", "#b97a00", "#7a4f00"]   # Na s3  panel

# Experimental anchors (spacing -> (label, colour, is_SC))
ANCHORS = [
    (5.5, "anhydrous\n(not SC)", C_SEC, False),
    (6.9, "1-layer hydrate\n(not SC)", C_SEC, False),
    (9.9, "2-layer hydrate\n(SC 4.5 K)", C_RED, True),
]


def use_house_style():
    """Serif / STIX text, mathtext for symbols, thin spines, ticks-in, no grid."""
    mpl.rcParams.update({
        "font.family": "serif",
        "font.serif": ["STIXGeneral", "STIX Two Text", "Times New Roman",
                       "Times", "DejaVu Serif"],
        "mathtext.fontset": "stix",
        "font.size": 8.5,
        "axes.titlesize": 9,
        "axes.labelsize": 8.5,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 7.2,
        "axes.linewidth": 0.7,
        "axes.grid": False,
        "axes.axisbelow": True,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.major.width": 0.7,
        "ytick.major.width": 0.7,
        "xtick.minor.width": 0.5,
        "ytick.minor.width": 0.5,
        "xtick.major.size": 3.2,
        "ytick.major.size": 3.2,
        "xtick.minor.size": 1.8,
        "ytick.minor.size": 1.8,
        "xtick.top": True,
        "ytick.right": True,
        "legend.frameon": False,
        "axes.edgecolor": C_SEC,
        "text.color": C_INK,
        "axes.labelcolor": C_INK,
        "xtick.color": C_SEC,
        "ytick.color": C_SEC,
        "figure.dpi": 150,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "pdf.fonttype": 42,   # embed as TrueType (editable vector text)
        "ps.fonttype": 42,
    })


def panel_label(ax, letter, x=-0.02, y=1.02, **kw):
    """Bold (a)/(b) panel letter, top-left, in axes coords."""
    ax.text(x, y, f"({letter})", transform=ax.transAxes, fontsize=10,
            fontweight="bold", va="bottom", ha="right", color=C_INK, **kw)


def thin_spines(ax, hide=("top", "right")):
    for s in hide:
        ax.spines[s].set_visible(False)


def save(fig, stem):
    """Vector PDF + 300 dpi PNG, both into paper/figures/."""
    fig.savefig(HERE / f"{stem}.pdf")
    fig.savefig(HERE / f"{stem}.png", dpi=300)
    plt.close(fig)


# ======================================================================
# Data loaders (all read-only)
# ======================================================================
def load_results():
    """results.json -> DataFrame keyed by (element, cell, c)."""
    rows = json.loads((V4 / "results.json").read_text())
    return pd.DataFrame(rows)


def load_well_fits():
    """Canonical quartic (alpha,beta) + sextic refit (a_q,b_q,g_q) + depth."""
    return pd.read_csv(THEORY_RES / "well_fits_v4.csv")


def load_quantum():
    """Quantized levels E0/E1/E2, hw_eff, d_rms per (element,cell,c,mass)."""
    return pd.read_csv(THEORY_RES / "quantum_wells_v4.csv")


def load_coupling():
    """lambda(c), Tc(c) dome + gamma-sensitivity band (Tc_K_g1..g2)."""
    return pd.read_csv(THEORY_RES / "coupling_tc_vs_c.csv")


def load_magnetization():
    """|m|(element,c,delta) from the LSDA spin scans."""
    return pd.read_csv(SPIN / "magnetization_v2.csv")


_E_RE = re.compile(r"^!\s+total energy\s+=\s+(-?[\d.]+)\s+Ry", re.M)
_JOB_RE = re.compile(r"^([A-Za-z]+)_(1x1|s3)_c([0-9.]+)_d([0-9.]+)$")


def load_scans():
    """Re-parse raw E(delta) totals from runpod jobs/*/pw.out.

    Returns {(element, cell, c): (delta_array, E_meV_array)} with energies
    referenced to the centred (delta=0) configuration, in meV.
    """
    groups: dict[tuple, list[tuple[float, float]]] = {}
    for d in sorted((V4 / "jobs").iterdir()):
        m = _JOB_RE.match(d.name)
        if not m:
            continue
        el, cell, c, delta = m.group(1), m.group(2), float(m.group(3)), float(m.group(4))
        pw = d / "pw.out"
        if not pw.exists():
            continue
        txt = pw.read_text()
        if "JOB DONE" not in txt:
            continue
        hits = _E_RE.findall(txt)
        if not hits:
            continue
        E = float(hits[-1]) * RY_EV
        groups.setdefault((el, cell, c), []).append((delta, E))
    scans = {}
    for key, pts in groups.items():
        pts.sort()
        dd = np.array([p[0] for p in pts])
        EE = np.array([p[1] for p in pts])
        EE = (EE - EE[dd.argmin()]) * 1e3   # reference to delta=0, -> meV
        scans[key] = (dd, EE)
    return scans


def vpoly(x, a, b, g=0.0):
    """Well potential V(delta) = a d^2 + b d^4 + g d^6  (eV, delta in A)."""
    return a * x**2 + b * x**4 + g * x**6

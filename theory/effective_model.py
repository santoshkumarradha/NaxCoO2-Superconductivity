# /// script
# requires-python = ">=3.10"
# dependencies = ["numpy", "scipy", "matplotlib", "pandas"]
# ///
"""
Effective model for gallery-controlled superconductivity in Na_xCoO2·yH2O
=========================================================================

End-to-end pipeline from the canonical DFT dataset in this repository
(runpod/results_v4: Quantum ESPRESSO, PBE+PAW, spin-polarized, Co-top
alkali site) to the "money plot".  Run with

    uv run theory/effective_model.py          (or)
    theory/.venv/bin/python theory/effective_model.py

Pipeline (default, results_v4)
------------------------------
1.  Load runpod/results_v4/{results.json, manifest.json} and re-parse the
    raw E(delta) scans from jobs/*/pw.out (energies in Ry).  Series:
    Li 1x1 (x=1), Na 1x1 (x=1), Na sqrt3 x sqrt3 (x=1/3), CoO2-CoO2
    spacings c = 5.5 / 6.9 / 8.4 / 9.9 A.
2.  Landau alpha(c), beta(c): the canonical quartic fits from
    results.json.  For quantization, wells whose quartic fit fails
    (beta < 0 at large c) or whose minimum is not bracketed by the
    sampled delta range are refit from the pw.out energies with
    alpha d^2 + beta d^4 + gamma6 d^6 (gamma6 > 0 constrained,
    minimum capped at delta = c/2 - 2.0 A: hard-core alkali-O repulsion).
3.  Interpolate alpha(c) (PCHIP); report c*_Na (bracketed between 5.5 and
    6.9 A) and c*_Li (upper bound <= 5.5 A only: alpha < 0 already at the
    smallest spacing in this dataset).
4.  Exact 1D Schroedinger quantization of each well (finite differences,
    parity-resolved) with the physical alkali mass.  For deep wells the
    reported frequency is ALWAYS hbar*omega_eff = hbar^2/(2M<d^2>_0),
    never E1-E0 (which is an exponentially small tunnel splitting there).
5.  Electron-phonon vertex: the SSH4 tight-binding model of Radha &
    Lambrecht, SciPost Phys. 10, 057 (2021), tau2 -> tau2(1+gamma d),
    tau4 -> tau4(1-gamma d).  The linear deformation potential vanishes
    at d = 0 by mirror parity (verified numerically); the leading vertex
    is quadratic, E_band(d) = E_edge + (1/2) chi d^2 with chi from exact
    second-order perturbation theory.  lambda = 2 N(0) g^2 / (E2-E0),
    g = (1/2)|chi| <0|d^2|2>, with the DFT N(0) from results_v4 dos.x
    runs and the on-site splitting delta(c) calibrated on the results_v4
    Bader charges.  Allen-Dynes Tc, mu* = 0.10.  Where lambda > 2 the
    Migdal treatment is invalid: reported as "SC killed (polaronic)",
    never as a divergent Tc.
6.  Figures: diagnostics in theory/figures/, the money plot in
    theory/money_plot.{pdf,png}; numeric tables in theory/results/.

Legacy mode
-----------
    effective_model.py --legacy
re-runs the original proof-of-concept analysis on the old relaxed-scan
data (phase_transition/{licoo2,nacoo2}/structures_done.txt, Questaal lmf)
and prints the old-model summary (c*, well depths, Bader onset) used for
the old-vs-new comparison in MODEL.md.  Legacy tables go to
theory/results/legacy/.

Every number is traceable to (a) the results_v4 DFT data, (b) the SciPost
paper's fitted parameters, or (c) an explicitly labeled ansatz with a
stated range.  See theory/MODEL.md.
"""

from __future__ import annotations

import argparse
import json
import pickle
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator
from scipy.linalg import eigh_tridiagonal
from scipy.optimize import brentq, curve_fit

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

# ----------------------------------------------------------------------
# Paths
# ----------------------------------------------------------------------
HERE = Path(__file__).resolve().parent          # theory/
REPO = HERE.parent
V4 = REPO / "runpod" / "results_v4"             # canonical DFT dataset
PT = REPO / "phase_transition"                  # legacy proof-of-concept data
FIGDIR = HERE / "figures"
RESDIR = HERE / "results"
FIGDIR.mkdir(exist_ok=True)
RESDIR.mkdir(exist_ok=True)

# ----------------------------------------------------------------------
# Physical constants (CODATA)
# ----------------------------------------------------------------------
RY_EV = 13.605693122994          # Ry -> eV
HBAR2_OVER_2AMU = 2.0901e-3      # hbar^2 / (2 * 1 amu) in eV * Angstrom^2
KB_EV = 8.617333e-5              # eV / K
MASS_AMU = {"Li": 6.941, "Na": 22.98977}

# ----------------------------------------------------------------------
# Model parameters
# ----------------------------------------------------------------------
# --- SciPost Phys. 10, 057 (2021), Radha & Lambrecht (fitted TB values) ---
TAU1 = 0.5      # eV  t_A^z         (on-site alkali sp_z pair hopping)
TAU3 = 0.5      # eV  t_CoO2^z      (hopping across the CoO2 slab)
TAU24 = 2.0     # eV  t_{A-CoO2}^z  (alkali <-> CoO2 interlayer hopping)
T_LI_XY = -0.6  # eV  in-plane alkali hopping
T_CO_XY = 0.09  # eV  in-plane CoO2 hopping
GAMMA_BANDS = 6 * abs(T_LI_XY) + 6 * abs(T_CO_XY)   # = 4.14 eV, band-overlap Gamma

# --- Explicit ansatz (stated range; see MODEL.md) ---
Z0_ANSATZ = 0.70          # Angstrom, hopping decay length t(z) ~ exp(-z/z0)
Z0_RANGE = (0.5, 1.0)     # sensitivity range -> gamma_ep = 1/z0 in [1, 2] 1/A
GAMMA_EP = 1.0 / Z0_ANSATZ
MU_STAR = 0.10            # Coulomb pseudopotential (task prescription)
LAMBDA_LOC = 2.0          # polaronic-localization criterion, ansatz range [1, 2]
Q_TOT = 1.0               # e/cell: full conversion of the alkali band (paper limit)
DELTA_FLOOR = 0.05        # eV, keeps the TB spectrum non-degenerate
N0_GATE = 0.10            # states/eV/cell/spin: "2DEG exists" gate (DFT N(0))
D_CONTACT = 2.0           # A, min alkali-to-O-plane distance (hard-core ansatz)

ANCHOR_SPACINGS = (5.5, 6.9, 9.9)   # anhydrous / monolayer / bilayer hydrate
TC_EXP = 4.5                        # K at 9.9 A (Takada et al. 2003)

# ----------------------------------------------------------------------
# Dataviz palette (validated reference palette from the dataviz skill)
# ----------------------------------------------------------------------
C_BLUE = "#2a78d6"     # series 1: Na 1x1 (primary)
C_AQUA = "#1baf7a"     # series 2: Li 1x1
C_YELL = "#eda100"     # series 3: Na sqrt3 (x=1/3)
C_VIOLET = "#4a3aa7"   # series 5
C_RED = "#e34948"      # reserved: SC anchor
C_GRAY = "#52514e"     # secondary text / neutral anchor
C_MUTED = "#9a9992"
CMAP_BLUE = LinearSegmentedColormap.from_list(
    "seq_blue", ["#cde2fb", "#5598e7", "#1c5cab", "#0d366b"])

plt.rcParams.update({
    "font.size": 9.5,
    "font.family": "sans-serif",
    "axes.linewidth": 0.8,
    "axes.grid": True,
    "grid.color": "#e6e5e0",
    "grid.linewidth": 0.6,
    "axes.axisbelow": True,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "legend.frameon": False,
    "figure.dpi": 150,
    "savefig.dpi": 300,
})

SERIES_STYLE = {  # (color, marker, label)
    ("Na", "1x1"): (C_BLUE, "o", "Na 1x1 (x=1)"),
    ("Li", "1x1"): (C_AQUA, "s", "Li 1x1 (x=1)"),
    ("Na", "s3"): (C_YELL, "^", r"Na $\sqrt{3}{\times}\sqrt{3}$ (x=1/3)"),
}


# ======================================================================
# 1. DATA EXTRACTION (results_v4)
# ======================================================================
E_RE = re.compile(r"^!\s+total energy\s+=\s+(-?[\d.]+)\s+Ry", re.M)


@dataclass
class Scan:
    """One constant-c scan: total energy vs alkali off-center displacement."""
    element: str
    cell: str                # "1x1" or "s3"
    c: float                 # CoO2-CoO2 spacing (A)
    d: np.ndarray            # off-center displacement delta >= 0 (A)
    E: np.ndarray            # energy relative to centered position (eV)


def load_v4() -> tuple[dict, list[Scan], dict]:
    """results.json records keyed by (element, cell, c), raw E(delta) scans
    re-parsed from jobs/*/pw.out, and the manifest (a_lat, zval, ...)."""
    manifest = json.loads((V4 / "manifest.json").read_text())
    records = {(r["element"], r["cell"], r["c"]): r
               for r in json.loads((V4 / "results.json").read_text())}
    groups: dict[tuple, list[tuple[float, float]]] = {}
    for j in manifest["jobs"]:
        if j["type"] != "scf":
            continue
        txt = (V4 / "jobs" / j["name"] / "pw.out").read_text()
        if "JOB DONE" not in txt:
            continue
        e = float(E_RE.findall(txt)[-1]) * RY_EV
        groups.setdefault((j["element"], j["cell"], j["c"]), []).append(
            (j["delta"], e))
    scans = []
    for key in sorted(groups):
        pts = sorted(groups[key])
        d = np.array([p[0] for p in pts])
        E = np.array([p[1] for p in pts])
        E -= E[d.argmin()]                       # reference: centered point
        scans.append(Scan(element=key[0], cell=key[1], c=key[2], d=d, E=E))
    return records, scans, manifest


# ======================================================================
# 2. WELL POTENTIALS
#    canonical Landau alpha/beta = quartic fits from results.json;
#    quantization potential = sextic refit where the quartic fails
# ======================================================================
def vpoly(x, a, b, g=0.0):
    return a * x**2 + b * x**4 + g * x**6


def well_minimum(a, b, g, dmax=5.0):
    dd = np.linspace(0.0, dmax, 5001)
    V = vpoly(dd, a, b, g)
    i = int(np.argmin(V))
    return float(dd[i]), float(-V[i])            # d0, depth (>0 if double well)


def refit_sextic(d, E, c) -> tuple[float, float, float, float]:
    """Constrained refit E(d) = E0 + a d^2 + b d^4 + g d^6 used only for
    quantization.  g is scanned (>= 0); (E0, a, b) solved by linear least
    squares at each g.  Constraints: the potential must be confining
    (V rises again beyond the well) and the minimum must respect the
    hard-core cap d0 <= c/2 - D_CONTACT (the alkali cannot penetrate the
    O plane; ansatz, stated in MODEL.md).  Returns (a, b, g, rms_eV)."""
    cap = c / 2.0 - D_CONTACT
    best = None
    for g in np.concatenate([[0.0], np.logspace(-3, 1.6, 120)]):
        A = np.stack([np.ones_like(d), d**2, d**4], axis=1)
        coef, *_ = np.linalg.lstsq(A, E - g * d**6, rcond=None)
        e0, a, b = coef
        dd = np.linspace(0.0, 5.0, 2001)
        V = vpoly(dd, a, b, g)
        if V[-1] < 0.5:                           # not confining
            continue
        if dd[int(np.argmin(V))] > cap + 1e-9:    # min beyond hard-core cap
            continue
        rms = float(np.sqrt(np.mean((e0 + vpoly(d, a, b, g) - E) ** 2)))
        if best is None or rms < best[3]:
            best = (float(a), float(b), float(g), rms)
    if best is None:
        raise RuntimeError(f"no admissible sextic fit at c = {c}")
    return best


def build_well_table(records: dict, scans: list[Scan]) -> pd.DataFrame:
    """One row per (element, cell, c): canonical alpha/beta (results.json)
    plus the quantization potential (a_q, b_q, g_q) and its provenance."""
    rows = []
    for s in scans:
        rec = records[(s.element, s.cell, s.c)]
        a_can, b_can = rec["alpha"], rec["beta"]
        dmax = float(s.d.max())
        d0_can, _ = well_minimum(a_can, b_can, 0.0)
        bracketed = int(np.argmin(s.E)) < len(s.E) - 1   # data min interior
        quartic_ok = (b_can > 0) and (d0_can <= dmax) and bracketed
        if quartic_ok:
            a_q, b_q, g_q = a_can, b_can, 0.0
            rms = float(np.sqrt(np.mean((vpoly(s.d, a_q, b_q) - s.E) ** 2)))
            refit = False
        else:
            a_q, b_q, g_q, rms = refit_sextic(s.d, s.E, s.c)
            refit = True
        d0, depth = well_minimum(a_q, b_q, g_q)
        rows.append({
            "element": s.element, "cell": s.cell, "c": s.c,
            "alpha": a_can, "beta": b_can,           # canonical (results.json)
            "a_q": a_q, "b_q": b_q, "g_q": g_q,      # quantization potential
            "refit_sextic": refit, "rms_meV": rms * 1e3,
            "d0_A": d0, "depth_eV": depth,
            "min_bracketed": bracketed,              # False -> depth is a
            "E_edge_meV": float(s.E[-1]) * 1e3,      # lower bound >= |E_edge|
            "d_max": dmax, "n_pts": len(s.d),
        })
    return pd.DataFrame(rows)


def critical_spacing(cs: np.ndarray, alpha: np.ndarray):
    """c* where PCHIP-interpolated alpha(c) crosses zero.  Returns
    (c*_pchip, c*_linear, (c_lo, c_hi) DFT bracket) or None if alpha < 0
    already at the smallest spacing (only an upper bound then)."""
    if alpha[0] <= 0:
        return None
    i = int(np.argmax(alpha < 0))                 # first negative
    p = PchipInterpolator(cs, alpha)
    c_pchip = float(brentq(p, cs[i - 1], cs[i]))
    c_lin = float(cs[i - 1] + (cs[i] - cs[i - 1])
                  * alpha[i - 1] / (alpha[i - 1] - alpha[i]))
    return c_pchip, c_lin, (float(cs[i - 1]), float(cs[i]))


# ======================================================================
# 3. QUANTUM WELL: parity-resolved 1D Schroedinger solver (exact levels)
# ======================================================================
@dataclass
class WellQuantum:
    E_even: np.ndarray       # eigenvalues of even-parity sector (eV)
    E_odd: np.ndarray        # eigenvalues of odd-parity sector (eV)
    d2_00: float             # <0|d^2|0>   (A^2)
    d2_02: float             # <0_e|d^2|1_e>  matrix element (A^2)
    omega_eff: float         # hbar*omega_eff = hbar^2/(2M<d^2>_0) (eV)
    split_01: float          # tunneling splitting E1-E0 (eV)
    gap_02: float            # even-sector gap E2-E0 (eV)


def solve_well(a: float, b: float, g: float, mass_amu: float,
               N: int = 2400, nstates: int = 4) -> WellQuantum:
    """Parity-resolved finite differences on the half line (staggered grid
    x_j = (j+1/2) h).  Even sector: psi'(0)=0 -> H_00 = t + V0; odd sector:
    psi(0)=0 -> H_00 = 3t + V0.  Separates the tunneling doublet exactly.
    The box size follows the well: L = max(2.6, 2*d0 + 1.5) A."""
    d0, _ = well_minimum(a, b, g)
    L = max(2.6, 2.0 * d0 + 1.5)
    h = L / N
    x = (np.arange(N) + 0.5) * h
    V = vpoly(x, a, b, g)
    t = HBAR2_OVER_2AMU / mass_amu / h**2
    off = -t * np.ones(N - 1)

    def solve(parity: str):
        dg = 2 * t + V
        dg[0] = (t if parity == "even" else 3 * t) + V[0]
        w, v = eigh_tridiagonal(dg, off, select="i",
                                select_range=(0, nstates - 1))
        v = v / np.sqrt(2 * h * np.sum(v**2, axis=0))    # full-line norm
        return w, v

    we, ve = solve("even")
    wo, _ = solve("odd")
    d2_00 = 2 * h * float(np.sum(x**2 * ve[:, 0] ** 2))
    d2_02 = 2 * h * float(np.sum(x**2 * ve[:, 0] * ve[:, 1]))
    return WellQuantum(
        E_even=we, E_odd=wo,
        d2_00=d2_00, d2_02=abs(d2_02),
        omega_eff=(HBAR2_OVER_2AMU / mass_amu) / d2_00,
        split_01=float(wo[0] - we[0]),
        gap_02=float(we[1] - we[0]),
    )


def quantum_table(wells: pd.DataFrame) -> pd.DataFrame:
    """Quantize every well with its physical alkali mass (plus the Li mass
    in the Na wells as an isotope-style cross-check)."""
    rows = []
    for _, r in wells.iterrows():
        for mlab in ({r["element"], "Li"} if r["element"] == "Na"
                     else {r["element"]}):
            q = solve_well(r["a_q"], r["b_q"], r["g_q"], MASS_AMU[mlab])
            rows.append({
                "element": r["element"], "cell": r["cell"], "c": r["c"],
                "mass": mlab,
                "E0_meV": q.E_even[0] * 1e3, "E1_meV": q.E_odd[0] * 1e3,
                "E2_meV": q.E_even[1] * 1e3,
                "split01_meV": q.split_01 * 1e3,
                "gap02_meV": q.gap_02 * 1e3,
                "d2_A2": q.d2_00, "d_rms_A": np.sqrt(q.d2_00),
                "d2_02_A2": q.d2_02,
                "hw_eff_meV": q.omega_eff * 1e3,
            })
    return pd.DataFrame(rows)


# ======================================================================
# 4. CARRIERS: DFT N(0) + Bader-returned gallery charge
# ======================================================================
CS_DFT = (5.5, 6.9, 8.4, 9.9)


@dataclass
class Carriers:
    cs: np.ndarray           # DFT spacings
    N0: dict                 # (el, cell) -> N(0) per spin per cell (1/eV)
    dq: dict                 # (el, cell) -> Bader charge RETURNED to alkali (e)
    n2d: dict                # (el, cell) -> dq / A_cell (e/cm^2)
    area: dict               # (el, cell) -> cell area (A^2)


def build_carriers(records: dict, manifest: dict) -> Carriers:
    """N(0) from the dos.x runs (per spin, per cell) and the Bader charge
    *returned* to the alkali relative to the closed gallery at c = 5.5 A:
    dq(c) = q_bader(c) - q_bader(5.5).  Rising dq = electrons coming back
    into the gallery = the 2DEG forming (the donated charge Z_val - q_bader
    falls with c)."""
    cs = np.array(CS_DFT)
    N0, dq, n2d, area = {}, {}, {}, {}
    for el, cell in SERIES_STYLE:
        a = manifest["a_lat"][el] * (np.sqrt(3) if cell == "s3" else 1.0)
        A = np.sqrt(3) / 2 * a**2
        n0 = np.array([records[(el, cell, c)]["N0"] for c in cs])
        q = np.array([records[(el, cell, c)]["q_bader"] for c in cs])
        N0[(el, cell)] = n0
        dq[(el, cell)] = q - q[0]
        n2d[(el, cell)] = (q - q[0]) / A * 1e16
        area[(el, cell)] = A
    return Carriers(cs=cs, N0=N0, dq=dq, n2d=n2d, area=area)


def delta_of_c(car: Carriers, key=("Na", "1x1")):
    """On-site splitting delta(c) for the TB vertex, calibrated on the
    results_v4 Bader charges via the SciPost band-overlap formula
    q_alk = Q_tot (Gamma - 2 delta)/(2 Gamma)  ->
    delta = (Gamma/2)(1 - 2 dq/Q_tot), floored at DELTA_FLOOR."""
    dlt = np.maximum(GAMMA_BANDS / 2 * (1 - 2 * car.dq[key] / Q_TOT),
                     DELTA_FLOOR)
    return PchipInterpolator(car.cs, dlt)


# ======================================================================
# 5. SSH4 TIGHT-BINDING VERTEX (Radha & Lambrecht) — parity-correct
# ======================================================================
def h4(kz: float, delta: float, d: float = 0.0,
       gamma: float = GAMMA_EP) -> np.ndarray:
    """4x4 Bloch Hamiltonian, basis (A1, A2, Co1, Co2) stacked along z.
    tau2/tau4 (alkali <-> CoO2) modulated antisymmetrically by d."""
    t2 = TAU24 * (1 + gamma * d)
    t4 = TAU24 * (1 - gamma * d)
    ph = np.exp(1j * kz)
    return np.array([
        [delta, TAU1, 0, t4 * np.conj(ph)],
        [TAU1, delta, t2, 0],
        [0, t2, -delta, TAU3],
        [t4 * ph, 0, TAU3, -delta]], dtype=complex)


def dh4_dd(kz: float, gamma: float = GAMMA_EP) -> np.ndarray:
    ph = np.exp(1j * kz)
    g = gamma * TAU24
    return np.array([
        [0, 0, 0, -g * np.conj(ph)],
        [0, 0, g, 0],
        [g, 0, 0, 0],
        [-g * ph, 0, 0, 0]], dtype=complex)


ALKALI_BAND = 2   # first alkali-dominated band above the charge-transfer gap


def chi_band(delta: float, gamma: float = GAMMA_EP, kz: float = 0.0,
             band: int = ALKALI_BAND) -> float:
    """Curvature chi = d^2 E_band / dd^2 at d = 0 from exact second-order
    perturbation theory (d^2H/dd^2 = 0 for linearized hoppings):
    chi_n = 2 sum_m |<n|dH/dd|m>|^2 / (E_n - E_m).  The FIRST derivative
    dE/dd vanishes identically at d = 0 by mirror parity (z -> -z maps
    d -> -d) — checked numerically in main()."""
    w, v = np.linalg.eigh(h4(kz, delta))
    dH = dh4_dd(kz, gamma)
    chi = 0.0
    for m in range(4):
        if m == band:
            continue
        me = v[:, band].conj() @ dH @ v[:, m]
        chi += 2 * abs(me) ** 2 / (w[band] - w[m])
    return float(chi)


# ======================================================================
# 6. ELECTRON-PHONON COUPLING AND Tc
# ======================================================================
def allen_dynes_tc(lam: float, hw_eV: float, mu_star: float = MU_STAR) -> float:
    """Allen-Dynes Tc (K) for a single Einstein mode (omega_log =
    <omega^2>^1/2 = omega, f2 = 1); strong-coupling factor f1 included."""
    den = lam - mu_star * (1 + 0.62 * lam)
    if den <= 0 or lam <= 0 or hw_eV <= 0:
        return 0.0
    lam1 = 2.46 * (1 + 3.8 * mu_star)
    f1 = (1 + (lam / lam1) ** 1.5) ** (1 / 3)
    return f1 * hw_eV / KB_EV / 1.2 * np.exp(-1.04 * (1 + lam) / den)


@dataclass
class CouplingResult:
    c: np.ndarray
    alpha: np.ndarray        # canonical Landau alpha (PCHIP)
    N0: np.ndarray           # DFT N(0), per spin per cell (1/eV)
    hw_eff: np.ndarray       # eV, hbar^2/(2M<d^2>_0)
    hw_02: np.ndarray        # eV, even-sector boson energy E2-E0
    split01: np.ndarray      # eV, tunneling splitting
    d2: np.ndarray           # <d^2>_0, A^2
    chi: np.ndarray          # eV/A^2
    lam: np.ndarray          # even-channel (two-phonon) lambda
    tc: np.ndarray           # K; 0 where gated off or polaronic
    status: np.ndarray       # "no-2DEG" | "SC" | "polaronic"


def coupling_vs_c(cgrid, wells_series: pd.DataFrame, car: Carriers,
                  dlt_of_c, mass_amu: float,
                  gamma: float = GAMMA_EP,
                  key=("Na", "1x1")) -> CouplingResult:
    """PCHIP-interpolate the quantization potential (a, b, g6) and the DFT
    N(0) between the four spacings, quantize on the grid, and evaluate the
    even-channel lambda and Allen-Dynes Tc.  Between grid points where the
    interpolated potential would be unbounded (b < 0, g6 ~ 0) a minimal
    g6 = 1e-3 guard keeps it confining."""
    ws = wells_series.sort_values("c")
    ai = PchipInterpolator(ws["c"], ws["a_q"])
    bi = PchipInterpolator(ws["c"], ws["b_q"])
    gi = PchipInterpolator(ws["c"], ws["g_q"])
    al_can = PchipInterpolator(ws["c"], ws["alpha"])
    n0i = PchipInterpolator(car.cs, car.N0[key])
    out = {k: np.zeros_like(cgrid) for k in
           ("hw_eff", "hw_02", "split01", "d2", "chi", "lam", "tc")}
    status = np.empty(len(cgrid), dtype=object)
    for i, c in enumerate(cgrid):
        a, b, g = float(ai(c)), float(bi(c)), max(float(gi(c)), 0.0)
        if b < 0 and g < 1e-3:
            g = 1e-3
        wq = solve_well(a, b, g, mass_amu)
        chi = abs(chi_band(max(float(dlt_of_c(c)), DELTA_FLOOR), gamma))
        n0 = float(n0i(c))
        gcp = 0.5 * chi * wq.d2_02               # quadratic-vertex coupling
        lam = 2 * n0 * gcp**2 / wq.gap_02 if wq.gap_02 > 0 else np.inf
        out["hw_eff"][i] = wq.omega_eff
        out["hw_02"][i] = wq.gap_02
        out["split01"][i] = max(wq.split_01, 0.0)
        out["d2"][i] = wq.d2_00
        out["chi"][i] = chi
        out["lam"][i] = lam
        if n0 < N0_GATE:
            status[i], out["tc"][i] = "no-2DEG", 0.0
        elif lam > LAMBDA_LOC:
            status[i], out["tc"][i] = "polaronic", 0.0
        else:
            status[i] = "SC"
            out["tc"][i] = allen_dynes_tc(lam, wq.gap_02)
    return CouplingResult(c=cgrid, alpha=al_can(cgrid), N0=n0i(cgrid),
                          status=status, **out)


# ======================================================================
# 7. FIGURES
# ======================================================================
def mark_anchors(ax, labels=False):
    x0, x1 = ax.get_xlim()
    for c, lab, col in ((5.5, "anhydrous (no SC)", C_GRAY),
                        (6.9, "1-layer hydrate (no SC)", C_GRAY),
                        (9.9, "2-layer hydrate (SC 4.5 K)", C_RED)):
        if x0 <= c <= x1:
            ax.axvline(c, color=col, ls="--", lw=0.9, alpha=0.6, zorder=0)
            if labels:
                ax.text(c, 0.03, " " + lab, transform=ax.get_xaxis_transform(),
                        fontsize=6.0, color=col, ha="right", va="bottom",
                        rotation=90)
    ax.set_xlim(x0, x1)


def fig_potential_fits(scans, wells):
    fig, axes = plt.subplots(1, 3, figsize=(11.4, 3.5), sharex=True)
    panel = {("Na", "1x1"): 0, ("Na", "s3"): 1, ("Li", "1x1"): 2}
    norm = plt.Normalize(5.5, 9.9)
    dd = np.linspace(-1.6, 1.6, 400)
    for s in scans:
        ax = axes[panel[(s.element, s.cell)]]
        r = wells[(wells.element == s.element) & (wells.cell == s.cell)
                  & (wells.c == s.c)].iloc[0]
        col = CMAP_BLUE(norm(s.c))
        ax.plot(np.concatenate([-s.d[::-1], s.d]),
                np.concatenate([s.E[::-1], s.E]) * 1e3,
                "o", ms=3.4, color=col, mew=0, alpha=0.8)
        ax.plot(dd, vpoly(dd, r["a_q"], r["b_q"], r["g_q"]) * 1e3,
                color=col, lw=1.2,
                ls="--" if r["refit_sextic"] else "-")
    for (el, cell), j in panel.items():
        axes[j].set_title(SERIES_STYLE[(el, cell)][2], fontsize=9)
        axes[j].set_xlabel(r"off-center displacement $\delta$ (Å)")
        axes[j].set_ylim(-650, 900)
    axes[0].set_ylabel(r"$E(\delta) - E(0)$ (meV)")
    sm = plt.cm.ScalarMappable(cmap=CMAP_BLUE, norm=norm)
    fig.colorbar(sm, ax=axes, label=r"CoO$_2$–CoO$_2$ spacing $c$ (Å)",
                 pad=0.015)
    fig.suptitle(r"results_v4 DFT $E(\delta)$ scans + quantization potentials"
                 "  (solid = results.json quartic, dashed = constrained "
                 "sextic refit)", fontsize=10)
    fig.savefig(FIGDIR / "fig1_potential_fits.png", bbox_inches="tight")
    plt.close(fig)


def fig_landau(wells, cstars):
    fig, axes = plt.subplots(1, 3, figsize=(11.0, 3.4))
    cg = np.linspace(5.5, 9.9, 200)
    for (el, cell), (col, mk, lab) in SERIES_STYLE.items():
        w = wells[(wells.element == el) & (wells.cell == cell)].sort_values("c")
        p = PchipInterpolator(w["c"], w["alpha"])
        axes[0].plot(cg, p(cg), color=col, lw=1.3)
        axes[0].plot(w["c"], w["alpha"], mk, ms=5, color=col, mfc="white",
                     mew=1.3, label=lab)
        axes[1].plot(w["c"], w["depth_eV"] * 1e3, mk + "-", ms=5, lw=1.2,
                     color=col, mfc="white", mew=1.3)
        for _, r in w.iterrows():
            if not r["min_bracketed"]:
                axes[1].annotate("", (r["c"], r["depth_eV"] * 1e3 * 1.25),
                                 (r["c"], r["depth_eV"] * 1e3),
                                 arrowprops=dict(arrowstyle="->", color=col,
                                                 lw=0.9))
        axes[2].plot(w["c"], w["d0_A"], mk + "-", ms=5, lw=1.2, color=col,
                     mfc="white", mew=1.3)
    for key, cst in cstars.items():
        if cst is not None:
            axes[0].axvline(cst[0], color=SERIES_STYLE[key][0], ls=":", lw=1)
    axes[0].axhline(0, color="k", lw=0.7)
    axes[0].set_ylabel(r"$\alpha$ (eV/Å$^2$)  [results.json quartic]")
    axes[0].legend(fontsize=7.5)
    axes[1].set_yscale("log")
    axes[1].set_ylabel("well depth (meV)   [↑ = lower bound, min not\n"
                       "bracketed by the sampled δ range]")
    axes[2].set_ylabel(r"well position $d_0$ (Å)")
    for ax in axes:
        ax.set_xlabel(r"$c$ (Å)")
        mark_anchors(ax)
    fig.suptitle("Landau parameters and well geometry (results_v4)",
                 fontsize=10.5)
    fig.tight_layout()
    fig.savefig(FIGDIR / "fig2_landau_parameters.png", bbox_inches="tight")
    plt.close(fig)


def fig_quantum(qt):
    fig, axes = plt.subplots(1, 3, figsize=(11.0, 3.4))
    for (el, cell), (col, mk, lab) in SERIES_STYLE.items():
        q = qt[(qt.element == el) & (qt.cell == cell)
               & (qt.mass == el)].sort_values("c")
        axes[0].semilogy(q["c"], q["hw_eff_meV"], mk + "-", ms=5, lw=1.2,
                         color=col, mfc="white", mew=1.3, label=lab)
        axes[1].semilogy(q["c"], q["gap02_meV"], mk + "-", ms=5, lw=1.2,
                         color=col, mfc="white", mew=1.3)
        axes[1].semilogy(q["c"], np.clip(q["split01_meV"], 1e-13, None),
                         mk + ":", ms=3.5, lw=0.9, color=col, alpha=0.7)
        axes[2].plot(q["c"], q["d_rms_A"], mk + "-", ms=5, lw=1.2, color=col,
                     mfc="white", mew=1.3)
    qn = qt[(qt.element == "Na") & (qt.cell == "1x1")
            & (qt.mass == "Li")].sort_values("c")
    axes[0].semilogy(qn["c"], qn["hw_eff_meV"], "s--", ms=3.5, lw=0.9,
                     color=C_VIOLET, label="Na 1x1 well, Li mass")
    axes[0].set_ylabel(r"$\hbar\omega_{\rm eff}=\hbar^2/2M\langle\delta^2"
                       r"\rangle_0$ (meV)")
    axes[0].legend(fontsize=7)
    axes[1].set_ylabel("levels (meV): solid $E_2{-}E_0$,\n"
                       "dotted $E_1{-}E_0$ (tunneling)")
    axes[1].set_ylim(1e-13, 3e2)
    axes[2].set_ylabel(r"$\sqrt{\langle\delta^2\rangle_0}$ (Å)")
    for ax in axes:
        ax.set_xlabel(r"$c$ (Å)")
        mark_anchors(ax)
    fig.suptitle("Quantized alkali wells (exact 1D Schrödinger, physical "
                 "masses)", fontsize=10.5)
    fig.tight_layout()
    fig.savefig(FIGDIR / "fig3_quantum_wells.png", bbox_inches="tight")
    plt.close(fig)


def fig_carriers(car, dlt):
    fig, axes = plt.subplots(1, 3, figsize=(11.0, 3.4))
    for (el, cell), (col, mk, lab) in SERIES_STYLE.items():
        axes[0].plot(car.cs, car.N0[(el, cell)], mk + "-", ms=5, lw=1.2,
                     color=col, mfc="white", mew=1.3, label=lab)
        axes[1].plot(car.cs, car.dq[(el, cell)], mk + "-", ms=5, lw=1.2,
                     color=col, mfc="white", mew=1.3, label=lab)
    axes[0].axhline(N0_GATE, color=C_RED, ls=":", lw=1,
                    label=f"2DEG gate {N0_GATE}")
    axes[0].set_ylabel(r"DFT $N(0)$ (states/eV/cell/spin)")
    axes[0].legend(fontsize=7)
    axes[1].set_ylabel(r"Bader charge returned to alkali $\Delta q$ (e)")
    axes[1].legend(fontsize=7)
    cg = np.linspace(5.5, 9.9, 200)
    axes[2].plot(cg, dlt(cg), color=C_BLUE, lw=1.4)
    axes[2].axhline(GAMMA_BANDS / 2, color=C_GRAY, ls="--", lw=1,
                    label=r"$\Gamma/2$ = 2.07 eV")
    axes[2].set_ylabel(r"mapped $\delta(c)$ (eV), Na 1x1 Bader")
    axes[2].legend(fontsize=7.5)
    for ax in axes:
        ax.set_xlabel(r"$c$ (Å)")
        mark_anchors(ax)
    fig.suptitle("Carrier turn-on: DFT N(0) and Bader gallery charge tell "
                 "the same story", fontsize=10.5)
    fig.tight_layout()
    fig.savefig(FIGDIR / "fig4_carriers.png", bbox_inches="tight")
    plt.close(fig)


def fig_coupling(res_mid, res_lo, res_hi, cstar):
    fig, axes = plt.subplots(1, 3, figsize=(11.0, 3.4))
    ax = axes[0]
    kz = np.linspace(0, np.pi, 60)
    for dlt, col in ((2.05, C_BLUE), (1.5, C_VIOLET), (0.4, C_YELL)):
        ax.plot(kz, [chi_band(dlt, GAMMA_EP, k) for k in kz], lw=1.4,
                color=col, label=rf"$\delta$ = {dlt:.2f} eV")
    ax.set_xlabel(r"$k_z$ (rad)")
    ax.set_ylabel(r"$\chi = \partial^2 E_{alk}/\partial\delta^2$ (eV/Å$^2$)")
    ax.set_title(rf"TB vertex ($\gamma$ = {GAMMA_EP:.2f}/Å)", fontsize=9.5)
    ax.legend(fontsize=7.5)
    ax = axes[1]
    ax.plot(res_mid.c, res_mid.lam, lw=1.7, color=C_BLUE,
            label=rf"$\lambda$, $\gamma$={GAMMA_EP:.2f}/Å")
    ax.fill_between(res_mid.c, res_lo.lam, res_hi.lam, color=C_BLUE,
                    alpha=0.18, lw=0, label=r"$\gamma \in [1,2]$ /Å")
    ax.axhline(LAMBDA_LOC, color=C_RED, ls=":", lw=1.2,
               label=rf"polaron criterion $\lambda$={LAMBDA_LOC:g}")
    ax.axvline(cstar, color=C_GRAY, ls=":", lw=1)
    ax.set_yscale("log")
    ax.set_ylim(1e-4, 3e3)
    ax.set_xlabel(r"$c$ (Å)")
    ax.set_ylabel(r"$\lambda$ (even channel)")
    ax.legend(fontsize=6.8)
    mark_anchors(ax)
    ax = axes[2]
    ax.plot(res_mid.c, res_mid.hw_02 * 1e3, lw=1.7, color=C_BLUE,
            label=r"$\hbar\Omega = E_2-E_0$")
    ax.plot(res_mid.c, res_mid.hw_eff * 1e3, ls="--", lw=1.2, color=C_VIOLET,
            label=r"$\hbar\omega_{\rm eff}$")
    ax.axvline(cstar, color=C_GRAY, ls=":", lw=1)
    ax.set_yscale("log")
    ax.set_xlabel(r"$c$ (Å)")
    ax.set_ylabel("mode energy (meV)")
    ax.legend(fontsize=7.5)
    mark_anchors(ax)
    fig.suptitle("Electron–phonon coupling on the Na 1x1 axis (DFT N(0), "
                 "Bader-calibrated vertex)", fontsize=10.5)
    fig.tight_layout()
    fig.savefig(FIGDIR / "fig5_coupling.png", bbox_inches="tight")
    plt.close(fig)


def money_plot(wells, cstars, car, res_mid, res_lo, res_hi):
    fig, axes = plt.subplots(3, 1, figsize=(6.8, 8.8), sharex=True,
                             gridspec_kw={"hspace": 0.09})
    C_ANH, C_MONO, C_BI = ANCHOR_SPACINGS
    cst = cstars[("Na", "1x1")]
    cs_na, brack = cst[0], cst[2]

    for ax in axes:
        ax.axvline(C_ANH, color=C_GRAY, ls="--", lw=0.9, alpha=0.7)
        ax.axvline(C_MONO, color=C_GRAY, ls="--", lw=1.1)
        ax.axvline(C_BI, color=C_RED, ls="--", lw=1.1)
        ax.axvspan(*brack, color=C_BLUE, alpha=0.05, lw=0)
        ax.axvline(cs_na, color=C_BLUE, ls=":", lw=1.1)

    # (i) alpha(c)
    ax = axes[0]
    ax.axhline(0, color="k", lw=0.7)
    cg = np.linspace(5.5, 9.9, 200)
    for (el, cell), (col, mk, lab) in SERIES_STYLE.items():
        w = wells[(wells.element == el) & (wells.cell == cell)].sort_values("c")
        ax.plot(cg, PchipInterpolator(w["c"], w["alpha"])(cg), color=col,
                lw=1.3, alpha=0.9)
        ax.plot(w["c"], w["alpha"], mk, ms=5.5, color=col, mfc="white",
                mew=1.3, label=lab)
    ax.annotate(rf"$c^*_{{\rm Na}}$ = {cs_na:.2f} Å"
                f"\n(DFT bracket {brack[0]:.1f}–{brack[1]:.1f} Å)",
                (cs_na + 0.1, 1.4), fontsize=8, color=C_BLUE)
    ax.annotate(r"$c^*_{\rm Li} \leq$ 5.5 Å ($\alpha<0$ already"
                "\nat the smallest spacing)",
                (7.35, 0.35), fontsize=7.5, color=C_AQUA)
    ax.set_ylabel(r"$\alpha(c)$  (eV/Å$^2$)")
    ax.legend(fontsize=7.5, loc="upper right")

    # (ii) carriers: DFT N(0) + Bader-returned gallery charge
    ax = axes[1]
    for (el, cell) in (("Na", "1x1"), ("Li", "1x1")):
        col, mk, lab = SERIES_STYLE[(el, cell)]
        p = PchipInterpolator(car.cs, car.N0[(el, cell)])
        ax.plot(cg, p(cg), color=col, lw=1.3)
        ax.plot(car.cs, car.N0[(el, cell)], mk, ms=5.5, color=col,
                mfc="white", mew=1.3, label=f"N(0), {lab}")
    ax.set_ylabel(r"DFT $N(0)$ (st/eV/cell/spin)")
    ax2 = ax.twinx()
    ax2.grid(False)
    for (el, cell) in (("Na", "1x1"), ("Li", "1x1")):
        col, mk, lab = SERIES_STYLE[(el, cell)]
        p = PchipInterpolator(car.cs, car.n2d[(el, cell)] / 1e13)
        ax2.plot(cg, p(cg), color=col, lw=1.1, ls="--", alpha=0.8)
        ax2.plot(car.cs, car.n2d[(el, cell)] / 1e13, mk, ms=4, color=col,
                 alpha=0.8, label=f"Bader, {lab}")
    ax2.set_ylabel(r"Bader $n_{\rm 2D}=\Delta q/A_{\rm cell}$"
                   "  ($10^{13}$ e/cm$^2$, dashed)")
    ax.legend(fontsize=7, loc="upper left")
    ax2.legend(fontsize=7, loc="lower right")
    ax.text(0.02, 0.55, "2DEG turns on between 5.5 and 6.9 Å\n"
            "in BOTH probes (band N(0) and real-space Bader)",
            transform=ax.transAxes, fontsize=7.5, color=C_GRAY)

    # (iii) Tc(c)
    ax = axes[2]
    def sc_only(res):
        return np.where(res.status == "SC", res.tc, np.nan)
    ax.plot(res_mid.c, sc_only(res_mid), lw=2.0, color=C_BLUE,
            solid_capstyle="round",
            label=rf"Allen–Dynes $T_c$, $\gamma$={GAMMA_EP:.2f}/Å "
                  rf"($z_0$={Z0_ANSATZ} Å), $\mu^*$={MU_STAR}")
    ax.plot(res_lo.c, sc_only(res_lo), lw=1.1, ls="-.", color=C_AQUA,
            label=r"$\gamma$=1.0/Å")
    ax.plot(res_hi.c, sc_only(res_hi), lw=1.1, ls=":", color=C_VIOLET,
            label=r"$\gamma$=2.0/Å")
    pol = res_mid.status == "polaronic"
    if pol.any():
        ax.axvspan(res_mid.c[pol].min(), res_mid.c[pol].max(),
                   color=C_MUTED, alpha=0.20, lw=0)
        ax.text(0.5 * (res_mid.c[pol].min() + min(res_mid.c[pol].max(), 9.9)),
                18.5, "SC killed (polaronic, $\\lambda>2$):\nNa frozen "
                "off-center, Migdal invalid", fontsize=7.5, color=C_GRAY,
                ha="center")
    ng = res_mid.status == "no-2DEG"
    if ng.any():
        ax.text(res_mid.c[ng].mean(), 18.5, "no\n2DEG", fontsize=7.5,
                color=C_GRAY, ha="center")
    ax.plot([C_BI], [TC_EXP], "*", ms=14, color=C_RED, zorder=6,
            label=rf"expt: bilayer hydrate $T_c$ = {TC_EXP} K")
    ax.set_ylabel(r"$T_c$  (K)")
    ax.set_xlabel(r"CoO$_2$–CoO$_2$ interlayer spacing $c$  (Å)")
    ax.set_xlim(5.4, 10.15)
    ax.set_ylim(0, 32)
    ax.legend(fontsize=6.6, loc="upper right")
    ax.annotate("model prediction at 9.9 Å: polaronic, no SC.\n"
                "Reconciling the observed 4.5 K requires the H$_2$O bilayer\n"
                "to soften/re-symmetrize the Na well back toward $c^*$\n"
                "(UNTESTED prediction — needs E($\\delta$) with explicit "
                "water)", (C_BI, 12.5), fontsize=7, color=C_RED, ha="right",
                xytext=(9.75, 8.5),
                arrowprops=dict(arrowstyle="->", color=C_RED, lw=0.9))

    # anchor labels on the top panel
    ymax = axes[0].get_ylim()[1]
    axes[0].text(C_ANH + 0.1, ymax * 1.05, "anhydrous\n(5.5 Å, no SC)",
                 fontsize=7.5, ha="left", color=C_GRAY)
    axes[0].text(C_MONO, ymax * 1.05, "monolayer hydrate\n(6.9 Å, no SC)",
                 fontsize=7.5, ha="center", color=C_GRAY)
    axes[0].text(C_BI, ymax * 1.05, "bilayer hydrate\n(9.9 Å, SC 4.5 K)",
                 fontsize=7.5, ha="center", color=C_RED)

    fig.suptitle("Gallery-controlled Na double well and interlayer-2DEG "
                 "superconductivity (results_v4)", y=0.995, fontsize=11)
    fig.text(0.12, 0.002,
             "Data: QE PBE+PAW, spin-polarized, Co-top alkali site, fixed "
             "$z_O$ (runpod/results_v4; k-mesh ±1.6 meV, spin treatment "
             "<7% on well drops,\n$z_O$±0.06 Å → well depths carry a ±50% "
             "systematic — $\\alpha$ SIGNS are robust). No water anywhere "
             "in the model.",
             fontsize=7, color=C_GRAY, va="bottom")
    fig.savefig(HERE / "money_plot.pdf", bbox_inches="tight")
    fig.savefig(HERE / "money_plot.png", bbox_inches="tight")
    plt.close(fig)


# ======================================================================
# 8. LEGACY MODE (old phase_transition proof-of-concept data)
# ======================================================================
def legacy_quartic(d, a, b):
    return a * d**2 + b * d**4


def legacy_load_scans(compound: str):
    """Parse phase_transition/<compound>/structures_done.txt (JSON-lines of
    pymatgen Structure dicts, energies in Ry)."""
    fname = PT / compound / "structures_done.txt"
    groups: dict[float, dict[float, float]] = {}
    a_lat = None
    for line in fname.read_text().splitlines():
        if not line.strip():
            continue
        r = json.loads(line)
        m = np.array(r["lattice"]["matrix"])
        z33 = m[2][2]
        a_lat = r["lattice"]["a"]
        d = r["sites"][1]["xyz"][2] - (r["sites"][0]["xyz"][2] - z33 / 2.0)
        groups.setdefault(round(z33, 3), {})[round(d, 4)] = r["energy"] * RY_EV
    scans = []
    for key in sorted(groups):
        g = groups[key]
        d0_key = min(g, key=lambda k: abs(k))
        d = np.array(sorted(g))
        E = np.array([g[k] for k in sorted(g)]) - g[d0_key]
        scans.append((key, d, E))
    return scans, a_lat


class _Stub:
    def __init__(self, *a, **k):
        self._args = a

    def __setstate__(self, state):
        if isinstance(state, dict):
            self.__dict__.update(state)
        else:
            self._state = state


class _StubUnpickler(pickle.Unpickler):
    _cache: dict[str, type] = {}

    def find_class(self, module, name):
        if module == "numpy" or module.startswith("numpy."):
            return super().find_class(module, name)
        key = f"{module}.{name}"
        if key not in self._cache:
            self._cache[key] = type(name, (_Stub,), {})
        return self._cache[key]


def legacy_load_bader() -> pd.DataFrame:
    with open(PT / "licoo2" / "charge_data", "rb") as f:
        objs = _StubUnpickler(f).load()
    rows = []
    for s in objs:
        m = np.array(s.__dict__["_lattice"].__dict__["_matrix"])
        rows.append({"c": m[2][2],
                     "q_Li": s.__dict__["_sites"][1].__dict__["charge"]})
    return pd.DataFrame(rows).sort_values("c").reset_index(drop=True)


def legacy_main():
    """Old proof-of-concept pipeline (Questaal lmf relaxed scans): quartic
    fits, c*, quantization, Bader onset.  Prints the summary used for the
    old-vs-new comparison in MODEL.md; tables in theory/results/legacy/."""
    legdir = RESDIR / "legacy"
    legdir.mkdir(exist_ok=True)
    print("=" * 72)
    print("LEGACY MODE: phase_transition/ proof-of-concept data (Questaal lmf)")
    summary = {}
    for comp, el in (("licoo2", "Li"), ("nacoo2", "Na")):
        scans, a_lat = legacy_load_scans(comp)
        rows = []
        for c, d, E in scans:
            ds = np.concatenate([d, -d])
            Es = np.concatenate([E, E])
            (a, b), _ = curve_fit(legacy_quartic, ds, Es, p0=[1.0, 1.0])
            d0 = float(np.sqrt(-a / (2 * b))) if (a < 0 and b > 0) else 0.0
            depth = float(a**2 / (4 * b)) if (a < 0 and b > 0) else 0.0
            q = solve_well(a, b, 0.0, MASS_AMU[el]) if b > 0 else None
            rows.append({"c": c, "alpha": a, "beta": b, "d0_A": d0,
                         "depth_eV": depth,
                         "hw_eff_meV": q.omega_eff * 1e3 if q else np.nan,
                         "gap02_meV": q.gap_02 * 1e3 if q else np.nan})
        fits = pd.DataFrame(rows)
        fits.to_csv(legdir / f"well_fits_{comp}.csv", index=False)
        cs = critical_spacing(fits["c"].values, fits["alpha"].values)
        summary[el] = (fits, cs)
        print(f"  {el}CoO2: {len(fits)} spacings c = {fits['c'].min():.2f}..."
              f"{fits['c'].max():.2f} A (a = {a_lat:.3f} A)")
        if cs:
            print(f"    c* = {cs[0]:.3f} A (pchip; linear {cs[1]:.3f}, "
                  f"bracket {cs[2][0]:.2f}-{cs[2][1]:.2f})")
        idx = (fits["c"] - 6.9).abs().idxmin()
        r69 = fits.loc[idx]
        print(f"    nearest to 6.9 A: c = {r69['c']:.2f}, depth = "
              f"{r69['depth_eV']*1e3:.0f} meV, d0 = {r69['d0_A']:.2f} A")
        print(f"    deepest well: {fits['depth_eV'].max()*1e3:.0f} meV "
              f"at c = {fits.loc[fits['depth_eV'].idxmax(), 'c']:.2f} A")
    try:
        bader = legacy_load_bader()
        i_min = int(bader["q_Li"].idxmin())
        print(f"  Bader (GPAW, LiCoO2): onset (q_Li minimum) at c = "
              f"{bader['c'][i_min]:.2f} A")
    except Exception as e:  # noqa: BLE001 -- legacy pickle is best-effort
        print(f"  Bader legacy data unavailable: {e}")
    print("  (legacy figures are no longer generated; the canonical "
          "pipeline is the default results_v4 mode)")
    print("done (legacy).")


# ======================================================================
# MAIN (results_v4)
# ======================================================================
def main_v4():
    print("=" * 72)
    print("1. DATA EXTRACTION (runpod/results_v4)")
    records, scans, manifest = load_v4()
    print(f"   {len(scans)} E(delta) scans: "
          + ", ".join(sorted({f'{s.element}-{s.cell}' for s in scans}))
          + f" at c = {sorted({s.c for s in scans})} A")
    print(f"   a_lat: Na {manifest['a_lat']['Na']} A, "
          f"Li {manifest['a_lat']['Li']} A; z_O = {manifest['z_O']} A "
          "(fixed, Co-top site)")

    print("=" * 72)
    print("2. WELL POTENTIALS (canonical quartic + constrained sextic refits)")
    wells = build_well_table(records, scans)
    wells.to_csv(RESDIR / "well_fits_v4.csv", index=False)
    for _, r in wells.iterrows():
        tag = "sextic refit" if r["refit_sextic"] else "results.json quartic"
        lb = "" if r["min_bracketed"] else \
            f"  [min NOT bracketed by data: depth >= " \
            f"{-r['E_edge_meV']:.0f} meV is extrapolated]"
        print(f"   {r['element']}-{r['cell']} c={r['c']:>4}: "
              f"alpha={r['alpha']:+.3f} eV/A^2, d0={r['d0_A']:.2f} A, "
              f"depth={r['depth_eV']*1e3:6.0f} meV ({tag}, "
              f"rms {r['rms_meV']:.1f} meV){lb}")

    cstars = {}
    for el, cell in SERIES_STYLE:
        w = wells[(wells.element == el) & (wells.cell == cell)].sort_values("c")
        cstars[(el, cell)] = critical_spacing(w["c"].values, w["alpha"].values)
    cst_na = cstars[("Na", "1x1")]
    print(f"   c*_Na (1x1)  = {cst_na[0]:.2f} A (PCHIP; linear "
          f"{cst_na[1]:.2f} A) -- bracketed by DFT points at "
          f"{cst_na[2][0]:.1f} and {cst_na[2][1]:.1f} A")
    cst_s3 = cstars[("Na", "s3")]
    print(f"   c*_Na (x=1/3) = {cst_s3[0]:.2f} A (PCHIP; linear "
          f"{cst_s3[1]:.2f} A)")
    assert cstars[("Li", "1x1")] is None
    print("   c*_Li: alpha(5.5 A) = "
          f"{wells[(wells.element == 'Li')]['alpha'].iloc[0]:+.3f} < 0 -- "
          "only an UPPER BOUND c*_Li <= 5.5 A in this dataset")

    print("=" * 72)
    print("3. QUANTUM WELLS (exact 1D Schroedinger, physical masses)")
    qt = quantum_table(wells)
    qt.to_csv(RESDIR / "quantum_wells_v4.csv", index=False)
    qn = qt[(qt.element == "Na") & (qt.cell == "1x1") & (qt.mass == "Na")]
    print("   Na 1x1 (Na mass): hw_eff = "
          + ", ".join(f"{r['hw_eff_meV']:.3g} meV @ {r['c']} A"
                      for _, r in qn.sort_values("c").iterrows()))
    print("   deep wells: omega_eff = hbar^2/(2M<d^2>_0) BY CONSTRUCTION "
          "(E1-E0 is a tunnel splitting there, "
          f"{qn['split01_meV'].min():.1e} meV at c=9.9)")

    print("=" * 72)
    print("4. CARRIERS: DFT N(0) vs Bader gallery charge")
    car = build_carriers(records, manifest)
    for key in (("Na", "1x1"), ("Li", "1x1")):
        n0, dq = car.N0[key], car.dq[key]
        r = np.corrcoef(n0, dq)[0, 1]
        print(f"   {key[0]}-{key[1]}: N(0) = "
              + "/".join(f"{v:.2f}" for v in n0)
              + " st/eV/cell/spin;  Bader dq = "
              + "/".join(f"{v:.2f}" for v in dq)
              + f" e  (corr {r:.2f}; both turn on between 5.5 and 6.9 A)")
    n0s3 = car.N0[("Na", "s3")]
    print(f"   Na-s3 (x=1/3): N(0) = {n0s3.min():.1f}-{n0s3.max():.1f} "
          "flat -- Co-d metallic background at x=1/3 masks the gallery "
          "turn-on in N(0); Bader dq_s3 = "
          + "/".join(f"{v:+.2f}" for v in car.dq[("Na", "s3")])
          + " e (weak/non-monotone)")
    dlt = delta_of_c(car)
    print("   Bader-mapped on-site splitting delta(c) = "
          + ", ".join(f"{float(dlt(c)):.2f} eV @ {c} A" for c in car.cs))

    print("=" * 72)
    print("5. SSH4 VERTEX (parity-correct) + LAMBDA + Tc")
    for kz in (0.0, 1.0, np.pi):
        w0 = np.linalg.eigvalsh(h4(kz, 2.0, d=1e-4))[ALKALI_BAND]
        w1 = np.linalg.eigvalsh(h4(kz, 2.0, d=-1e-4))[ALKALI_BAND]
        assert abs(w0 - w1) / 2e-4 < 1e-6, "linear DP should vanish by parity"
    print("   linear deformation potential dE/dd|_{d=0} = 0 (parity) -- "
          "verified; leading vertex is quadratic (chi)")
    print("   NOTE: results.json's 'I' is a linear slope fitted over the "
          "one-sided delta scan (finite-delta secant), not the parity-"
          "vanishing derivative at 0; it is superseded here.")
    wna = wells[(wells.element == "Na") & (wells.cell == "1x1")]
    cgrid = np.linspace(5.5, 9.9, 221)
    res_mid = coupling_vs_c(cgrid, wna, car, dlt, MASS_AMU["Na"], GAMMA_EP)
    res_lo = coupling_vs_c(cgrid, wna, car, dlt, MASS_AMU["Na"],
                           1.0 / Z0_RANGE[1])
    res_hi = coupling_vs_c(cgrid, wna, car, dlt, MASS_AMU["Na"],
                           1.0 / Z0_RANGE[0])
    pd.DataFrame({
        "c_A": res_mid.c, "alpha": res_mid.alpha, "N0": res_mid.N0,
        "hw_eff_meV": res_mid.hw_eff * 1e3, "hw02_meV": res_mid.hw_02 * 1e3,
        "split01_meV": res_mid.split01 * 1e3, "d2_A2": res_mid.d2,
        "chi_eV_A2": res_mid.chi, "lambda": res_mid.lam,
        "Tc_K": res_mid.tc, "status": res_mid.status,
        "lambda_g1": res_lo.lam, "Tc_K_g1": res_lo.tc,
        "lambda_g2": res_hi.lam, "Tc_K_g2": res_hi.tc,
    }).to_csv(RESDIR / "coupling_tc_vs_c.csv", index=False)

    for gname, res in ((f"gamma={GAMMA_EP:.2f}", res_mid),
                       ("gamma=1.00", res_lo), ("gamma=2.00", res_hi)):
        ok = res.status == "SC"
        if ok.any() and np.nanmax(res.tc[ok]) > 0.01:
            i_pk = int(np.argmax(np.where(ok, res.tc, -1)))
            print(f"   [{gname}] SC window c = {res.c[ok].min():.2f}..."
                  f"{res.c[ok].max():.2f} A, peak Tc = {res.tc[i_pk]:.1f} K "
                  f"at c = {res.c[i_pk]:.2f} A (lambda = {res.lam[i_pk]:.2f},"
                  f" hbar*Omega = {res.hw_02[i_pk]*1e3:.1f} meV)")
        else:
            print(f"   [{gname}] no SC window")
    for target in ANCHOR_SPACINGS:
        j = int(np.argmin(np.abs(res_mid.c - target)))
        st = res_mid.status[j]
        lamtxt = (f"lambda = {res_mid.lam[j]:.2f}" if res_mid.lam[j] <= 5
                  else "lambda >> 2")
        print(f"   at c = {target} A: {lamtxt}, N(0) = {res_mid.N0[j]:.2f} "
              f"-> {st}" + ("  [expt: SC 4.5 K -- reconciliation requires "
                            "water-softened well]" if target == 9.9 else ""))

    print("=" * 72)
    print("6. FIGURES")
    fig_potential_fits(scans, wells)
    fig_landau(wells, cstars)
    fig_quantum(qt)
    fig_carriers(car, dlt)
    fig_coupling(res_mid, res_lo, res_hi, cst_na[0])
    money_plot(wells, cstars, car, res_mid, res_lo, res_hi)
    print(f"   figures -> {FIGDIR}")
    print(f"   money plot -> {HERE / 'money_plot.pdf'} (+ .png)")
    print("done.")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--legacy", action="store_true",
                    help="re-run the old phase_transition/ proof-of-concept "
                         "analysis (comparison numbers only)")
    args = ap.parse_args()
    if args.legacy:
        legacy_main()
    else:
        main_v4()


if __name__ == "__main__":
    main()

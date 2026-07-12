# /// script
# requires-python = ">=3.10"
# dependencies = ["numpy", "scipy", "matplotlib", "pandas"]
# ///
"""
Effective model for gallery-controlled superconductivity in Na_xCoO2·yH2O
=========================================================================

End-to-end pipeline, from the raw DFT data in this repository to the
"money plot".  Run with

    uv run theory/effective_model.py          (or)
    theory/.venv/bin/python theory/effective_model.py

Pipeline
--------
1.  Parse the lmf (Questaal) total-energy scans E(d; c) for LiCoO2 and
    NaCoO2 from phase_transition/{licoo2,nacoo2}/structures_done.txt
    (JSON-lines of pymatgen Structure dicts, energies in Ry).
2.  Fit V(d; c) = alpha(c) d^2 + beta(c) d^4 per interlayer spacing c,
    extract the critical spacing c* where alpha crosses zero.
3.  Solve the 1D Schroedinger equation in each well (finite differences,
    parity-resolved) for Li and Na masses: levels, tunneling splitting,
    <d^2>, omega_eff = hbar/(2M<d^2>_0).
4.  Wire d into the SSH4 tight-binding model of Radha & Lambrecht,
    SciPost Phys. 10, 057 (2021): tau2 -> tau2(1+gamma d),
    tau4 -> tau4(1-gamma d).  The linear deformation potential vanishes
    by parity; the leading electron-lattice vertex is quadratic,
    E_band(d) = E_edge + (1/2) chi d^2, with chi from exact second-order
    perturbation theory in the 4x4 Bloch Hamiltonian.
5.  lambda(c) and Allen-Dynes Tc(c) for the single anharmonic Na mode,
    with the 2DEG turn-on gate calibrated on the repo's Bader-charge data.
6.  Publication figures: intermediate diagnostics in theory/figures/,
    the money plot in theory/money_plot.{pdf,png}; numeric tables in
    theory/results/.

Every number is traceable to (a) the repo DFT data, (b) the SciPost
paper's fitted parameters, or (c) an explicitly labeled ansatz with a
stated range.  See theory/MODEL.md.
"""

from __future__ import annotations

import json
import pickle
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.linalg import eigh_tridiagonal
from scipy.optimize import curve_fit

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

# ----------------------------------------------------------------------
# Paths
# ----------------------------------------------------------------------
HERE = Path(__file__).resolve().parent          # theory/
REPO = HERE.parent
PT = REPO / "phase_transition"
FIGDIR = HERE / "figures"
RESDIR = HERE / "results"
FIGDIR.mkdir(exist_ok=True)
RESDIR.mkdir(exist_ok=True)

# ----------------------------------------------------------------------
# Physical constants (CODATA)
# ----------------------------------------------------------------------
RY_EV = 13.605693122994          # Ry -> eV
HBAR2_OVER_2AMU = 2.0901e-3      # hbar^2 / (2 * 1 amu) in eV * Angstrom^2
HBAR2_OVER_2ME = 3.80998         # hbar^2 / (2 m_e)     in eV * Angstrom^2
KB_EV = 8.617333e-5              # eV / K
M_LI = 6.941                     # amu
M_NA = 22.98977                  # amu

# ----------------------------------------------------------------------
# Model parameters
# ----------------------------------------------------------------------
# --- SciPost Phys. 10, 057 (2021), Radha & Lambrecht (fitted TB values) ---
TAU1 = 0.5      # eV  t_Li^z        (on-site alkali sp_z pair hopping)
TAU3 = 0.5      # eV  t_CoO2^z      (hopping across the CoO2 slab)
TAU24 = 2.0     # eV  t_{Li-CoO2}^z (alkali <-> CoO2 interlayer hopping)
T_LI_XY = -0.6  # eV  in-plane alkali hopping
T_CO_XY = 0.09  # eV  in-plane CoO2 hopping
GAMMA_BANDS = 6 * abs(T_LI_XY) + 6 * abs(T_CO_XY)   # = 4.14 eV, band-overlap Gamma

# --- Explicit ansatz (stated range; see MODEL.md) ---
Z0_ANSATZ = 0.70          # Angstrom, hopping decay length t(z) ~ exp(-z/z0)
Z0_RANGE = (0.5, 1.0)     # sensitivity range -> gamma = 1/z0 in [1, 2] 1/Angstrom
GAMMA_EP = 1.0 / Z0_ANSATZ
MU_STAR = 0.10            # Coulomb pseudopotential (task prescription)
LAMBDA_LOC = 2.0          # polaronic-localization criterion, ansatz range [1, 2]
Q_TOT = 1.0               # e/cell: full conversion of the alkali band (paper limit)

# ----------------------------------------------------------------------
# Dataviz palette (validated reference palette from the dataviz skill)
# ----------------------------------------------------------------------
C_BLUE = "#2a78d6"     # series 1: primary model / NaCoO2
C_AQUA = "#1baf7a"     # series 2: LiCoO2 / secondary
C_YELL = "#eda100"     # series 3
C_VIOLET = "#4a3aa7"   # series 5
C_RED = "#e34948"      # reserved: SC anchor
C_GRAY = "#52514e"     # secondary text / neutral anchor
C_MUTED = "#9a9992"
CMAP_BLUE = LinearSegmentedColormap.from_list(
    "seq_blue", ["#cde2fb", "#5598e7", "#1c5cab", "#0d366b"])
CMAP_AQUA = LinearSegmentedColormap.from_list(
    "seq_aqua", ["#c8ecdd", "#4cc39a", "#0f7a55", "#064732"])

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


# ======================================================================
# 1. DATA EXTRACTION
# ======================================================================
@dataclass
class WellScan:
    """One constant-c DFT scan: total energy vs alkali off-center displacement."""
    c: float                 # CoO2-CoO2 interlayer spacing along z (Angstrom)
    c_lat: float             # |third lattice vector| (the notebooks' "c")
    d: np.ndarray            # off-center displacement (Angstrom), one-sided (<= 0)
    E: np.ndarray            # energy relative to centered position (eV)


def load_scans(compound: str) -> tuple[list[WellScan], float]:
    """Parse structures_done.txt (JSON-lines, energies in Ry -> eV).

    The alkali off-center displacement is computed geometrically:
    d = z_alkali - (z_Co - z33/2), where z33 is the z-component of the
    third lattice vector (= the CoO2-CoO2 layer spacing).  Every c-group
    contains the centered point d = 0, which is used as energy reference.
    Returns the scans and the in-plane lattice constant a.
    """
    fname = PT / compound / "structures_done.txt"
    groups: dict[float, dict[float, float]] = {}
    latmeta = {}
    a_lat = None
    for line in fname.read_text().splitlines():
        if not line.strip():
            continue
        r = json.loads(line)
        m = np.array(r["lattice"]["matrix"])
        z33 = m[2][2]
        a_lat = r["lattice"]["a"]
        z_alk = r["sites"][1]["xyz"][2]
        z_co = r["sites"][0]["xyz"][2]
        d = z_alk - (z_co - z33 / 2.0)
        key = round(z33, 3)
        groups.setdefault(key, {})[round(d, 4)] = r["energy"] * RY_EV
        latmeta[key] = r["lattice"]["c"]
    scans = []
    for key in sorted(groups):
        g = groups[key]
        d0_key = min(g, key=lambda k: abs(k))     # centered point
        d = np.array(sorted(g))
        E = np.array([g[k] for k in sorted(g)]) - g[d0_key]
        scans.append(WellScan(c=key, c_lat=latmeta[key], d=d, E=E))
    return scans, a_lat


class _Stub:
    """Generic stand-in for pymatgen classes when unpickling charge_data."""
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


def load_bader() -> pd.DataFrame:
    """GPAW/Bader valence charges for the minimum-energy LiCoO2 structures
    (phase_transition/licoo2/charge_data: pickled pymatgen structures with
    per-site Bader electron counts; unpickled with stub classes so the
    pipeline has no pymatgen dependency)."""
    with open(PT / "licoo2" / "charge_data", "rb") as f:
        objs = _StubUnpickler(f).load()
    rows = []
    for s in objs:
        m = np.array(s.__dict__["_lattice"].__dict__["_matrix"])
        sites = s.__dict__["_sites"]
        rows.append({
            "c": m[2][2],
            "q_Co": sites[0].__dict__["charge"],
            "q_Li": sites[1].__dict__["charge"],
            "q_O1": sites[2].__dict__["charge"],
            "q_O2": sites[3].__dict__["charge"],
        })
    return pd.DataFrame(rows).sort_values("c").reset_index(drop=True)


# ======================================================================
# 2. QUARTIC FITS  V(d) = alpha d^2 + beta d^4
# ======================================================================
def quartic(d, a, b):
    return a * d**2 + b * d**4


def fit_wells(scans: list[WellScan]) -> pd.DataFrame:
    rows = []
    for s in scans:
        ds = np.concatenate([s.d, -s.d])
        Es = np.concatenate([s.E, s.E])
        (a, b), _ = curve_fit(quartic, ds, Es, p0=[1.0, 1.0])
        dmin = float(np.sqrt(-a / (2 * b))) if (a < 0 and b > 0) else 0.0
        depth = float(a**2 / (4 * b)) if (a < 0 and b > 0) else 0.0
        rms = float(np.sqrt(np.mean((quartic(s.d, a, b) - s.E) ** 2)))
        rows.append({"c": s.c, "c_lat": s.c_lat, "alpha": a, "beta": b,
                     "d_min": dmin, "depth_eV": depth, "rms_eV": rms,
                     "n_pts": len(s.d), "d_range": float(-s.d.min())})
    return pd.DataFrame(rows)


def critical_spacing(fits: pd.DataFrame, nfit: int = 4) -> tuple[float, float]:
    """c* where alpha(c) crosses zero: local linear fit alpha = a0*(c - c*)
    through the nfit points nearest the crossing."""
    c = fits["c"].values
    al = fits["alpha"].values
    i = int(np.argmin(np.abs(al)))
    lo = max(0, i - nfit // 2)
    hi = min(len(c), lo + nfit)
    lo = max(0, hi - nfit)
    p = np.polyfit(c[lo:hi], al[lo:hi], 1)
    cstar = -p[1] / p[0]
    return float(cstar), float(p[0])   # c*, slope a0 (eV/A^2 per A)


# ======================================================================
# 3. QUANTUM WELL: parity-resolved 1D Schroedinger solver
# ======================================================================
@dataclass
class WellQuantum:
    E_even: np.ndarray       # eigenvalues of even-parity sector (eV)
    E_odd: np.ndarray        # eigenvalues of odd-parity sector (eV)
    d2_00: float             # <0|d^2|0>   (A^2)
    d2_02: float             # <0_e|d^2|1_e>  matrix element (A^2)
    omega_eff: float         # hbar/(2M<d^2>_0) as energy hbar*omega_eff (eV)
    split_01: float          # tunneling splitting E1-E0 (eV)
    gap_02: float            # even-sector gap E2-E0 (eV)


def solve_well(alpha: float, beta: float, mass_amu: float,
               L: float = 2.6, N: int = 1300, nstates: int = 4) -> WellQuantum:
    """Parity-resolved finite differences on the half line (staggered grid
    x_j = (j+1/2) h, j = 0..N-1).  Even sector: psi(-x)=psi(x) ->
    H_00 = t + V0; odd sector: psi(-x)=-psi(x) -> H_00 = 3t + V0.
    This separates the tunneling doublet exactly (no degeneracy issues)."""
    h = L / N
    x = (np.arange(N) + 0.5) * h
    V = quartic(x, alpha, beta)
    t = HBAR2_OVER_2AMU / mass_amu / h**2
    diag = 2 * t + V
    off = -t * np.ones(N - 1)

    def solve(parity: str):
        dg = diag.copy()
        dg[0] = (t if parity == "even" else 3 * t) + V[0]
        w, v = eigh_tridiagonal(dg, off, select="i",
                                select_range=(0, nstates - 1))
        v = v / np.sqrt(2 * h * np.sum(v**2, axis=0))    # full-line norm
        return w, v

    we, ve = solve("even")
    wo, vo = solve("odd")
    d2_00 = 2 * h * float(np.sum(x**2 * ve[:, 0] ** 2))
    d2_02 = 2 * h * float(np.sum(x**2 * ve[:, 0] * ve[:, 1]))
    return WellQuantum(
        E_even=we, E_odd=wo,
        d2_00=d2_00, d2_02=abs(d2_02),
        omega_eff=(HBAR2_OVER_2AMU / mass_amu) / d2_00,   # hbar^2/(2M<d^2>) [eV]
        split_01=float(wo[0] - we[0]),
        gap_02=float(we[1] - we[0]),
    )


def quantum_table(fits: pd.DataFrame, mass_amu: float, label: str) -> pd.DataFrame:
    rows = []
    for _, r in fits.iterrows():
        q = solve_well(r["alpha"], r["beta"], mass_amu)
        rows.append({
            "c": r["c"], "alpha": r["alpha"], "beta": r["beta"],
            "E0_meV": q.E_even[0] * 1e3, "E1_meV": q.E_odd[0] * 1e3,
            "E2_meV": q.E_even[1] * 1e3, "E3_meV": q.E_odd[1] * 1e3,
            "split01_meV": q.split_01 * 1e3, "gap02_meV": q.gap_02 * 1e3,
            "d2_A2": q.d2_00, "d_rms_A": np.sqrt(q.d2_00),
            "hw_eff_meV": q.omega_eff * 1e3,
            "mass": label,
        })
    return pd.DataFrame(rows)


# ======================================================================
# 4. SSH4 TIGHT-BINDING MODEL (Radha & Lambrecht) + displacement coupling
# ======================================================================
def h4(kz: float, delta: float, d: float = 0.0, gamma: float = GAMMA_EP,
       fL: float = 0.0, fC: float = 0.0) -> np.ndarray:
    """4x4 Bloch Hamiltonian, basis (Li1, Li2, Co1, Co2) stacked along z.
    tau1 within the alkali sp_z pair, tau2/tau4 alkali<->CoO2 (modulated
    antisymmetrically by the off-center displacement d), tau3 across the
    CoO2 slab.  On-site +/-delta (alkali/CoO2 electronegativity splitting),
    in-plane dispersions fL, fC added on the respective diagonals."""
    t2 = TAU24 * (1 + gamma * d)
    t4 = TAU24 * (1 - gamma * d)
    ph = np.exp(1j * kz)
    return np.array([
        [delta + fL, TAU1, 0, t4 * np.conj(ph)],
        [TAU1, delta + fL, t2, 0],
        [0, t2, -delta + fC, TAU3],
        [t4 * ph, 0, TAU3, -delta + fC]], dtype=complex)


def dh4_dd(kz: float, gamma: float = GAMMA_EP) -> np.ndarray:
    ph = np.exp(1j * kz)
    g = gamma * TAU24
    return np.array([
        [0, 0, 0, -g * np.conj(ph)],
        [0, 0, g, 0],
        [g, 0, 0, 0],
        [-g * ph, 0, 0, 0]], dtype=complex)


ALKALI_BAND = 2   # first Li-dominated band above the charge-transfer gap


def chi_band(delta: float, gamma: float = GAMMA_EP, kz: float = 0.0,
             band: int = ALKALI_BAND) -> float:
    """Curvature chi = d^2 E_band / dd^2 at d=0 from exact second-order
    perturbation theory (dH/dd is the full derivative; d^2H/dd^2 = 0 for
    the linearized hoppings):  chi_n = 2 sum_m |<n|dH/dd|m>|^2/(E_n-E_m).
    The FIRST derivative dE/dd vanishes identically at d=0 by mirror
    parity (z -> -z maps d -> -d) -- checked numerically in main()."""
    w, v = np.linalg.eigh(h4(kz, delta))
    dH = dh4_dd(kz, gamma)
    chi = 0.0
    for m in range(4):
        if m == band:
            continue
        me = v[:, band].conj() @ dH @ v[:, m]
        chi += 2 * abs(me) ** 2 / (w[band] - w[m])
    return float(chi)


def inplane_f(kx: float, ky: float, t: float, a: float) -> float:
    """f_i(k) = 2 t_i sum_delta cos(k . delta) on the triangular lattice."""
    return 2 * t * (np.cos(kx * a)
                    + 2 * np.cos(kx * a / 2) * np.cos(np.sqrt(3) * ky * a / 2))


def alkali_mass_and_dos(delta: float, a: float) -> tuple[float, float, float]:
    """In-plane effective mass of the alkali band bottom (k|| = Gamma,
    kz = 0 -- verified to be the band minimum for t_Li^xy < 0) and the
    2D step DOS N(0) = m*/(2 pi hbar^2), per spin, per unit area and per
    unit cell (A_cell = sqrt(3)/2 a^2)."""
    def Eb(kx):
        fL = inplane_f(kx, 0.0, T_LI_XY, a)
        fC = inplane_f(kx, 0.0, T_CO_XY, a)
        w = np.linalg.eigvalsh(h4(0.0, delta, fL=fL, fC=fC))
        return w[ALKALI_BAND]
    k = 0.02
    curv = (Eb(k) - 2 * Eb(0.0) + Eb(-k)) / k**2       # = 2 * hbar^2/2m*
    mstar_me = HBAR2_OVER_2ME / (curv / 2.0)
    n0_area = mstar_me / (2 * np.pi * HBAR2_OVER_2ME)  # 1/(eV A^2) per spin
    A_cell = np.sqrt(3) / 2 * a**2
    return mstar_me, n0_area, n0_area * A_cell


# ----------------------------------------------------------------------
# 4a. Gallery-opening mapping c -> (delta, q_alkali), calibrated on Bader
# ----------------------------------------------------------------------
@dataclass
class GalleryMap:
    """q_alkali/q_CoO2 = (Gamma - 2 delta)/(Gamma + 2 delta) for 2delta <
    Gamma, else 0 [SciPost paper, 2D band-overlap argument].  With
    q_alkali + q_CoO2 = Q_TOT this gives q_alkali = Q_TOT (Gamma - 2
    delta)/(2 Gamma).  We take delta(c) linear, crossing Gamma/2 at the
    DFT charge-onset spacing c_on with width w calibrated so the model
    slope dq/dc reproduces the Bader slope of the repo data."""
    c_on: float      # onset spacing (A) -- Bader minimum, coincides with c*
    w: float         # A; delta(c) = Gamma/2 (1 - (c - c_on)/w)
    A_cell: float    # A^2

    def delta(self, c):
        return np.maximum(GAMMA_BANDS / 2 * (1 - (np.asarray(c, float) - self.c_on) / self.w),
                          0.05)   # floor keeps TB spectrum non-degenerate

    def q(self, c):   # e / cell
        return np.clip(Q_TOT * (np.asarray(c, float) - self.c_on) / (2 * self.w),
                       0.0, Q_TOT)

    def n_cm2(self, c):
        return self.q(c) / self.A_cell * 1e16   # e / cm^2


def calibrate_gallery(bader: pd.DataFrame, cstar_li: float, cstar_na: float,
                      a_na: float) -> tuple[GalleryMap, float, float, np.ndarray]:
    """Onset = Bader minimum of q_Li(c).  The pre-onset points define a
    linear baseline (chemical background, decreasing with c); the gallery
    charge is Delta q = q_Li - baseline.  The model slope dq/dc is fit to
    the post-onset Delta q.  Transferred to the Na axis by the shift
    c*_Na - c*_Li."""
    i_min = int(bader["q_Li"].idxmin())
    c_on_li = float(bader["c"][i_min])
    pre = bader.iloc[: i_min + 1]
    base = np.polyfit(pre["c"], pre["q_Li"], 1)
    dq = np.maximum(bader["q_Li"].values - np.polyval(base, bader["c"].values),
                    0.0)
    post = bader.iloc[i_min:]
    slope = float(np.polyfit(post["c"], dq[i_min:], 1)[0])  # e/A
    w = Q_TOT / (2 * slope)
    c_on_na = c_on_li + (cstar_na - cstar_li)
    A_cell = np.sqrt(3) / 2 * a_na**2
    return GalleryMap(c_on=c_on_na, w=w, A_cell=A_cell), c_on_li, slope, dq


# ======================================================================
# 5. ELECTRON-PHONON COUPLING AND Tc
# ======================================================================
def allen_dynes_tc(lam: float, hw_eV: float, mu_star: float = MU_STAR) -> float:
    """Allen-Dynes Tc (K) for a single Einstein mode (omega_log = <omega^2>^1/2
    = omega, so f2 = 1); strong-coupling factor f1 included."""
    den = lam - mu_star * (1 + 0.62 * lam)
    if den <= 0 or lam <= 0 or hw_eV <= 0:
        return 0.0
    lam1 = 2.46 * (1 + 3.8 * mu_star)
    f1 = (1 + (lam / lam1) ** 1.5) ** (1 / 3)
    return f1 * hw_eV / KB_EV / 1.2 * np.exp(-1.04 * (1 + lam) / den)


@dataclass
class CouplingResult:
    c: np.ndarray
    alpha: np.ndarray
    beta: np.ndarray
    q: np.ndarray            # e/cell in alkali band
    n_cm2: np.ndarray
    hw_eff: np.ndarray       # eV, hbar*omega_eff = hbar^2/(2M<d^2>_0)
    hw_02: np.ndarray        # eV, even-sector boson energy E2-E0
    split01: np.ndarray      # eV, tunneling splitting
    d2: np.ndarray           # <d^2>_0, A^2
    chi: np.ndarray          # eV/A^2
    lam: np.ndarray          # primary: two-phonon/even-channel lambda
    lam_task: np.ndarray     # cross-check: N0 I^2/(M w_eff^2), I = chi*d_rms
    tc: np.ndarray           # K (primary lambda, gated by 2DEG existence)


def coupling_vs_c(cgrid: np.ndarray, fits: pd.DataFrame, gal: GalleryMap,
                  mass_amu: float, n0_cell: float,
                  gamma: float = GAMMA_EP) -> CouplingResult:
    """Interpolate alpha(c), beta(c) (linear inside the DFT range, constant
    outside = conservative rigid extrapolation), solve the quantum well and
    evaluate lambda and Tc on the grid."""
    al = np.interp(cgrid, fits["c"], fits["alpha"])
    be = np.interp(cgrid, fits["c"], fits["beta"])
    out = {k: np.zeros_like(cgrid) for k in
           ("hw_eff", "hw_02", "split01", "d2", "chi", "lam", "lam_task", "tc")}
    q = gal.q(cgrid)
    for i, c in enumerate(cgrid):
        wq = solve_well(al[i], be[i], mass_amu)
        chi = abs(chi_band(float(gal.delta(c)), gamma))
        out["hw_eff"][i] = wq.omega_eff
        out["hw_02"][i] = wq.gap_02
        out["split01"][i] = max(wq.split_01, 0.0)
        out["d2"][i] = wq.d2_00
        out["chi"][i] = chi
        # -- primary: quadratic (two-phonon / even-channel) coupling.
        #    H_ep = (1/2) chi d^2 n_el; boson = |0_e> -> |1_e| transition:
        #    g = (1/2) chi <0|d^2|2>,  lambda = 2 N(0) g^2 / (E2 - E0)
        g = 0.5 * chi * wq.d2_02
        lam = 2 * n0_cell * g**2 / wq.gap_02 if wq.gap_02 > 0 else 0.0
        # -- task prescription cross-check: lambda = N0 I^2 / (M w_eff^2),
        #    I = chi * sqrt(<d^2>) (rms Hellmann-Feynman derivative; the
        #    linear deformation potential at d = 0 vanishes by parity)
        Mw2 = wq.omega_eff**2 / (2 * HBAR2_OVER_2AMU / mass_amu)  # eV/A^2
        lam_task = n0_cell * (chi**2 * wq.d2_00) / Mw2 if Mw2 > 0 else 0.0
        out["lam"][i] = lam
        out["lam_task"][i] = lam_task
        out["tc"][i] = allen_dynes_tc(lam, wq.gap_02) if q[i] > 0 else 0.0
    return CouplingResult(c=cgrid, alpha=al, beta=be, q=q,
                          n_cm2=gal.n_cm2(cgrid), **out)


# ======================================================================
# 6. FIGURES
# ======================================================================
def fig_potential_fits(scans_li, fits_li, scans_na, fits_na):
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.6), sharey=True)
    for ax, scans, fits, cmap, name in (
            (axes[0], scans_li, fits_li, CMAP_AQUA, "LiCoO$_2$"),
            (axes[1], scans_na, fits_na, CMAP_BLUE, "NaCoO$_2$")):
        cvals = fits["c"].values
        norm = plt.Normalize(cvals.min(), cvals.max())
        dd = np.linspace(-1.35, 1.35, 300)
        for s, (_, r) in zip(scans, fits.iterrows()):
            col = cmap(norm(r["c"]))
            ax.plot(np.concatenate([s.d, -s.d[::-1]]),
                    np.concatenate([s.E, s.E[::-1]]) * 1e3,
                    "o", ms=2.6, color=col, alpha=0.55, mew=0)
            ax.plot(dd, quartic(dd, r["alpha"], r["beta"]) * 1e3,
                    color=col, lw=1.1, alpha=0.9)
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        fig.colorbar(sm, ax=ax, label=r"CoO$_2$–CoO$_2$ spacing $c$ (Å)")
        ax.set_xlabel(r"alkali off-center displacement $d$ (Å)")
        ax.set_title(f"{name}: DFT $E(d)$ + quartic fits", fontsize=10)
        ax.set_ylim(-550, 900)
    axes[0].set_ylabel(r"$E(d) - E(0)$ (meV)")
    fig.tight_layout()
    fig.savefig(FIGDIR / "fig1_potential_fits.png", bbox_inches="tight")
    plt.close(fig)


def fig_landau_parameters(fits_li, cs_li, fits_na, cs_na):
    fig, axes = plt.subplots(2, 2, figsize=(8.6, 5.6), sharex=True)
    for fits, cstar, col, name in ((fits_li, cs_li, C_AQUA, "LiCoO$_2$"),
                                   (fits_na, cs_na, C_BLUE, "NaCoO$_2$")):
        axes[0, 0].plot(fits["c"], fits["alpha"], "o-", color=col, ms=3.5,
                        lw=1.4, label=f"{name}, $c^*$={cstar:.2f} Å")
        axes[0, 0].axvline(cstar, color=col, ls=":", lw=1)
        axes[0, 1].plot(fits["c"], fits["beta"], "o-", color=col, ms=3.5, lw=1.4)
        axes[1, 0].plot(fits["c"], fits["d_min"], "o-", color=col, ms=3.5, lw=1.4)
        axes[1, 1].plot(fits["c"], fits["depth_eV"] * 1e3, "o-", color=col,
                        ms=3.5, lw=1.4)
    axes[0, 0].axhline(0, color="k", lw=0.7)
    axes[0, 0].set_ylabel(r"$\alpha$ (eV/Å$^2$)")
    axes[0, 0].legend(fontsize=8)
    axes[0, 1].set_ylabel(r"$\beta$ (eV/Å$^4$)")
    axes[1, 0].set_ylabel(r"well position $d_0$ (Å)")
    axes[1, 1].set_ylabel("well depth (meV)")
    for ax in axes[1]:
        ax.set_xlabel(r"CoO$_2$–CoO$_2$ spacing $c$ (Å)")
    fig.suptitle(r"Landau parameters of $V(d;c)=\alpha d^2+\beta d^4$ (repo DFT)",
                 fontsize=10.5)
    fig.tight_layout()
    fig.savefig(FIGDIR / "fig2_landau_parameters.png", bbox_inches="tight")
    plt.close(fig)


def fig_quantum_wells(fits_na, qt_na_na, qt_na_li, qt_li_li, cs_na, cs_li):
    fig = plt.figure(figsize=(9.6, 6.4))
    gs = fig.add_gridspec(2, 3, height_ratios=[1, 1.15])
    # top row: three representative wells with levels (NaCoO2, Na mass)
    reps = [np.argmin(np.abs(fits_na["c"].values - t))
            for t in (5.85, 6.15, 6.9)]
    dd = np.linspace(-1.6, 1.6, 400)
    for j, irep in enumerate(reps):
        ax = fig.add_subplot(gs[0, j])
        r = fits_na.iloc[irep]
        V = quartic(dd, r["alpha"], r["beta"]) * 1e3
        ax.plot(dd, V, color=C_BLUE, lw=1.6)
        wq = solve_well(r["alpha"], r["beta"], M_NA)
        levels = sorted(list(wq.E_even[:2]) + list(wq.E_odd[:2]))
        for E in levels:
            ax.axhline(E * 1e3, color=C_GRAY, lw=0.8, alpha=0.8,
                       xmin=0.28, xmax=0.72)
        ax.set_title(f"$c$ = {r['c']:.2f} Å", fontsize=9.5)
        ax.set_xlabel("$d$ (Å)")
        lo = min(V.min(), levels[0] * 1e3)
        ax.set_ylim(lo * 1.25 - 8, max(60.0, levels[-1] * 1e3 * 1.6))
        if j == 0:
            ax.set_ylabel("$V$, levels (meV)")
    # bottom left: omega_eff softening
    ax = fig.add_subplot(gs[1, 0])
    ax.plot(qt_na_na["c"], qt_na_na["hw_eff_meV"], "o-", ms=3.5, lw=1.4,
            color=C_BLUE, label="NaCoO$_2$ well, Na mass")
    ax.plot(qt_na_li["c"], qt_na_li["hw_eff_meV"], "s--", ms=3, lw=1.1,
            color=C_VIOLET, label="NaCoO$_2$ well, Li mass")
    ax.plot(qt_li_li["c"] + (cs_na - cs_li), qt_li_li["hw_eff_meV"], "^--",
            ms=3, lw=1.1, color=C_AQUA,
            label="LiCoO$_2$ well, Li mass\n($c$ shifted by $c^*_{Na}-c^*_{Li}$)")
    ax.axvline(cs_na, color=C_GRAY, ls=":", lw=1)
    ax.set_yscale("log")
    ax.set_xlabel(r"$c$ (Å)")
    ax.set_ylabel(r"$\hbar\omega_{\rm eff}=\hbar^2/2M\langle d^2\rangle_0$ (meV)")
    ax.legend(fontsize=7)
    # bottom middle: level gaps
    ax = fig.add_subplot(gs[1, 1])
    ax.semilogy(qt_na_na["c"], qt_na_na["gap02_meV"], "o-", ms=3.5, lw=1.4,
                color=C_BLUE, label=r"$E_2-E_0$ (even channel)")
    ax.semilogy(qt_na_na["c"],
                np.clip(qt_na_na["split01_meV"], 1e-9, None), "s-", ms=3,
                lw=1.1, color=C_YELL, label=r"$E_1-E_0$ (tunneling)")
    ax.axvline(cs_na, color=C_GRAY, ls=":", lw=1)
    ax.set_ylim(1e-6, 3e2)
    ax.set_xlabel(r"$c$ (Å)")
    ax.set_ylabel("level spacing (meV)")
    ax.legend(fontsize=7.5)
    ax.set_title("Na mass, NaCoO$_2$ wells", fontsize=9)
    # bottom right: <d^2>
    ax = fig.add_subplot(gs[1, 2])
    ax.plot(qt_na_na["c"], qt_na_na["d2_A2"], "o-", ms=3.5, lw=1.4,
            color=C_BLUE, label="Na mass")
    ax.plot(qt_na_li["c"], qt_na_li["d2_A2"], "s--", ms=3, lw=1.1,
            color=C_VIOLET, label="Li mass")
    ax.axvline(cs_na, color=C_GRAY, ls=":", lw=1)
    ax.set_xlabel(r"$c$ (Å)")
    ax.set_ylabel(r"$\langle d^2\rangle_0$ (Å$^2$)")
    ax.legend(fontsize=7.5)
    fig.suptitle("Quantized alkali double-well mode (finite differences, exact parity)",
                 fontsize=10.5)
    fig.tight_layout()
    fig.savefig(FIGDIR / "fig3_quantum_wells.png", bbox_inches="tight")
    plt.close(fig)


def fig_charge_map(bader, gal, c_on_li, cs_li, cs_na, dq):
    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.3))
    ax = axes[0]
    ax.plot(bader["c"], bader["q_Li"] * 1e3, "o-", ms=4, lw=1.4, color=C_AQUA)
    ax.axvline(c_on_li, color=C_GRAY, ls="--", lw=1,
               label=f"Bader onset {c_on_li:.2f} Å")
    ax.axvline(cs_li, color=C_AQUA, ls=":", lw=1.2,
               label=f"$c^*_{{Li}}$ = {cs_li:.2f} Å")
    ax.set_xlabel(r"$c$ (Å)  [LiCoO$_2$]")
    ax.set_ylabel(r"Bader $q_{\rm Li}$ (10$^{-3}$ e)")
    ax.set_title("DFT Bader charge on Li", fontsize=9.5)
    ax.legend(fontsize=7.5)
    ax = axes[1]
    cg = np.linspace(5.5, 10.2, 200)
    ax.plot(cg, gal.delta(cg), color=C_BLUE, lw=1.6)
    ax.axhline(GAMMA_BANDS / 2, color=C_GRAY, ls="--", lw=1,
               label=r"$\Gamma/2$ = 2.07 eV")
    ax.axvline(gal.c_on, color=C_GRAY, ls=":", lw=1)
    ax.set_xlabel(r"$c$ (Å)  [Na axis]")
    ax.set_ylabel(r"$\delta(c)$ (eV)")
    ax.set_title("Mapped on-site splitting", fontsize=9.5)
    ax.legend(fontsize=7.5)
    ax = axes[2]
    ax.plot(cg, gal.n_cm2(cg), color=C_BLUE, lw=1.6, label="model $q(c)/A_{cell}$")
    ax.plot(bader["c"] + (cs_na - cs_li), dq / gal.A_cell * 1e16, "o", ms=4,
            color=C_AQUA, label=r"LiCoO$_2$ Bader $\Delta q$ (shifted $c$)")
    ax.set_xlabel(r"$c$ (Å)  [Na axis]")
    ax.set_ylabel(r"$n_{\rm 2DEG}$ (e/cm$^2$)")
    ax.set_title("2DEG density turn-on", fontsize=9.5)
    ax.legend(fontsize=7.5)
    fig.tight_layout()
    fig.savefig(FIGDIR / "fig4_charge_turnon.png", bbox_inches="tight")
    plt.close(fig)


def fig_coupling(res_mid, res_lo, res_hi, cs_na, cmax_dft):
    fig, axes = plt.subplots(1, 3, figsize=(10.8, 3.4))
    ax = axes[0]
    kz = np.linspace(0, np.pi, 60)
    for dlt, col in ((2.05, C_BLUE), (1.5, C_VIOLET), (1.0, C_YELL)):
        ax.plot(kz, [chi_band(dlt, GAMMA_EP, k) for k in kz], lw=1.4, color=col,
                label=rf"$\delta$ = {dlt:.2f} eV")
    ax.set_xlabel(r"$k_z$ (rad)")
    ax.set_ylabel(r"$\chi = \partial^2 E_{alk}/\partial d^2$ (eV/Å$^2$)")
    ax.set_title(rf"TB vertex curvature ($\gamma$ = {GAMMA_EP:.2f}/Å)", fontsize=9.5)
    ax.legend(fontsize=7.5)
    ax = axes[1]
    ax.plot(res_mid.c, res_mid.lam, lw=1.7, color=C_BLUE,
            label=rf"$\lambda$ (even channel), $\gamma$={GAMMA_EP:.2f}")
    ax.fill_between(res_mid.c, res_lo.lam, res_hi.lam, color=C_BLUE, alpha=0.18,
                    lw=0, label=r"$\gamma \in [1,2]$ /Å")
    ax.plot(res_mid.c, res_mid.lam_task, ls="--", lw=1.2, color=C_YELL,
            label=r"$\lambda$ (task rms recipe)")
    ax.axhline(LAMBDA_LOC, color=C_RED, ls=":", lw=1.2,
               label=rf"polaron criterion $\lambda$={LAMBDA_LOC:g}")
    ax.axvline(cs_na, color=C_GRAY, ls=":", lw=1)
    ax.set_yscale("log")
    ax.set_ylim(1e-3, 3e3)
    ax.set_xlabel(r"$c$ (Å)")
    ax.set_ylabel(r"$\lambda$")
    ax.legend(fontsize=6.8)
    ax = axes[2]
    ax.plot(res_mid.c, res_mid.hw_02 * 1e3, lw=1.7, color=C_BLUE,
            label=r"$\hbar\Omega = E_2-E_0$")
    ax.plot(res_mid.c, res_mid.hw_eff * 1e3, ls="--", lw=1.2, color=C_VIOLET,
            label=r"$\hbar\omega_{\rm eff}$")
    ax.axvline(cs_na, color=C_GRAY, ls=":", lw=1)
    ax.axvspan(cmax_dft, res_mid.c.max(), color=C_MUTED, alpha=0.12, lw=0)
    ax.set_yscale("log")
    ax.set_xlabel(r"$c$ (Å)")
    ax.set_ylabel("mode energy (meV)")
    ax.legend(fontsize=7.5)
    fig.tight_layout()
    fig.savefig(FIGDIR / "fig5_coupling.png", bbox_inches="tight")
    plt.close(fig)


def money_plot(fits_na, cs_na, gal, res_mid, res_lo, res_hi,
               bader, cs_li, cmax_dft, dq):
    fig, axes = plt.subplots(3, 1, figsize=(6.8, 8.6), sharex=True,
                             gridspec_kw={"hspace": 0.08})
    C_ANH, C_MONO, C_BI = 5.55, 6.9, 9.9   # experimental anchors (spacings, A)

    for ax in axes:
        ax.axvspan(cmax_dft, 10.35, color=C_MUTED, alpha=0.10, lw=0)
        ax.axvline(C_ANH, color=C_GRAY, ls="--", lw=0.9, alpha=0.7)
        ax.axvline(C_MONO, color=C_GRAY, ls="--", lw=1.1)
        ax.axvline(C_BI, color=C_RED, ls="--", lw=1.1)
        ax.axvline(cs_na, color=C_BLUE, ls=":", lw=1.1)

    # (i) alpha(c)
    ax = axes[0]
    ax.axhline(0, color="k", lw=0.7)
    ax.plot(res_mid.c, res_mid.alpha, lw=1.6, color=C_BLUE)
    ax.plot(fits_na["c"], fits_na["alpha"], "o", ms=4.5, color=C_BLUE,
            mfc="white", mew=1.2, label=r"NaCoO$_2$ DFT (this repo)")
    ax.text(cs_na + 0.08, 0.72, rf"$c^*$ = {cs_na:.2f} Å", fontsize=8.5,
            color=C_BLUE)
    ax.set_ylabel(r"$\alpha(c)$  (eV/Å$^2$)")
    ax.legend(fontsize=8, loc="upper right")
    ax.text(0.985, 0.06,
            "double well opens:\nsoft anharmonic Na mode",
            transform=ax.transAxes, ha="right", fontsize=7.5, color=C_GRAY)

    # (ii) n_2DEG(c)
    ax = axes[1]
    ax.plot(res_mid.c, res_mid.n_cm2 / 1e13, lw=1.6, color=C_BLUE,
            label=r"model: $q_{alk}=Q\,(\Gamma-2\delta)/2\Gamma$, $\delta(c)$ from Bader")
    ax.plot(bader["c"] + (cs_na - cs_li), dq / gal.A_cell * 1e16 / 1e13,
            "o", ms=4.5, color=C_AQUA, mfc="white", mew=1.2,
            label=r"LiCoO$_2$ Bader $\Delta q_{\rm Li}$ (axis shifted by $c^*_{Na}-c^*_{Li}$)")
    ax.set_ylabel(r"$n_{\rm 2DEG}$  ($10^{13}$ e/cm$^2$)")
    ax.legend(fontsize=7.5, loc="upper left")
    ax.text(0.985, 0.06, "interlayer 2DEG turns on at $c^*$",
            transform=ax.transAxes, ha="right", fontsize=7.5, color=C_GRAY)

    # (iii) Tc(c)
    ax = axes[2]
    def valid_tc(res):
        ok = (res.lam <= LAMBDA_LOC) & (res.q > 0)
        return np.where(ok, res.tc, np.nan)
    ax.plot(res_mid.c, valid_tc(res_mid), lw=2.0, color=C_BLUE, solid_capstyle="round",
            label=rf"Allen–Dynes $T_c$, $\gamma$={GAMMA_EP:.2f} Å$^{{-1}}$ ($z_0$={Z0_ANSATZ} Å)")
    ax.plot(res_lo.c, valid_tc(res_lo), lw=1.1, ls="-.", color=C_AQUA,
            label=r"$\gamma$=1.0 Å$^{-1}$ ($z_0$=1.0 Å)")
    ax.plot(res_hi.c, valid_tc(res_hi), lw=1.1, ls=":", color=C_VIOLET,
            label=r"$\gamma$=2.0 Å$^{-1}$ ($z_0$=0.5 Å)")
    pol = np.where((res_mid.lam > LAMBDA_LOC) & (res_mid.q > 0), res_mid.tc,
                   np.nan)
    ax.plot(res_mid.c, pol, lw=1.0, ls="--", color=C_YELL, alpha=0.85,
            label=rf"raw AD $T_c$ where $\lambda>{LAMBDA_LOC:g}$: polaronic,"
                  "\nMigdal invalid — SC expected destroyed")
    ax.plot([C_BI], [4.5], "*", ms=14, color=C_RED, zorder=6,
            label=r"expt: bilayer hydrate $T_c$ = 4.5 K")
    ax.set_ylabel(r"$T_c$  (K)")
    ax.set_xlabel(r"CoO$_2$–CoO$_2$ interlayer spacing $c$  (Å)")
    ax.set_xlim(5.45, 10.35)
    ax.set_ylim(0, 48)
    ax.legend(fontsize=6.8, loc="upper right")
    ax.annotate("Na freezes off-center:\nmode spectral weight collapses",
                (7.6, 6), fontsize=7.5, color=C_GRAY, ha="center")

    # anchor labels on the top panel
    ymax = axes[0].get_ylim()[1]
    axes[0].text(C_ANH + 0.12, ymax * 1.06, "anhydrous\n(5.5 Å, no SC)",
                 fontsize=7.5, ha="left", color=C_GRAY)
    axes[0].text(C_MONO, ymax * 1.06, "monolayer hydrate\n(6.9 Å, no SC)",
                 fontsize=7.5, ha="center", color=C_GRAY)
    axes[0].text(C_BI, ymax * 1.06, "bilayer hydrate\n(9.9 Å, SC 4.5 K)",
                 fontsize=7.5, ha="center", color=C_RED)

    fig.suptitle("Gallery-controlled Na double well and interlayer-2DEG "
                 "superconductivity", y=0.99, fontsize=11)
    fig.text(0.12, 0.002,
             "Model caveat: $c$ is the CoO$_2$–CoO$_2$ spacing of the anhydrous Na(Li)CoO$_2$ DFT series (this repo);\n"
             "the hydrate $c$ values are placed on this axis as a rigid extrapolation (shaded beyond the DFT range).\n"
             "The hydrate galleries contain H$_2$O, absent from the model; water screening is expected to shallow\n"
             "the Na well and pull the system back toward the critical point (see MODEL.md).",
             fontsize=7, color=C_GRAY, va="bottom")
    fig.savefig(HERE / "money_plot.pdf", bbox_inches="tight")
    fig.savefig(HERE / "money_plot.png", bbox_inches="tight")
    plt.close(fig)


# ======================================================================
# MAIN
# ======================================================================
def main():
    print("=" * 72)
    print("1. DATA EXTRACTION")
    scans_li, a_li = load_scans("licoo2")
    scans_na, a_na = load_scans("nacoo2")
    bader = load_bader()
    print(f"   LiCoO2: {len(scans_li)} spacings c = {scans_li[0].c:.3f}..."
          f"{scans_li[-1].c:.3f} A, a = {a_li:.4f} A")
    print(f"   NaCoO2: {len(scans_na)} spacings c = {scans_na[0].c:.3f}..."
          f"{scans_na[-1].c:.3f} A, a = {a_na:.4f} A")
    print(f"   Bader (LiCoO2, GPAW): {len(bader)} points, q_Li "
          f"{bader['q_Li'].min():.4f}...{bader['q_Li'].max():.4f} e")

    print("=" * 72)
    print("2. QUARTIC FITS")
    fits_li = fit_wells(scans_li)
    fits_na = fit_wells(scans_na)
    cs_li, a0_li = critical_spacing(fits_li)
    cs_na, a0_na = critical_spacing(fits_na)
    fits_li.to_csv(RESDIR / "well_fits_licoo2.csv", index=False)
    fits_na.to_csv(RESDIR / "well_fits_nacoo2.csv", index=False)
    print(f"   c*_Li = {cs_li:.3f} A (slope {a0_li:.2f} eV/A^2/A), "
          f"c*_Na = {cs_na:.3f} A (slope {a0_na:.2f})")
    print(f"   NaCoO2 deepest well: {fits_na['depth_eV'].max()*1e3:.0f} meV at "
          f"d0 = {fits_na['d_min'].iloc[-1]:.2f} A (c = {fits_na['c'].iloc[-1]:.2f} A)")

    print("=" * 72)
    print("3. QUANTUM WELLS")
    qt_li_li = quantum_table(fits_li, M_LI, "Li")
    qt_li_na = quantum_table(fits_li, M_NA, "Na")
    qt_na_na = quantum_table(fits_na, M_NA, "Na")
    qt_na_li = quantum_table(fits_na, M_LI, "Li")
    pd.concat([qt_li_li, qt_li_na]).to_csv(RESDIR / "quantum_wells_licoo2.csv",
                                           index=False)
    pd.concat([qt_na_na, qt_na_li]).to_csv(RESDIR / "quantum_wells_nacoo2.csv",
                                           index=False)
    for lbl, qt in (("NaCoO2/Na", qt_na_na),):
        print(f"   {lbl}: hbar*w_eff {qt['hw_eff_meV'].max():.1f} -> "
              f"{qt['hw_eff_meV'].min():.3f} meV; tunneling splitting at "
              f"largest c: {qt['split01_meV'].iloc[-1]:.2e} meV")

    print("=" * 72)
    print("4. SSH4 MODEL")
    # parity check: linear deformation potential at d=0
    for kz in (0.0, 1.0, np.pi):
        w0 = np.linalg.eigvalsh(h4(kz, 2.0, d=1e-4))[ALKALI_BAND]
        w1 = np.linalg.eigvalsh(h4(kz, 2.0, d=-1e-4))[ALKALI_BAND]
        assert abs(w0 - w1) / 2e-4 < 1e-6, "linear DP should vanish by parity"
    print("   linear deformation potential dE/dd|_{d=0} = 0 (parity) -- verified")
    gal, c_on_li, slope, dq = calibrate_gallery(bader, cs_li, cs_na, a_na)
    mstar, n0_area, n0_cell = alkali_mass_and_dos(float(gal.delta(cs_na)), a_na)
    print(f"   Bader onset (Li axis): {c_on_li:.3f} A  [c*_Li = {cs_li:.3f} A]"
          f"  slope dq/dc = {slope:.4f} e/A -> w = {gal.w:.1f} A")
    print(f"   alkali band: m* = {mstar:.3f} m_e, N(0) = {n0_area:.4f} /(eV A^2 spin)"
          f" = {n0_cell:.3f} /(eV cell spin)")
    print(f"   band splitting of alkali sp_z pair: 2*tau1 = {2*TAU1:.1f} eV; "
          f"chi(c*) = {chi_band(float(gal.delta(cs_na))):.2f} eV/A^2 "
          f"(gamma = {GAMMA_EP:.2f} /A)")

    print("=" * 72)
    print("5. COUPLING AND Tc")
    cmax_dft = float(fits_na["c"].max())
    cgrid = np.linspace(5.55, 10.3, 240)
    res_mid = coupling_vs_c(cgrid, fits_na, gal, M_NA, n0_cell, GAMMA_EP)
    res_lo = coupling_vs_c(cgrid, fits_na, gal, M_NA, n0_cell, 1.0 / Z0_RANGE[1])
    res_hi = coupling_vs_c(cgrid, fits_na, gal, M_NA, n0_cell, 1.0 / Z0_RANGE[0])
    pd.DataFrame({
        "c_A": res_mid.c, "alpha": res_mid.alpha, "beta": res_mid.beta,
        "q_e_per_cell": res_mid.q, "n_cm2": res_mid.n_cm2,
        "hw_eff_meV": res_mid.hw_eff * 1e3, "hw02_meV": res_mid.hw_02 * 1e3,
        "split01_meV": res_mid.split01 * 1e3, "d2_A2": res_mid.d2,
        "chi_eV_A2": res_mid.chi, "lambda": res_mid.lam,
        "lambda_task": res_mid.lam_task, "Tc_K": res_mid.tc,
        "lambda_g1": res_lo.lam, "Tc_K_g1": res_lo.tc,
        "lambda_g2": res_hi.lam, "Tc_K_g2": res_hi.tc,
    }).to_csv(RESDIR / "coupling_tc_vs_c.csv", index=False)

    sc_ok = (res_mid.lam <= LAMBDA_LOC) & (res_mid.q > 0)
    if sc_ok.any():
        i_pk = int(np.nanargmax(np.where(sc_ok, res_mid.tc, np.nan)))
        print(f"   SC window (q>0 & lambda<={LAMBDA_LOC:g}): c = "
              f"{res_mid.c[sc_ok].min():.2f}...{res_mid.c[sc_ok].max():.2f} A")
        print(f"   peak Tc = {res_mid.tc[i_pk]:.2f} K at c = {res_mid.c[i_pk]:.2f} A "
              f"(lambda = {res_mid.lam[i_pk]:.2f}, hbar*Omega = "
              f"{res_mid.hw_02[i_pk]*1e3:.1f} meV)")
    else:
        print("   NO SC window within the model -- reported as negative result")
    for gname, res in (("gamma=1.43", res_mid), ("gamma=1.0", res_lo),
                       ("gamma=2.0", res_hi)):
        ok = (res.lam <= LAMBDA_LOC) & (res.q > 0)
        if ok.any():
            i_pk = int(np.nanargmax(np.where(ok, res.tc, np.nan)))
            print(f"   [{gname}] SC window {res.c[ok].min():.2f}..."
                  f"{res.c[ok].max():.2f} A, peak Tc = {res.tc[i_pk]:.1f} K "
                  f"at {res.c[i_pk]:.2f} A")
        else:
            print(f"   [{gname}] no SC window")
    for target in (6.9, 9.9):
        j = int(np.argmin(np.abs(res_mid.c - target)))
        print(f"   at c = {target} A: lambda = {res_mid.lam[j]:.2f} "
              f"(gamma range {res_lo.lam[j]:.2f}..{res_hi.lam[j]:.2f}), "
              f"raw AD Tc = {res_mid.tc[j]:.2f} K, "
              f"{'POLARONIC (lambda>2, SC destroyed)' if res_mid.lam[j] > LAMBDA_LOC else 'SC allowed'}")

    print("=" * 72)
    print("6. FIGURES")
    fig_potential_fits(scans_li, fits_li, scans_na, fits_na)
    fig_landau_parameters(fits_li, cs_li, fits_na, cs_na)
    fig_quantum_wells(fits_na, qt_na_na, qt_na_li, qt_li_li, cs_na, cs_li)
    fig_charge_map(bader, gal, c_on_li, cs_li, cs_na, dq)
    fig_coupling(res_mid, res_lo, res_hi, cs_na, cmax_dft)
    money_plot(fits_na, cs_na, gal, res_mid, res_lo, res_hi, bader, cs_li,
               cmax_dft, dq)
    print(f"   figures -> {FIGDIR}")
    print(f"   money plot -> {HERE / 'money_plot.pdf'} (+ .png)")
    print("done.")


if __name__ == "__main__":
    main()

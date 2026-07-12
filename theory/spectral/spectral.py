#!/usr/bin/env python
"""
ARPES-comparable spectral analysis for the Na_xCoO2 band chains.

Read-only parsers over runpod/results_bands (band + projwfc "filproj" output)
and runpod/results_v4/jobs/nscf_* (dense uniform-grid eigenvalues).  Nothing
outside theory/spectral/ and paper/figures/{fig7,fig8}.* is written.

filproj format (pw.proj.projwfc_{up,down}), one run per spin:
  header: nr1x.. nat ntyp / ibrav celldm / vol nelec nkstot .. / ntyp species
          lines / nat atom lines / (natomwfc nkstot nbnd) / (noncolin lspinorb)
  then, per atomic wavefunction:
      "<iw> <atom#> <symbol> <orblabel> <n> <l> <m>"
      followed by nkstot*nbnd lines "<ik> <ibnd> <|proj|^2>"
So the Na-orbital projected weight per (band,k) is the sum over the atomic
wavefunctions carrying symbol "Na" of |<phi_Na|psi_{k,ibnd}>|^2.
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
BANDS = REPO / "runpod" / "results_bands"
NSCF = REPO / "runpod" / "results_v4" / "jobs"

# k-path Gamma-M-K-Gamma, crystal_b, 60 pts, per-vertex counts 22/12/25/1.
# Vertices therefore land at cumulative indices 0, 22, 34, 59.
KPATH_LABELS = [r"$\Gamma$", "M", "K", r"$\Gamma$"]
KPATH_IDX = [0, 22, 34, 59]


# ----------------------------------------------------------------------
def fermi_energy(job: str) -> float:
    """Fermi energy (eV) from the bands-run pw.out."""
    txt = (BANDS / job / "pw.out").read_text()
    m = re.search(r"the Fermi energy is\s+(-?[\d.]+)\s*ev", txt)
    if not m:
        raise ValueError(f"no Fermi energy in {job}")
    return float(m.group(1))


def total_magnetization(job: str) -> float:
    """Last 'total magnetization' (Bohr mag/cell) from the bands-run pw.out."""
    txt = (BANDS / job / "pw.out").read_text()
    hits = re.findall(r"total magnetization\s+=\s+(-?[\d.]+)\s+Bohr", txt)
    return float(hits[-1]) if hits else float("nan")


def load_bands_gnu(path: Path):
    """Parse a pw.bands.*.dat.gnu file.

    Returns (kdist[nk], E[nband, nk]) in eV.  Blocks (one per band) are
    separated by blank lines; every block lists the same nk k-distances.
    """
    blocks, cur = [], []
    for ln in path.read_text().splitlines():
        if ln.strip() == "":
            if cur:
                blocks.append(cur)
                cur = []
            continue
        k, e = ln.split()
        cur.append((float(k), float(e)))
    if cur:
        blocks.append(cur)
    nk = len(blocks[0])
    kdist = np.array([p[0] for p in blocks[0]])
    E = np.array([[p[1] for p in blk] for blk in blocks])
    assert all(len(b) == nk for b in blocks), "ragged band blocks"
    return kdist, E


_HDR = re.compile(
    r"^\s*\d+\s+\d+\s+([A-Za-z]+)\s+(\w+)\s+\d+\s+(\d+)\s+\d+\s*$")


def load_projection(path: Path, nband: int, nk: int):  # noqa: C901
    """Parse a filproj pw.proj.projwfc_* file into per-species weights.

    Returns dict of arrays shape (nband, nk):
      'Na', 'Co_d', 'O_p', 'tot' (sum over all atomic wfc, ~ orbital-resolved
      completeness, <= 1 per state).
    """
    W = {k: np.zeros((nband, nk)) for k in ("Na", "Co_d", "O_p", "tot")}
    sym = orb_l = None
    target = None
    for ln in path.read_text().splitlines():
        m = _HDR.match(ln)
        if m:
            sym, orb, orb_l = m.group(1), m.group(2), int(m.group(3))
            # which accumulators does this atomic wfc feed?
            target = []
            if sym == "Na":
                target.append("Na")
            if sym == "Co" and orb_l == 2:
                target.append("Co_d")
            if sym == "O" and orb_l == 1:
                target.append("O_p")
            continue
        p = ln.split()
        if len(p) != 3 or sym is None:
            continue
        try:
            ik, ib = int(p[0]), int(p[1])
            val = float(p[2])
        except ValueError:
            continue
        # nspin=2 filproj numbers spin-down k-points 61..120; fold to 0..nk-1
        ki = (ik - 1) % nk
        if not (0 <= ib - 1 < nband):
            continue
        W["tot"][ib - 1, ki] += val
        for t in target:
            W[t][ib - 1, ki] += val
    return W


def load_spinsummed(job: str):
    """Band chain summed over spin, aligned to E_F.

    Returns dict:
      kdist[nk], ticks (positions of Gamma,M,K,Gamma),
      E[2*nband, nk]  (eV, already E - E_F, both spins stacked),
      Na, Co_d, O_p   (matching weights, spin-summed per state),
      ef, mag
    """
    ef = fermi_energy(job)
    d = BANDS / job
    ku, Eu = load_bands_gnu(d / "pw.bands.up.dat.gnu")
    kd, Ed = load_bands_gnu(d / "pw.bands.dw.dat.gnu")
    nb, nk = Eu.shape
    Wu = load_projection(d / "pw.proj.projwfc_up", nb, nk)
    Wd = load_projection(d / "pw.proj.projwfc_down", nb, nk)
    # stack the two spins as independent (band,k) sheets
    E = np.vstack([Eu, Ed]) - ef
    Na = np.vstack([Wu["Na"], Wd["Na"]])
    Cod = np.vstack([Wu["Co_d"], Wd["Co_d"]])
    Op = np.vstack([Wu["O_p"], Wd["O_p"]])
    ticks = [ku[i] for i in KPATH_IDX]
    return dict(kdist=ku, ticks=ticks, E=E, Na=Na, Co_d=Cod, O_p=Op,
                ef=ef, mag=total_magnetization(job))


def spectral_image(kdist, E, weight, kfine, omega, gamma):
    """ARPES-like intensity I(omega, k) on a (len(omega), len(kfine)) grid.

    Each band is linearly interpolated in k onto kfine, then Lorentzian
    broadened in omega with HWHM gamma and multiplied by its weight.
    """
    img = np.zeros((omega.size, kfine.size))
    g2 = gamma * gamma
    for b in range(E.shape[0]):
        eb = np.interp(kfine, kdist, E[b])
        wb = np.interp(kfine, kdist, weight[b])
        if wb.max() <= 0 and weight is not None and np.all(weight[b] == 0):
            # still contributes to total intensity if weight is all-ones proxy
            pass
        # Lorentzian in omega for each fine-k column
        dE = omega[:, None] - eb[None, :]           # (nomega, nkfine)
        lor = (gamma / np.pi) / (dE * dE + g2)
        img += lor * wb[None, :]
    return img


# ======================================================================
# Dense NSCF parsing for the Fermi-surface map
# ======================================================================
_KLINE = re.compile(r"^\s*k =\s*(-?[\d.]+)\s+(-?[\d.]+)\s+(-?[\d.]+)")


def load_nscf_kz0(job: str, kz_tol=1e-3):
    """Parse nscf pw.out; return kz~0 k-points (2pi/alat cart) + eigenvalues.

    Returns (kxy[N,2], Ebands[N, nband] in eV, ef).  Spin-up and spin-down
    blocks are both read and treated as independent sheets (concatenated),
    so the map is spin-summed.
    """
    txt = (NSCF / job / "pw.out").read_text()
    ef = float(re.search(r"the Fermi energy is\s+(-?[\d.]+)\s*ev", txt).group(1))
    lines = txt.splitlines()
    kxy, Eb = [], []
    i = 0
    while i < len(lines):
        m = _KLINE.match(lines[i])
        if m and "bands (ev)" in lines[i]:
            kx, ky, kz = float(m.group(1)), float(m.group(2)), float(m.group(3))
            # collect eigenvalue floats until a blank line
            vals, j = [], i + 2
            while j < len(lines) and lines[j].strip() != "":
                vals += [float(x) for x in lines[j].split()]
                j += 1
            if abs(kz) < kz_tol:
                kxy.append((kx, ky))
                Eb.append(vals)
            i = j
        else:
            i += 1
    return np.array(kxy), np.array(Eb), ef


# hexagonal reciprocal lattice (2pi/alat cart) from the nscf pw.out header
B1 = np.array([1.0, 0.577350])
B2 = np.array([0.0, 1.154701])


def hex_symmetrize(kxy, Eb):
    """Apply the 12 C6v operations (in cartesian kx-ky) to fill the BZ plane.

    The irreducible kz=0 wedge is expanded to the full zone so the Fermi
    surface reads as the hexagon ARPES measures.  Returns (K[M,2], E[M,nb]).
    """
    th = np.pi / 3
    rots = [np.array([[np.cos(n * th), -np.sin(n * th)],
                      [np.sin(n * th), np.cos(n * th)]]) for n in range(6)]
    mir = np.array([[1.0, 0.0], [0.0, -1.0]])          # a mirror
    ops = rots + [r @ mir for r in rots]
    Ks, Es = [], []
    for op in ops:
        Ks.append(kxy @ op.T)
        Es.append(Eb)
    K = np.vstack(Ks)
    E = np.vstack(Es)
    # de-duplicate near-coincident images
    key = np.round(K / 1e-4).astype(np.int64)
    _, idx = np.unique(key, axis=0, return_index=True)
    return K[idx], E[idx]


if __name__ == "__main__":
    for job in ("bands_Na_1x1_c5.5", "bands_Na_1x1_c6.9", "bands_Na_1x1_c9.9"):
        d = load_spinsummed(job)
        print(f"{job}: E_F={d['ef']:.3f} eV  mag={d['mag']:+.2f}  "
              f"E{d['E'].shape}  Na max={d['Na'].max():.3f}  "
              f"tot mean={d['E'].shape}")
        # Na weight of states within +/-0.5 eV of E_F
        near = np.abs(d["E"]) < 0.5
        print(f"    states within 0.5 eV of E_F: {near.sum()}, "
              f"max Na weight among them = "
              f"{(d['Na'][near].max() if near.any() else 0):.3f}")

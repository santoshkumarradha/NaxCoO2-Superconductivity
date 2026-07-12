#!/usr/bin/env python3
"""Generate Quantum ESPRESSO pw.x inputs for the Na_xCoO2 interlayer-2DEG /
anharmonic-rattler study.

Default stage writes all SCF jobs (sets A and B) into jobs/<name>/pw.in and a
manifest.json.  After the SCF phase has run, `--stage nscf` parses the SCF
energies, finds the E(delta) minimum per (element, cell, c) and writes the
dense-k NSCF + dos.x jobs (set C), appending them to the manifest.

`--stage extra` writes the follow-up job sets E-H into jobs_extra/<name>/pw.in
and manifest_extra.json (run with run_extra.sh):
  E: K_xCoO2 1x1 E(delta) scans (mirror of set A, element K)      24 SCF
  F: gated LiCoO2, tot_charge jellium proxy for IL gating          12 SCF
  G: Na_1/3CoO2.yH2O explicit-water bilayer + vacuum reference     10 SCF
  H: Na2CoSe2O gallery-band check (SCF + dense-k NSCF)              2 jobs

`--stage bands` writes the referee-response job set I into jobs_bands/ and
manifest_bands.json (run with run_bands.sh):
  I1: Na fatband chains, 1x1 c = 5.5/6.9/9.9 + s3 c = 5.5/9.9 (delta = 0):
      pw scf -> pw 'bands' (Gamma-M-K-Gamma, 60 pts) -> bands.x (up/dw)
      -> projwfc.x (all-atom projections: Na s/p AND Co d)      5 chains
  I2: hydrate relaxed-water check, c = 9.9, delta = 0.00/0.30:
      BFGS relax, Na z + all H2O free, CoO2 frozen (if_pos)      2 relax
  I3: hydrate vdW spot-check, c = 9.9, delta = 0.00/0.30/0.50:
      SCF + vdw_corr = 'grimme-d3', compare E(delta) vs set G    3 SCF

`--stage mobile` writes the manuscript-TODO job set J into jobs_mobile/ and
manifest_mobile.json (run with run_mobile.sh): the mobile-water ADIABATIC
E(delta).  All jobs are hydrate BFGS relaxations (c = 9.9, set-G ansatz):
  J1: Na fully frozen at delta = 0.00-1.00 A (6 pts), CoO2 frozen,
      all 12 H2O atoms free -> E_mobile(delta)                   6 relax
  J2: J1 twin + vdw_corr = 'grimme-d3' (dispersion sensitivity)  6 relax
  J3: joint-minimum search at delta = 0.30: waters free AND Na z
      free (x,y pinned) -- is the joint minimum off-center?      1 relax

`--stage zscan` writes set K into jobs_zscan/ and manifest_zscan.json (run
with run_zscan.sh): vacuum sqrt3 x sqrt3 Na_1/3CoO2 E(delta) scans at fixed
c = 9.9 A and z_O = 0.90/1.02 A                              10 SCF

stdlib + numpy only.
"""
import argparse
import json
import os
import re
import sys

import numpy as np

BOHR = 0.529177210903  # Angstrom

# ---------------------------------------------------------------- physics ---
# K in-plane a = 2.90 A is an ANSATZ: K_xCoO2 (P2/P3 phases) has a ~ 2.84-2.87
# A experimentally; 2.90 A is used here as a round-number upper estimate for
# the larger K+ ion, mirroring the Na scans (task spec).
A_LAT = {"Na": 2.888, "Li": 2.82, "K": 2.90}  # in-plane lattice constant, Ang
C_LIST = [5.5, 6.9, 8.4, 9.9]         # CoO2 plane-to-plane spacing, Angstrom
DELTAS_A = [0.0, 0.15, 0.30, 0.50, 0.75, 1.00]  # off-center displacement, Ang
DELTAS_B = [0.0, 0.25, 0.5, 0.75]
Z_O = 0.96                             # O height above Co plane, Angstrom

PSEUDOS = {
    # NOTE: Co.pbe-spn-kjpaw_psl.1.0.0.UPF is NOT hosted on
    # pseudopotentials.quantum-espresso.org; the closest available pslibrary
    # PAW variant (same 3s3p semicore spn flavour) is used instead.
    "Co": "Co.pbe-spn-kjpaw_psl.0.3.1.UPF",
    "O": "O.pbe-n-kjpaw_psl.1.0.0.UPF",
    "Na": "Na.pbe-spn-kjpaw_psl.1.0.0.UPF",
    "Li": "Li.pbe-s-kjpaw_psl.1.0.0.UPF",
    # Extra-stage species (sets E-H). Filenames verified to exist upstream
    # (HTTP 200 on pseudopotentials.quantum-espresso.org/upf_files, 2026-07-12).
    # Semicore flavours chosen to match the Co.spn / Na.spn convention:
    #   K  spn PAW, z_val = 9  (3s 3p 4s)
    #   Se dn  PAW, z_val = 16 (3d 4s 4p)
    #   H      PAW, z_val = 1
    "K": "K.pbe-spn-kjpaw_psl.1.0.0.UPF",
    "Se": "Se.pbe-dn-kjpaw_psl.1.0.0.UPF",
    "H": "H.pbe-kjpaw_psl.1.0.0.UPF",
}
MASS = {"Co": 58.933, "O": 15.999, "Na": 22.98977, "Li": 6.941,
        "K": 39.0983, "Se": 78.971, "H": 1.008}
ZVAL = {"Co": 17, "O": 6, "Na": 9, "Li": 3,
        "K": 9, "Se": 16, "H": 1}  # PAW valence electrons

KPTS_1X1 = (8, 8, 3)
KPTS_S3 = (6, 6, 3)
KPTS_NSCF = (24, 24, 8)
KPTS_G = (6, 6, 2)                     # set G: gamma-centered 6x6x2


def nbnd_for(nelec):
    """Filled bands + ~8-16 empty ones so the alkali band is always printed."""
    return nelec // 2 + max(8, nelec // 6)


# ------------------------------------------------------------- structures ---
def atoms_1x1(element, c, delta):
    zo = Z_O / c
    za = (c / 2.0 + delta) / c
    return [("Co", 0.0, 0.0, 0.0),
            ("O", 1 / 3, 2 / 3, zo),
            ("O", 2 / 3, 1 / 3, 1.0 - zo),
            (element, 0.0, 0.0, za)]  # Co-top site (lowest-energy, SciPost Fig.5 pos.1)


def atoms_s3(element, c, delta, zo=Z_O):
    """sqrt3 x sqrt3 R30 supercell, 1 alkali per 3 Co (x = 1/3), 10 atoms.

    Primitive fractional (u,v) maps to supercell fractional
    ((u+v)/3, (2v-u)/3); Co sublattice images are (0,0), (1/3,2/3), (2/3,1/3).
    """
    zo = zo / c
    za = (c / 2.0 + delta) / c
    offs = [(0.0, 0.0), (1 / 3, 2 / 3), (2 / 3, 1 / 3)]
    atoms = []
    for ox, oy in offs:
        atoms.append(("Co", ox % 1.0, oy % 1.0, 0.0))
    for base, z in (((1 / 3, 1 / 3), zo), ((1 / 3, 0.0), 1.0 - zo)):
        for ox, oy in offs:
            atoms.append(("O", (base[0] + ox) % 1.0, (base[1] + oy) % 1.0, z))
    atoms.append((element, 0.0, 0.0, za))  # Co-top site
    return atoms


# ------------------------------------------------------------ input files ---
def pw_input(calc, element, a, c, atoms, kpts, prefix="pw", outdir="./out",
             species=None, extra_system=None, notes=None, extra_control=None,
             ions=False, if_pos=None, kpath=None):
    """Hexagonal (ibrav=4) pw.x input.

    species      : explicit ATOMIC_SPECIES order (default Co, O, element)
    extra_system : extra lines appended inside &SYSTEM (e.g. tot_charge)
    notes        : '!'-comment lines placed before ATOMIC_POSITIONS
                   (pw.x card parser skips lines starting with '!'/'#')
    extra_control: extra lines appended inside &CONTROL (e.g. nstep)
    ions         : append an &IONS namelist (BFGS) for calc = 'relax'
    if_pos       : per-atom (ix,iy,iz) 0/1 movability flags (relax only)
    kpath        : list of (kx, ky, kz, npts) -> K_POINTS crystal_b card
                   (band-structure path) instead of an automatic grid; the
                   `kpts` argument is ignored (pass None)
    """
    if species is None:
        species = ["Co", "O", element]
    nelec = sum(ZVAL[s] for s, *_ in atoms)
    nbnd = nbnd_for(nelec)
    lines = [
        "&CONTROL",
        f"  calculation = '{calc}'",
        f"  prefix = '{prefix}'",
        f"  outdir = '{outdir}'",
        "  pseudo_dir = '../../pseudo'",
        "  verbosity = 'high'",
        "  tprnfor = .true.",
    ]
    if extra_control:
        lines.extend(extra_control)
    lines += [
        "/",
        "&SYSTEM",
        "  ibrav = 4",
        f"  celldm(1) = {a / BOHR:.8f}",
        f"  celldm(3) = {c / a:.8f}",
        f"  nat = {len(atoms)}",
        f"  ntyp = {len(species)}",
        f"  nbnd = {nbnd}",
        "  ecutwfc = 60.0",
        "  ecutrho = 480.0",
        "  occupations = 'smearing'",
        "  smearing = 'mv'",
        "  degauss = 0.02",
        "  nspin = 2",
        "  starting_magnetization(1) = 0.3",
    ]
    if extra_system:
        lines.extend(extra_system)
    lines += [
        "/",
        "&ELECTRONS",
        "  conv_thr = 1.0d-7",
        "  mixing_beta = 0.3",
        "  mixing_mode = 'local-TF'",
        "  electron_maxstep = 300",
        "/",
    ]
    if ions:
        lines += ["&IONS", "  ion_dynamics = 'bfgs'", "/"]
    lines.append("ATOMIC_SPECIES")
    for s in species:
        lines.append(f"  {s}  {MASS[s]:.5f}  {PSEUDOS[s]}")
    if notes:
        lines.extend(notes)
    lines.append("ATOMIC_POSITIONS crystal")
    for i, (s, x, y, z) in enumerate(atoms):
        flags = f"  {if_pos[i][0]} {if_pos[i][1]} {if_pos[i][2]}" \
            if if_pos else ""
        lines.append(f"  {s}  {x:.10f}  {y:.10f}  {z:.10f}{flags}")
    if kpath:
        lines.append("K_POINTS crystal_b")
        lines.append(f"  {len(kpath)}")
        for kx, ky, kz, n in kpath:
            lines.append(f"  {kx:.10f}  {ky:.10f}  {kz:.10f}  {n}")
    else:
        lines.append("K_POINTS automatic")
        lines.append(f"  {kpts[0]} {kpts[1]} {kpts[2]} 0 0 0")
    lines.append("")
    return "\n".join(lines)


def dos_input(outdir):
    return "\n".join([
        "&DOS",
        "  prefix = 'pw'",
        f"  outdir = '{outdir}'",
        "  degauss = 0.02",
        "  DeltaE = 0.01",
        "  fildos = 'pw.dos'",
        "/",
        "",
    ])


def write_job(root, name, text, extra=None, jobs_dir="jobs"):
    d = os.path.join(root, jobs_dir, name)
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, "pw.in"), "w") as f:
        f.write(text)
    if extra:
        for fname, ftext in extra.items():
            with open(os.path.join(d, fname), "w") as f:
                f.write(ftext)


# ------------------------------------------------------------- SCF stage ----
def generate_scf(root):
    jobs = []
    # Set A: 1x1 E(delta) scans, Na and Li
    for el in ("Na", "Li"):
        a = A_LAT[el]
        for c in C_LIST:
            for d in DELTAS_A:
                name = f"{el}_1x1_c{c:.1f}_d{d:.2f}"
                text = pw_input("scf", el, a, c, atoms_1x1(el, c, d), KPTS_1X1)
                write_job(root, name, text)
                jobs.append(dict(name=name, element=el, cell="1x1", c=c,
                                 delta=d, type="scf", kpts=list(KPTS_1X1),
                                 nat=4))
    # Set B: sqrt3 supercell, Na only, x = 1/3
    el = "Na"
    a = A_LAT[el] * np.sqrt(3.0)
    for c in C_LIST:
        for d in DELTAS_B:
            name = f"Na_s3_c{c:.1f}_d{d:.2f}"
            text = pw_input("scf", el, a, c, atoms_s3(el, c, d), KPTS_S3)
            write_job(root, name, text)
            jobs.append(dict(name=name, element=el, cell="s3", c=c,
                             delta=d, type="scf", kpts=list(KPTS_S3), nat=10))
    manifest = dict(pseudos=PSEUDOS, zval=ZVAL, a_lat=A_LAT, z_O=Z_O,
                    jobs=jobs)
    with open(os.path.join(root, "manifest.json"), "w") as f:
        json.dump(manifest, f, indent=2)
    print(f"wrote {len(jobs)} SCF jobs + manifest.json")


# ------------------------------------------------------------ NSCF stage ----
ENERGY_RE = re.compile(r"^!\s+total energy\s+=\s+(-?[\d.]+)\s+Ry", re.M)


def final_energy(pw_out):
    try:
        with open(pw_out) as f:
            txt = f.read()
    except OSError:
        return None
    if "JOB DONE" not in txt:
        return None
    m = ENERGY_RE.findall(txt)
    return float(m[-1]) if m else None


def generate_nscf(root):
    with open(os.path.join(root, "manifest.json")) as f:
        manifest = json.load(f)
    jobs = manifest["jobs"]
    existing = {j["name"] for j in jobs}
    # group SCF jobs by (element, cell, c)
    groups = {}
    for j in jobs:
        if j["type"] != "scf":
            continue
        groups.setdefault((j["element"], j["cell"], j["c"]), []).append(j)
    # set C: Na 1x1, Li 1x1, Na s3 -> 12 NSCF+DOS jobs
    n_new = 0
    for (el, cell, c), grp in sorted(groups.items()):
        energies = []
        for j in grp:
            e = final_energy(os.path.join(root, "jobs", j["name"], "pw.out"))
            if e is not None:
                energies.append((e, j))
        if not energies:
            print(f"skip nscf for {el}/{cell}/c={c}: no finished SCF")
            continue
        e_min, parent = min(energies, key=lambda t: t[0])
        name = f"nscf_{el}_{cell}_c{c:.1f}"
        if name in existing:
            continue
        a = A_LAT[el] * (np.sqrt(3.0) if cell == "s3" else 1.0)
        atoms = (atoms_s3 if cell == "s3" else atoms_1x1)(el, c,
                                                          parent["delta"])
        outdir = f"../{parent['name']}/out"  # reuse parent charge density
        text = pw_input("nscf", el, a, c, atoms, KPTS_NSCF, outdir=outdir)
        write_job(root, name, text, extra={"dos.in": dos_input(outdir)})
        jobs.append(dict(name=name, element=el, cell=cell, c=c,
                         delta=parent["delta"], type="nscf",
                         parent=parent["name"], kpts=list(KPTS_NSCF),
                         nat=len(atoms)))
        n_new += 1
        print(f"nscf job {name}: parent {parent['name']} "
              f"(delta_min = {parent['delta']:.2f} A, E = {e_min:.6f} Ry)")
    with open(os.path.join(root, "manifest.json"), "w") as f:
        json.dump(manifest, f, indent=2)
    print(f"wrote {n_new} NSCF jobs, manifest updated")


# ------------------------------------------------------------ extra stage ---
# Sets E-H, written to jobs_extra/ + manifest_extra.json, run via run_extra.sh.
EXTRA_DIR = "jobs_extra"
DELTAS_F = [0.0, 0.15, 0.30, 0.50]
CHARGES_F = [-0.1, -0.2, -0.3]          # tot_charge: added e- + jellium
DELTAS_G = [0.0, 0.15, 0.30, 0.50, 0.75]
D_OW = 1.2                              # water-O height above/below Na plane, A
R_OH = 0.9572                           # gas-phase water O-H bond length, A
ANG_HOH = 104.52                        # gas-phase water H-O-H angle, deg


def _cart_hex(u, v, a):
    """Fractional (u,v) -> Cartesian (x,y), a1 = (a,0), a2 = (-a/2, a*sqrt3/2)."""
    return u * a - v * a / 2.0, v * a * np.sqrt(3.0) / 2.0


def _frac_hex(x, y, a):
    """Cartesian (x,y) -> fractional (u,v), wrapped into [0,1)."""
    v = 2.0 * y / (a * np.sqrt(3.0))
    u = x / a + v / 2.0
    return u % 1.0, v % 1.0


# Water orientation angles (local frame per O site; see water_at docstring).
# Any two same-layer water columns in this cell are min-image 2.888 A apart,
# so coplanar H-H axes always produce ~1.7-2.0 A H..H image clashes; tilting
# the molecular planes removes them.  These angles were found by maximizing
# the minimum intermolecular contact distance over the whole delta grid; the
# resulting floor, 2.05 A, is the delta=0 Na-Ow distance itself (i.e. the
# H orientation is no longer the limiting contact).
PSI_W = np.radians(50.0)                       # bisector tilt off +-z
PHI_W = {-1: np.radians(80.0), +1: np.radians(60.0)}  # H-H axis azimuth


def water_at(site_uv, a, c, z_ow, sgn):
    """One rigid H2O: O at fractional (site_uv, z_ow/c), sgn = -1/+1 for the
    lower/upper water layer.  Orientation in the local (radial r, tangential
    t, z) frame of the O site relative to the Na column at (0,0):
      bisector  b = sin(psi) r + sgn cos(psi) z
      H-H axis  e = cos(phi) t + sin(phi) (cos(psi) r - sgn sin(psi) z)
    i.e. the H side points away from Na both vertically (toward the nearest
    CoO2 plane) and radially; the layer z-dipoles cancel exactly (2 down,
    2 up)."""
    x0, y0 = _cart_hex(site_uv[0], site_uv[1], a)
    o = np.array([x0, y0, z_ow])
    rho = np.arctan2(y0, x0)
    r = np.array([np.cos(rho), np.sin(rho), 0.0])
    t = np.array([-np.sin(rho), np.cos(rho), 0.0])
    zhat = np.array([0.0, 0.0, 1.0])
    psi, phi = PSI_W, PHI_W[sgn]
    b = np.sin(psi) * r + sgn * np.cos(psi) * zhat
    e = (np.cos(phi) * t
         + np.sin(phi) * (np.cos(psi) * r - sgn * np.sin(psi) * zhat))
    half = np.radians(ANG_HOH / 2.0)
    lb, lt = R_OH * np.cos(half), R_OH * np.sin(half)
    atoms = [("O", *_frac_hex(x0, y0, a), z_ow / c)]
    for s in (1.0, -1.0):
        p = o + lb * b + s * lt * e
        atoms.append(("H", *_frac_hex(p[0], p[1], a), p[2] / c))
    return atoms


# Water in-plane sites (sqrt3 cell fractional) and which layer each occupies.
# Na sits on the (0,0) Co-top site; the 4 waters take the two OTHER Co-top
# images and two O-top images (of the lower framework-O column).  Layer
# assignment is steric: the close O-top sites (in-plane 1.67 A from Na) go in
# the LOWER layer that Na moves away from as delta > 0 (Na-Ow 2.05-2.57 A over
# the delta grid); the farther Co-top images (in-plane 2.89 A) go in the UPPER
# layer Na moves toward (Na-Ow 2.92-3.13 A).  Same-layer O-O = 2.89 A.
WATER_SITES = [
    ((1 / 3, 0.0), -1),     # lower-O-top image, lower water layer (H down)
    ((0.0, 1 / 3), -1),     # lower-O-top image, lower water layer (H down)
    ((1 / 3, 2 / 3), +1),   # Co-top image,      upper water layer (H up)
    ((2 / 3, 1 / 3), +1),   # Co-top image,      upper water layer (H up)
]

G_NOTES = [
    "! Set G water-placement ANSATZ (crystallographically-informed guess",
    "! after Jorgensen et al., PRB 68, 214517 (2003): in Na_x CoO2 . 1.3 H2O",
    "! the Na plane at z = c/2 is sandwiched by two water layers ~1.2 A",
    "! below/above it, with Na octahedrally coordinated by water O):",
    "!   - Na on the (0,0) Co-top site at z = c/2 + delta;",
    "!   - 4 rigid gas-phase H2O (O-H 0.9572 A, H-O-H 104.52 deg), water O at",
    "!     z = c/2 -/+ 1.2 A, on in-plane sites NOT occupied by Na: the two",
    "!     O-top images (1/3,0), (0,1/3) in the LOWER layer (Na moves away",
    "!     from them, Na-Ow 2.05-2.57 A over the delta grid) and the two",
    "!     other Co-top images (1/3,2/3), (2/3,1/3) in the UPPER layer (Na",
    "!     moves toward them, Na-Ow 2.92-3.13 A);",
    "!   - H point away from Na, toward the nearest CoO2 plane and radially",
    "!     outward (bisector tilted 50 deg off +-z); the layer z-dipoles",
    "!     cancel exactly (2 down, 2 up).  Orientation angles maximize the",
    "!     minimum intermolecular contact: ALL water-water, water-framework",
    "!     and water-Na contacts stay >= 2.05 A over the whole delta grid;",
    "!   - waters stay FIXED while Na scans delta, so E(delta) probes the",
    "!     double well inside a static solvation cage. The twin *_vac job",
    "!     differs ONLY by deleting the 12 H2O atoms -> water effect on",
    "!     E(delta) is the clean difference of the two scans.",
]


def atoms_s3_hydrate(c, delta, with_water=True):
    """sqrt3 x sqrt3 Na_1/3 CoO2 cell +- 4 H2O (22 atoms with water, 10 without)."""
    atoms = list(atoms_s3("Na", c, delta))
    if with_water:
        a = A_LAT["Na"] * np.sqrt(3.0)
        for site, sgn in WATER_SITES:
            atoms += water_at(site, a, c, c / 2.0 + sgn * D_OW, sgn)
    return atoms


# --- Set H: Na2CoSe2O structure -- VERIFIED (not a skeleton) -----------------
# Source: Y. et al., J. Am. Chem. Soc. 146, 5908 (2024), DOI 10.1021/jacs.3c11968
# (arXiv:2310.17464); single-crystal XRD refinement in the Supporting
# Information (ja3c11968_si_001.pdf, Tables S1-S2), checked 2026-07-12:
#   space group R-3m (No. 166), hexagonal setting, Z = 3,
#   a = 3.5161(9) A, c = 28.745(11) A;
#   published sites: Se (1/3,2/3,0.37783), Co (2/3,1/3,1/3), O (0,0,1/2),
#   Na (2/3,1/3,0.4584); reduced to the standard origin these are
#   Co 3a (0,0,0), O 3b (0,0,1/2), Se 6c (0,0,0.28884), Na 6c (0,0,0.1251).
# Stacking: edge-sharing CoSe6 [CoSe2] layers alternate with anti-CdCl2-like
# [Na2O] layers along c -> the [Na2O] slab is the "gallery" whose Na-derived
# band we test for weight near E_F.
# All sites sit on (0,0,z_hex), so in the RHOMBOHEDRAL PRIMITIVE cell
# (r1+r2+r3 = c_hex) their crystal coordinates are exactly (z,z,z).
NA2COSE2O = dict(a_hex=3.5161, c_hex=28.745, z_se=0.28884, z_na=0.1251)
KPTS_H_SCF = (8, 8, 8)
KPTS_H_NSCF = (16, 16, 16)   # dense k: eigenvalue dump for the gallery band

H_NOTES = [
    "! Set H: Na2CoSe2O primitive rhombohedral cell (1 f.u., 6 atoms).",
    "! Structure VERIFIED from JACS 146, 5908 (2024) SI Tables S1-S2",
    "! (arXiv:2310.17464): R-3m (166), a_hex = 3.5161 A, c_hex = 28.745 A;",
    "! Co 3a (0,0,0), O 3b (0,0,1/2), Se 6c z = 0.28884, Na 6c z = 0.1251",
    "! (standard-origin hexagonal setting; on the 3-fold axis the",
    "! rhombohedral crystal coordinates are exactly (z,z,z)).",
    "! Purpose: does a Na-derived [Na2O]-gallery band sit near E_F?",
    "! Inspect the printed eigenvalues/occupations (verbosity high) and the",
    "! dos.x output of the dense-k NSCF twin.",
]


def pw_input_h(calc, kpts, outdir="./out"):
    """pw.x input for the Na2CoSe2O rhombohedral primitive cell (ibrav=5)."""
    s = NA2COSE2O
    a2, c2 = s["a_hex"] ** 2, s["c_hex"] ** 2
    a_rh = np.sqrt(a2 / 3.0 + c2 / 9.0)              # 9.7944 A
    cos_alpha = (2.0 * c2 - 3.0 * a2) / (2.0 * (c2 + 3.0 * a2))
    t, u = s["z_se"], s["z_na"]
    atoms = [("Co", 0.0, 0.0, 0.0),
             ("O", 0.5, 0.5, 0.5),
             ("Se", t, t, t),
             ("Se", 1.0 - t, 1.0 - t, 1.0 - t),
             ("Na", u, u, u),
             ("Na", 1.0 - u, 1.0 - u, 1.0 - u)]
    species = ["Na", "Co", "Se", "O"]
    nelec = sum(ZVAL[sp] for sp, *_ in atoms)
    nbnd = nbnd_for(nelec)
    lines = [
        "&CONTROL",
        f"  calculation = '{calc}'",
        "  prefix = 'pw'",
        f"  outdir = '{outdir}'",
        "  pseudo_dir = '../../pseudo'",
        "  verbosity = 'high'",
        "  tprnfor = .true.",
        "/",
        "&SYSTEM",
        "  ibrav = 5",
        f"  celldm(1) = {a_rh / BOHR:.8f}",
        f"  celldm(4) = {cos_alpha:.8f}",
        f"  nat = {len(atoms)}",
        f"  ntyp = {len(species)}",
        f"  nbnd = {nbnd}",
        "  ecutwfc = 60.0",
        "  ecutrho = 480.0",
        "  occupations = 'smearing'",
        "  smearing = 'mv'",
        "  degauss = 0.02",
        "  nspin = 2",
        f"  starting_magnetization({species.index('Co') + 1}) = 0.3",
        "/",
        "&ELECTRONS",
        "  conv_thr = 1.0d-7",
        "  mixing_beta = 0.3",
        "  mixing_mode = 'local-TF'",
        "  electron_maxstep = 300",
        "/",
        "ATOMIC_SPECIES",
    ]
    for sp in species:
        lines.append(f"  {sp}  {MASS[sp]:.5f}  {PSEUDOS[sp]}")
    lines.extend(H_NOTES)
    lines.append("ATOMIC_POSITIONS crystal")
    for sp, x, y, z in atoms:
        lines.append(f"  {sp}  {x:.10f}  {y:.10f}  {z:.10f}")
    lines.append("K_POINTS automatic")
    lines.append(f"  {kpts[0]} {kpts[1]} {kpts[2]} 0 0 0")
    lines.append("")
    return "\n".join(lines)


def generate_set_h(root, add):
    name_scf = "Na2CoSe2O_scf"
    write_job(root, name_scf, pw_input_h("scf", KPTS_H_SCF),
              jobs_dir=EXTRA_DIR)
    add(dict(name=name_scf, set="H", element="Na2CoSe2O", cell="r3m-prim",
             type="scf", kpts=list(KPTS_H_SCF), nat=6))
    name_nscf = "Na2CoSe2O_nscf"
    outdir = f"../{name_scf}/out"
    write_job(root, name_nscf, pw_input_h("nscf", KPTS_H_NSCF, outdir=outdir),
              extra={"dos.in": dos_input(outdir)}, jobs_dir=EXTRA_DIR)
    add(dict(name=name_nscf, set="H", element="Na2CoSe2O", cell="r3m-prim",
             type="nscf", parent=name_scf, kpts=list(KPTS_H_NSCF), nat=6))


def generate_extra(root):
    jobs = []
    counts = {}

    def add(job):
        jobs.append(job)
        counts[job["set"]] = counts.get(job["set"], 0) + 1

    # --- Set E: K_x CoO2 1x1 E(delta) scans (exact mirror of set A, K) ------
    el, a = "K", A_LAT["K"]
    e_note = ["! Set E: in-plane a = 2.90 A is an ANSATZ for K_xCoO2 (see",
              "! A_LAT note in generate_inputs.py); mirrors the Na 1x1 scans."]
    for c in C_LIST:
        for d in DELTAS_A:
            name = f"K_1x1_c{c:.1f}_d{d:.2f}"
            text = pw_input("scf", el, a, c, atoms_1x1(el, c, d), KPTS_1X1,
                            notes=e_note)
            write_job(root, name, text, jobs_dir=EXTRA_DIR)
            add(dict(name=name, set="E", element=el, cell="1x1", c=c, delta=d,
                     type="scf", kpts=list(KPTS_1X1), nat=4))

    # --- Set F: gated LiCoO2 (ionic-liquid-gating proxy) --------------------
    el, a, c = "Li", A_LAT["Li"], 5.5
    for d in DELTAS_F:
        for q in CHARGES_F:
            extra = [
                "  ! Set F: tot_charge < 0 ADDS |tot_charge| electrons with a",
                "  ! compensating uniform jellium background -- an electrostatic",
                "  ! proxy for ionic-liquid gating of LiCoO2. Check via the",
                "  ! printed eigenvalues/occupations (verbosity high) whether",
                "  ! the added charge occupies the interlayer gallery band and",
                "  ! whether the delta well stays quantum-paraelectric.",
                f"  tot_charge = {q:.1f}",
            ]
            name = f"Li_gate_c{c:.1f}_d{d:.2f}_q{abs(q):.1f}"
            text = pw_input("scf", el, a, c, atoms_1x1(el, c, d), KPTS_1X1,
                            extra_system=extra)
            write_job(root, name, text, jobs_dir=EXTRA_DIR)
            add(dict(name=name, set="F", element=el, cell="1x1", c=c, delta=d,
                     tot_charge=q, type="scf", kpts=list(KPTS_1X1), nat=4))

    # --- Set G: explicit-water bilayer hydrate + vacuum reference -----------
    el, c = "Na", 9.9
    a = A_LAT[el] * np.sqrt(3.0)
    for with_water, tag in ((True, "hyd"), (False, "vac")):
        for d in DELTAS_G:
            atoms = atoms_s3_hydrate(c, d, with_water=with_water)
            species = ["Co", "O", "Na", "H"] if with_water else None
            notes = list(G_NOTES) if with_water else [
                "! Set G vacuum reference: identical to the *_hyd twin with",
                "! the 12 H2O atoms deleted (see G_NOTES ansatz there)."]
            name = f"Na_s3{tag}_c{c:.1f}_d{d:.2f}"
            text = pw_input("scf", el, a, c, atoms, KPTS_G,
                            species=species, notes=notes)
            write_job(root, name, text, jobs_dir=EXTRA_DIR)
            add(dict(name=name, set="G", element=el, cell="s3", c=c, delta=d,
                     water=with_water, type="scf", kpts=list(KPTS_G),
                     nat=len(atoms)))

    # --- Set H: Na2CoSe2O gallery-band check ---------------------------------
    generate_set_h(root, add)

    manifest = dict(pseudos=PSEUDOS, zval=ZVAL, a_lat=A_LAT, z_O=Z_O,
                    jobs=jobs)
    with open(os.path.join(root, "manifest_extra.json"), "w") as f:
        json.dump(manifest, f, indent=2)
    per_set = ", ".join(f"{k}: {v}" for k, v in sorted(counts.items()))
    print(f"wrote {len(jobs)} extra jobs ({per_set}) + manifest_extra.json")


# ------------------------------------------------------------ bands stage ---
# Set I (referee response), written to jobs_bands/ + manifest_bands.json,
# run via run_bands.sh.
BANDS_DIR = "jobs_bands"

# Hexagonal-BZ band path Gamma-M-K-Gamma in crystal (reciprocal-lattice
# fractional) coordinates: Gamma (0,0,0), M (1/2,0,0), K (1/3,1/3,0).
# crystal_b semantics: the integer after each point is the number of k-points
# from that point up to (excluding) the next one; pw.x appends the final
# point, so totals are 22 + 12 + 25 + 1 = 60 k-points.  Segment counts are
# proportional to segment lengths |GM| : |MK| : |KG| = 1 : 1/sqrt3 : 2/sqrt3.
KPATH_HEX = [
    (0.0, 0.0, 0.0, 22),        # Gamma
    (0.5, 0.0, 0.0, 12),        # M
    (1 / 3, 1 / 3, 0.0, 25),    # K
    (0.0, 0.0, 0.0, 1),         # Gamma (final point, count ignored)
]

DELTAS_I2 = [0.0, 0.30]
DELTAS_I3 = [0.0, 0.30, 0.50]

I1_NOTES = [
    "! Set I1 (referee response): band structure with fatband projections.",
    "! Chain per geometry: pw scf -> pw 'bands' (Gamma-M-K-Gamma, 60 pts,",
    "! K_POINTS crystal_b) -> bands.x (spin up/dw) -> projwfc.x.  projwfc",
    "! projects onto ALL atoms, so one run yields both the Na s/p fatbands",
    "! (interlayer-gallery weight, the reviewer figure) and the Co d",
    "! fatbands (e_g' pocket check).",
]
I1_NOTES_S3 = I1_NOTES + [
    "! sqrt3 x sqrt3 R30 cell: the path is Gamma-M-K-Gamma of the SUPERCELL",
    "! BZ (rotated 30 deg, 1/3 the area of the 1x1 BZ).  Zone folding: the",
    "! 1x1 K points fold onto the supercell Gamma and the 1x1 M points onto",
    "! the supercell M, so the 1x1 K-point physics appears at Gamma here.",
]
I2_NOTES = [
    "! Set I2 (referee response): relaxed-water check.  Same hydrate ansatz",
    "! as set G (see G_NOTES in generate_inputs.py) but calculation='relax'",
    "! (BFGS): Na z and ALL 12 H2O atoms are free; the CoO2 framework (3 Co",
    "! + 6 lattice O) is frozen and Na is pinned to its Co-top column",
    "! (x,y fixed) via the if_pos 0/1 triplets on each position line.",
    "! Purpose: does the bounded off-center Na minimum survive when the",
    "! water cage is allowed to respond?  Compare the final Na z for the",
    "! delta = 0.00 vs 0.30 A starting points.",
]
I3_NOTES = [
    "! Set I3 (referee response): vdW spot-check.  Geometry identical to the",
    "! set G hydrate SCF at the same delta, plus vdw_corr = 'grimme-d3'",
    "! (DFT-D3, supported by QE 7.3.1).  Compare E(delta) at 0.00/0.30/0.50",
    "! against the non-vdW set G energies: does D3 change the well shape?",
]


def bands_x_input(spin, outdir="./out"):
    """bands.x input; nspin=2 runs need one file per spin_component."""
    tag = {1: "up", 2: "dw"}[spin]
    return "\n".join([
        "&BANDS",
        "  prefix = 'pw'",
        f"  outdir = '{outdir}'",
        f"  spin_component = {spin}",
        f"  filband = 'pw.bands.{tag}.dat'",
        "  lsym = .true.",
        "/",
        "",
    ])


def projwfc_input(outdir="./out"):
    return "\n".join([
        "&PROJWFC",
        "  prefix = 'pw'",
        f"  outdir = '{outdir}'",
        "  ngauss = 0",
        "  degauss = 0.02",
        "  DeltaE = 0.01",
        "  filpdos = 'pw.pdos'",
        "  filproj = 'pw.proj'",
        "/",
        "! projwfc.x projects onto ALL atomic wavefunctions by default (no",
        "! atom selection): Na s/p for the gallery fatbands AND Co d for the",
        "! e_g' pocket check come from the same filproj output.",
        "",
    ])


def generate_bands(root):
    jobs = []
    counts = {}

    def add(job):
        jobs.append(job)
        counts[job["set"]] = counts.get(job["set"], 0) + 1

    el = "Na"

    # --- I1: fatband chains (scf -> bands -> bands.x up/dw -> projwfc.x) ----
    for cell, c in (("1x1", 5.5), ("1x1", 6.9), ("1x1", 9.9),
                    ("s3", 5.5), ("s3", 9.9)):
        s3 = cell == "s3"
        a = A_LAT[el] * (np.sqrt(3.0) if s3 else 1.0)
        atoms = (atoms_s3 if s3 else atoms_1x1)(el, c, 0.0)
        kpts = KPTS_S3 if s3 else KPTS_1X1
        notes = I1_NOTES_S3 if s3 else I1_NOTES
        name = f"bands_{el}_{cell}_c{c:.1f}"
        scf = pw_input("scf", el, a, c, atoms, kpts, notes=notes)
        nscf = pw_input("bands", el, a, c, atoms, None, notes=notes,
                        kpath=KPATH_HEX)
        write_job(root, name, scf, jobs_dir=BANDS_DIR, extra={
            "pw_bands.in": nscf,
            "bands_up.in": bands_x_input(1),
            "bands_dw.in": bands_x_input(2),
            "projwfc.in": projwfc_input(),
        })
        add(dict(name=name, set="I1", element=el, cell=cell, c=c, delta=0.0,
                 type="bands_chain", kpts=list(kpts),
                 kpath="G-M-K-G_60pts", nat=len(atoms)))

    # --- I2: hydrate relaxed-water check (BFGS, CoO2 frozen) -----------------
    c = 9.9
    a = A_LAT[el] * np.sqrt(3.0)
    for d in DELTAS_I2:
        atoms = atoms_s3_hydrate(c, d, with_water=True)   # 22 atoms
        # atoms_s3 order: 3 Co (0-2), 6 lattice O (3-8), Na (9), then 4x
        # (O,H,H) water (10-21).  Freeze CoO2, free Na z only, free water.
        if_pos = ([(0, 0, 0)] * 9) + [(0, 0, 1)] + ([(1, 1, 1)] * 12)
        name = f"relax_{el}_s3hyd_c{c:.1f}_d{d:.2f}"
        text = pw_input("relax", el, a, c, atoms, KPTS_G,
                        species=["Co", "O", "Na", "H"], ions=True,
                        if_pos=if_pos, notes=I2_NOTES,
                        extra_control=["  forc_conv_thr = 1.0d-3",
                                       "  nstep = 60"])
        write_job(root, name, text, jobs_dir=BANDS_DIR)
        add(dict(name=name, set="I2", element=el, cell="s3", c=c, delta=d,
                 water=True, type="relax", kpts=list(KPTS_G),
                 nat=len(atoms)))

    # --- I3: hydrate vdW (grimme-d3) spot-check ------------------------------
    for d in DELTAS_I3:
        atoms = atoms_s3_hydrate(c, d, with_water=True)
        name = f"vdw_{el}_s3hyd_c{c:.1f}_d{d:.2f}"
        text = pw_input("scf", el, a, c, atoms, KPTS_G,
                        species=["Co", "O", "Na", "H"], notes=I3_NOTES,
                        extra_system=["  vdw_corr = 'grimme-d3'"])
        write_job(root, name, text, jobs_dir=BANDS_DIR)
        add(dict(name=name, set="I3", element=el, cell="s3", c=c, delta=d,
                 water=True, vdw="grimme-d3", type="scf", kpts=list(KPTS_G),
                 nat=len(atoms)))

    manifest = dict(pseudos=PSEUDOS, zval=ZVAL, a_lat=A_LAT, z_O=Z_O,
                    kpath=[list(p) for p in KPATH_HEX], jobs=jobs)
    with open(os.path.join(root, "manifest_bands.json"), "w") as f:
        json.dump(manifest, f, indent=2)
    per_set = ", ".join(f"{k}: {v}" for k, v in sorted(counts.items()))
    print(f"wrote {len(jobs)} bands-set jobs ({per_set}) "
          f"+ manifest_bands.json")


# ----------------------------------------------------------- mobile stage ---
# Set J (the manuscript's remaining TODO), written to jobs_mobile/ +
# manifest_mobile.json, run via run_mobile.sh.
#
# Physics: the RIGID water cage (set G) bounds the Na runaway seen in vacuum
# but still leaves a ~198 meV off-center well.  The manuscript claim under
# test is that letting the water respond ADIABATICALLY -- i.e. fully relaxing
# all 12 H2O atoms at each PINNED Na displacement delta -- flattens that
# well.  E_mobile(delta) is the relaxed total energy at each pinned delta;
# E_mobile(delta) - E_mobile(0) is the mobile-water potential, to be compared
# against the rigid-cage curve from the set G hydrate SCF scan
# (results_extra, Na_s3hyd_c9.9_d*).
MOBILE_DIR = "jobs_mobile"
DELTAS_J = [0.0, 0.15, 0.30, 0.50, 0.75, 1.00]
DELTA_J3 = 0.30

J1_NOTES = [
    "! Set J (mobile-water adiabatic E(delta), manuscript TODO): same",
    "! hydrate ansatz as set G (see G_NOTES in generate_inputs.py), but",
    "! calculation='relax' (BFGS) with Na FULLY FROZEN (if_pos 0 0 0) at",
    "! z = c/2 + delta while ALL 12 H2O atoms are free (if_pos 1 1 1); the",
    "! CoO2 framework (3 Co + 6 lattice O) is frozen as in set I2.",
    "! E_mobile(delta) = relaxed total energy at this pinned delta.  The",
    "! rigid cage (set G) bounds the Na runaway but leaves a ~198 meV",
    "! off-center well; the claim under test is that adiabatically",
    "! relaxing the water at each Na position flattens that well.",
]
J2_NOTES = J1_NOTES + [
    "! Set J2 twin: identical geometry and constraints, plus",
    "! vdw_corr = 'grimme-d3' (DFT-D3) -> dispersion sensitivity of the",
    "! adiabatic E_mobile(delta) curve.",
]
J3_NOTES = [
    "! Set J3 (joint-minimum search): same hydrate relax as the set J1",
    "! delta = 0.30 job (waters free, CoO2 frozen), but Na z is ALSO free",
    "! (if_pos 0 0 1; x,y pinned to the Co-top column), starting",
    "! off-center at delta = 0.30 A.  Tests whether the JOINT (Na + water)",
    "! minimum is off-center at all: if Na slides back to z = c/2 the",
    "! adiabatic well is not merely flattened but centered.",
]

# atoms_s3_hydrate order: 3 Co (0-2), 6 lattice O (3-8), Na (9), then 4x
# (O,H,H) water (10-21).  Freeze CoO2 AND Na, free all 12 water atoms.
IF_POS_J_PINNED = ([(0, 0, 0)] * 9) + [(0, 0, 0)] + ([(1, 1, 1)] * 12)
# J3: as above but Na z free (x,y still pinned to the Co-top column).
IF_POS_J_NAFREE = ([(0, 0, 0)] * 9) + [(0, 0, 1)] + ([(1, 1, 1)] * 12)

J_RELAX_CONTROL = ["  forc_conv_thr = 1.0d-3", "  nstep = 100"]


def generate_mobile(root):
    jobs = []
    counts = {}

    def add(job):
        jobs.append(job)
        counts[job["set"]] = counts.get(job["set"], 0) + 1

    el, c = "Na", 9.9
    a = A_LAT[el] * np.sqrt(3.0)
    species = ["Co", "O", "Na", "H"]

    # --- J1/J2: adiabatic scan, Na pinned, waters relax (+- grimme-d3) ------
    for vdw, tag, jset, notes in ((False, "", "J1", J1_NOTES),
                                  (True, "vdw_", "J2", J2_NOTES)):
        for d in DELTAS_J:
            atoms = atoms_s3_hydrate(c, d, with_water=True)   # 22 atoms
            name = f"mobile_{tag}{el}_s3hyd_c{c:.1f}_d{d:.2f}"
            extra_system = ["  vdw_corr = 'grimme-d3'"] if vdw else None
            text = pw_input("relax", el, a, c, atoms, KPTS_G,
                            species=species, ions=True,
                            if_pos=IF_POS_J_PINNED, notes=notes,
                            extra_system=extra_system,
                            extra_control=J_RELAX_CONTROL)
            write_job(root, name, text, jobs_dir=MOBILE_DIR)
            job = dict(name=name, set=jset, element=el, cell="s3", c=c,
                       delta=d, water=True, na_constraint="frozen",
                       type="relax", kpts=list(KPTS_G), nat=len(atoms))
            if vdw:
                job["vdw"] = "grimme-d3"
            add(job)

    # --- J3: joint minimum search, Na z free, off-center start --------------
    d = DELTA_J3
    atoms = atoms_s3_hydrate(c, d, with_water=True)
    name = f"mobile_freeNa_{el}_s3hyd_c{c:.1f}_d{d:.2f}"
    text = pw_input("relax", el, a, c, atoms, KPTS_G,
                    species=species, ions=True,
                    if_pos=IF_POS_J_NAFREE, notes=J3_NOTES,
                    extra_control=J_RELAX_CONTROL)
    write_job(root, name, text, jobs_dir=MOBILE_DIR)
    add(dict(name=name, set="J3", element=el, cell="s3", c=c, delta=d,
             water=True, na_constraint="z-free", type="relax",
             kpts=list(KPTS_G), nat=len(atoms)))

    manifest = dict(
        pseudos=PSEUDOS, zval=ZVAL, a_lat=A_LAT, z_O=Z_O,
        analyze=("E_mobile(delta) = final '!  total energy' of each set-J1 "
                 "relax; E_mobile(delta) - E_mobile(0) is the mobile-water "
                 "potential.  Compare against the rigid-cage curve "
                 "E_hyd(delta) - E_hyd(0) from the set G hydrate SCF scan in "
                 "results_extra (Na_s3hyd_c9.9_d*).  J2 (grimme-d3 twin) "
                 "gives the dispersion sensitivity of the adiabatic curve; "
                 "J3 (Na z free, off-center start) tests whether the joint "
                 "Na+water minimum is off-center at all."),
        jobs=jobs)
    with open(os.path.join(root, "manifest_mobile.json"), "w") as f:
        json.dump(manifest, f, indent=2)
    per_set = ", ".join(f"{k}: {v}" for k, v in sorted(counts.items()))
    print(f"wrote {len(jobs)} mobile-set jobs ({per_set}) "
          f"+ manifest_mobile.json")


# ------------------------------------------------------------ zscan stage ---
# Set K, written to jobs_zscan/ + manifest_zscan.json, run via run_zscan.sh.
# At fixed c, vary the CoO2 layer thickness via z_O and fit the Na E(delta)
# well shape.  The z_O = 0.96 A baseline is the existing set-G vacuum scan.
ZSCAN_DIR = "jobs_zscan"
ZOS_K = [0.90, 1.02]
DELTAS_K = [0.0, 0.15, 0.30, 0.50, 0.75]

K_NOTES = [
    "! Set K: vacuum sqrt3 x sqrt3 Na_1/3CoO2 z_O sensitivity scan.",
    "! At fixed c = 9.9 A, vary the CoO2 layer thickness (O height z_O)",
    "! and scan Na displacement delta to extract the well-shape parameter",
    "! alpha(z_O).  The z_O = 0.96 A baseline is the existing set-G vacuum",
    "! scan in results_extra (Na_s3vac_c9.9_d*); it is not regenerated.",
]


def generate_zscan(root):
    jobs = []
    el, c = "Na", 9.9
    a = A_LAT[el] * np.sqrt(3.0)

    for zo in ZOS_K:
        for d in DELTAS_K:
            atoms = atoms_s3(el, c, d, zo=zo)
            name = f"zscan_{el}_s3_c{c:.1f}_zO{zo:.2f}_d{d:.2f}"
            text = pw_input("scf", el, a, c, atoms, KPTS_G, notes=K_NOTES)
            write_job(root, name, text, jobs_dir=ZSCAN_DIR)
            jobs.append(dict(name=name, set="K", element=el, cell="s3", c=c,
                             z_O=zo, delta=d, type="scf",
                             kpts=list(KPTS_G), nat=len(atoms)))

    manifest = dict(
        pseudos=PSEUDOS, zval=ZVAL, a_lat=A_LAT, z_O=Z_O,
        analyze=("Fit alpha(z_O), the well-shape/curvature parameter, at fixed "
                 "c = 9.9 A across the delta scan for each z_O, including the "
                 "existing z_O = 0.96 A baseline curve from set G in "
                 "results_extra (Na_s3vac_c9.9_d*), which is not regenerated. "
                 "Then compare the resulting slope d(alpha)/d(z_O) against "
                 "the existing set D single-point sensitivity result of "
                 "roughly +/-100 meV per -/+0.06 A in z_O."),
        jobs=jobs)
    with open(os.path.join(root, "manifest_zscan.json"), "w") as f:
        json.dump(manifest, f, indent=2)
    print(f"wrote {len(jobs)} zscan-set jobs (K: {len(jobs)}) "
          f"+ manifest_zscan.json")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--stage", choices=["scf", "nscf", "extra", "bands",
                                       "mobile", "zscan"],
                   default="scf")
    p.add_argument("--root", default=os.path.dirname(os.path.abspath(__file__)))
    args = p.parse_args()
    if args.stage == "scf":
        generate_scf(args.root)
    elif args.stage == "extra":
        generate_extra(args.root)
    elif args.stage == "bands":
        generate_bands(args.root)
    elif args.stage == "mobile":
        generate_mobile(args.root)
    elif args.stage == "zscan":
        generate_zscan(args.root)
    else:
        generate_nscf(args.root)


if __name__ == "__main__":
    sys.exit(main())

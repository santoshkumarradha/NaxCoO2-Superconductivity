#!/usr/bin/env python3
"""Generate Quantum ESPRESSO pw.x inputs for the Na_xCoO2 interlayer-2DEG /
anharmonic-rattler study.

Default stage writes all SCF jobs (sets A and B) into jobs/<name>/pw.in and a
manifest.json.  After the SCF phase has run, `--stage nscf` parses the SCF
energies, finds the E(delta) minimum per (element, cell, c) and writes the
dense-k NSCF + dos.x jobs (set C), appending them to the manifest.

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
A_LAT = {"Na": 2.888, "Li": 2.82}     # in-plane lattice constant, Angstrom
C_LIST = [5.5, 6.9, 8.4, 9.9]         # CoO2 plane-to-plane spacing, Angstrom
DELTAS_A = [0.0, 0.15, 0.30, 0.50, 0.75, 1.00]  # off-center displacement, Ang
DELTAS_B = [0.0, 0.50]
Z_O = 0.96                             # O height above Co plane, Angstrom

PSEUDOS = {
    # NOTE: Co.pbe-spn-kjpaw_psl.1.0.0.UPF is NOT hosted on
    # pseudopotentials.quantum-espresso.org; the closest available pslibrary
    # PAW variant (same 3s3p semicore spn flavour) is used instead.
    "Co": "Co.pbe-spn-kjpaw_psl.0.3.1.UPF",
    "O": "O.pbe-n-kjpaw_psl.1.0.0.UPF",
    "Na": "Na.pbe-spn-kjpaw_psl.1.0.0.UPF",
    "Li": "Li.pbe-s-kjpaw_psl.1.0.0.UPF",
}
MASS = {"Co": 58.933, "O": 15.999, "Na": 22.98977, "Li": 6.941}
ZVAL = {"Co": 17, "O": 6, "Na": 9, "Li": 3}  # PAW valence electrons

KPTS_1X1 = (8, 8, 3)
KPTS_S3 = (6, 6, 3)
KPTS_NSCF = (24, 24, 8)


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
            (element, 1 / 3, 2 / 3, za)]


def atoms_s3(element, c, delta):
    """sqrt3 x sqrt3 R30 supercell, 1 alkali per 3 Co (x = 1/3), 10 atoms.

    Primitive fractional (u,v) maps to supercell fractional
    ((u+v)/3, (2v-u)/3); Co sublattice images are (0,0), (1/3,2/3), (2/3,1/3).
    """
    zo = Z_O / c
    za = (c / 2.0 + delta) / c
    offs = [(0.0, 0.0), (1 / 3, 2 / 3), (2 / 3, 1 / 3)]
    atoms = []
    for ox, oy in offs:
        atoms.append(("Co", ox % 1.0, oy % 1.0, 0.0))
    for base, z in (((1 / 3, 1 / 3), zo), ((1 / 3, 0.0), 1.0 - zo)):
        for ox, oy in offs:
            atoms.append(("O", (base[0] + ox) % 1.0, (base[1] + oy) % 1.0, z))
    atoms.append((element, 1 / 3, 1 / 3, za))  # image of primitive (1/3,2/3)
    return atoms


# ------------------------------------------------------------ input files ---
def pw_input(calc, element, a, c, atoms, kpts, prefix="pw", outdir="./out"):
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
        "/",
        "&ELECTRONS",
        "  conv_thr = 1.0d-8",
        "  mixing_beta = 0.3",
        "  electron_maxstep = 300",
        "/",
        "ATOMIC_SPECIES",
    ]
    for s in species:
        lines.append(f"  {s}  {MASS[s]:.5f}  {PSEUDOS[s]}")
    lines.append("ATOMIC_POSITIONS crystal")
    for s, x, y, z in atoms:
        lines.append(f"  {s}  {x:.10f}  {y:.10f}  {z:.10f}")
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


def write_job(root, name, text, extra=None):
    d = os.path.join(root, "jobs", name)
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


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--stage", choices=["scf", "nscf"], default="scf")
    p.add_argument("--root", default=os.path.dirname(os.path.abspath(__file__)))
    args = p.parse_args()
    if args.stage == "scf":
        generate_scf(args.root)
    else:
        generate_nscf(args.root)


if __name__ == "__main__":
    sys.exit(main())

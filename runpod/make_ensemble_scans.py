#!/usr/bin/env python3
"""Pod-side stage M generator for the "cage-ensemble" set L/M.

Run from runpod/ AFTER jobs_ensemble/md_cage/pw.out has finished (or made
enough progress -- run_ensemble.sh gates this on >= 3000 MD steps):

    python3 make_ensemble_scans.py

Parses the set-L Born-Oppenheimer MD trajectory (jobs_ensemble/md_cage/pw.out,
water-only dynamics at Na delta = 0, see generate_inputs.py generate_ensemble
docstring), discards the first 1000 steps as equilibration, takes 10 snapshots
evenly spaced over the remainder, and for each snapshot writes 7 SCF inputs
with the water atoms FROZEN at the snapshot coordinates and Na displaced to
delta in {-0.45, -0.30, -0.15, 0.00, +0.15, +0.30, +0.45} A along c from the
midplane (the CoO2 framework is taken from the original, static, input
geometry). These 70 SCF jobs use the EXACT set-G electronic settings (K_POINTS
6 6 2, conv_thr 1e-7, same 60/480 Ry cutoffs, mv 0.02 smearing, nspin=2 with
the set-G starting_magnetization) for direct comparability with the rigid-cage
(set G) and adiabatic (set J) results. Job dirs:
jobs_ensemble/snap<NN>_d<+/-0.NN>/pw.in.

Pure stdlib -- no numpy, no third-party deps (deliberately, even though numpy
is available in the NGC image: the geometry math needed here is trivial
fractional-coordinate arithmetic, so avoiding the import keeps this script
robust to any pod/image quirk).
"""
import json
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
MD_DIR = os.path.join(HERE, "jobs_ensemble", "md_cage")
MD_IN = os.path.join(MD_DIR, "pw.in")
MD_OUT = os.path.join(MD_DIR, "pw.out")
ENSEMBLE_DIR = os.path.join(HERE, "jobs_ensemble")
MANIFEST = os.path.join(HERE, "manifest_ensemble.json")

BOHR = 0.529177210903  # Angstrom

EQUIL_STEPS = 1000       # discard as equilibration
N_SNAPSHOTS = 10          # evenly spaced over the remaining trajectory
DELTAS_M = [-0.45, -0.30, -0.15, 0.00, 0.15, 0.30, 0.45]  # Angstrom
MIN_STEPS_REQUIRED = 3000  # must match run_ensemble.sh's phase-1 gate

# Same electronic-structure constants as generate_inputs.py (set G), kept
# self-contained here so this script has zero third-party dependencies.
PSEUDOS = {
    "Co": "Co.pbe-spn-kjpaw_psl.0.3.1.UPF",
    "O": "O.pbe-n-kjpaw_psl.1.0.0.UPF",
    "Na": "Na.pbe-spn-kjpaw_psl.1.0.0.UPF",
    "H": "H.pbe-kjpaw_psl.1.0.0.UPF",
}
MASS = {"Co": 58.933, "O": 15.999, "Na": 22.98977, "H": 1.008}
ZVAL = {"Co": 17, "O": 6, "Na": 9, "H": 1}
KPTS_M = (6, 6, 2)  # exact set-G K-mesh


# ------------------------------------------------------------- pw.in parse --
ATOM_LINE_RE = re.compile(
    r"^\s*([A-Za-z][A-Za-z0-9]*)\s+(-?[0-9.eEdD+-]+)\s+(-?[0-9.eEdD+-]+)"
    r"\s+(-?[0-9.eEdD+-]+)(?:\s+([01])\s+([01])\s+([01]))?\s*$"
)


def parse_md_input(path):
    """Return (a_ang, c_ang, species_order, atoms, if_pos) from the set-L
    pw.in, where atoms = [(species, x, y, z), ...] in file order and if_pos
    is a parallel list of (ix,iy,iz) ints (None entries default to (1,1,1),
    but the set-L input always writes explicit flags)."""
    with open(path) as f:
        txt = f.read()
    c1 = float(re.search(r"celldm\(1\)\s*=\s*([\d.eEdD+-]+)", txt).group(1))
    c3 = float(re.search(r"celldm\(3\)\s*=\s*([\d.eEdD+-]+)", txt).group(1))
    a = c1 * BOHR
    c = c3 * a

    m = re.search(r"ATOMIC_SPECIES\s*\n(.*?)(?=\n\s*(?:ATOMIC_POSITIONS|!))",
                  txt, re.S)
    species_order = []
    for line in m.group(1).splitlines():
        line = line.strip()
        if not line:
            continue
        species_order.append(line.split()[0])

    m = re.search(r"ATOMIC_POSITIONS\s+crystal\s*\n(.*?)(?=\nK_POINTS)",
                  txt, re.S)
    atoms, if_pos = [], []
    for line in m.group(1).splitlines():
        line = line.strip()
        if not line or line.startswith(("!", "#")):
            continue
        mm = ATOM_LINE_RE.match(line)
        if not mm:
            continue
        sp, x, y, z = mm.group(1), float(mm.group(2)), float(mm.group(3)), \
            float(mm.group(4))
        atoms.append((sp, x, y, z))
        if mm.group(5) is not None:
            if_pos.append((int(mm.group(5)), int(mm.group(6)),
                           int(mm.group(7))))
        else:
            if_pos.append((1, 1, 1))
    return a, c, species_order, atoms, if_pos


# ------------------------------------------------------------ pw.out parse --
def parse_md_trajectory(path, nat):
    """Return a list of snapshots, each a list of nat (species,x,y,z) tuples
    in crystal (fractional) coordinates, in trajectory order (one entry per
    ATOMIC_POSITIONS block found, including the initial geometry)."""
    with open(path) as f:
        txt = f.read()
    snapshots = []
    for m in re.finditer(r"ATOMIC_POSITIONS[^\n]*\n", txt):
        start = m.end()
        lines = txt[start:start + 4000].splitlines()  # generous slice
        block = []
        for line in lines:
            mm = ATOM_LINE_RE.match(line)
            if not mm:
                if block:
                    break
                continue
            block.append((mm.group(1), float(mm.group(2)),
                         float(mm.group(3)), float(mm.group(4))))
            if len(block) == nat:
                break
        if len(block) == nat:
            snapshots.append(block)
    return snapshots


def count_md_steps(path):
    """Cheap step count: number of ATOMIC_POSITIONS markers (~= nstep + 1
    for the initial geometry); matches run_ensemble.sh's phase-1 gate logic
    so this script and the bash gate agree."""
    if not os.path.exists(path):
        return 0
    with open(path) as f:
        txt = f.read()
    return len(re.findall(r"ATOMIC_POSITIONS", txt))


# -------------------------------------------------------------- pw.x input --
def nbnd_for(nelec):
    return nelec // 2 + max(8, nelec // 6)


def pw_input_scf(a, c, species_order, atoms, delta, notes):
    """SCF pw.in text, EXACT set-G electronic settings. `atoms` already has
    Na at z = c/2 + delta and water at the frozen snapshot coordinates."""
    nelec = sum(ZVAL[sp] for sp, *_ in atoms)
    nbnd = nbnd_for(nelec)
    lines = [
        "&CONTROL",
        "  calculation = 'scf'",
        "  prefix = 'pw'",
        "  outdir = './out'",
        "  pseudo_dir = '../../pseudo'",
        "  verbosity = 'high'",
        "  tprnfor = .true.",
        "/",
        "&SYSTEM",
        "  ibrav = 4",
        f"  celldm(1) = {a / BOHR:.8f}",
        f"  celldm(3) = {c / a:.8f}",
        f"  nat = {len(atoms)}",
        f"  ntyp = {len(species_order)}",
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
        "  conv_thr = 1.0d-7",
        "  mixing_beta = 0.3",
        "  mixing_mode = 'local-TF'",
        "  electron_maxstep = 300",
        "/",
        "ATOMIC_SPECIES",
    ]
    for sp in species_order:
        lines.append(f"  {sp}  {MASS[sp]:.5f}  {PSEUDOS[sp]}")
    lines.extend(notes)
    lines.append("ATOMIC_POSITIONS crystal")
    for sp, x, y, z in atoms:
        lines.append(f"  {sp}  {x:.10f}  {y:.10f}  {z:.10f}")
    lines.append("K_POINTS automatic")
    lines.append(f"  {KPTS_M[0]} {KPTS_M[1]} {KPTS_M[2]} 0 0 0")
    lines.append("")
    return "\n".join(lines)


def snapshot_name(idx):
    return f"snap{idx:02d}"


def delta_tag(d):
    if d > 0:
        return f"+{d:.2f}"
    if d < 0:
        return f"-{abs(d):.2f}"
    return f"{0.0:.2f}"


def main():
    if not os.path.exists(MD_IN):
        print(f"ERROR: {MD_IN} not found (run generate_inputs.py --stage "
              "ensemble first)", file=sys.stderr)
        return 1
    if not os.path.exists(MD_OUT):
        print(f"ERROR: {MD_OUT} not found (MD has not run yet)",
              file=sys.stderr)
        return 1

    n_steps = count_md_steps(MD_OUT)
    print(f"MD trajectory: {n_steps} ATOMIC_POSITIONS blocks found in "
          f"{MD_OUT}")
    if n_steps < MIN_STEPS_REQUIRED:
        print(f"ERROR: only {n_steps} steps (< {MIN_STEPS_REQUIRED} "
              "required) -- MD did not progress far enough", file=sys.stderr)
        return 1

    a, c, species_order, atoms0, if_pos = parse_md_input(MD_IN)
    nat = len(atoms0)
    mobile_idx = [i for i, fp in enumerate(if_pos) if fp == (1, 1, 1)]
    na_idx = [i for i, (sp, *_) in enumerate(atoms0) if sp == "Na"]
    if len(na_idx) != 1:
        print(f"ERROR: expected exactly 1 Na atom, found {len(na_idx)}",
              file=sys.stderr)
        return 1
    na_idx = na_idx[0]
    na_x, na_y = atoms0[na_idx][1], atoms0[na_idx][2]
    print(f"parsed md_cage/pw.in: nat={nat}, {len(mobile_idx)} mobile "
          f"(water) atoms, Na at index {na_idx} (x={na_x:.6f}, "
          f"y={na_y:.6f}), a={a:.4f} A, c={c:.4f} A")

    snapshots = parse_md_trajectory(MD_OUT, nat)
    print(f"parsed {len(snapshots)} full-nat trajectory snapshots from "
          f"{MD_OUT}")
    if len(snapshots) <= EQUIL_STEPS:
        print(f"ERROR: only {len(snapshots)} parsed snapshots, <= "
              f"{EQUIL_STEPS} equilibration steps to discard", file=sys.stderr)
        return 1

    kept = snapshots[EQUIL_STEPS:]
    n_take = min(N_SNAPSHOTS, len(kept))
    if n_take < N_SNAPSHOTS:
        print(f"WARNING: only {len(kept)} post-equilibration snapshots "
              f"available; taking {n_take} instead of {N_SNAPSHOTS}")
    if n_take == 1:
        pick_local = [0]
    else:
        pick_local = [round(i * (len(kept) - 1) / (n_take - 1))
                      for i in range(n_take)]
    picks = [(EQUIL_STEPS + j, kept[j]) for j in pick_local]

    notes_template = [
        "! Set M (cage-ensemble frozen-snapshot Na scan): water frozen at a",
        "! thermal snapshot from the set-L 290 K Born-Oppenheimer MD "
        "(jobs_ensemble/md_cage), Na rigidly displaced along c; CoO2 "
        "framework from the original static set-G geometry. Exact set-G",
        "! electronic settings (K_POINTS 6 6 2, conv_thr 1e-7) for direct",
        "! comparability with the rigid-cage (set G) and adiabatic",
        "! (set J) E(delta) curves.",
    ]

    jobs = []
    for snap_i, (step_idx, snap_atoms) in enumerate(picks, start=1):
        # sanity: species order in the trajectory block must match pw.in
        mismatches = sum(1 for i in range(nat)
                         if snap_atoms[i][0] != atoms0[i][0])
        if mismatches:
            print(f"WARNING: snapshot {snap_i} (step {step_idx}) has "
                  f"{mismatches} species-order mismatches vs pw.in -- "
                  "skipping", file=sys.stderr)
            continue
        for d in DELTAS_M:
            new_atoms = []
            for i in range(nat):
                sp = atoms0[i][0]
                if i == na_idx:
                    z = (c / 2.0 + d) / c
                    new_atoms.append((sp, na_x, na_y, z))
                elif i in mobile_idx:
                    _, x, y, z = snap_atoms[i]
                    new_atoms.append((sp, x, y, z))
                else:
                    new_atoms.append(atoms0[i])
            notes = list(notes_template) + [
                f"! snapshot {snap_i:02d}: MD trajectory index {step_idx} "
                f"(step {step_idx} of {n_steps} parsed steps), delta = "
                f"{d:+.2f} A."
            ]
            text = pw_input_scf(a, c, species_order, new_atoms, d, notes)
            name = f"{snapshot_name(snap_i)}_d{delta_tag(d)}"
            job_dir = os.path.join(ENSEMBLE_DIR, name)
            os.makedirs(job_dir, exist_ok=True)
            with open(os.path.join(job_dir, "pw.in"), "w") as f:
                f.write(text)
            jobs.append(dict(name=name, set="M", element="Na", cell="s3",
                             c=c, delta=d, snapshot=snap_i,
                             md_step_index=step_idx, water=True,
                             type="scf", kpts=list(KPTS_M), nat=nat))

    print(f"wrote {len(jobs)} set-M SCF jobs "
          f"({len(picks)} snapshots x {len(DELTAS_M)} deltas)")

    # update manifest_ensemble.json (append, keep the set-L md_cage entry)
    if os.path.exists(MANIFEST):
        with open(MANIFEST) as f:
            manifest = json.load(f)
    else:
        manifest = dict(pseudos=PSEUDOS, zval=ZVAL, jobs=[])
    existing_names = {j["name"] for j in manifest.get("jobs", [])}
    n_new = 0
    for j in jobs:
        if j["name"] not in existing_names:
            manifest.setdefault("jobs", []).append(j)
            n_new += 1
    manifest["ensemble_md_steps_parsed"] = n_steps
    manifest["ensemble_snapshots_used"] = len(picks)
    with open(MANIFEST, "w") as f:
        json.dump(manifest, f, indent=2)
    print(f"manifest_ensemble.json updated (+{n_new} new job entries)")
    return 0


if __name__ == "__main__":
    sys.exit(main())

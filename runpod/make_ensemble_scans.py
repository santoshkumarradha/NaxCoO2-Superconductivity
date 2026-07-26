#!/usr/bin/env python3
"""Pod-side stage M generator for the "cage-ensemble" set L/M.

Run from runpod/ AFTER the set-L walker MDs (jobs_ensemble/md_w1..md_w4)
have finished or progressed far enough (run_ensemble.sh gates this on
>= 1000 MD steps per walker; walkers below the gate are skipped --
survivors-only degradation):

    python3 make_ensemble_scans.py

For each SURVIVING walker, parses its Born-Oppenheimer MD trajectory
(jobs_ensemble/md_w<K>/pw.out, water-only dynamics at Na delta = 0, see
generate_inputs.py generate_ensemble), discards the first 400 steps as
equilibration, takes 3 snapshots evenly spaced over the remainder (~steps
600/900/1200 of a full 1200-step run), and for each snapshot writes 7 SCF
inputs with the water atoms FROZEN at the snapshot coordinates and Na
displaced to delta in {-0.45, -0.30, -0.15, 0.00, +0.15, +0.30, +0.45} A
along c from the midplane (the CoO2 framework is taken from the original,
static, input geometry). Full production run: 4 walkers x 3 snapshots x 7
deltas = 84 SCF jobs, dirs jobs_ensemble/w<K>s<J>_d<+/-0.NN>/pw.in.

These SCFs use the EXACT set-G electronic settings (K_POINTS 6 6 2,
conv_thr 1e-7, same 60/480 Ry cutoffs, mv 0.02 smearing, nspin=2 with the
set-G starting_magnetization) for direct comparability with the rigid-cage
(set G) and adiabatic (set J) results. Warm restart: the delta != 0 inputs
carry startingwfc = 'file' / startingpot = 'file'; run_ensemble.sh phase 3
runs each snapshot unit on one GPU, delta = 0.00 first, then copies its
./out (charge density + wavefunctions) into the remaining jobs of the unit
in |delta|-increasing order. (pw.x falls back to the standard atomic guess
with a warning if the files are absent, so a missing delta=0 out is a
slowdown, not a failure.)

Pure stdlib -- no numpy, no third-party deps (deliberately, even though
numpy is available in the NGC image: the geometry math needed here is
trivial fractional-coordinate arithmetic, so avoiding the import keeps this
script robust to any pod/image quirk).
"""
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ENSEMBLE_DIR = os.path.join(HERE, "jobs_ensemble")
MANIFEST = os.path.join(HERE, "manifest_ensemble.json")

BOHR = 0.529177210903  # Angstrom

EQUIL_STEPS = 400         # discard per walker as equilibration
N_SNAP_PER_WALKER = 3     # evenly spaced over the remaining trajectory
DELTAS_M = [-0.45, -0.30, -0.15, 0.00, 0.15, 0.30, 0.45]  # Angstrom
MIN_WALKER_STEPS = 1000   # must match run_ensemble.sh's phase-1 gate

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
    """Return (a_ang, c_ang, species_order, atoms, if_pos) from a set-L
    walker pw.in, where atoms = [(species, x, y, z), ...] in file order and
    if_pos is a parallel list of (ix,iy,iz) ints (None entries default to
    (1,1,1), but the set-L inputs always write explicit flags)."""
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
    ATOMIC_POSITIONS block found)."""
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
    """Cheap step count: number of ATOMIC_POSITIONS markers; matches
    run_ensemble.sh's phase-1 gate logic so this script and the bash gate
    agree."""
    if not os.path.exists(path):
        return 0
    with open(path) as f:
        txt = f.read()
    return len(re.findall(r"ATOMIC_POSITIONS", txt))


# -------------------------------------------------------------- pw.x input --
def nbnd_for(nelec):
    return nelec // 2 + max(8, nelec // 6)


def pw_input_scf(a, c, species_order, atoms, notes, warm=False):
    """SCF pw.in text, EXACT set-G electronic settings. `atoms` already has
    Na at z = c/2 + delta and water at the frozen snapshot coordinates.
    warm=True adds startingwfc/startingpot='file' (delta != 0 jobs; the
    delta=0 ./out is copied in by run_ensemble.sh phase 3 before the run)."""
    nelec = sum(ZVAL[sp] for sp, *_ in atoms)
    nbnd = nbnd_for(nelec)
    electrons = [
        "&ELECTRONS",
        "  conv_thr = 1.0d-7",
        "  mixing_beta = 0.3",
        "  mixing_mode = 'local-TF'",
        "  electron_maxstep = 300",
    ]
    if warm:
        electrons += ["  startingwfc = 'file'", "  startingpot = 'file'"]
    electrons += ["/"]
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
    ] + electrons + ["ATOMIC_SPECIES"]
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


def delta_tag(d):
    if d > 0:
        return f"+{d:.2f}"
    if d < 0:
        return f"-{abs(d):.2f}"
    return f"{0.0:.2f}"


def process_walker(k, jobs):
    """Extract snapshots from walker k and write its set-M SCF inputs.
    Returns the number of snapshots written (0 if the walker is skipped)."""
    md_dir = os.path.join(ENSEMBLE_DIR, f"md_w{k}")
    md_in = os.path.join(md_dir, "pw.in")
    md_out = os.path.join(md_dir, "pw.out")
    if not os.path.exists(md_in):
        print(f"walker {k}: no pw.in, skip")
        return 0
    n_steps = count_md_steps(md_out)
    if n_steps < MIN_WALKER_STEPS:
        print(f"walker {k}: only {n_steps} steps (< {MIN_WALKER_STEPS}), "
              "skip (survivors-only)")
        return 0

    a, c, species_order, atoms0, if_pos = parse_md_input(md_in)
    nat = len(atoms0)
    mobile_idx = [i for i, fp in enumerate(if_pos) if fp == (1, 1, 1)]
    na_idx = [i for i, (sp, *_) in enumerate(atoms0) if sp == "Na"]
    if len(na_idx) != 1:
        print(f"walker {k}: ERROR expected exactly 1 Na, found "
              f"{len(na_idx)}, skip", file=sys.stderr)
        return 0
    na_idx = na_idx[0]
    na_x, na_y = atoms0[na_idx][1], atoms0[na_idx][2]

    snapshots = parse_md_trajectory(md_out, nat)
    if len(snapshots) <= EQUIL_STEPS + 1:
        print(f"walker {k}: only {len(snapshots)} parsed snapshots "
              f"(<= {EQUIL_STEPS} equilibration), skip", file=sys.stderr)
        return 0
    last = len(snapshots) - 1
    # 3 snapshots at ~1/3, 2/3, 3/3 of the post-equilibration trajectory
    # (~steps 600/900/1200 of a full 1200-step run)
    pick = sorted({EQUIL_STEPS + round((last - EQUIL_STEPS) * f)
                   for f in (1.0 / 3.0, 2.0 / 3.0, 1.0)})
    print(f"walker {k}: {n_steps} steps, {len(snapshots)} parsed blocks, "
          f"snapshot indices {pick}")

    n_written = 0
    for sj, step_idx in enumerate(pick, start=1):
        snap_atoms = snapshots[step_idx]
        mismatches = sum(1 for i in range(nat)
                         if snap_atoms[i][0] != atoms0[i][0])
        if mismatches:
            print(f"walker {k} snap {sj} (step {step_idx}): {mismatches} "
                  "species-order mismatches vs pw.in, skip", file=sys.stderr)
            continue
        unit = f"w{k}s{sj}"
        for d in DELTAS_M:
            new_atoms = []
            for i in range(nat):
                sp = atoms0[i][0]
                if i == na_idx:
                    new_atoms.append((sp, na_x, na_y, (c / 2.0 + d) / c))
                elif i in mobile_idx:
                    _, x, y, z = snap_atoms[i]
                    new_atoms.append((sp, x, y, z))
                else:
                    new_atoms.append(atoms0[i])
            notes = [
                "! Set M (cage-ensemble frozen-snapshot Na scan): water",
                f"! frozen at walker-{k} MD step {step_idx} (of {n_steps} "
                "parsed; 290 K svr BOMD,",
                "! jobs_ensemble/md_w*), Na rigidly displaced along c; CoO2",
                "! framework from the original static set-G geometry. Exact",
                "! set-G electronic settings (K_POINTS 6 6 2, conv_thr 1e-7)",
                "! for direct comparability with the rigid-cage (set G) and",
                "! adiabatic (set J) E(delta) curves. delta != 0 jobs warm-",
                "! start from the unit's delta=0 charge density (out copied",
                "! in by run_ensemble.sh phase 3).",
                f"! unit {unit}: delta = {d:+.2f} A.",
            ]
            text = pw_input_scf(a, c, species_order, new_atoms, notes,
                                warm=(d != 0.0))
            name = f"{unit}_d{delta_tag(d)}"
            job_dir = os.path.join(ENSEMBLE_DIR, name)
            os.makedirs(job_dir, exist_ok=True)
            with open(os.path.join(job_dir, "pw.in"), "w") as f:
                f.write(text)
            jobs.append(dict(name=name, set="M", element="Na", cell="s3",
                             c=c, delta=d, walker=k, snapshot=sj, unit=unit,
                             md_step_index=step_idx, water=True,
                             type="scf", kpts=list(KPTS_M), nat=nat))
        n_written += 1
    return n_written


def main():
    walker_ids = sorted(
        int(m.group(1)) for m in
        (re.match(r"md_w(\d+)$", n) for n in os.listdir(ENSEMBLE_DIR))
        if m) if os.path.isdir(ENSEMBLE_DIR) else []
    if not walker_ids:
        print("ERROR: no jobs_ensemble/md_w* walkers found (run "
              "generate_inputs.py --stage ensemble first)", file=sys.stderr)
        return 1

    jobs = []
    n_snaps = 0
    survivors = 0
    for k in walker_ids:
        n = process_walker(k, jobs)
        if n > 0:
            survivors += 1
            n_snaps += n
    if survivors == 0:
        print("ERROR: no surviving walkers (all < "
              f"{MIN_WALKER_STEPS} steps)", file=sys.stderr)
        return 1

    print(f"wrote {len(jobs)} set-M SCF jobs ({survivors} surviving "
          f"walkers, {n_snaps} snapshots x {len(DELTAS_M)} deltas)")

    # update manifest_ensemble.json (append, keep the set-L walker entries)
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
    manifest["ensemble_surviving_walkers"] = survivors
    manifest["ensemble_snapshots_used"] = n_snaps
    with open(MANIFEST, "w") as f:
        json.dump(manifest, f, indent=2)
    print(f"manifest_ensemble.json updated (+{n_new} new job entries)")
    return 0


if __name__ == "__main__":
    sys.exit(main())

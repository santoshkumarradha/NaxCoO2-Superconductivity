#!/usr/bin/env python3
"""
parse_magnetization.py -- harvest spin data from the RunPod QE batch.

Points at runpod/jobs/ (directories named  <El>_<cell>_c<c>_d<delta>, each
containing a spin-polarized pw.x output  pw.out) and extracts

  * total magnetization      (Bohr mag/cell, last printed = converged)
  * absolute magnetization   (Bohr mag/cell)
  * per-atom moments from the 'Magnetic moment per site' block (new QE
    format 'atom   1 (R=...) charge= ... magn= ...' and the older
    'atom:    1    charge: ...    magn: ...' format are both supported),
    aggregated per element,

and tabulates moment vs (element, c, delta) as a CSV -- ready to run the
moment the DFT batch syncs back, to test the prediction that the gallery
polarization is strongest just past carrier turn-on and dies where SC lives
(c = 9.9 A).

Usage:
  python3 parse_magnetization.py                      # default ../runpod/jobs
  python3 parse_magnetization.py /path/to/jobs -o mag.csv
  python3 parse_magnetization.py --selftest           # fabricated pw.out check

No dependencies beyond the standard library.
"""

import argparse
import csv
import re
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

# --------------------------------------------------------------------------
# Regexes
# --------------------------------------------------------------------------
RE_DIR = re.compile(
    r"^(?P<el>[A-Z][a-z]?)_(?P<cell>.+)_c(?P<c>[\d.]+)_d(?P<delta>[\d.]+)$")
RE_TOTAL = re.compile(
    r"total magnetization\s*=\s*(-?[\d.]+)\s*Bohr mag/cell")
RE_ABS = re.compile(
    r"absolute magnetization\s*=\s*(-?[\d.]+)\s*Bohr mag/cell")
# new-style per-site line:  atom   1 (R=0.413)  charge= 10.2327  magn=  1.4581
RE_SITE_NEW = re.compile(
    r"atom\s+(\d+)\s*\(R=\s*[\d.]+\)\s*charge=\s*(-?[\d.]+)\s*magn=\s*(-?[\d.]+)")
# old-style:  atom:    1    charge:   10.2327    magn:    1.4581    constr: ...
RE_SITE_OLD = re.compile(
    r"atom:\s*(\d+)\s+charge:\s*(-?[\d.]+)\s+magn:\s*(-?[\d.]+)")
# species table in pw.out:      1           Co  tau(   1) = ( ... )
RE_TAU = re.compile(
    r"^\s*(\d+)\s+([A-Z][a-z]?\d*)\s+tau\(\s*\d+\)\s*=", re.MULTILINE)
# fallback: ATOMIC_POSITIONS block in pw.in
RE_POS_LINE = re.compile(r"^\s*([A-Z][a-z]?\d*)\s+(-?[\d.]+)\s+(-?[\d.]+)")


def species_from_pwout(text):
    """Atom-index -> element map from the pw.out coordinate table."""
    out = {}
    for m in RE_TAU.finditer(text):
        idx, el = int(m.group(1)), re.sub(r"\d+$", "", m.group(2))
        out.setdefault(idx, el)
    return out


def species_from_pwin(path):
    out = {}
    if not path.exists():
        return out
    lines = path.read_text(errors="replace").splitlines()
    in_block = False
    idx = 0
    for ln in lines:
        if ln.strip().upper().startswith("ATOMIC_POSITIONS"):
            in_block = True
            continue
        if in_block:
            m = RE_POS_LINE.match(ln)
            if not m:
                break
            idx += 1
            out[idx] = re.sub(r"\d+$", "", m.group(1))
    return out


def parse_pwout(path):
    """Return dict with converged (last) magnetizations and per-atom moments,
    or None if the file has no spin output at all."""
    text = path.read_text(errors="replace")
    totals = RE_TOTAL.findall(text)
    absol = RE_ABS.findall(text)
    if not totals and not absol:
        return None

    # per-site moments: take the LAST block (converged); both formats
    sites = RE_SITE_NEW.findall(text) or RE_SITE_OLD.findall(text)
    per_atom = {}
    for idx_s, charge_s, magn_s in sites:  # later blocks overwrite earlier
        per_atom[int(idx_s)] = (float(charge_s), float(magn_s))

    spec = species_from_pwout(text) or species_from_pwin(
        path.with_name("pw.in"))

    by_el = defaultdict(list)
    for idx, (_q, m) in sorted(per_atom.items()):
        by_el[spec.get(idx, f"atom{idx}")].append(m)

    return {
        "total_mag": float(totals[-1]) if totals else float("nan"),
        "abs_mag": float(absol[-1]) if absol else float("nan"),
        "moments_by_element": dict(by_el),
        "job_done": "JOB DONE" in text,
    }


def scan_jobs(jobs_dir):
    rows = []
    elements_seen = set()
    jobs_dir = Path(jobs_dir)
    for d in sorted(p for p in jobs_dir.iterdir() if p.is_dir()):
        m = RE_DIR.match(d.name)
        if not m:
            continue
        out = d / "pw.out"
        if not out.exists():
            continue
        parsed = parse_pwout(out)
        if parsed is None:
            continue
        row = {
            "job": d.name,
            "alkali": m.group("el"),
            "cell": m.group("cell"),
            "c_A": float(m.group("c")),
            "delta_A": float(m.group("delta")),
            "total_mag_muB": parsed["total_mag"],
            "abs_mag_muB": parsed["abs_mag"],
            "converged": parsed["job_done"],
        }
        for el, moms in parsed["moments_by_element"].items():
            elements_seen.add(el)
            row[f"m_{el}_sum_muB"] = sum(moms)
            row[f"m_{el}_max_muB"] = max(moms, key=abs)
        rows.append(row)
    return rows, elements_seen


def write_csv(rows, elements_seen, out_path):
    base = ["job", "alkali", "cell", "c_A", "delta_A",
            "total_mag_muB", "abs_mag_muB", "converged"]
    extra = sorted(k for el in elements_seen
                   for k in (f"m_{el}_sum_muB", f"m_{el}_max_muB"))
    fields = base + extra
    with open(out_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, restval="")
        w.writeheader()
        for r in sorted(rows, key=lambda r: (r["alkali"], r["cell"],
                                             r["c_A"], r["delta_A"])):
            w.writerow(r)


def print_table(rows):
    if not rows:
        print("No spin-polarized pw.out files found (batch not synced yet?).")
        return
    hdr = f"{'job':28s}{'c(A)':>7s}{'d(A)':>7s}{'M_tot':>8s}{'M_abs':>8s}" \
          f"{'m_Co(max)':>11s}{'ok':>4s}"
    print(hdr)
    print("-" * len(hdr))
    for r in sorted(rows, key=lambda r: (r["alkali"], r["cell"],
                                         r["c_A"], r["delta_A"])):
        mco = r.get("m_Co_max_muB", "")
        mco = f"{mco:11.3f}" if mco != "" else f"{'--':>11s}"
        print(f"{r['job']:28s}{r['c_A']:7.2f}{r['delta_A']:7.2f}"
              f"{r['total_mag_muB']:8.3f}{r['abs_mag_muB']:8.3f}{mco}"
              f"{'  y' if r['converged'] else '  n':>4s}")


# --------------------------------------------------------------------------
# Self test on a fabricated pw.out snippet
# --------------------------------------------------------------------------
FAKE_PWOUT_NEW = """\
     Program PWSCF v.7.3.1 starts ...
   Cartesian axes

     site n.     atom                  positions (alat units)
         1           Co  tau(   1) = (   0.0000000   0.0000000   0.0000000  )
         2           O   tau(   2) = (   0.3333333   0.6666667   0.1391304  )
         3           O   tau(   3) = (   0.6666667   0.3333333   0.8608696  )
         4           Na  tau(   4) = (   0.3333333   0.6666667   0.5724638  )

     iteration #  1     ecut=    60.00 Ry
     total magnetization       =     1.20 Bohr mag/cell
     absolute magnetization    =     1.90 Bohr mag/cell

     Magnetic moment per site  (integrated on atomic sphere of radius R)
     atom   1 (R=0.357)  charge= 14.5843  magn=  0.6100
     atom   2 (R=0.339)  charge=  6.0110  magn=  0.0300
     atom   3 (R=0.339)  charge=  6.0110  magn=  0.0250
     atom   4 (R=0.443)  charge=  8.5000  magn= -0.0450

     iteration # 12     ecut=    60.00 Ry
     total magnetization       =     0.51 Bohr mag/cell
     absolute magnetization    =     0.71 Bohr mag/cell

     Magnetic moment per site  (integrated on atomic sphere of radius R)
     atom   1 (R=0.357)  charge= 14.5843  magn=  0.4210
     atom   2 (R=0.339)  charge=  6.0110  magn=  0.0150
     atom   3 (R=0.339)  charge=  6.0110  magn=  0.0140
     atom   4 (R=0.443)  charge=  8.5000  magn= -0.0320

     convergence has been achieved in  12 iterations

     JOB DONE.
"""

FAKE_PWOUT_OLD = """\
     Program PWSCF v.6.4 starts ...
     site n.     atom                  positions (alat units)
         1           Co  tau(   1) = (   0.0000000   0.0000000   0.0000000  )
         2           O   tau(   2) = (   0.3333333   0.6666667   0.1391304  )
         3           O   tau(   3) = (   0.6666667   0.3333333   0.8608696  )
         4           Li  tau(   4) = (   0.3333333   0.6666667   0.5724638  )

     total magnetization       =     0.02 Bohr mag/cell
     absolute magnetization    =     0.05 Bohr mag/cell
 atom:    1    charge:   14.2000    magn:    0.0300    constr:    0.0000
 atom:    2    charge:    6.0000    magn:    0.0100    constr:    0.0000
 atom:    3    charge:    6.0000    magn:    0.0100    constr:    0.0000
 atom:    4    charge:    2.1000    magn:   -0.0050    constr:    0.0000

     JOB DONE.
"""


def selftest():
    ok = True
    with tempfile.TemporaryDirectory() as td:
        jobs = Path(td) / "jobs"
        for name, blob in (("Na_1x1_c9.9_d0.50", FAKE_PWOUT_NEW),
                           ("Li_1x1_c6.9_d0.15", FAKE_PWOUT_OLD)):
            d = jobs / name
            d.mkdir(parents=True)
            (d / "pw.out").write_text(blob)
        # a job dir without output must be skipped silently
        (jobs / "Na_1x1_c5.5_d0.00").mkdir()

        rows, elements = scan_jobs(jobs)
        checks = [
            ("two jobs parsed", len(rows) == 2),
        ]
        na = next(r for r in rows if r["job"] == "Na_1x1_c9.9_d0.50")
        li = next(r for r in rows if r["job"] == "Li_1x1_c6.9_d0.15")
        checks += [
            ("last (converged) total mag kept", na["total_mag_muB"] == 0.51),
            ("last absolute mag kept", na["abs_mag_muB"] == 0.71),
            ("dir name -> c", na["c_A"] == 9.9 and li["c_A"] == 6.9),
            ("dir name -> delta", na["delta_A"] == 0.50),
            ("new-format Co moment", abs(na["m_Co_max_muB"] - 0.4210) < 1e-9),
            ("new-format Na moment (2DEG, antiparallel)",
             abs(na["m_Na_sum_muB"] - (-0.0320)) < 1e-9),
            ("new-format O sum", abs(na["m_O_sum_muB"] - 0.029) < 1e-9),
            ("old-format Co moment", abs(li["m_Co_max_muB"] - 0.03) < 1e-9),
            ("old-format Li moment", abs(li["m_Li_sum_muB"] + 0.005) < 1e-9),
            ("JOB DONE detected", na["converged"] and li["converged"]),
        ]
        out = Path(td) / "mag.csv"
        write_csv(rows, elements, out)
        content = out.read_text().splitlines()
        checks.append(("csv has header + 2 rows", len(content) == 3))
        checks.append(("csv header has element columns",
                       "m_Co_max_muB" in content[0]))

        for name, passed in checks:
            print(f"  [{'PASS' if passed else 'FAIL'}] {name}")
            ok &= passed
    print("selftest:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def main():
    ap = argparse.ArgumentParser(
        description="Extract magnetization vs (element, c, delta) from QE "
                    "pw.out files in runpod/jobs/")
    default_jobs = Path(__file__).resolve().parent.parent / "runpod" / "jobs"
    ap.add_argument("jobs_dir", nargs="?", default=str(default_jobs))
    ap.add_argument("-o", "--output", default="magnetization_table.csv",
                    help="CSV output path (default: ./magnetization_table.csv)")
    ap.add_argument("--selftest", action="store_true",
                    help="run parser checks on fabricated pw.out snippets")
    args = ap.parse_args()

    if args.selftest:
        sys.exit(selftest())

    jobs_dir = Path(args.jobs_dir)
    if not jobs_dir.is_dir():
        sys.exit(f"jobs dir not found: {jobs_dir}")
    rows, elements = scan_jobs(jobs_dir)
    print_table(rows)
    if rows:
        write_csv(rows, elements, args.output)
        print(f"\nwrote {args.output}  ({len(rows)} jobs)")
        print("Prediction to test: absolute/Co+gallery moments peak just past"
              " the carrier\nturn-on (intermediate c ~ 8.4 A) and fall at"
              " c = 9.9 A where SC lives.")


if __name__ == "__main__":
    main()

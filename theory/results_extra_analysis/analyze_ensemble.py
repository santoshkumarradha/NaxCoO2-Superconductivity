# /// script
# requires-python = ">=3.10"
# dependencies = ["numpy", "scipy", "pandas", "matplotlib"]
# ///
"""
Set L/M -- the "cage-ensemble" thermal Na E(delta) at c = 9.9 A.

Physics recap: set J showed that a single ADIABATIC (fully-relaxed) water
cage is the wrong statistical object for the Na potential -- the 4-water
H-bond network has multiple local minima spanning >= 0.4 eV at fixed Na, and
the Na mode (~200 fs) is fast against H-bond reorganization (~ps, Jalarvo
QENS). The physical Na potential is set by the INSTANTANEOUS cage
configuration, drawn from a thermal ensemble, not the single deepest
adiabatic minimum. Set L/M tests this directly:

  SET L  runpod/results_ensemble/md_cage/pw.out          Born-Oppenheimer MD
         of the set-G hydrate cell at 290 K (svr thermostat, dt=10 Ry a.u.,
         nstep=4000, ~1.9 ps). Water free, Na + CoO2 frozen at delta = 0.
  SET M  runpod/results_ensemble/snap<NN>_d<+/-0.NN>/pw.out   10 snapshots
         (first 1000 MD steps discarded as equilibration, evenly spaced over
         the remainder) x 7 Na displacements delta in
         {-0.45,-0.30,-0.15,0.00,+0.15,+0.30,+0.45} A, water frozen at each
         snapshot's coordinates, exact set-G electronic settings (K_POINTS
         6x6x2, conv_thr 1e-7) for direct comparability.

This script parses every snap*_d* pw.out, builds a per-snapshot E(delta)
curve (referenced to that snapshot's own delta=0, mirroring the set-J
convention: the well SHAPE, not the large per-configuration energy offset,
is the physically meaningful comparison), fits each to the quartic Landau
form E(delta) = E0 + alpha*delta^2 + beta*delta^4, and reports per-snapshot
alpha, well depth, and minimum position. The ensemble MEAN curve (averaging
E(delta) - E(0) across all snapshots that have a delta=0 point, at each
common delta) is the physical, statistically-averaged Na potential; it is
compared against:
  - the rigid-cage set-G curve (runpod/results_extra/jobs_extra/
    Na_s3hyd_c9.9_d*): frozen gas-phase cage, single configuration.
  - the adiabatic set-J minimum (runpod/results_mobile/
    mobile_Na_s3hyd_c9.9_d*): fully-relaxed water at each pinned Na.

Robust to missing jobs (partial harvests): any snapshot/delta pw.out that
is absent, unconverged, or missing its own delta=0 anchor is skipped with a
warning rather than raising.

Reuses theory/effective_model.py (RY_EV, MASS_AMU, HBAR2_OVER_2AMU) and the
parsing conventions of analyze_setJ.py. Every number traceable to the
"!    total energy" lines in each pw.out.

Run with:
    uv run theory/results_extra_analysis/analyze_ensemble.py   (or)
    theory/.venv/bin/python theory/results_extra_analysis/analyze_ensemble.py
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
THEORY = REPO / "theory"
EXTRA = REPO / "runpod" / "results_extra" / "jobs_extra"
MOBILE = REPO / "runpod" / "results_mobile"
ENSEMBLE = REPO / "runpod" / "results_ensemble"
sys.path.insert(0, str(THEORY))

import effective_model as em  # noqa: E402

RY_EV = em.RY_EV
H2M = em.HBAR2_OVER_2AMU            # hbar^2/(2 amu), eV*A^2
M_NA = dict(em.MASS_AMU)["Na"]

E_RE = re.compile(r"^!\s+total energy\s+=\s+(-?[\d.]+)\s+Ry", re.M)
DELTAS_M = [-0.45, -0.30, -0.15, 0.00, 0.15, 0.30, 0.45]

out = {}
lines = []


def log(s=""):
    print(s)
    lines.append(s)


def read_out(path: Path):
    """Return (E_eV, converged) from a pw.out, or None if absent/unparsable."""
    if not path.exists():
        return None
    txt = path.read_text()
    hits = E_RE.findall(txt)
    if not hits:
        return None
    return float(hits[-1]) * RY_EV, "JOB DONE" in txt


def delta_tag(d):
    if d > 0:
        return f"+{d:.2f}"
    if d < 0:
        return f"-{abs(d):.2f}"
    return f"{0.0:.2f}"


def discover_snapshots(base: Path):
    """Snapshot indices present under results_ensemble/, from dir names
    snap<NN>_d*.  Robust to any subset actually harvested."""
    if not base.exists():
        return []
    idxs = set()
    for p in base.iterdir():
        m = re.match(r"snap(\d+)_d", p.name)
        if m:
            idxs.add(int(m.group(1)))
    return sorted(idxs)


def quartic_fit(d, E_meV):
    """Fit E(delta) [meV] = E0 + alpha*delta^2 + beta*delta^4 [meV/A^2,
    meV/A^4]. Returns (E0, alpha, beta, d_min, well_depth_meV) or None if
    fewer than 4 points. d_min from a dense scan of the fitted quartic over
    the sampled delta range (avoids spurious extrapolated minima)."""
    if len(d) < 4:
        return None
    A = np.column_stack([np.ones_like(d), d ** 2, d ** 4])
    coef, *_ = np.linalg.lstsq(A, E_meV, rcond=None)
    E0, alpha, beta = coef
    dd = np.linspace(d.min(), d.max(), 2001)
    EE = E0 + alpha * dd ** 2 + beta * dd ** 4
    i = int(np.argmin(EE))
    return dict(E0_meV=float(E0), alpha_meV_per_A2=float(alpha),
               beta_meV_per_A4=float(beta), d_min_A=float(dd[i]),
               well_depth_meV=float(EE[i] - E0))


log("=" * 78)
log("SET L/M -- cage-ensemble thermal Na E(delta), sqrt3xsqrt3, c = 9.9 A")
log("=" * 78)

if not ENSEMBLE.exists():
    log(f"\nWARNING: {ENSEMBLE} not found -- nothing harvested yet. "
        "Run the rsync harvest command first (see run_ensemble.sh / the "
        "campaign launch report), then re-run this script.")

snap_idxs = discover_snapshots(ENSEMBLE)
log(f"\nfound {len(snap_idxs)} snapshot indices under {ENSEMBLE}: {snap_idxs}")

per_snapshot = {}
for si in snap_idxs:
    d_ok, E_abs = [], []
    for d in DELTAS_M:
        p = ENSEMBLE / f"snap{si:02d}_d{delta_tag(d)}" / "pw.out"
        r = read_out(p)
        if r is None:
            log(f"  snap{si:02d} d={d:+.2f}: missing/unparsable pw.out, skip")
            continue
        E, done = r
        if not done:
            log(f"  snap{si:02d} d={d:+.2f}: pw.out present but not JOB DONE, skip")
            continue
        d_ok.append(d)
        E_abs.append(E)
    if not d_ok:
        log(f"  snap{si:02d}: no converged deltas, dropping snapshot")
        continue
    d_ok = np.array(d_ok)
    E_abs = np.array(E_abs)
    if 0.0 not in d_ok:
        log(f"  snap{si:02d}: missing its own delta=0.00 anchor, dropping "
            "(no valid own-reference)")
        continue
    i0 = int(np.where(np.isclose(d_ok, 0.0))[0][0])
    E_meV = (E_abs - E_abs[i0]) * 1e3
    fit = quartic_fit(d_ok, E_meV)
    per_snapshot[si] = dict(d=d_ok.tolist(), E_meV=E_meV.tolist(),
                             E_abs_eV=E_abs.tolist(), n_deltas=len(d_ok),
                             fit=fit)

log(f"\n{len(per_snapshot)}/{len(snap_idxs)} snapshots usable "
    "(have delta=0 anchor + >=1 other point).")

if per_snapshot:
    log("")
    log(f"{'snap':>6} {'n_d':>4} {'alpha(meV/A^2)':>16} {'beta(meV/A^4)':>15} "
        f"{'d_min(A)':>9} {'depth(meV)':>11}")
    for si, rec in sorted(per_snapshot.items()):
        f = rec["fit"]
        if f is None:
            log(f"{si:6d} {rec['n_deltas']:4d} "
                f"{'(< 4 pts, no fit)':>16}")
            continue
        log(f"{si:6d} {rec['n_deltas']:4d} {f['alpha_meV_per_A2']:16.1f} "
            f"{f['beta_meV_per_A4']:15.1f} {f['d_min_A']:9.3f} "
            f"{f['well_depth_meV']:11.1f}")

# --- ensemble mean curve: average E(delta)-E(0) across snapshots at each ---
# common delta (only snapshots that report that delta contribute).
log("")
log("Ensemble mean curve (average over available snapshots at each delta):")
mean_curve = {}
for d in DELTAS_M:
    vals = []
    for si, rec in per_snapshot.items():
        if d in rec["d"]:
            vals.append(rec["E_meV"][rec["d"].index(d)])
    if vals:
        mean_curve[d] = dict(mean_meV=float(np.mean(vals)),
                             std_meV=float(np.std(vals)), n=len(vals))
if mean_curve:
    log(f"{'delta(A)':>9} {'mean(meV)':>10} {'std(meV)':>9} {'n':>3}")
    for d in sorted(mean_curve):
        r = mean_curve[d]
        log(f"{d:9.2f} {r['mean_meV']:10.1f} {r['std_meV']:9.1f} {r['n']:3d}")
    d_arr = np.array(sorted(mean_curve))
    E_arr = np.array([mean_curve[d]["mean_meV"] for d in d_arr])
    mean_fit = quartic_fit(d_arr, E_arr)
    if mean_fit:
        log("")
        log(f"  ensemble mean quartic fit: alpha = "
            f"{mean_fit['alpha_meV_per_A2']:.1f} meV/A^2, "
            f"beta = {mean_fit['beta_meV_per_A4']:.1f} meV/A^4, "
            f"well depth = {mean_fit['well_depth_meV']:.1f} meV at "
            f"delta_min = {mean_fit['d_min_A']:.3f} A.")
else:
    mean_fit = None
    log("  (no common deltas across usable snapshots -- nothing to average)")

# --- compare against set G (rigid, single config) and set J (adiabatic) ---
log("")
log("=" * 78)
log("Comparison: rigid-cage (set G) vs adiabatic-minimum (set J) vs "
    "cage-ensemble mean (set L/M)")
log("=" * 78)


def scan_family(base: Path, stem: str, ds):
    d_ok, E_ok = [], []
    for d in ds:
        r = read_out(base / f"{stem}_d{d:.2f}" / "pw.out")
        if r is None:
            continue
        E, done = r
        if not done:
            continue
        d_ok.append(d)
        E_ok.append(E)
    if not d_ok:
        return np.array([]), np.array([])
    d_ok = np.array(d_ok)
    E_ok = np.array(E_ok)
    return d_ok, (E_ok - E_ok[0]) * 1e3


DS_STD = [0.00, 0.15, 0.30, 0.50, 0.75, 1.00]
dR, ER = scan_family(EXTRA, "Na_s3hyd_c9.9", DS_STD)
dM, EM = scan_family(MOBILE, "mobile_Na_s3hyd_c9.9", DS_STD)

if len(dR):
    iR = int(np.argmin(ER))
    log(f"  RIGID (set G)    : min {ER[iR]:.0f} meV at delta={dR[iR]:.2f} A "
        "(gas-phase cage, single frozen configuration).")
else:
    log("  RIGID (set G)    : not found (harvest runpod/results_extra first)")

if len(dM):
    iM = int(np.argmin(EM))
    log(f"  ADIABATIC (set J): min {EM[iM]:.0f} meV at delta={dM[iM]:.2f} A "
        "(fully-relaxed water, single deepest configuration).")
else:
    log("  ADIABATIC (set J): not found (harvest runpod/results_mobile first)")

if mean_fit:
    log(f"  ENSEMBLE (set L/M): mean-curve depth "
        f"{mean_fit['well_depth_meV']:.0f} meV at "
        f"delta={mean_fit['d_min_A']:.2f} A ({len(per_snapshot)} thermal "
        "snapshots, 290 K).")
    alphas = [rec["fit"]["alpha_meV_per_A2"] for rec in per_snapshot.values()
              if rec["fit"]]
    if alphas:
        log(f"  per-snapshot alpha spread: {min(alphas):.0f} to "
            f"{max(alphas):.0f} meV/A^2 (mean {np.mean(alphas):.0f}, std "
            f"{np.std(alphas):.0f}) -- this spread IS the statistical-"
            "dynamical softening/broadening set J's timescale argument "
            "predicts (Na is fast vs ~ps H-bond reorganization, so it "
            "samples an ensemble of instantaneous cage curvatures rather "
            "than one adiabatic well).")
else:
    log("  ENSEMBLE (set L/M): not enough harvested snapshots for a mean "
        "curve yet.")

# --- write outputs ------------------------------------------------------
out = dict(
    n_snapshots_found=len(snap_idxs),
    n_snapshots_usable=len(per_snapshot),
    per_snapshot={str(k): v for k, v in per_snapshot.items()},
    mean_curve=mean_curve,
    mean_fit=mean_fit,
    rigid=dict(d=dR.tolist(), E_meV=ER.tolist()),
    adiabatic=dict(d=dM.tolist(), E_meV=EM.tolist()),
)
(HERE / "summary_ensemble.json").write_text(json.dumps(out, indent=2))
(HERE / "summary_ensemble.txt").write_text("\n".join(lines))
log("")
log(f"wrote {HERE/'summary_ensemble.json'}")
log(f"wrote {HERE/'summary_ensemble.txt'}")

#!/usr/bin/env python3
"""Local post-processing for the NaxCoO2 interlayer-2DEG / anharmonic-rattler
job set (run AFTER syncing jobs/ + manifest.json back from RunPod).

Per (element, cell, c):
  * fit E(delta) = E0 + alpha*delta^2 + beta*delta^4 to the SCF scan,
  * solve the 1D Schroedinger equation in the fitted potential (alkali mass,
    finite differences) -> level splitting E10 = hbar*omega_eff,
  * extract the alkali-band deformation potential I = d eps_Gamma / d delta
    from the Fermi-referenced Gamma eigenvalue blocks (verbosity='high'),
  * read N(0) from dos.x output (per spin), n_2DEG from Bader ACF.dat,
  * lambda = N(0) I^2 / (M omega_eff^2),  Allen-Dynes Tc with a single-mode
    alpha^2F and mu* = 0.10,
and produce triple_plot.pdf: alpha(c), n_2DEG(c), Tc(c) with the experimental
anchors (c = 6.9 A non-SC, c = 9.9 A SC with Tc = 4.5 K).

Self-test (no data needed):  python3 postprocess.py --selftest
"""
import argparse
import json
import os
import re
import sys

import numpy as np

RY_EV = 13.605693122994
K_AMU = 2.09008e-3        # hbar^2 / (2 m) for m = 1 amu, in eV * Angstrom^2
EV_K = 11604.518           # 1 eV in Kelvin
MASS = {"Na": 22.98977, "Li": 6.941}
ZVAL = {"Na": 9, "Li": 3}
A_LAT = {"Na": 2.888, "Li": 2.82}

# ------------------------------------------------------------------ parsing -
E_RE = re.compile(r"^!\s+total energy\s+=\s+(-?[\d.]+)\s+Ry", re.M)
EF_RE = re.compile(r"the Fermi energy is\s+(-?[\d.]+)\s+ev")
K_RE = re.compile(r"k\s*=\s*(-?[\d.]+)\s+(-?[\d.]+)\s+(-?[\d.]+).*bands \(ev\)")
NUM_RE = re.compile(r"-?\d+\.\d+")


def parse_pw_out(path):
    """Return dict(E_ry, ef, bands={spin: [eigs at Gamma]}) or None."""
    try:
        with open(path) as f:
            lines = f.read().splitlines()
    except OSError:
        return None
    txt = "\n".join(lines)
    if "JOB DONE" not in txt:
        return None
    e = E_RE.findall(txt)
    ef = EF_RE.findall(txt)
    out = dict(E_ry=float(e[-1]) if e else None,
               ef=float(ef[-1]) if ef else None, bands={})
    spin = "up"
    i = 0
    while i < len(lines):
        ln = lines[i]
        if "SPIN UP" in ln:
            spin = "up"
        elif "SPIN DOWN" in ln:
            spin = "down"
        else:
            m = K_RE.search(ln)
            if m and all(abs(float(g)) < 1e-4 for g in m.groups()):
                eigs, j = [], i + 1
                while j < len(lines):
                    lj = lines[j].strip()
                    if lj == "" and eigs:
                        break
                    if "occupation numbers" in lj:
                        break
                    eigs += [float(x) for x in NUM_RE.findall(lj)]
                    j += 1
                if eigs:
                    out["bands"][spin] = eigs  # keep LAST Gamma block per spin
                i = j
                continue
        i += 1
    return out


def parse_dos(path):
    """dos.x fildos -> N(0) per spin per cell (states/eV) at E_F."""
    try:
        with open(path) as f:
            header = f.readline()
            data = np.loadtxt(f)
    except OSError:
        return None
    m = re.search(r"EFermi\s*=\s*(-?[\d.]+)", header)
    if m is None or data.size == 0:
        return None
    ef = float(m.group(1))
    e = data[:, 0]
    dos_tot = data[:, 1] + data[:, 2] if data.shape[1] >= 4 else data[:, 1]
    return float(np.interp(ef, e, dos_tot)) / 2.0  # per spin


def parse_acf(path, alkali_index):
    """Bader ACF.dat -> charge on the alkali atom (1-based index)."""
    try:
        with open(path) as f:
            lines = f.readlines()
    except OSError:
        return None
    for ln in lines:
        parts = ln.split()
        if len(parts) >= 5 and parts[0] == str(alkali_index):
            return float(parts[4])
    return None


# ------------------------------------------------------------------ physics -
def fit_even_poly(deltas, energies_ev):
    """E(d) = E0 + a d^2 + b d^4 (b fixed to 0 with < 4 points)."""
    d = np.asarray(deltas, float)
    e = np.asarray(energies_ev, float)
    quartic = len(d) >= 4
    cols = [np.ones_like(d), d ** 2] + ([d ** 4] if quartic else [])
    coef, *_ = np.linalg.lstsq(np.stack(cols, axis=1), e, rcond=None)
    e0, alpha = coef[0], coef[1]
    beta = coef[2] if quartic else 0.0
    return float(e0), float(alpha), float(beta)


def solve_schrodinger(alpha, beta, mass_amu, nlev=4, n=2001):
    """1D FD Schroedinger in V(x) = alpha x^2 + beta x^4 (eV, Angstrom).
    Returns (levels[nlev] in eV, <x^2> of the ground state in A^2), or None
    if the potential is unbound.  hbar*omega_eff = hbar^2/(2 M <x^2>) equals
    E1-E0 for a harmonic well but stays finite and physical in a deep double
    well, where E1-E0 is an exponentially small tunnel splitting."""
    k = K_AMU / mass_amu
    if beta > 0:
        x0 = np.sqrt(max(-alpha, 0.0) / (2 * beta)) if alpha < 0 else 0.0
        xmax = max(2.5, 3.0 * x0, (1.0 / beta) ** 0.25 + x0)
    elif alpha > 0:
        xmax = 0.9 * np.sqrt(alpha / (2 * abs(beta))) if beta < 0 else \
            max(2.5, 3.0 * np.sqrt(0.5 / alpha))
    else:
        return None  # unbound
    x = np.linspace(-xmax, xmax, n)
    h = x[1] - x[0]
    v = alpha * x ** 2 + beta * x ** 4
    diag = 2 * k / h ** 2 + v
    off = np.full(n - 1, -k / h ** 2)
    try:
        from scipy.linalg import eigh_tridiagonal
        w, vec = eigh_tridiagonal(diag, off, select="i",
                                  select_range=(0, nlev - 1))
    except ImportError:
        m = np.diag(diag) + np.diag(off, 1) + np.diag(off, -1)
        w, vec = np.linalg.eigh(m)
        w, vec = w[:nlev], vec[:, :nlev]
    psi0 = vec[:, 0] / np.sqrt(np.sum(vec[:, 0] ** 2))
    x2 = float(np.sum(psi0 ** 2 * x ** 2))
    return [float(wi) for wi in w], x2


def deformation_potential(deltas, parsed_list, window=(-1.0, 4.0)):
    """Max |d(eps - E_F)/d(delta)| over bands near E_F at Gamma (eV/Ang)."""
    best = None
    for spin in ("up", "down"):
        series = [(d, p) for d, p in zip(deltas, parsed_list)
                  if p and p.get("ef") is not None and spin in p["bands"]]
        if len(series) < 2:
            continue
        nb = min(len(p["bands"][spin]) for _, p in series)
        d = np.array([s[0] for s in series])
        eps = np.array([[p["bands"][spin][b] - p["ef"] for b in range(nb)]
                        for _, p in series])
        for b in range(nb):
            if not (window[0] <= eps[:, b].mean() <= window[1]):
                continue
            slope = np.polyfit(d, eps[:, b], 1)[0]
            if best is None or abs(slope) > abs(best[0]):
                best = (slope, spin, b)
    return best  # (slope eV/A, spin, band index) or None


def allen_dynes_tc(lam, wlog_K, mu=0.10):
    if lam <= mu * (1 + 0.62 * lam) or lam <= 0:
        return 0.0
    return (wlog_K / 1.2) * np.exp(-1.04 * (1 + lam)
                                   / (lam - mu * (1 + 0.62 * lam)))


# ----------------------------------------------------------------- pipeline -
def analyze(root, mu):
    with open(os.path.join(root, "manifest.json")) as f:
        manifest = json.load(f)
    jobs = manifest["jobs"]
    scf = {}
    for j in jobs:
        if j["type"] == "scf":
            scf.setdefault((j["element"], j["cell"], j["c"]), []).append(j)
    nscf = {(j["element"], j["cell"], j["c"]): j
            for j in jobs if j["type"] == "nscf"}

    results = []
    for key in sorted(scf):
        el, cell, c = key
        grp = sorted(scf[key], key=lambda j: j["delta"])
        parsed = [parse_pw_out(os.path.join(root, "jobs", j["name"], "pw.out"))
                  for j in grp]
        pts = [(j["delta"], p) for j, p in zip(grp, parsed)
               if p and p["E_ry"] is not None]
        r = dict(element=el, cell=cell, c=c, n_scf=len(pts))
        if len(pts) >= 2:
            d = [x[0] for x in pts]
            e = [x[1]["E_ry"] * RY_EV for x in pts]
            e0, alpha, beta = fit_even_poly(d, e)
            r.update(alpha=alpha, beta=beta,
                     delta_min=float(d[int(np.argmin(e))]))
            sol = solve_schrodinger(alpha, beta, MASS[el])
            if sol:
                lev, x2 = sol
                hw_eff = (K_AMU / MASS[el]) / x2  # hbar^2 / (2 M <x^2>), eV
                r["levels"] = lev
                r["E10"] = lev[1] - lev[0]        # tunnel/level splitting
                r["x2_gs"] = x2
                r["hw_eff"] = hw_eff
                r["omega_meV"] = 1e3 * hw_eff
            dp = deformation_potential([x[0] for x in pts],
                                       [x[1] for x in pts])
            if dp:
                r.update(I=abs(dp[0]), I_spin=dp[1], I_band=dp[2])
        # N(0) from the dense-k dos.x run
        nj = nscf.get(key)
        if nj:
            n0 = parse_dos(os.path.join(root, "jobs", nj["name"], "pw.dos"))
            if n0:
                r["N0"] = n0  # states/eV/cell per spin
        # Bader charge of the alkali at the E(delta) minimum geometry
        ref = nscf.get(key, {}).get("parent")
        if ref is None and len(pts) >= 1:
            ref = min(zip(grp, parsed),
                      key=lambda t: t[1]["E_ry"] if t[1] and t[1]["E_ry"]
                      is not None else np.inf)[0]["name"]
        if ref:
            q = parse_acf(os.path.join(root, "jobs", ref, "ACF.dat"),
                          alkali_index=grp[0]["nat"])  # alkali is last atom
            if q is not None:
                donated = ZVAL[el] - q
                a = A_LAT[el] * (np.sqrt(3) if cell == "s3" else 1.0)
                area = np.sqrt(3) / 2 * a ** 2
                r["q_bader"] = q
                r["n2deg_cm2"] = donated / area * 1e16
        # lambda = N(0) I^2 / (M omega_eff^2) and Allen-Dynes Tc
        if all(k in r for k in ("N0", "I", "hw_eff")) and r["hw_eff"] > 0:
            m_w2 = r["hw_eff"] ** 2 / (2 * K_AMU / MASS[el])  # eV/A^2
            lam = r["N0"] * r["I"] ** 2 / m_w2
            r["lambda"] = lam
            r["Tc"] = allen_dynes_tc(lam, r["hw_eff"] * EV_K, mu)
        results.append(r)
    return results


def make_plot(results, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    series = {}
    for r in results:
        series.setdefault((r["element"], r["cell"]), []).append(r)
    fig, axes = plt.subplots(3, 1, figsize=(6, 10), sharex=True)
    style = {("Na", "1x1"): ("o-", "tab:blue", "Na 1x1"),
             ("Li", "1x1"): ("s-", "tab:green", "Li 1x1"),
             ("Na", "s3"): ("^--", "tab:orange", "Na x=1/3")}
    panels = [("alpha", r"$\alpha$ (eV/$\mathrm{\AA}^2$)"),
              ("n2deg_cm2", r"$n_{\rm 2DEG}$ (cm$^{-2}$)"),
              ("Tc", r"$T_c$ (K)")]
    for ax, (key, ylab) in zip(axes, panels):
        for skey, rs in sorted(series.items()):
            rs = sorted((r for r in rs if key in r), key=lambda r: r["c"])
            if not rs:
                continue
            mk, col, lab = style.get(skey, ("x-", "gray", str(skey)))
            ax.plot([r["c"] for r in rs], [r[key] for r in rs], mk,
                    color=col, label=lab)
        ax.set_ylabel(ylab)
        ax.axvline(6.9, color="crimson", ls=":", lw=1)
        ax.axvline(9.9, color="teal", ls=":", lw=1)
        ax.grid(alpha=0.3)
    axes[0].axhline(0, color="k", lw=0.8)
    axes[0].set_title(r"anchors: $c=6.9\,\mathrm{\AA}$ non-SC (red), "
                      r"$c=9.9\,\mathrm{\AA}$ SC (teal)", fontsize=9)
    axes[2].plot([9.9], [4.5], "*", ms=16, color="teal",
                 label="exp: $T_c=4.5$ K")
    axes[2].plot([6.9], [0.0], "x", ms=10, color="crimson",
                 label="exp: non-SC")
    axes[2].set_xlabel(r"$c$ = CoO$_2$ plane spacing ($\mathrm{\AA}$)")
    for ax in axes:
        ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path)
    print(f"wrote {path}")


# ----------------------------------------------------------------- selftest -
def selftest():
    ok = True

    def check(name, cond, detail=""):
        nonlocal ok
        print(f"  [{'PASS' if cond else 'FAIL'}] {name} {detail}")
        ok = ok and cond

    print("selftest: quartic fit recovery")
    rng = np.random.default_rng(0)
    a_true, b_true = -0.8, 1.2
    d = np.array([0.0, 0.15, 0.30, 0.50, 0.75, 1.00])
    e = -1000.0 + a_true * d ** 2 + b_true * d ** 4 \
        + rng.normal(0, 1e-6, d.size)
    _, a, b = fit_even_poly(d, e)
    check("alpha", abs(a - a_true) < 0.01 * abs(a_true), f"{a:.4f}")
    check("beta", abs(b - b_true) < 0.01 * abs(b_true), f"{b:.4f}")

    print("selftest: harmonic oscillator (Na mass)")
    alpha = 0.5  # eV/A^2 -> hbar*omega = 2 sqrt(alpha * K)
    lev, x2 = solve_schrodinger(alpha, 0.0, MASS["Na"])
    hw_exact = 2 * np.sqrt(alpha * K_AMU / MASS["Na"])
    e10 = lev[1] - lev[0]
    check("E1-E0 = hbar*omega", abs(e10 - hw_exact) < 0.01 * hw_exact,
          f"num {e10*1e3:.4f} meV vs exact {hw_exact*1e3:.4f} meV")
    e21 = lev[2] - lev[1]
    check("equal spacing", abs(e21 - e10) < 0.01 * hw_exact)
    hw_eff = (K_AMU / MASS["Na"]) / x2
    check("hw_eff = hbar^2/(2M<x^2>) = hbar*omega in harmonic limit",
          abs(hw_eff - hw_exact) < 0.01 * hw_exact, f"{hw_eff*1e3:.4f} meV")

    # light 1 amu mass: with Na's mass this test well's tunnel splitting is
    # ~1e-20 eV, below double precision -- physical, but not testable.
    print("selftest: double well (alpha=-0.5, beta=0.5, depth 0.125 eV, m=1)")
    lev, x2 = solve_schrodinger(-0.5, 0.5, 1.0)
    split = lev[1] - lev[0]
    gap = lev[2] - lev[1]
    check("ordered levels", lev[0] < lev[1] < lev[2] < lev[3])
    check("tunnel splitting < intra-well gap", 0 < split < gap,
          f"split {split*1e3:.4f} meV, gap {gap*1e3:.4f} meV")
    check("ground state above well bottom", lev[0] > -0.125)
    check("<x^2> ~ well position^2 (x0=0.707)", 0.2 < x2 < 1.0,
          f"<x^2> = {x2:.3f} A^2")

    print("selftest: Allen-Dynes sanity")
    tc = allen_dynes_tc(1.0, 100.0, 0.10)
    check("Tc(lambda=1, wlog=100K) in (1, 20) K", 1 < tc < 20, f"{tc:.2f} K")
    check("Tc = 0 below threshold", allen_dynes_tc(0.1, 100.0, 0.10) == 0.0)

    print("selftest:", "ALL PASS" if ok else "FAILURES")
    return 0 if ok else 1


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root",
                   default=os.path.dirname(os.path.abspath(__file__)))
    p.add_argument("--mu", type=float, default=0.10, help="mu* (default 0.10)")
    p.add_argument("--out", default="triple_plot.pdf")
    p.add_argument("--selftest", action="store_true")
    args = p.parse_args()
    if args.selftest:
        return selftest()

    results = analyze(args.root, args.mu)
    out_json = os.path.join(args.root, "results.json")
    with open(out_json, "w") as f:
        json.dump(results, f, indent=2)
    print(f"wrote {out_json}")
    hdr = (f"{'series':>10} {'c':>5} {'alpha':>8} {'beta':>7} {'w(meV)':>7} "
           f"{'I(eV/A)':>8} {'N0':>6} {'n2D(cm-2)':>10} {'lam':>6} {'Tc(K)':>7}")
    print(hdr)
    for r in results:
        print(f"{r['element'] + '-' + r['cell']:>10} {r['c']:>5.1f} "
              f"{r.get('alpha', float('nan')):>8.3f} "
              f"{r.get('beta', float('nan')):>7.3f} "
              f"{r.get('omega_meV', float('nan')):>7.2f} "
              f"{r.get('I', float('nan')):>8.3f} "
              f"{r.get('N0', float('nan')):>6.2f} "
              f"{r.get('n2deg_cm2', float('nan')):>10.3e} "
              f"{r.get('lambda', float('nan')):>6.3f} "
              f"{r.get('Tc', float('nan')):>7.3f}")
    plottable = [r for r in results
                 if any(k in r for k in ("alpha", "n2deg_cm2", "Tc"))]
    if plottable:
        try:
            make_plot(results, os.path.join(args.root, args.out))
        except ImportError:
            print("matplotlib not installed -- skipping triple_plot.pdf "
                  "(pip install matplotlib)")
    else:
        print("no plottable data yet (sync pw.out / pw.dos / ACF.dat first)")
    return 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""
stoner_model.py -- Stoner phase diagram of the interlayer (gallery) 2DEG in
Na_xCoO2.yH2O, built on the SSH4 surface-charge model of Radha & Lambrecht,
SciPost Phys. 10, 057 (2021).

Physics
-------
The gallery 2DEG carrier fraction per alkali follows the paper's two-surface-
band overlap formula

    q_A / q_CoO2 = (Gamma - 2*delta) / (Gamma + 2*delta),   q_A + q_CoO2 = 1
    =>  q_A = (Gamma - 2*delta) / (2*Gamma)   if Gamma > 2*delta, else 0
    Gamma = Gamma_A + Gamma_Co = 6(|t_A^xy| + |t_CoO2^xy|)

where Gamma_A, Gamma_Co are the half-widths of the alkali and CoO2 surface
bands and 2*delta is the splitting of their centers (paper, Sec. 3 + App. E).

Two knobs move Gamma_A:
  (1) alkali dilution x: the alkali-alkali spacing is d(x) = a/sqrt(x); the
      in-plane hopping decays exponentially, t(d) = t0*exp(-(d-a)/lambda).
      (paper, App. C: half coverage in rows -> bandwidth ~3x smaller and the
      2DEG empties; here modelled by the exponential distance law -- ANSATZ.)
  (2) gallery opening c: opening the CoO2-CoO2 spacing decompresses the
      alkali sp_z orbital and develops the double well, growing the
      full-coverage bandwidth Gamma0(c) -- ANSATZ, linear, calibrated so the
      x=0.35 carrier turn-on falls between the monolayer hydrate (c=6.9 A,
      not SC) and the bilayer hydrate (c=9.9 A, SC), the one experimental
      anchor pair we are allowed.

Stoner criterion on the gallery band (2D parabolic band bottom):
    N(0) = nu / Gamma_A(x,c)  per spin per alkali site,  nu = sqrt(3)/(2*pi)
    polarized  <=>  I_s * N(0) > 1
so polarization is strongest just past carrier turn-on (narrow band, high
N(0)) and dies as the band widens.  Stoner enhancement S = 1/(1 - I_s N(0))
defines a paramagnon-boosted strip next to the polarized one.

Every parameter is tagged PAPER / DERIVED / ANSATZ(range) in PARAMS below.

Outputs: printed tables (Li-vs-Na threshold inequality, sensitivity scan,
region boundaries) and spin_phase_diagram.pdf.
"""

import argparse
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D

# ----------------------------------------------------------------------------
# Parameters.  source: PAPER = Radha & Lambrecht SciPost Phys. 10, 057 (2021);
# DERIVED = algebra on PAPER numbers; ANSATZ = assumed here, with range.
# ----------------------------------------------------------------------------
PARAMS = {
    # in-plane hoppings of the two surface Wannier bands (eV)
    "t_Li_xy": (0.60, "PAPER: |t_Li^xy| = 0.6 eV (negative, pi-type sp_z)"),
    "t_Co_xy": (0.09, "PAPER: t_CoO2^xy = 0.09 eV"),
    "t_Na_xy": (0.90, "ANSATZ [0.75, 1.10] eV: > t_Li because Na 3s/3p is more"
                      " extended (<r>_3s/<r>_2s ~ 1.25); PAPER App. C only says"
                      " 'larger lateral hopping between Na', no number"),
    # on-site half-splitting delta (eV)
    "delta_Li": (0.414, "DERIVED: q_Li = 0.4e (PAPER Bader) with"
                        " Gamma = 6*(0.6+0.09) = 4.14 eV =>"
                        " 2*delta = Gamma*(1-2q) = 0.828 eV"),
    "delta_Na_over_Li": (1.00, "ANSATZ [1.00, 1.15]: electronegativity (Li 0.98"
                               " vs Na 0.93) and IP (5.39 vs 5.14 eV) place the"
                               " Na level slightly HIGHER, i.e. delta_Na >="
                               " delta_Li.  NOTE: this is the honest direction;"
                               " it works AGAINST Na, and is overwhelmed by the"
                               " Gamma advantage."),
    # in-plane lattice constants (Angstrom)
    "a_Li": (2.82, "experimental LiCoO2 in-plane a"),
    "a_Na": (2.84, "experimental Na_xCoO2 in-plane a (hydrate ~2.82-2.84)"),
    # hopping decay lengths (Angstrom)
    "lam_Li": (0.840, "ANSATZ [0.72, 0.96]: lambda = 1/kappa,"
                      " kappa = sqrt(2*IP), IP_Li = 5.39 eV"),
    "lam_Na": (0.861, "ANSATZ [0.74, 0.98]: same rule, IP_Na = 5.14 eV"),
    # Stoner
    "nu": (np.sqrt(3.0) / (2.0 * np.pi),
           "DERIVED: 2D triangular band bottom, m* = hbar^2/(3 t d^2) =>"
           " N(0) = nu/Gamma_A per spin per alkali site, nu = 0.2757"),
    "I_s": (1.3, "ANSATZ [1.0, 1.6] eV: effective Stoner I of the gallery band"
                 " (alkali sp_z with Co-d admixture; bare Co-d I ~ 0.9 eV)."
                 " Strip exists only if I_s > (2*delta - Gamma_Co)/nu"
                 " = 1.05 eV -- printed below."),
    "S_star": (3.0, "ANSATZ [2, 5]: Stoner enhancement above which we call the"
                    " metal 'paramagnon-boosted'"),
    # gallery-opening calibration Gamma0(c) = A + B*c  (eV, c in Angstrom)
    "c_on": (8.0, "ANSATZ [7.2, 9.2] A: c at which the x=0.35 gallery band"
                  " turns on; only constrained to lie in (6.9, 9.9) by the"
                  " monolayer/bilayer hydrate on/off anchor"),
    "c_full": (9.9, "bilayer hydrate spacing; Gamma0(c_full) = 6|t_Na^xy|"
                    " (fully developed orbital) -- ANSATZ"),
}

P = {k: v[0] for k, v in PARAMS.items()}


# ----------------------------------------------------------------------------
# Model functions
# ----------------------------------------------------------------------------
def f_dilution(x, a, lam):
    """Exponential in-plane hopping suppression at coverage x (ANSATZ).
    d(x) = a/sqrt(x); f = exp(-(d-a)/lam); f(1) = 1."""
    d = a / np.sqrt(np.asarray(x, dtype=float))
    return np.exp(-(d - a) / lam)


def gamma_alkali_halfwidth(x, t_xy, a, lam, scale=1.0):
    """Half-width of the alkali surface band, Gamma_A = 6|t^xy| f(x) * scale."""
    return 6.0 * t_xy * f_dilution(x, a, lam) * scale


def q_carrier(gamma_alk, gamma_co, two_delta):
    """Paper's overlap formula: q per alkali; 0 below threshold."""
    gam = gamma_alk + gamma_co
    q = (gam - two_delta) / (2.0 * gam)
    return np.where(gam > two_delta, q, 0.0)


def x_star(t_xy, a, lam, two_delta, gamma_co):
    """Coverage threshold: Gamma_A(x*) = 2*delta - Gamma_Co (q -> 0)."""
    rhs = two_delta - gamma_co
    if rhs <= 0:
        return 0.0  # threshold cleared at any coverage
    f_needed = rhs / (6.0 * t_xy)
    if f_needed >= 1.0:
        return np.inf  # never clears, even at full coverage
    return (1.0 + (lam / a) * np.log(1.0 / f_needed)) ** -2


def gamma0_of_c(c):
    """Gallery-opening law Gamma0(c) = A + B*c (ANSATZ, calibrated).
    Constraints: Gamma_A(x=0.35, c_on) = 2*delta  and
                 Gamma0(c_full) = 6|t_Na^xy|."""
    two_delta = 2.0 * P["delta_Li"] * P["delta_Na_over_Li"]
    f35 = f_dilution(0.35, P["a_Na"], P["lam_Na"])
    g_on = (two_delta - 6.0 * P["t_Co_xy"]) / f35   # Gamma0 needed at turn-on
    g_full = 6.0 * P["t_Na_xy"]
    B = (g_full - g_on) / (P["c_full"] - P["c_on"])
    A = g_full - B * P["c_full"]
    return np.maximum(A + B * np.asarray(c, dtype=float), 0.0)


def classify(x, c):
    """Region code on a (x, c) grid.
    0 no 2DEG | 1 polarized 2DEG | 2 paramagnon-boosted PM | 3 plain PM."""
    X, C = np.meshgrid(x, c)
    w_alk = gamma0_of_c(C) * f_dilution(X, P["a_Na"], P["lam_Na"])
    two_delta = 2.0 * P["delta_Li"] * P["delta_Na_over_Li"]
    gamma_co = 6.0 * P["t_Co_xy"]
    q = q_carrier(w_alk, gamma_co, two_delta)
    with np.errstate(divide="ignore"):
        F = P["I_s"] * P["nu"] / w_alk          # Stoner factor I_s*N(0)
        S = np.where(F < 1.0, 1.0 / (1.0 - F), np.inf)
    region = np.full(X.shape, 3, dtype=int)
    region[(q > 0) & (F < 1.0) & (S >= P["S_star"])] = 2
    region[(q > 0) & (F >= 1.0)] = 1
    region[q == 0.0] = 0
    return X, C, q, F, S, region, w_alk


# ----------------------------------------------------------------------------
# Reporting
# ----------------------------------------------------------------------------
def li_vs_na_report():
    gamma_co = 6.0 * P["t_Co_xy"]
    rows = []
    for el, t, a, lam, dscale in (
        ("Li", P["t_Li_xy"], P["a_Li"], P["lam_Li"], 1.0),
        ("Na", P["t_Na_xy"], P["a_Na"], P["lam_Na"], P["delta_Na_over_Li"]),
    ):
        two_delta = 2.0 * P["delta_Li"] * dscale
        xs = x_star(t, a, lam, two_delta, gamma_co)
        w35 = gamma_alkali_halfwidth(0.35, t, a, lam)
        q35 = float(q_carrier(w35, gamma_co, two_delta))
        rows.append((el, 6 * t, two_delta, xs, w35 + gamma_co, q35))

    print("=" * 78)
    print("Li vs Na carrier threshold  (q > 0  <=>  Gamma(x) > 2*delta)")
    print("  Gamma(x) = 6|t_Co^xy| + 6|t_A^xy| exp(-(a/sqrt(x)-a)/lambda_A)")
    print("  x*_A = [1 + (lambda_A/a) ln(6|t_A^xy| / (2 delta_A"
          " - 6|t_Co^xy|))]^-2")
    print("-" * 78)
    print(f"{'':3s}{'6|t_xy| (eV)':>14s}{'2*delta (eV)':>14s}{'x* (turn-on)':>14s}"
          f"{'Gamma(0.35)':>13s}{'q(x=0.35)':>12s}")
    for el, g6, td, xs, g35, q35 in rows:
        print(f"{el:3s}{g6:14.2f}{td:14.3f}{xs:14.3f}{g35:13.3f}{q35:12.3f}")
    print("-" * 78)
    qli, qna = rows[0][5], rows[1][5]
    print(f"Baseline verdict: x*_Li = {rows[0][3]:.3f} vs x*_Na = {rows[1][3]:.3f};"
          f"  q_Na/q_Li at x=0.35 = {qna / max(qli, 1e-12):.1f}"
          if qli > 0 else
          f"Baseline verdict: Li below threshold at x=0.35, Na above.")
    print("HONESTY NOTE: with the baseline ansatz Li sits *marginally above*")
    print("its threshold at x=0.35 (q_Li ~ %.3f, i.e. within model error of 0);"
          % qli)
    print("the robust statements are the inequality x*_Li > x*_Na and the")
    print("factor ~%.0f carrier ratio.  The sensitivity scan below shows how"
          % (qna / max(qli, 1e-12)))
    print("often q_Li(0.35) = 0 exactly across the stated ansatz ranges.")
    return rows


def sensitivity_scan(n=20000, seed=7):
    rng = np.random.default_rng(seed)
    gamma_co = 6.0 * P["t_Co_xy"]
    lam_li = rng.uniform(0.72, 0.96, n)
    lam_na = rng.uniform(0.74, 0.98, n)
    t_na = rng.uniform(0.75, 1.10, n)
    d_li = rng.uniform(0.38, 0.45, n)          # delta_Li around DERIVED 0.414
    d_ratio = rng.uniform(1.00, 1.15, n)       # delta_Na / delta_Li

    w_li = gamma_alkali_halfwidth(0.35, P["t_Li_xy"], P["a_Li"], lam_li)
    w_na = gamma_alkali_halfwidth(0.35, t_na, P["a_Na"], lam_na)
    q_li = q_carrier(w_li, gamma_co, 2.0 * d_li)
    q_na = q_carrier(w_na, gamma_co, 2.0 * d_li * d_ratio)

    p_li0 = np.mean(q_li == 0.0)
    p_na0 = np.mean(q_na == 0.0)
    both = (q_li == 0.0) & (q_na > 0.0)
    ratio = np.median(q_na[(q_li > 0)] / q_li[q_li > 0])
    print("=" * 78)
    print(f"Sensitivity scan over ansatz ranges (n = {n}):")
    print(f"  P[q_Li(0.35) = 0]                    = {p_li0:5.1%}")
    print(f"  P[q_Na(0.35) = 0]                    = {p_na0:5.1%}")
    print(f"  P[Li dead AND Na alive]              = {np.mean(both):5.1%}")
    print(f"  median q_Na/q_Li where both nonzero  = {ratio:.1f}")
    print(f"Reading: Li straddles its threshold at x = 0.35 -- {p_li0:.0%} of the"
          f" ansatz volume\nkills it outright and the rest leaves a ~{ratio:.0f}x"
          f" weaker 2DEG -- while Na is\nabove threshold in"
          f" {1 - p_na0:.0%} of the same volume.")


def boundaries_at_x(x=0.35):
    """Critical c values along the x cut."""
    two_delta = 2.0 * P["delta_Li"] * P["delta_Na_over_Li"]
    gamma_co = 6.0 * P["t_Co_xy"]
    f = float(f_dilution(x, P["a_Na"], P["lam_Na"]))
    crit = {
        "carrier turn-on (q>0)": (two_delta - gamma_co) / f,
        "Stoner boundary (I_s N(0)=1)": P["I_s"] * P["nu"] / f,
        "paramagnon edge (S=S*)": P["I_s"] * P["nu"]
        * P["S_star"] / (P["S_star"] - 1.0) / f,
    }
    # invert gamma0_of_c (linear)
    g1, g2 = gamma0_of_c(np.array([8.0, 9.0]))
    B = g2 - g1
    A = g1 - 8.0 * B
    print("=" * 78)
    print(f"Region boundaries along x = {x} (c in A; Gamma0(c) linear ANSATZ):")
    out = {}
    for name, g0 in crit.items():
        c_val = (g0 - A) / B
        out[name] = c_val
        print(f"  {name:34s} c = {c_val:5.2f}")
    for c_ref, label in ((6.9, "monolayer hydrate (no SC)"),
                         (9.9, "bilayer hydrate (SC, Tc = 4.5 K)")):
        w = gamma0_of_c(c_ref) * f
        q = float(q_carrier(w, gamma_co, two_delta))
        if q == 0.0:
            print(f"  at c = {c_ref} ({label}): q = 0 (no carriers; Stoner"
                  f" factor moot)")
            continue
        F = P["I_s"] * P["nu"] / w
        S = 1.0 / (1.0 - F) if F < 1 else np.inf
        print(f"  at c = {c_ref} ({label}): q = {q:.3f}, I_s N(0) = {F:.2f}, "
              f"S = {'inf (ordered)' if not np.isfinite(S) else f'{S:.1f}'}")
    return out


# ----------------------------------------------------------------------------
# Figure (palette: dataviz reference palette, light mode)
# ----------------------------------------------------------------------------
INK = "#0b0b0b"
INK2 = "#52514e"
FILL = {0: "#ebebe8",   # no 2DEG      (neutral)
        1: "#f5c4c4",   # polarized    (red tint)
        2: "#f8e2b4",   # paramagnon   (yellow tint)
        3: "#cde2fb"}   # plain PM     (blue tint, sequential step 100)
EDGE = {"turnon": "#3d3c39", "stoner": "#b23230", "sstar": "#8a6100"}
BLUE = "#2a78d6"


def make_figure(outfile):
    x = np.linspace(0.12, 0.80, 400)
    c = np.linspace(5.5, 10.6, 400)
    X, C, q, F, S, region, w_alk = classify(x, c)

    fig = plt.figure(figsize=(11.2, 4.9))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.55, 1.0], hspace=0.12,
                          wspace=0.28)
    ax = fig.add_subplot(gs[:, 0])
    axq = fig.add_subplot(gs[0, 1])
    axs = fig.add_subplot(gs[1, 1], sharex=axq)

    # --- (a) phase diagram ---------------------------------------------------
    cmap = ListedColormap([FILL[i] for i in range(4)])
    ax.pcolormesh(X, C, region, cmap=cmap, vmin=-0.5, vmax=3.5, shading="auto",
                  rasterized=True)

    two_delta = 2.0 * P["delta_Li"] * P["delta_Na_over_Li"]
    gamma_co = 6.0 * P["t_Co_xy"]
    lv_on = two_delta - gamma_co
    lv_st = P["I_s"] * P["nu"]
    lv_ss = lv_st * P["S_star"] / (P["S_star"] - 1.0)
    for lv, key, ls in ((lv_on, "turnon", "-"), (lv_st, "stoner", "-"),
                        (lv_ss, "sstar", "--")):
        ax.contour(X, C, w_alk, levels=[lv], colors=EDGE[key],
                   linewidths=1.4, linestyles=ls)

    # Stoner-enhancement contours in the plain-PM region
    S_plot = np.where(region == 3, S, np.nan)
    cs = ax.contour(X, C, S_plot, levels=[1.5, 2.0], colors=["#9ec5f4",
                                                             "#5598e7"],
                    linewidths=1.0)
    ax.clabel(cs, fmt=lambda v: f"S={v:g}", fontsize=7.5, colors=INK2)

    # region labels (direct labels, text ink -- identity never color-alone)
    ax.text(0.24, 6.1, "no 2DEG\n(gallery band empty)", color=INK2,
            ha="center", fontsize=9)
    ax.annotate("polarized 2DEG: magnetic order,\nkills singlet SC;"
                " equal-spin\nf-wave triplet allowed",
                xy=(0.585, 6.55), xytext=(0.50, 5.82), fontsize=8.2,
                color="#7c1f1e", ha="center",
                arrowprops=dict(arrowstyle="->", color="#7c1f1e", lw=1.0))
    ax.annotate("paramagnon-boosted PM (S > S*)",
                xy=(0.70, 6.62), xytext=(0.62, 7.35), fontsize=8.2,
                color="#6d4c00", ha="center",
                arrowprops=dict(arrowstyle="->", color="#6d4c00", lw=1.0))
    ax.text(0.60, 8.55, "paramagnetic 2DEG:\nphonon singlet + paramagnon"
            " channel", color="#184f95", ha="center", fontsize=8.2)

    # anchors
    ax.axvline(0.35, color=INK2, lw=0.7, ls=":")
    ax.plot([0.35], [6.9], marker="x", ms=8, mew=2, color=INK)
    ax.annotate("monolayer hydrate\nc = 6.9 A, no SC", (0.35, 6.9),
                xytext=(0.14, 6.85), fontsize=8, color=INK)
    ax.plot([0.35], [9.9], marker="o", ms=7, mfc="#104281", mec="white")
    ax.annotate("bilayer hydrate\nc = 9.9 A, Tc = 4.5 K", (0.35, 9.9),
                xytext=(0.14, 10.05), fontsize=8, color=INK)
    ax.annotate("", xytext=(0.345, 9.8), xy=(0.31, 9.2),
                arrowprops=dict(arrowstyle="->", color="#7c1f1e", lw=1.1))
    ax.text(0.135, 8.75, "dehydration / valence shift:\nmagnetic phase"
            " adjacent to the\nSC dome (Sakurai 2015)", fontsize=7.6,
            color="#7c1f1e")

    ax.set_xlabel("Na content x", color=INK)
    ax.set_ylabel("CoO$_2$ plane spacing c (A)   [gallery opening]", color=INK)
    ax.set_title("(a) Gallery-2DEG spin phase diagram (model, ansatz-labelled)",
                 fontsize=10, color=INK, loc="left")

    ax2 = ax.twinx()  # secondary axis: same variable c mapped to Gamma0 (not a
    ax2.set_ylim(ax.get_ylim())                       # second data axis)
    ticks = np.arange(6, 11)
    ax2.set_yticks(ticks)
    ax2.set_yticklabels([f"{g:.1f}" for g in gamma0_of_c(ticks)])
    ax2.set_ylabel(r"$\Gamma_0(c)$ = full-coverage gallery half-width (eV)",
                   fontsize=8.5, color=INK2)
    ax2.tick_params(colors=INK2, labelsize=8)

    # --- (b) cuts at x = 0.35 ------------------------------------------------
    f35 = float(f_dilution(0.35, P["a_Na"], P["lam_Na"]))
    cc = np.linspace(5.5, 10.6, 500)
    w35 = gamma0_of_c(cc) * f35
    q35 = q_carrier(w35, gamma_co, two_delta)
    with np.errstate(divide="ignore"):
        F35 = P["I_s"] * P["nu"] / w35
        S35 = np.where(F35 < 1, 1.0 / (1.0 - F35), np.nan)

    # shade regions on the cuts
    bounds = boundaries_cut(f35)
    for a_ in (axq, axs):
        a_.axvspan(5.5, bounds["on"], color=FILL[0])
        a_.axvspan(bounds["on"], bounds["st"], color=FILL[1])
        a_.axvspan(bounds["st"], bounds["ss"], color=FILL[2])
        a_.axvspan(bounds["ss"], 10.6, color=FILL[3])
        a_.set_xlim(5.5, 10.6)

    axq.plot(cc, q35, color=BLUE, lw=2)
    axq.set_ylabel("q per Na", fontsize=9)
    axq.set_title("(b) cut at x = 0.35", fontsize=10, loc="left", color=INK)
    axq.tick_params(labelbottom=False)

    axs.plot(cc, S35, color=BLUE, lw=2)
    axs.set_ylim(0, 8)
    axs.set_ylabel("Stoner S", fontsize=9)
    axs.set_xlabel("c (A)", fontsize=9)
    for a_ in (axq, axs):
        for c_ref, mk in ((6.9, "x"), (9.9, "o")):
            yv = np.interp(c_ref, cc, q35 if a_ is axq else
                           np.nan_to_num(S35, nan=0.0))
            a_.plot([c_ref], [yv], marker=mk, color=INK, ms=6,
                    mfc="none" if mk == "o" else INK)
        a_.grid(color="#e4e3df", lw=0.6)
        a_.set_axisbelow(True)
        for s in ("top", "right"):
            a_.spines[s].set_visible(False)

    # legend for boundary strokes
    handles = [Line2D([], [], color=EDGE["turnon"], lw=1.4,
                      label="carrier turn-on  $\\Gamma$ = 2$\\delta$"),
               Line2D([], [], color=EDGE["stoner"], lw=1.4,
                      label="Stoner  $I_s N(0)$ = 1"),
               Line2D([], [], color=EDGE["sstar"], lw=1.4, ls="--",
                      label=f"S = S* = {P['S_star']:g}")]
    ax.legend(handles=handles, loc="upper right", fontsize=7.5, frameon=False)

    fig.suptitle("", fontsize=1)
    fig.savefig(outfile, bbox_inches="tight")
    print(f"wrote {outfile}")


def boundaries_cut(f35):
    two_delta = 2.0 * P["delta_Li"] * P["delta_Na_over_Li"]
    gamma_co = 6.0 * P["t_Co_xy"]
    g1, g2 = gamma0_of_c(np.array([8.0, 9.0]))
    B = g2 - g1
    A = g1 - 8.0 * B
    inv = lambda g0: (g0 - A) / B
    return {"on": inv((two_delta - gamma_co) / f35),
            "st": inv(P["I_s"] * P["nu"] / f35),
            "ss": inv(P["I_s"] * P["nu"] * P["S_star"]
                      / (P["S_star"] - 1.0) / f35)}


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("-o", "--output", default="spin_phase_diagram.pdf")
    args = ap.parse_args()

    print("Parameters (source-tagged):")
    for k, (v, src) in PARAMS.items():
        print(f"  {k:18s} = {v:8.4f}   {src}")
    i_min = (2.0 * P["delta_Li"] - 6.0 * P["t_Co_xy"]) / P["nu"]
    print(f"\n  [check] polarized strip exists iff I_s > (2 delta - Gamma_Co)/nu"
          f" = {i_min:.2f} eV; baseline I_s = {P['I_s']:.2f} eV")

    li_vs_na_report()
    sensitivity_scan()
    boundaries_at_x(0.35)
    make_figure(args.output)


if __name__ == "__main__":
    main()

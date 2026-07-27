#!/usr/bin/env python
"""
Figure 4 - DFT+U summary (PBE+U, U = 4 eV on Co 3d unless noted).

(a) Magnetic configuration energies per Co, each composition referred to its
    own ground state.  x = 1/2 orders as an intra-row antiferromagnet, x = 1/3
    as a ferrimagnet; in both cases the nonmagnetic solution is the expensive
    one, and at x = 1/3 it is expensive by an order of magnitude more than
    anything else in the panel, so the axis is broken rather than clipped.
(b) Li-motif energetics at x = 1/2, referred to the row order, at U = 4 eV
    (filled) and U = 0 (open).  Row and zigzag straddle zero - their splitting
    is smaller than kB x 240 K and even changes sign with U - while the
    clumped motif costs several times that thermal scale at either U.
(c) The two single-number fingerprints of each phase: the single-particle gap
    (left axis) and the moment on Co(4+) (right axis).  The nonmagnetic
    solutions are metallic and carry no moment, so both of their bars are
    zero by construction; that is stated in the panel rather than hidden.

Input: results/dftu_summary.json, the distilled output of the three DFT+U
campaigns.  This script computes nothing except differences of numbers that
are already in that file, and it renders no other data source: there is no
synthetic-data path.

Usage
-----
    python fig4_dftu.py
    python fig4_dftu.py --png-dir DIR      # extra 400 dpi inspection PNG
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, Patch

import _style as S

HERE = Path(__file__).resolve().parent
DEFAULT_JSON = HERE.parents[1] / "results" / "dftu_summary.json"
PNG = ("/private/tmp/claude-501/-Users-santoshkumar/"
       "57bc6bbc-f3db-4267-a42f-5531905aedd5/scratchpad")

# ----------------------------------------------------------------------
# Colour grammar for this figure.  Configuration identity is carried by hue
# AND by a direct text label on every level / bar / marker, never by hue
# alone.  Hues are the validated colourblind-safe set in _style.
# ----------------------------------------------------------------------
C_GS = S.C_BLUE          # the magnetic ground state (AFM at x=1/2, ferri at 1/3)
C_FM = S.C_RED           # ferromagnetic alignment
C_NM = S.C_AQUA          # nonmagnetic (spin-restricted) solution
C_MOT = S.C_BLUE         # Li-motif energies, panel (b)
INK_MOT = S.INK_BLUE
C_GAP = S.C_YELL         # gap, left axis of panel (c)
INK_GAP = S.INK_AMBER
C_MOM = S.C_BLUE         # moment, right axis of panel (c)
INK_MOM = S.INK_BLUE


# ----------------------------------------------------------------------
# Input
# ----------------------------------------------------------------------
def load(path: Path) -> dict:
    """Read results/dftu_summary.json and verify every field this figure
    consumes, so a schema drift fails here and not halfway through a draw."""
    if not path.exists():
        raise SystemExit(f"fig4_dftu: input not found: {path}")
    d = json.loads(path.read_text())
    need = {
        "row": ("energies_meV_per_Co_rel_NM", "gap_eV_ground_state_afm_intra",
                "moment_co4_muB", "nm_state"),
        "order_contrast": ("dE_meV_per_fu_rel_row_U4",
                           "dE_meV_per_fu_rel_row_U0", "kBT240K_meV"),
        "x13": ("energies_meV_per_Co_rel_ferri", "gap_eV_ground_state_ferri",
                "moments_muB", "nm_state"),
    }
    for block, keys in need.items():
        if block not in d:
            raise SystemExit(f"fig4_dftu: block '{block}' missing from {path}")
        missing = [k for k in keys if k not in d[block]]
        if missing:
            raise SystemExit(f"fig4_dftu: {block} is missing {missing}")
    for block in ("row", "x13"):
        if "metallic" not in d[block]["nm_state"]:
            raise SystemExit(
                f"fig4_dftu: {block} nm_state is '{d[block]['nm_state']}'; "
                f"panel (c) draws the nonmagnetic solution as gapless.")
    return d


def levels(data):
    """Configuration energies per Co, each composition referred to its own
    ground state.  The JSON stores x=1/2 relative to NM and x=1/3 relative to
    the ferrimagnet, so both are re-zeroed here on their minimum."""
    e12 = data["row"]["energies_meV_per_Co_rel_NM"]
    raw12 = {"gs": e12["AFM_intra_B1"], "FM": e12["FM_B1"], "NM": e12["NM"]}
    e13 = data["x13"]["energies_meV_per_Co_rel_ferri"]
    raw13 = {"gs": e13["ferri"], "FM": e13["FM"], "NM": e13["NM"]}
    out = []
    for raw, gs_name in ((raw12, "AFM (intra-row)"), (raw13, "ferri")):
        z = min(raw.values())
        out.append(({k: v - z for k, v in raw.items()}, gs_name))
    return out


def fmt(v):
    """Level values: an exact zero prints as 0, everything else signed."""
    if abs(v) < 5e-2:
        return "0"
    return f"{v:+.1f}".rstrip("0").rstrip(".") if abs(v) < 100 else f"{v:+.0f}"


# ----------------------------------------------------------------------
# (a) configuration energies, broken axis
# ----------------------------------------------------------------------
def panel_a(ax_hi, ax_lo, data):
    groups = levels(data)
    xt = [0.0, 1.45]
    half_l, half_r = 0.30, 0.02      # level mark spans [c-0.30, c+0.02]
    order = [("NM", C_NM, "NM"), ("FM", C_FM, "FM"), ("gs", C_GS, None)]

    # Split: everything except the x = 1/3 nonmagnetic level lives below.
    # Each level is drawn on exactly the one segment that contains it, so no
    # label can leak out of its own axes.
    span = {ax_lo: (-7.0, 78.0), ax_hi: (222.0, 238.0)}
    for ax, (ymin, ymax) in span.items():
        for c, (vals, gs_name) in zip(xt, groups):
            for key, col, tag in order:
                v = vals[key]
                if not (ymin <= v <= ymax):
                    continue
                name = gs_name if tag is None else tag
                ax.plot([c - half_l, c + half_r], [v, v], color=col, lw=2.6,
                        solid_capstyle="butt", zorder=4)
                ax.text(c + half_r + 0.05, v, f"{name}   {fmt(v)}",
                        fontsize=7.0, color=col, ha="left", va="center",
                        zorder=6, clip_on=False)
        ax.set_ylim(ymin, ymax)
    ax_lo.set_yticks([0, 20, 40, 60])
    ax_hi.set_yticks([230])

    for ax in (ax_hi, ax_lo):
        ax.set_xlim(-0.60, 2.20)
    ax_hi.spines["bottom"].set_visible(False)
    ax_lo.spines["top"].set_visible(False)
    ax_hi.tick_params(axis="x", bottom=False, labelbottom=False, top=False)
    ax_lo.tick_params(axis="x", top=False, length=0)
    ax_lo.set_xticks(xt)
    ax_lo.set_xticklabels([r"$x=1/2$", r"$x=1/3$"])
    ax_hi.set_xticks([])

    # Break marks on both broken spines: the gap in the axis is drawn, never
    # implied.
    brk = dict(marker=[(-1.0, -0.55), (1.0, 0.55)], markersize=5.5,
               linestyle="none", color=S.C_SEC, mec=S.C_SEC, mew=0.9,
               clip_on=False)
    ax_hi.plot([0, 1], [0, 0], transform=ax_hi.transAxes, **brk)
    ax_lo.plot([0, 1], [1, 1], transform=ax_lo.transAxes, **brk)

    ax_lo.set_ylabel(r"$E-E_{\mathrm{g.s.}}$   (meV per Co)")
    ax_lo.yaxis.set_label_coords(-0.145, 0.66)
    S.panel_label(ax_hi, "a")


# ----------------------------------------------------------------------
# (b) Li-motif energetics at x = 1/2
# ----------------------------------------------------------------------
def panel_b(ax, data):
    oc = data["order_contrast"]
    motifs = ["zigzag", "row", "clumped"]
    u4 = np.array([float(oc["dE_meV_per_fu_rel_row_U4"][m]) for m in motifs])
    u0 = np.array([float(oc["dE_meV_per_fu_rel_row_U0"][m]) for m in motifs])
    kbt = float(oc["kBT240K_meV"])
    xs = np.arange(len(motifs), dtype=float)
    dx = 0.145

    ax.set_xlim(-0.62, 2.78)
    ax.set_ylim(-30.0, 116.0)

    # The thermal scale of the transition, as a band about zero: anything
    # inside it is degenerate as far as the 240 K transformation is concerned.
    ax.axhspan(-kbt, kbt, color=S.FILL_NONE, zorder=0)
    ax.axhline(0.0, color=S.C_MUT, lw=0.7, zorder=1)
    ax.axhline(kbt, color=S.C_SEC, lw=0.9, ls=(0, (3.2, 2.0)), zorder=2)
    ax.text(-0.55, kbt + 3.0,
            rf"$k_{{\mathrm{{B}}}}\times 240\,$K$\,=\,{kbt:.1f}$ meV",
            fontsize=6.8, color=S.C_SEC, ha="left", va="bottom", zorder=6)

    for x, a, b in zip(xs, u4, u0):
        ax.plot([x - dx, x + dx], [a, b], color=S.C_MUT, lw=0.7, zorder=3)
    ax.plot(xs - dx, u4, ls="none", marker="o", ms=5.2, mfc=C_MOT,
            mec=INK_MOT, mew=0.9, zorder=5)
    ax.plot(xs + dx, u0, ls="none", marker="o", ms=5.2, mfc="white",
            mec=C_MOT, mew=1.2, zorder=5)
    # A value of exactly zero is lifted clear of the zero rule, which would
    # otherwise strike through the glyph.
    for x, v, ha in [(x - 0.08, v, "right") for x, v in zip(xs - dx, u4)] + \
                    [(x + 0.08, v, "left") for x, v in zip(xs + dx, u0)]:
        zero = abs(v) < 5e-2
        ax.text(x, v + (3.2 if zero else 0.0), "0" if zero else f"{v:+.1f}",
                fontsize=7.0, color=INK_MOT, ha=ha,
                va="bottom" if zero else "center", zorder=6)

    # The clumping penalty, measured against the same thermal scale.
    ax.add_patch(FancyArrowPatch((2.0, 0.0), (2.0, float(u4[2])),
                                 arrowstyle="<|-|>", mutation_scale=6,
                                 lw=0.9, color=S.C_INK, shrinkA=0, shrinkB=0,
                                 zorder=5))
    ax.text(2.06, 0.5 * float(u4[2]),
            rf"${u4[2] / kbt:.1f}\,k_{{\mathrm{{B}}}}T$",
            fontsize=7.0, color=S.C_INK, ha="left", va="center", zorder=6)

    ax.text(0.62, 66.0,
            "row and zigzag straddle zero:\n"
            "near-degenerate stripe manifold",
            fontsize=6.8, color=S.C_SEC, ha="center", va="bottom",
            linespacing=1.35, zorder=6)

    ax.set_xticks(xs)
    ax.set_xticklabels(motifs)
    ax.tick_params(axis="x", length=0, top=False)
    ax.set_xlabel(r"Li ordering motif at $x=1/2$")
    ax.set_ylabel(r"$E-E_{\mathrm{row}}$   (meV per f.u.)")

    handles = [Line2D([], [], ls="none", marker="o", ms=5.2, mfc=C_MOT,
                      mec=INK_MOT, mew=0.9, label=r"$U=4$ eV"),
               Line2D([], [], ls="none", marker="o", ms=5.2, mfc="white",
                      mec=C_MOT, mew=1.2, label=r"$U=0$")]
    leg = ax.legend(handles=handles, loc="upper left",
                    bbox_to_anchor=(0.135, 1.015), ncol=2, fontsize=7.0,
                    handletextpad=0.35, columnspacing=1.1,
                    borderaxespad=0.0)
    leg.set_zorder(9)
    S.panel_label(ax, "b")


# ----------------------------------------------------------------------
# (c) gaps and moments
# ----------------------------------------------------------------------
def panel_c(ax, data):
    axr = ax.twinx()
    gaps = np.array([float(data["row"]["gap_eV_ground_state_afm_intra"]),
                     float(data["x13"]["gap_eV_ground_state_ferri"]),
                     0.0])
    moms = np.array([float(data["row"]["moment_co4_muB"]),
                     float(data["x13"]["moments_muB"]["Co4_a"]),
                     0.0])
    xs = np.arange(3, dtype=float)
    w = 0.30

    ax.set_xlim(-0.62, 2.62)
    ax.set_ylim(0.0, 1.85)          # eV,  left
    axr.set_ylim(0.0, 1.65)         # muB, right
    ax.set_yticks([0.0, 0.5, 1.0, 1.5])
    axr.set_yticks([0.0, 0.5, 1.0, 1.5])

    ax.bar(xs[:2] - 0.17, gaps[:2], width=w, color=S.FILL_AMBER,
           edgecolor=INK_GAP, lw=0.9, zorder=3)
    axr.bar(xs[:2] + 0.17, moms[:2], width=w, color=S.FILL_BLUE,
            edgecolor=C_MOM, lw=0.9, zorder=3)
    # Each value is printed above its own bar in that bar's own ink, and the
    # two scales differ (1.5 eV and 1.5 muB sit at different heights), so a
    # bar can never be read off the wrong axis.
    for x, g in zip(xs[:2] - 0.17, gaps[:2]):
        ax.text(x, g + 0.05, f"{g:.3f}", fontsize=7.0, color=INK_GAP,
                ha="center", va="bottom", zorder=6)
    for x, m in zip(xs[:2] + 0.17, moms[:2]):
        axr.text(x, m + 0.045, f"{m:.2f}", fontsize=7.0, color=INK_MOM,
                 ha="center", va="bottom", zorder=6)

    # Nonmagnetic column: both quantities are zero, and that is written out.
    ax.plot([2 - 0.17 - w / 2, 2 + 0.17 + w / 2], [0, 0], color=S.C_SEC,
            lw=1.6, solid_capstyle="butt", zorder=4)
    ax.text(2.0, 0.11, "metallic:\ngap $=0$,  $m=0$", fontsize=6.8,
            color=S.C_SEC, ha="center", va="bottom", linespacing=1.35,
            zorder=6)

    ax.set_xticks(xs)
    ax.set_xticklabels([r"$x=1/2$" "\n" "AFM", r"$x=1/3$" "\n" "ferri",
                        "nonmagnetic" "\n" r"(either $x$)"])
    ax.tick_params(axis="x", length=0, top=False)
    ax.set_ylabel(r"band gap   (eV)", color=INK_GAP)
    axr.set_ylabel(r"Co$^{4+}$ moment   ($\mu_{\mathrm{B}}$)", color=INK_MOM)
    ax.tick_params(axis="y", colors=INK_GAP)
    axr.tick_params(axis="y", colors=INK_MOM, direction="in", width=0.8,
                    length=3.0)
    ax.spines["left"].set_color(INK_GAP)
    axr.spines["left"].set_color(INK_GAP)
    axr.spines["right"].set_color(INK_MOM)
    axr.spines["top"].set_visible(False)

    handles = [Patch(facecolor=S.FILL_AMBER, edgecolor=INK_GAP, lw=0.9,
                     label="gap (eV)"),
               Patch(facecolor=S.FILL_BLUE, edgecolor=C_MOM, lw=0.9,
                     label=r"$m$ ($\mu_{\mathrm{B}}$)")]
    leg = ax.legend(handles=handles, loc="upper right",
                    bbox_to_anchor=(1.005, 1.02), fontsize=7.0,
                    handlelength=1.1, handletextpad=0.4, borderaxespad=0.0)
    leg.set_zorder(9)
    S.panel_label(ax, "c")


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--json", type=Path, default=DEFAULT_JSON,
                    help=f"summary JSON (default: {DEFAULT_JSON})")
    ap.add_argument("--png-dir", default=PNG,
                    help="directory for the 400 dpi inspection PNG")
    args = ap.parse_args()

    data = load(args.json)
    S.use_house_style()

    fig = plt.figure(figsize=(3.375, 4.80))
    gs = GridSpec(3, 1, figure=fig, height_ratios=[1.30, 1.30, 1.26],
                  left=0.163, right=0.845, top=0.978, bottom=0.088,
                  hspace=0.46)
    inner = GridSpecFromSubplotSpec(2, 1, subplot_spec=gs[0],
                                    height_ratios=[1.0, 4.6], hspace=0.10)
    ax_hi = fig.add_subplot(inner[0])
    ax_lo = fig.add_subplot(inner[1])
    panel_a(ax_hi, ax_lo, data)
    panel_b(fig.add_subplot(gs[1]), data)
    panel_c(fig.add_subplot(gs[2]), data)

    out = S.save(fig, "fig4_dftu", png_dir=args.png_dir)
    print(f"wrote {out}  (from {args.json})")


if __name__ == "__main__":
    main()

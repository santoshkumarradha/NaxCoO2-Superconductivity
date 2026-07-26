# Publication figures — Na$_x$CoO$_2\cdot y$H$_2$O superconductivity mechanism

Figures for the mechanism paper (PRL/PRB register). Everything here is generated
from the real DFT / theory data under `runpod/`, `theory/`, and `spin_analysis/`
— nothing is fabricated, and none of those source directories are modified.

## Folder layout (sources vs rendered artefacts, kept separate)

- `paper/figure_scripts/` — **sources**: every `figN.py`, the shared `_style.py`,
  the two TikZ `.tex` figures, and this README.
- `paper/figures/` — **rendered artefacts only**: the `figN.pdf`/`figN.png` that
  `main.tex` consumes (`\graphicspath{{figures/}}`). `_style.save()` always writes
  here (`HERE.parent/"figures"`), so running a script from `figure_scripts/` emits
  its PDF+PNG into `figures/` and never clutters the source folder.

## Build

```sh
# from this directory (paper/figure_scripts/)
PY=../../theory/.venv/bin/python
for n in 1 2 5 6 7 8 9 10; do $PY fig$n.py; done               # -> ../figures/figN.{pdf,png}
for s in dome thresholds zrncl zscan; do $PY figS_$s.py; done  # -> ../figures/figS_*.{pdf,png}
# LaTeX pieces (need pdflatex on PATH, e.g. /Library/TeX/texbin); copy the PDF to
# ../figures/ and rasterize a PNG preview there, keeping this folder sources-only:
for t in fig0_schematic fig3_model; do
  pdflatex -interaction=nonstopmode $t.tex
  cp $t.pdf ../figures/$t.pdf
  pdftoppm -png -r 150 $t.pdf ../figures/$t && mv ../figures/$t-1.png ../figures/$t.png
  rm -f $t.aux $t.log $t.pdf
done
```

Superseded generators (old `fig3.py`, `fig4.py`, `fig3_scoreboard.tex`) live in
`../archive/superseded_figure_scripts/`, with their rendered outputs in
`../archive/superseded_figures/`. `fig0_hero.{pdf,png}` in `../figures/` is the
AI-illustrated hero (Fig. 1); `fig0_schematic` is its deterministic TikZ fallback.

Each figure is emitted as a **vector PDF** (for the manuscript) and a **300 dpi PNG**
(for previews). `_style.py` holds the shared house style and the read-only data
loaders imported by every `figN.py`. Figs 7–9 additionally import the read-only
spectral parsers in `theory/spectral/spectral.py`.

## House style (shared, `_style.py`)

- **Serif / STIX** text, mathtext for symbols, ticks-in, thin spines, no gridlines;
  bold `(a)/(b)` panel letters. Single-column (3.4 in) for figs 1/4/5, full-width
  for figs 2/3.
- **Fixed, colourblind-safe series → colour mapping, identical across all figures**
  (validated with the dataviz palette validator: worst adjacent CVD $\Delta E = 21.6$,
  well above the 12 target; the aqua/yellow contrast WARN is relieved by a distinct
  marker per series + legends/direct labels, so identity is never colour-alone):

  | series | colour | marker |
  |---|---|---|
  | Na $1{\times}1$ ($x{=}1$) | blue `#2a78d6` | ● |
  | Li $1{\times}1$ ($x{=}1$) | aqua `#1baf7a` | ■ |
  | Na $\sqrt3{\times}\sqrt3$ ($x{=}1/3$) | yellow `#eda100` | ▲ |
  | **experimental SC anchor** (9.9 Å) | red `#e34948` ★ | reserved |

  Experimental numbers are always visually distinct from theory (open red star /
  grey–red anchor lines + on-figure citation).

---

## Figure 0 — Mechanism schematic  (`fig0_schematic.tex` → `.pdf`/`.png`)

**Provenance:** illustrative TikZ schematic; the well shapes / level ordering follow
`theory/results/well_fits_v4.csv` and `quantum_wells_v4.csv`. Colours here are
didactic and deliberately different from the data-plot series colours.

**Caption draft.** *Mechanism schematic.* As the CoO$_2$–CoO$_2$ gallery spacing $c$
grows with hydration, the alkali potential along the interlayer coordinate $\delta$
evolves from a stiff single well (anhydrous, 5.5 Å; no gallery 2DEG) through a soft
double well ($c\approx6.4$ Å) where the mode goes soft ($\alpha<0$,
$\hbar\omega_{\rm eff}\to0$) and a two-dimensional electron gas turns on in the
CoO$_2$ planes — the superconducting window — to a deep, polaronic double well at the
bilayer spacing (9.9 Å). The bilayer superconductor therefore sits in the
bare-geometry polaronic regime; reconciling it needs water to soften the Na well back
into the window.

## Figure 1 — The transition  (`fig1.py`)

**Provenance:** raw $E(\delta)$ totals re-parsed from `runpod/results_v4/jobs/*/pw.out`
(`!  total energy`, Ry → meV, referenced to $\delta{=}0$); fit coefficients
($a_q,b_q,g_q$), depths, and the `min_bracketed` flag from
`theory/results/well_fits_v4.csv`.

**Caption draft.** *The single→double-well transition.* DFT total energies $E(\delta)$
versus the alkali off-centre displacement for (a) Na $1{\times}1$ and (b) Na
$\sqrt3{\times}\sqrt3$ at the four CoO$_2$ spacings; points are DFT, curves the
quartic/sextic Landau fits (dashed = constrained sextic refit). The centred single
well ($\alpha>0$) at 5.5 Å gives way to a symmetric double well ($\alpha<0$) by
6.9 Å. For $c\ge8.4$ Å the scan is **unterminated** (open markers, down-arrows): the
alkali is still relaxing outward at the largest sampled $\delta$, so those well depths
are **lower bounds** (adsorption regime), shown honestly rather than extrapolated.

## Figure 2 — The money plot  (`fig2.py`)  ★ centrepiece

**Provenance:** (a) `alpha` from `theory/results/well_fits_v4.csv` (PCHIP through the
four DFT points); (b) `N0`, `n2deg_cm2` from `runpod/results_v4/results.json`;
(c) `Tc_K`, `Tc_K_g1`, `Tc_K_g2`, `lambda`, `status` from
`theory/results/coupling_tc_vs_c.csv`.

**Caption draft.** (a) Landau $\alpha(c)$: the interlayer mode softens ($\alpha$ crosses
zero) at $c^\ast_{\rm Na}\approx6.40$ Å ($6.57$ Å for the $\sqrt3$ cell), while Li is
already unstable at $\le5.5$ Å — no soft-mode window is accessible to Li. (b) The
gallery 2DEG switches on across the same window: DFT $N(0)$ rises through the
metallicity gate for the $x{=}1$ cells. (c) Allen–Dynes $T_c(c)$ forms a narrow dome
peaking at $\approx14$ K (γ-sensitivity band $11$–$23$ K over the electron–phonon range
$z_0=0.5$–$1.0$ Å); for $c\gtrsim6.44$ Å the coupling exceeds $\lambda=2$ and the
Migdal treatment is invalid (polaronic self-trapping — SC killed, hatched). The
experimental anchors (5.5/6.9 Å not SC, grey; 9.9 Å SC 4.5 K, red star) place the
bilayer superconductor **inside** the bare-geometry polaronic region; reconciling it
requires water-induced softening of the Na well (prediction, untested).

*Note on carrier density:* the Bader areal density $n_{\rm 2D}$ (`n2deg_cm2`,
$10^{14}$ cm$^{-2}$) is **1.16/0.96/0.61/0.35** (Na $1{\times}1$, $c=5.5/6.9/8.4/9.9$ Å)
— it *falls* with dilution and so was left off panel (b) to keep the "turn-on"
message legible; $N(0)$ is the direct 2DEG measure and $n_{\rm 2D}$ merely corroborates.

## Figure 3 — Theory vs experiment scoreboard  (`fig3_scoreboard.tex` — native LaTeX table)

**Delivered as a real LaTeX table**, not an image. The standalone `fig3_scoreboard.tex`
(booktabs + tabularx, `standalone` class) is kept for preview and renders
`fig3_scoreboard.pdf/.png`. **In the manuscript it is now Table~I** (`\label{tab:score}`),
built directly in `main.tex` as a `table*` float with a plain **booktabs** `tabular` and
fixed `p{}` column widths — *not* tabularx. This is deliberate: `revtex4-2` in the current
TeX build conflicts with the `array`/`tabularx` packages ("Extra \or" in the column
preamble), so those are not loaded in the manuscript; only `booktabs` + `pifont` are.
Rules are `\toprule`/`\midrule`/`\bottomrule` (no vertical rules); verdict glyphs carry a
colour *and* a word (green check "explained" / amber triangle "prediction"). (`fig3.py`
renders an equivalent raster version if a bitmap is ever needed, but the table is canonical.)

**Provenance:** verdicts summarise the mechanism results; $\lambda\approx20$ at 6.9 Å is
the canonical `theory/MODEL.md` value; anchors cite Foo *et al.* (2003) and
Takada *et al.* (2003).

**Caption draft.** *Scoreboard.* The soft-mode / charge-bifurcation mechanism against
the six established facts: five are explained; the bilayer-hydrate superconductivity at
9.9 Å is a standing **prediction** requiring water to soften the Na well back into the
SC window.

## Figure 4 — Spin turn-on  (`fig4.py`)

**Provenance:** `spin_analysis/magnetization_v2.csv` (`abs_mag_muB`; per $(c)$ the
marker is the mean over the $\delta$ scan and the band is min–max).

**Caption draft.** *Spin turn-on.* LSDA cell moment $\langle|m|\rangle$ vs $c$ for Li and
Na $1{\times}1$: the local moment switches on precisely at the $\alpha<0$ well transition,
so the magnetic phase borders the SC dome, and the alkali gallery moment is antiparallel
to Co. The Na $\sqrt3$ cell is over-polarised in LSDA ($|m|\approx$ const $\approx2.3\,
\mu_B$) and is shown only as a shaded caveat band — an LSDA systematic, not a physical
turn-on.

## Figure 5 — Mode softening  (`fig5.py`)

**Provenance:** `theory/results/quantum_wells_v4.csv` (`hw_eff_meV`, physical-mass rows;
`E0/E1/E2_meV`); inset well shape from `well_fits_v4.csv` (Na $1{\times}1$, 6.9 Å).

**Caption draft.** *Mode softening.* The effective zero-point energy
$\hbar\omega_{\rm eff}=\hbar^2/2M\langle\delta^2\rangle_0$ (log scale) collapses from
$\approx30$ meV (stiff single well) to $\approx10\,\mu$eV (soft double well) as the gallery
opens, entering a quantum-paraelectric window ($\hbar\omega_{\rm eff}<1$ meV). Inset: the
quantized levels in the $c=6.9$ Å Na double well, with $E_0,E_1$ a near-degenerate
tunnelling doublet below the $E_2$ excited state.

---

## Honesty / systematics summary (carried on-figure and in captions)

- **$c\ge8.4$ Å well depths are lower bounds** — the $E(\delta)$ scans are unterminated
  (adsorption regime); flagged with open markers + down-arrows (Fig. 1).
- **$T_c$ band** spans the electron–phonon ansatz $z_0=0.5$–$1.0$ Å ($\gamma=1/z_0$);
  the central curve uses $z_0=0.7$ Å (Fig. 2c).
- **$\lambda>2$ regions are not reliable BCS $T_c$** — they are polaronic (SC killed),
  shown hatched rather than as a dome (Fig. 2c).
- **Na $\sqrt3$ moment** is an LSDA over-polarisation artefact, excluded from the
  turn-on and marked as a caveat band (Fig. 4).
- **Bilayer (9.9 Å) SC is a prediction, not a postdiction** — the bare geometry is
  polaronic there; SC needs water softening of the Na well (Figs. 2c, 3, 0).

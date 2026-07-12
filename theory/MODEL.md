# Effective model: gallery-controlled Na double well + interlayer 2DEG superconductivity in Na_xCoO2·yH2O

**TL;DR.** The canonical DFT dataset of this repo (`runpod/results_v4`: Quantum ESPRESSO,
PBE+PAW, spin-polarized, Co-top alkali site, c = 5.5/6.9/8.4/9.9 Å) shows that opening the
alkali gallery does two things between 5.5 and 6.9 Å: (i) the centered Na position becomes
unstable (Landau α(c) crosses zero — interpolated c\*_Na ≈ 6.4 Å, rigorously bracketed by
the DFT points at 5.5 and 6.9 Å), and (ii) the interlayer 2DEG turns on — seen *twice
independently*, in the band N(0) (0.03 → 1.06 states/eV/cell/spin) and in the Bader charge
flowing back onto the alkali (Δq: 0 → 0.14 e). Quantizing the DFT wells exactly (Na mass)
and wiring the displacement into the SSH4 tight-binding model of Radha & Lambrecht with the
parity-correct quadratic vertex gives a narrow superconducting window pinned just above c\*
(c ≈ 6.35–6.45 Å, peak Tc ≈ 11–23 K over the vertex-ansatz range), bounded below by "no
2DEG" and above by "Na frozen off-center" (λ > 2, polaronic — SC killed). Both anhydrous
(5.5 Å) and monolayer-hydrate (6.9 Å) anchors correctly come out non-superconducting. The
bilayer-hydrate anchor (9.9 Å, Tc = 4.5 K) does **not** land in the window under this
water-free model: there Na is deeply bound toward one CoO2 layer. The model's reconciling
prediction — that the H2O bilayer softens/re-symmetrizes the Na well back toward
criticality — is stated on the money plot as an **untested prediction**.

Everything here is produced end-to-end by `theory/effective_model.py` from
`runpod/results_v4`.

## Reproduce

```bash
uv run theory/effective_model.py            # PEP 723 inline deps (numpy/scipy/matplotlib/pandas)
# or: theory/.venv/bin/python theory/effective_model.py
theory/.venv/bin/python theory/effective_model.py --legacy   # old proof-of-concept numbers
```

Outputs: `theory/results/*.csv` (all tables), `theory/figures/fig1–fig5` (diagnostics),
`theory/money_plot.{pdf,png}`. The `--legacy` flag re-runs the original
`phase_transition/` analysis (Questaal lmf relaxed scans) and prints the comparison
summary quoted in §8; legacy tables go to `theory/results/legacy/`.

---

## 1. Input data and provenance (canonical: runpod/results_v4)

| Quantity | Value / file | Provenance |
|---|---|---|
| E(δ; c) scans: Li 1×1 (x=1), Na 1×1 (x=1), Na √3×√3 (x=1/3); c = 5.5/6.9/8.4/9.9 Å; δ = 0…1.0 Å (1×1, 6 pts) / 0…0.75 Å (√3, 4 pts) | `runpod/results_v4/jobs/*/pw.out` (energies re-parsed in Ry) | QE `pw.x`, PBE, PAW (kjpaw psl), ecutwfc/rho = 60/480 Ry, nspin=2, MV smearing 0.02 Ry, 8×8×3 k, Co-top alkali site, fixed z_O = 0.96 Å |
| Quartic fits α(c), β(c) | `runpod/results_v4/results.json` | same dataset (canonical Landau parameters) |
| N(0) per spin per cell | `results.json` (`N0`, from dense-k `dos.x` runs `jobs/nscf_*/pw.dos`) | same dataset |
| Bader charge on the alkali q(c) | `results.json` (`q_bader`, from `ACF.dat`) | same dataset (Z_val: Na 9, Li 3 — semicore PPs) |
| a_lat: Na 2.888 Å, Li 2.82 Å; cell areas 7.22 / 6.89 / 21.67 Å² (Na 1×1 / Li 1×1 / Na √3) | `manifest.json` | same dataset |
| τ1 = t_A^z = 0.5 eV, τ3 = t_CoO2^z = 0.5 eV, τ2 = τ4 = 2.0 eV; t_A^xy = −0.6, t_CoO2^xy = +0.09 eV → Γ = 4.14 eV; q_alk = Q(Γ−2δ)/(2Γ) | SSH4 fit + band-overlap formula | Radha & Lambrecht, SciPost Phys. 10, 057 (2021) |
| γ_ep = 1/z0, z0 = 0.7 Å (range 0.5–1.0 Å → γ_ep ∈ [1, 2] Å⁻¹) | **ansatz** | hopping-decay scale, full range shown |
| μ\* = 0.10; λ_loc = 2 (polaron criterion, range 1–2); Q_tot = 1 e/cell | **ansatz** / paper limit | task prescription / bipolaron literature / SciPost paper |
| Hard-core cap: well minimum at δ ≤ c/2 − 2.0 Å | **ansatz** (alkali–O-plane contact) | used only for the sextic refits, §2 |
| Anchors: 5.5 Å anhydrous (no SC), 6.9 Å monolayer hydrate (no SC), 9.9 Å bilayer hydrate (Tc = 4.5 K) | experimental | Takada et al., Nature 422, 53 (2003) |

### Verification (runpod/results_v4/checks/, all runs converged, JOB DONE)

* **k-mesh**: 12×12×4 vs 8×8×3 at the c = 6.9 Å, δ = 0.5 geometry: total energy shifts by
  **−1.6 meV** — negligible on the ~100 meV wells.
* **Spin treatment**: nspin=1 vs nspin=2 well drops E(δ)−E(0) are identical to within
  6–7 %: at c = 6.9 Å, δ = 0.5: **−89.1 vs −88.9 meV** (0.2 %); at c = 9.9 Å, δ = 0.75:
  −121.5 vs −113.8 meV (6.3 %). Magnetism is irrelevant to the well physics at this level.
* **z_O sensitivity**: moving the O planes by ±0.06 Å (z_O = 0.90 / 1.02 Å) shifts the
  single-point energy at (c = 6.9, δ = 0.5) by **+117 / −21 meV** — of the same order as
  the well depths themselves. Consequence: **well DEPTHS carry a ±50 % systematic** from
  the frozen (unrelaxed) z_O; the α *sign flips* (single ↔ double well) are robust
  (see §6 systematics).
* **ecut / smearing**: ecut = 80 Ry and degauss = 0.01 Ry single points exist at the same
  geometry (absolute shifts −260 / −36 meV, as expected for PAW total energies and
  smearing free energies); the *relative* (well-drop) checks at higher ecut/lower smearing
  are **pending** — they need a second δ point each.

## 2. Landau parameters and well potentials

Canonical α(c), β(c) are the quartic fits stored in `results.json`. For quantization the
script keeps the quartic wherever it is valid (β > 0 **and** the E(δ) minimum is bracketed
by the sampled δ range) and otherwise refits the raw pw.out energies with
E0 + αδ² + βδ⁴ + γ₆δ⁶, γ₆ > 0 constrained (γ₆ scanned, (E0, α, β) by linear least squares),
subject to a confining potential and the hard-core cap δ_min ≤ c/2 − 2.0 Å. Refit α agree
with the canonical ones to ≤ 2 %. Tables: `results/well_fits_v4.csv`; figure `fig1`.

| series | c (Å) | α (eV/Å²) | well d₀ (Å) | depth (meV) | fit |
|---|---|---|---|---|---|
| Na 1×1 | 5.5 | **+2.455** | 0 (single well) | — | quartic |
| Na 1×1 | 6.9 | **−0.506** | 0.69 | 119 | quartic (min bracketed) |
| Na 1×1 | 8.4 | −0.421 | 1.75 | 758 (≥ 379 sampled) | sextic, extrapolated |
| Na 1×1 | 9.9 | −0.193 | 2.87 | 1531 (≥ 210 sampled) | sextic, extrapolated |
| Na √3 (x=1/3) | 5.5 | +3.049 | 0 | — | quartic |
| Na √3 | 6.9 | −0.346 | 0.56 | 54 | quartic (min bracketed) |
| Li 1×1 | 5.5 | −0.155 | 0.23 | 4 | quartic (min bracketed) |
| Li 1×1 | 6.9 | −0.890 | 1.29 | 831 (≥ 679 sampled) | sextic, extrapolated |

* **c\*_Na (1×1) = 6.40 Å** (PCHIP interpolation of α(c); linear interpolation gives
  6.66 Å). The rigorous statement is the DFT **bracket: 5.5 Å < c\*_Na < 6.9 Å** — with
  only four spacings the interpolated value is scheme-dependent at the ±0.3 Å level.
  Na √3×√3 (x=1/3): c\* = 6.57 Å (PCHIP; linear 6.76) — the physical stoichiometry shifts
  c\* by only ≈ +0.2 Å, the mechanism is not an x = 1 artifact.
* **c\*_Li ≤ 5.5 Å — an upper bound only.** α is already (weakly) negative at the smallest
  spacing in this dataset, so the Li crossing is not bracketed from below. Note the Li well
  at 5.5 Å is only 4 meV deep with ħω_eff = 8.8 meV ≫ depth: Li is a *quantum
  paraelectric* there — dynamically centered despite α < 0.
* **Beyond ~7 Å the "double well" is really surface adsorption.** At c ≥ 8.4 Å every E(δ)
  scan is still decreasing at the last sampled point (δ = 1.0 / 0.75 Å): the alkali is
  sliding toward a bound site on one CoO2 layer, and the sampled range does not bracket
  the minimum. Quoted depths there are extrapolations of the constrained sextic (the
  sampled drop is the firm lower bound, flagged in the table and in fig2 with arrows).
  Consistently, the fitted d₀ ≈ c/2 − 2.7…2.9 Å tracks a constant alkali–layer distance.
* α(c) is **non-monotonic**: |α| shrinks again from 6.9 → 9.9 Å (−0.51 → −0.19 eV/Å²),
  i.e. the *curvature at the midpoint* re-flattens as the two CoO2 layers decouple — but
  the well *depth* keeps growing because the minimum runs away toward the layer. Softening
  of α alone does not re-symmetrize the well.

## 3. Quantized wells (exact 1D Schrödinger, physical masses)

Parity-resolved finite differences (staggered half-line grid, even/odd sectors separate the
tunneling doublet exactly; N = 2400, box follows the well; levels converged ≪ 1 %).
For deep wells the mode scale is **always** ħω_eff = ħ²/(2M⟨δ²⟩₀) — never E1−E0, which is
an exponentially small tunnel splitting there (< 10⁻¹² meV at c ≥ 6.9 for Na).
`results/quantum_wells_v4.csv`, figure `fig3`.

Na 1×1 wells, Na mass:

| c (Å) | ħω_eff (meV) | E1−E0 (meV) | E2−E0 (meV) | √⟨δ²⟩₀ (Å) |
|---|---|---|---|---|
| 5.5 | 30.2 | 30.2 (harmonic) | 60.6 | 0.055 |
| 6.9 | 0.196 | < 10⁻¹² | 18.6 | 0.68 |
| 8.4 | 0.030 | < 10⁻¹² | 21.7 | 1.74 |
| 9.9 | 0.011 | < 10⁻¹² | 23.0 | 2.87 |

ω_eff collapses by three orders of magnitude across c\* (30 → 0.2 → 0.01 meV) while the
intra-well phonon E2−E0 re-hardens to ≈ 20 meV: past c\* the Na spectral weight at low
energy is static displacement, not dynamics — the polaronic signature. The Li-mass curves
(fig3) are the √(23/6.9) ≈ 1.8× stiffer isotope cross-check.

## 4. Carriers: DFT N(0) and Bader gallery charge tell the same story

Two independent probes of the interlayer 2DEG, both straight from results_v4
(`fig4`, money-plot middle panel):

* **Band N(0)** (dos.x, per spin per cell): Na 1×1: 0.03 / 1.06 / 1.32 / 1.36 at
  5.5/6.9/8.4/9.9 Å; Li 1×1: 0.03 / 0.96 / 1.34 / 1.37. Insulating at 5.5 Å, metallic
  from 6.9 Å on, saturating.
* **Bader charge returned to the alkali** Δq(c) = q_bader(c) − q_bader(5.5): Na 1×1:
  0 / 0.14 / 0.40 / 0.58 e; Li 1×1: 0 / 0.07 / 0.47 / 0.67 e (the *donated* charge
  Z_val − q_bader falls with c — electrons come back into the gallery as it opens).
  As a 2D density n₂D = Δq/A_cell: 2.0 / 5.5 / 8.0 ×10¹³ e/cm² (Na 1×1 at
  6.9/8.4/9.9 Å) — sane 2DEG scales.
* Agreement: both probes turn on between 5.5 and 6.9 Å and grow monotonically
  (correlation N(0) vs Δq = 0.85 Na, 0.83 Li over the four spacings). N(0) saturates
  earlier than Δq — the band edge crosses E_F once, then the *occupation* keeps growing;
  no contradiction. Caveat: in the √3 cell (x=1/3) N(0) is 3.4–4.0 and flat — the Co-d
  metallic background of the mixed-valent x=1/3 plane masks the gallery turn-on there
  (its Bader Δq is also weak/non-monotone, ±0.2 e); the turn-on statement rests on the
  1×1 series.
* The 2DEG gate in the model is N(0) ≥ 0.1 states/eV/cell/spin (anything between 0.04
  and 0.9 gives the same picture).

The Bader Δq also calibrates the TB on-site splitting via the paper's band-overlap
formula, δ(c) = (Γ/2)(1 − 2Δq/Q_tot), floored at 0.05 eV: δ = 2.07 / 1.48 / 0.42 / 0.05 eV
at 5.5/6.9/8.4/9.9 Å.

## 5. Vertex, λ(c) and Tc(c)

**Parity.** The SSH-type coupling (τ2 → τ2(1+γδ), τ4 → τ4(1−γδ)) is odd under z → −z, so
the linear deformation potential of the alkali band vanishes identically at δ = 0
(verified numerically in the script at several k_z). The `I` stored in results.json is a
linear *secant* fitted over the one-sided δ scan — a finite-displacement average, not the
(zero) derivative at the symmetric point; the λ/Tc columns of results.json built on it are
superseded here. The leading vertex is quadratic:

E_band(δ) = E_edge + ½χδ²,  χ = 2Σ_m |⟨alk|∂H/∂δ|m⟩|²/(E_alk−E_m)   (exact 2nd-order PT)

χ = −6.1 eV/Å² at (δ = 2.07 eV, γ_ep = 1.43 Å⁻¹), growing to −15.3 eV/Å² as δ(c) collapses
at large c; χ ∝ γ_ep². Electrons couple to the even-parity excitation |0⟩ → |2⟩
(⟨0|δ²|1⟩ = 0 by parity):

g = ½|χ|⟨0|δ²|2⟩,  ħΩ = E2−E0,  λ = 2 N(0) g²/ħΩ

with the **DFT N(0)** (PCHIP in c) and exact well matrix elements. Allen–Dynes Tc
(single Einstein mode, f1 factor, μ\* = 0.10), gated by N(0) ≥ 0.1; where λ > λ_loc = 2
the Migdal treatment is invalid and the physical expectation with ħΩ ~ 20 meV, λ ≫ 1 is
polaronic self-trapping: those regions are reported as **"SC killed (polaronic)"** — raw
Allen–Dynes numbers there are meaningless and are not reported.

| c (Å) | status (γ_ep = 1.43 Å⁻¹) | numbers |
|---|---|---|
| 5.5 (anhydrous) | **no 2DEG** — no SC ✓ | N(0) = 0.03, λ = 4×10⁻⁴ |
| 6.35–6.45 | **SC window** | peak Tc = 14 K at 6.42 Å (λ = 1.3, ħΩ = 11 meV) |
| 6.9 (monolayer hydrate) | **polaronic, SC killed ✓** | λ ≈ 20 ≫ 2, Na frozen (d₀ = 0.69 Å) |
| 9.9 (bilayer hydrate) | **polaronic, SC killed ✗ vs expt 4.5 K** | λ ≫ 2, Na bound to one layer |

Sensitivity (λ ∝ γ_ep⁴ through χ²): γ_ep = 1.0 Å⁻¹ → window 6.4–6.5 Å, peak 11 K;
γ_ep = 2.0 Å⁻¹ → window ends 6.36 Å, peak 23 K. The window location (just above c\*) is
stable; its width (≈ 0.1 Å) and peak Tc (11–23 K) are not precision statements. The
interpolation of the well shape between the four DFT spacings adds comparable, unquantified
width uncertainty. What is robust: **the dome exists, sits between the 2DEG turn-on and the
polaron crossover, and contains neither the 6.9 Å nor the 9.9 Å anchor.**

Reconciling the observed SC at 9.9 Å therefore requires the H2O bilayer (78 % of that
gallery) to screen the Na–CoO2 attraction and pull the Na potential back toward the
symmetric, near-critical shape — the money plot carries this as an explicit **untested
prediction**; the decisive calculation is E(δ) for Na between explicit water layers.

## 6. Systematics (in one paragraph)

The dominant systematic is the **frozen O-plane height**: z_O = 0.96 Å for all c and δ, and
the checks show ±0.06 Å moves single-point energies by up to ~117 meV — the size of the
6.9 Å well itself. We therefore treat all **well depths as ±50 %** (and the c ≥ 8.4 depths
additionally as sextic extrapolations beyond the sampled δ, firm only as lower bounds);
because λ ≫ 2 by an order of magnitude in the frozen regime, no conclusion flips within
that band. **α sign changes are robust**: the 5.5 Å single well (α = +2.5 eV/Å²) and the
6.9 Å double well (well drop −89 meV reproduced within 0.2 % without spin polarization,
−1.6 meV under k-mesh refinement) cannot be undone by ±50 % scaling. Spin treatment
matters ≤ 6–7 % for well drops (§1), i.e. not at all for this model. The **ecut and
smearing relative checks are pending** (single points only so far). Not varied at all:
exchange-correlation functional, vdW corrections, in-plane lattice constant, alkali site
registry (Co-top only), x for the 1×1 series, and — most importantly — water.

## 7. Old vs new (proof-of-concept `phase_transition/` → canonical `results_v4`)

Reproduce the old numbers with `--legacy`.

* **Datasets.** Old: Questaal lmf, 19 closely spaced c values (Na: 5.55–7.15 Å), relaxed
  scan protocol, GPAW Bader on LiCoO2 only. New: QE PBE+PAW, 4 spacings but 3 series
  (Li/Na 1×1, Na √3 at the physical x = 1/3), spin-polarized, fixed z_O, plus N(0) and
  Bader for every series, plus convergence checks.
* **c\*_Na: consistent.** Old (dense scan): 6.07 Å. New: bracket (5.5, 6.9) Å, PCHIP
  6.40 Å. The old value lies inside the new bracket; the new interpolated value differs by
  +0.33 Å, within the resolution of a 4-point α(c).
* **c\*_Li.** Old: 5.38 Å. New: only the upper bound ≤ 5.5 Å (consistent; the new dataset
  has no point below the crossing).
* **Well depths: the new wells are ~3× shallower where comparable.** At c ≈ 6.9 Å:
  old 354 meV (d₀ = 1.04 Å) vs new 119 meV (d₀ = 0.69 Å). Old deepest (c = 7.15 Å) 499 meV
  vs new sextic-refit 758 meV at 8.4 Å but only ≥ 379 meV firmly sampled. Plausible origins:
  relaxed vs frozen z_O (a ±50 % lever by itself, §6), different codes/protocols, and site
  registry. The qualitative sequence (stiff single well → ~0.1 eV double well → deep
  frozen well) is unchanged, so every model conclusion survives; quantitatively the new,
  shallower 6.9 Å well moves the polaron crossover slightly up in c.
* **Carriers.** Old: single Bader series (LiCoO2, GPAW), onset at c = 5.29 Å, transferred
  to the Na axis by the c\* shift, and the TB step-DOS N(0) = 0.36 /eV/cell/spin ansatz.
  New: per-series DFT N(0) *and* Bader, no axis transfer, no DOS ansatz. The old N(0)
  guess was ~3× smaller than the DFT 1.06–1.36 — one reason the new λ values in the frozen
  regime are even larger.
* **Story change at 9.9 Å.** The old model reached 9.9 Å only by rigid extrapolation from
  7.15 Å; the new dataset actually samples it and shows something the extrapolation could
  not: E(δ) still decreasing at δ = 1 Å — Na headed for a *surface site* on one CoO2
  layer, while α itself re-flattens. The conclusion (polaronic, non-SC without water) is
  the same, but its mechanism at large c is now "adsorption onto one layer", not "deeper
  symmetric double well".

## 8. Assumptions, in one place

1. **1D decoupled alkali mode**, one alkali per cell, no intersite coupling, no in-plane
   motion, clamped CoO2 layers (frozen z_O — see §6).
2. **Even-polynomial truncation** of V(δ) (quartic, or capped sextic where the quartic
   fails); mirror symmetry of the gallery assumed (scans are one-sided).
3. **Hard-core cap** δ_min ≤ c/2 − 2.0 Å for the extrapolated minima (ansatz; varying the
   contact distance 1.8–2.5 Å changes frozen-regime depths, not any conclusion).
4. **Coefficient interpolation** (PCHIP in c) of (α, β, γ₆), N(0), δ(c) between four
   spacings — the resolution limit of this dataset.
5. **TB transferability**: SciPost SSH4 parameters (fitted for LiCoO2 slabs) used
   unchanged; γ_ep = 1/z0 exponential-hopping ansatz, z0 ∈ [0.5, 1] Å, λ ∝ γ_ep⁴ —
   dominant model uncertainty, full range shown.
6. **δ(c) from Bader via the paper's band-overlap formula**, Q_tot = 1 e, floor 0.05 eV.
7. **Pairing as single-boson Migdal–Eliashberg (Allen–Dynes)** with the explicit
   λ_loc = 2 validity cut instead of a real bipolaron theory; μ\* = 0.10.
8. **No water anywhere** — the single biggest caveat for the 9.9 Å anchor (see §5).

## 9. What would falsify this

* **No soft, strongly anharmonic Na c-axis mode** (meV-scale) near criticality in the
  bilayer hydrate (inelastic neutron/Raman), or Na seen *statically* off-center (split
  site, no dynamics) in the SC phase — kills the mechanism.
* **No interlayer/alkali-derived carriers** at E_F in the bilayer hydrate
  (ARPES/quantum oscillations/XAS). Known caution: published ARPES sees only Co-derived
  bands — already uncomfortable for the model.
* **Isotope/alkali-mass null result** at fixed structure (Li↔Na, H2O↔D2O).
* **c-axis knob**: the model predicts SC (or strong quantum-paraelectric-like dielectric
  response) whenever carriers are present and the system is held near c\* ≈ 6.4 Å —
  pressure on the monolayer hydrate (reducing c toward c\*) should *induce* SC before the
  gallery closes; its absence would falsify the window placement.
* **Water irrelevance test**: if AIMD/explicit-water E(δ) at 9.9 Å still shows Na deeply
  bound to one CoO2 layer, the model cannot explain the 4.5 K SC and is wrong or
  incomplete. This is the decisive next calculation.
* **Internal falsifiers already faced and reported**: linear deformation potential vanishes
  by parity (model rebuilt on the quadratic vertex); results.json's secant-`I` λ/Tc
  superseded; quartic fits fail at large c (refit with constrained sextic, depths labeled
  extrapolations); the 9.9 Å anchor fails under the water-free model — reported, not
  hidden.

## 10. Honest bottom line

results_v4 upgrades the two-sided story from a proof of concept to a verified, multi-series
statement: one critical gallery spacing bracketed by DFT points (5.5 < c\*_Na < 6.9 Å,
≈ 6.4 Å interpolated, only weakly x-dependent), at which a 2DEG turns on in two independent
probes while the Na mode goes critical, and an electron–phonon model built strictly from
that data plus published TB parameters produces a superconducting dome pinned to c\* with
Tc of the right order (peak 11–23 K vs experimental 4.5 K). The anhydrous and
monolayer-hydrate anchors come out correctly non-SC. What the model still does not do is
place the 9.9 Å bilayer hydrate inside the window: without water the DFT says Na simply
adsorbs onto one CoO2 layer there. The mechanism therefore stands or falls with the
untested prediction that interlayer water re-symmetrizes the Na potential — the next
calculation for this repo is E(δ) with explicit water.

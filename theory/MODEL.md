# Effective model: gallery-controlled Na double well + interlayer 2DEG superconductivity in Na_xCoO2·yH2O

**TL;DR.** The repo's DFT data show that opening the alkali gallery of A-CoO2 (A = Li, Na)
does two things *at the same critical spacing c\**: (i) the centered alkali position becomes
unstable (symmetric double well, Landau coefficient α(c) crossing zero), and (ii) Bader
charge starts flowing back onto the alkali (interlayer-2DEG turn-on). Quantizing the DFT
double well for the Na mass and wiring the displacement into the SSH4 tight-binding model of
Radha & Lambrecht gives an electron–phonon problem whose coupling turns on abruptly at c\*
and immediately runs into the strong-coupling/polaronic regime. The resulting Tc(c) is a
narrow sliver just above c\* — the superconducting window is bounded below by "no 2DEG" and
above by "Na frozen off-center", exactly the quantum-paraelectric-like optimum hypothesized —
but under a rigid (water-free) extrapolation the bilayer-hydrate anchor at 9.9 Å lands deep
in the frozen regime, not in the window. The honest conclusion is at the bottom.

Everything here is produced end-to-end by `theory/effective_model.py` from raw repo data.

## Reproduce

```bash
uv run theory/effective_model.py            # PEP 723 inline deps (numpy/scipy/matplotlib/pandas)
# or: theory/.venv/bin/python theory/effective_model.py
```

Outputs: `theory/results/*.csv` (all tables), `theory/figures/fig1–fig5` (diagnostics),
`theory/money_plot.{pdf,png}`.

---

## 1. Input data and provenance

| Quantity | Value / file | Provenance |
|---|---|---|
| E(d; c), LiCoO2, 19 spacings c = 4.72–5.92 Å | `phase_transition/licoo2/structures_done.txt` | repo DFT (Questaal lmf, energies in Ry, structures as pymatgen JSON) |
| E(d; c), NaCoO2, 19 spacings c = 5.55–7.15 Å | `phase_transition/nacoo2/structures_done.txt` | repo DFT (same protocol) |
| Bader charge on Li vs c (19 points) | `phase_transition/licoo2/charge_data` (pickle; read with a stub unpickler, no pymatgen needed) | repo DFT (GPAW + Bader on the minimum-energy structure per c) |
| In-plane lattice constants a(Li) = 2.843 Å, a(Na) = 2.920 Å | same JSON files | repo DFT |
| τ1 = t_Li^z = 0.5 eV, τ3 = t_CoO2^z = 0.5 eV, τ2 = τ4 = t_{A–CoO2}^z = 2.0 eV | SSH4 fit | Radha & Lambrecht, SciPost Phys. 10, 057 (2021) |
| t_Li^xy = −0.6 eV, t_CoO2^xy = +0.09 eV → Γ = 6·0.6 + 6·0.09 = 4.14 eV | in-plane fit | same paper |
| q_A/q_CoO2 = (Γ−2δ)/(Γ+2δ) for 2δ < Γ, else 0 | band-overlap formula | same paper (2D steplike-DOS argument) |
| γ = 1/z0 with t(z) ∝ e^{−z/z0}, z0 = 0.7 Å (range 0.5–1.0 Å → γ ∈ [1, 2] Å⁻¹) | **ansatz** | stated range; standard hopping-decay scale |
| μ\* = 0.1 | **ansatz** (task prescription) | standard value |
| λ_loc = 2 (polaronic-localization criterion, range 1–2) | **ansatz** | Migdal/bipolaron crossover literature scale |
| Q_tot = 1 e/cell (fully converted alkali band) | paper limit | SciPost paper (one alkali sp_z electron per cell) |
| Anchors: 5.5 Å anhydrous (no SC), 6.9 Å monolayer hydrate (no SC), 9.9 Å bilayer hydrate (Tc = 4.5 K) | experimental CoO2–CoO2 spacings | Takada et al. Nature 422, 53 (2003) and follow-ups (as given in the task) |

Notes on the data:

* The coordinate **c is the CoO2–CoO2 layer spacing along z** (z-component of the third
  lattice vector), not |v3|. This is the quantity the hydration anchors (5.5/6.9/9.9 Å)
  refer to. The notebooks' "c" (|v3|) is kept in the tables as `c_lat`. Reassuringly, the
  smallest NaCoO2 spacing in the data set, 5.55 Å, is exactly the anhydrous experimental
  spacing.
* The displacement d is computed geometrically per structure:
  d = z_alkali − (z_Co − c/2), i.e. offset from the gallery midplane. Each c-group contains
  the centered point d = 0 (energy reference) and 5–14 displaced points to |d| = 0.5–1.2 Å.
* Energies in `structures_done.txt` come from `lmf.py:get_potential_energy` →
  `save.<ctrl>` files, which are in **Ry** (confirmed against the notebooks' ×13.6057
  conversion).
* The notebooks (`interlayer Superconductivity.ipynb`) also reference plain-text scans
  `Li_energy/5.46.txt`, `Na_energy/5.99.txt`, … (bohr/Ry columns). **Those files are not in
  the repo**; the JSON records above cover the same physics (Li: c_lat 4.99–6.14 Å,
  Na: c_lat 5.80–7.35 Å) and are used instead. Nothing was reconstructed from memory.

## 2. Landau fits V(d; c) = α(c) d² + β(c) d⁴

Fits per c-group on the parity-symmetrized data (rms residuals 11–30 meV — a d⁶ term is
visibly missing at the largest |d| but irrelevant below ~0.4 eV). Full tables:
`results/well_fits_{licoo2,nacoo2}.csv`, figures `fig1`, `fig2`.

* **LiCoO2**: α: +1.45 → −0.80 eV/Å², zero crossing **c\*_Li = 5.38 Å** (spacing;
  = 5.62 Å in the notebooks' |v3| convention), slope dα/dc = −1.98 eV/Å²/Å.
  Deepest well 348 meV at d0 = 0.93 Å (c = 5.92 Å). This is the "single well → double well
  past c\*, depth ~0.35 eV" behavior described for the repo data.
* **NaCoO2**: α: +1.28 → −0.68 eV/Å², zero crossing **c\*_Na = 6.07 Å**, slope
  −1.96 eV/Å²/Å. Well position grows to d0 = 1.21 Å, depth to 499 meV at c = 7.15 Å.
  At the monolayer-hydrate spacing 6.9 Å (within the DFT range!): depth ≈ 350 meV,
  d0 ≈ 1.03 Å — the Na well is already deep at the *non-superconducting* hydrate spacing.
* β(c) is smooth and positive throughout (0.83 → 0.23 eV/Å⁴ for Na).

**Key DFT coincidence** (this is the empirical linchpin of the whole story): the Bader
charge on Li has a minimum at c = 5.29 Å and rises steadily beyond it (0.0133 → 0.0251 e,
i.e. ~+90 % above the post-onset minimum, +150 % relative to the extrapolated decreasing
baseline at the largest c) — the charge turn-on and the α = 0 instability occur at the
*same* spacing within one c-step (5.29 vs 5.38 Å). Gallery opening simultaneously creates
the soft Na mode and the interlayer 2DEG.

## 3. Quantized well (1D Schrödinger, finite differences)

Parity-resolved finite-difference solver (staggered half-line grid, Neumann/Dirichlet at
d = 0 for even/odd sectors — separates the tunneling doublet exactly; N = 1300, L = 2.6 Å;
levels converged to ≪ 1 %). Both Li (6.941 u) and Na (22.99 u) masses, both compounds:
`results/quantum_wells_*.csv`, figure `fig3`.

NaCoO2 wells, Na mass (the physical combination):

| c (Å) | ħω_eff = ħ²/2M⟨d²⟩₀ (meV) | E1−E0 tunneling (meV) | E2−E0 even gap (meV) | √⟨d²⟩₀ (Å) |
|---|---|---|---|---|
| 5.55 | 21.8 | 21.8 | 43.7 | 0.065 |
| 6.05 (≈c\*) | 5.95 | 5.93 | 13.2 | 0.124 |
| 6.15 | 1.26 | 0.77 | 6.6 | 0.269 |
| 6.25 | 0.51 | 7.6×10⁻³ | 11.3 | 0.423 |
| 6.45 | 0.23 | < 10⁻⁶ | 16.8 | 0.635 |
| 6.92 | 0.086 | < 10⁻¹² (below numerical resolution) | 21.7 | 1.03 |
| 7.15 | 0.062 | ≪ 10⁻¹² | 22.1 | 1.21 |

The mode softening near c\* is dramatic: ħω_eff drops by a factor ~20 between the anhydrous
spacing and just past c\*, and the even-sector gap E2−E0 (the excitation electrons actually
couple to, see below) has a *minimum of 6.6 meV at c = 6.15 Å* before re-hardening to the
within-well frequency ≈ 22 meV. Past c ≈ 6.3 Å the tunneling splitting collapses
exponentially — Na is frozen off-center on any experimental timescale. Li-mass curves
(dashed in fig3) are ~√(23/6.9) ≈ 1.8× stiffer; the LiCoO2 wells with the Li mass behave
nearly identically when shifted by c\*_Na − c\*_Li = 0.70 Å (fig3, bottom-left), which is
what justifies transferring the LiCoO2 Bader calibration to the Na axis.

## 4. Wiring d into the SSH4 model

4×4 Bloch Hamiltonian in the z-ordered basis (A1, A2, Co1, Co2), on-site ±δ, hoppings
τ1 (alkali sp_z pair), τ2/τ4 (alkali↔CoO2, modulated τ2 → τ2(1+γd), τ4 → τ4(1−γd)),
τ3 (across the CoO2 slab); in-plane triangular-lattice dispersions f_A, f_Co on the
diagonal. Paper parameters throughout. The alkali-derived band is the first band above the
charge-transfer gap (85 % alkali weight at δ ≈ 2 eV); its sp_z-pair partner sits 2τ1 = 1.0 eV
higher (the "surface/interlayer band splitting").

**Deformation potential — an honest structural point.** The SSH-type coupling is odd under
the mirror z → −z, so the *linear* deformation potential of the alkali band vanishes
identically at d = 0, at every k (verified numerically in the script; parity + time
reversal). The task's I = dE_band/dd is therefore zero at the symmetric position. The leading
electron–lattice vertex is quadratic:

E_band(d) = E_edge + ½ χ d²,  χ = 2 Σ_m |⟨alk|∂H/∂d|m⟩|²/(E_alk − E_m)

evaluated exactly (second-order perturbation theory is exact here since ∂²H/∂d² = 0 for
linearized hoppings). At the band bottom (k∥ = Γ, k_z = 0): χ = −6.1 eV/Å² for
γ = 1.43 Å⁻¹, δ = 2.05 eV, scaling as γ² (χ = −3.0 / −12.0 eV/Å² at γ = 1 / 2 Å⁻¹) and
growing ~2× from k_z = 0 to π (fig5, left). χ < 0: off-centering pulls the alkali band
*down* — the same sign as the Bader charge increase past c\*, an internal consistency check
(the DFT double well itself already contains this electronic energy gain; we use the DFT
V(d) directly, so there is no double counting — the TB is used only for the carrier-phonon
vertex).

**(a) Carrier density n(c) and the δ(c), Γ mapping.** The paper's band-overlap formula
q_alk = Q_tot (Γ−2δ)/2Γ (2δ < Γ) is calibrated on the repo's Bader data: the gallery charge
Δq_Li (Bader minus its pre-onset linear baseline) turns on at c_on = 5.29 Å with slope
dq/dc = 0.0207 e/Å. Taking δ(c) linear and crossing Γ/2 at c_on — i.e. gallery opening drives
the alkali level down through the CoO2 band top — fixes δ(c) = (Γ/2)(1 − (c−c_on)/w) with
w = Q_tot/(2·dq/dc) = 24.1 Å, transferred to the Na axis by the shift c\*_Na − c\*_Li.
This gives n_2DEG = q/A_cell ≈ 2×10¹² e/cm² at c\*+0.05 Å, 2.6×10¹³ at 6.9 Å, and (rigid
extrapolation) 1.1×10¹⁴ at 9.9 Å — sane 2DEG scales. The justification for "δ decreases as
the gallery opens" is precisely the Bader turn-on: charge on the alkali means the alkali
band edge has descended to the Fermi level, which in the paper's formula is 2δ → Γ.
In-plane mass of the alkali band bottom from the TB: m\* = 1.17 m_e →
N(0) = m\*/2πħ² = 0.049 eV⁻¹Å⁻² per spin = 0.36 states/(eV·cell·spin) (2D step DOS; N(0)
is c-independent once the band is occupied — occupancy is the gate).

**(b) Electron–phonon coupling.** Because the linear vertex vanishes, the electrons couple
to d² — i.e. to the *even-parity* excitation |0⟩ → |2⟩ of the anharmonic well (⟨0|d²|1⟩ = 0
by parity; the soft odd/tunneling channel decouples from intraband pairing). Single-boson
coupling:

g(c) = ½|χ(c)|·⟨0|d²|2⟩,  ħΩ(c) = E2 − E0,  λ(c) = 2 N(0) g²/ħΩ

with exact matrix elements from the solver. In the frozen regime this correctly reduces to
the ordered-phase result g → |χ|d0·ℓ_well, ħΩ → ħω_well (couple to fluctuations around ±d0),
and below c\* to the weak two-phonon vertex. The task's rms prescription
λ_task = N(0)χ²⟨d²⟩₀/(Mω_eff²) is tabulated alongside (`results/coupling_tc_vs_c.csv`);
it agrees with λ near and below c\* but diverges as ⟨d²⟩³ in the frozen regime because it
counts the *static* displacement d0 as dynamical coupling — we regard the even-channel λ as
the defensible one.

**(c) Allen–Dynes.** Tc = f1·(ħΩ/1.2k_B)·exp[−1.04(1+λ)/(λ−μ\*(1+0.62λ))], μ\* = 0.1,
single Einstein mode (ω_log = ⟨ω²⟩^½ = Ω, f2 = 1), strong-coupling factor f1 included.
Tc is set to zero where the 2DEG is absent (q = 0). Where λ > λ_loc = 2 the
Migdal–Eliashberg treatment is not trusted: with ħΩ ~ 10–20 meV against E_F ≲ 50 meV and
λ ≫ 1, the physical expectation is polaronic self-trapping / bipolaronic insulation, i.e.
**SC destroyed, not enhanced** — this region is drawn dashed and treated as non-SC.

## 5. Results (γ = 1.43 Å⁻¹ central value)

| c (Å) | λ (even channel) | ħΩ (meV) | raw Allen–Dynes Tc (K) | verdict |
|---|---|---|---|---|
| 6.07 (= c\*) | 0.29 | 11.9 | 0.05 | SC, tiny |
| 6.13 | 1.60 | 7.3 | **11.2** | SC — peak of the window |
| 6.5 | 3.6 | 18.0 | (51) | λ > 2: polaronic, non-SC |
| 6.9 (monolayer hydrate) | 5.9 | 21.6 | (83) | λ > 2: polaronic, non-SC ✓ matches "no SC" anchor |
| 9.9 (bilayer hydrate, rigid extrapolation) | 11.2 | 22.1 | (121) | λ > 2: polaronic — **disagrees with the SC anchor** |

* **SC window (γ = 1.43): c = 6.01–6.13 Å, peak Tc ≈ 11 K.** Bounded below by the 2DEG
  turn-on and above by λ crossing 2 as the well deepens — the hypothesized
  quantum-paraelectric-like optimum emerges, and its peak Tc is the right order of magnitude
  (experiment: 4.5 K).
* **Sensitivity (the honest error bar).** λ ∝ γ⁴ through χ², so z0 matters enormously:
  γ = 2 Å⁻¹ narrows the window to 6.01–6.09 Å (peak 18 K); γ = 1 Å⁻¹ widens it to
  6.01–7.32 Å with peak ≈ 42 K, and then the *monolayer* hydrate at 6.9 Å would be
  superconducting, contradicting experiment. Within this model the anchors actually
  *calibrate* γ: requiring 6.9 Å non-SC and a window near c\* wants γ ≳ 1.2 Å⁻¹
  (z0 ≲ 0.8 Å) — comfortably inside the physical range.
* **The 9.9 Å anchor does not land in the window** under the rigid (water-free)
  extrapolation: the model there has a ~0.5+ eV double well, frozen Na, λ ≈ 3–43
  (γ range) — a polaronic insulator, not a 4.5 K superconductor. For the model to be right,
  the H2O bilayer must qualitatively change the Na potential — plausible (in the bilayer
  hydrate Na is octahedrally solvated by water and sits in a screened, symmetric cage;
  ~78 % of the gallery is water, and the naked A–CoO2 electrostatics that create the double
  well are screened), but **not demonstrated here**. Equivalently: on the rigid-c axis the
  window sits at 6.0–6.1 Å where no stable anhydrous phase exists; hydration is the
  experimental knob that could re-place the system near criticality, and the bilayer→
  monolayer contrast (SC → no SC as water is removed and Na re-binds to the CoO2 layers)
  is at least directionally consistent with "SC lives near the symmetric, dynamically
  fluctuating Na state".

Diagnostics: `fig4` (Bader turn-on, δ(c), n(c)), `fig5` (χ(k_z; δ), λ(c) with γ band and
both prescriptions, mode energies). Money plot: `theory/money_plot.pdf`.

## 6. Assumptions, in one place

1. **1D decoupled Na mode**, one Na per cell, no Na–Na intersite coupling, no in-plane Na
   motion, classical c (clamped layers). DFT-cell constraint: one formula unit — collective
   Na ordering (known 2D Na-vacancy orderings in Na_xCoO2) is outside the model.
2. **Quartic truncation** of V(d): rms 11–30 meV; fine for the low levels used.
3. **Na stoichiometry x = 1** in the DFT scans, vs x ≈ 0.3 in the superconducting hydrate.
   Both the well and Q_tot would rescale with x; not modeled.
4. **Rigid extrapolation beyond c = 7.15 Å** (α, β frozen at their last DFT values;
   conservative, trends suggest slightly deeper wells) — and **no water** anywhere in the
   model, the single biggest caveat for the 9.9 Å anchor.
5. **TB transferability**: SciPost parameters were fitted for LiCoO2 slabs; used unchanged
   for NaCoO2 galleries (the paper argues alkali-species insensitivity of the sp_z band).
6. **γ from an exponential hopping ansatz**, z0 ∈ [0.5, 1] Å; λ ∝ γ⁴ — dominant
   uncertainty, full range shown everywhere.
7. **δ(c) linear, calibrated on the LiCoO2 Bader slope and transferred by the c\* shift**;
   Q_tot = 1 e. The 2D step-DOS N(0) and the q-formula come from the paper.
8. **Pairing treated as single-boson Migdal–Eliashberg (Allen–Dynes)** despite
   Ω/E_F ~ 0.1–0.5 (nonadiabatic) and λ up to ≫ 1; patched by the explicit λ_loc = 2
   validity criterion rather than by a real bipolaron theory.
9. Bader turn-on measured on **LiCoO2** (the only Bader series in the repo); NaCoO2 lmf
   screen charges exist (`nacoo2/data`) but were judged unreliable in the original notebook
   ("Not working out. Need to look at bader analysis") and are not used.

## 7. What would falsify this

* **No soft Na mode in the hydrates.** The model requires a low-energy (meV-scale),
  strongly anharmonic Na c-axis mode in the bilayer hydrate near criticality. Inelastic
  neutron/Raman spectroscopy seeing only a hard (≳ 20 meV), harmonic Na mode — or
  diffraction/EXAFS/NMR seeing Na *statically* off-center (frozen, split site along c with
  no dynamics) in the superconducting phase — kills the mechanism.
* **No interlayer 2DEG.** If ARPES/quantum oscillations/XAS in the bilayer hydrate show no
  alkali/interlayer-derived carriers at E_F (all weight in the Co t2g a1g/e'g sheets), the
  pairing glue may survive but the paired carriers proposed here do not. (Known caution:
  most ARPES on Na_xCoO2·yH2O sees Co-derived bands; an interlayer band has not been
  reported — this is already uncomfortable for the model.)
* **Isotope/mass test.** Tc should respond to the alkali-site dynamics: Li↔Na substitution
  at fixed structure, or H2O↔D2O (which shifts the gallery dynamics and screening), should
  shift Tc via ω and λ; a null isotope/alkali-mass effect with unchanged structure argues
  against an alkali-phonon glue.
* **c-axis knob.** The model predicts SC (or at least strong Tc enhancement and a
  quantum-paraelectric-like dielectric response ε(T) saturation) whenever the system is held
  near c\* ≈ 6.05–6.15 Å spacing with carriers present — e.g. by intercalating other
  neutral spacers, pressure on the monolayer hydrate (which *reduces* c toward c\*: the
  model predicts pressure-induced SC in the monolayer hydrate before it closes the gallery),
  or electrochemical gating. Observing the monolayer hydrate remain non-SC under modest
  compression toward 6.1 Å would falsify the placement of the window.
* **Water irrelevance test.** If the bilayer hydrate's Na potential is shown (by AIMD with
  explicit water — the natural next calculation for this repo) to be just as deep and
  double-welled as the rigid extrapolation says, the model cannot explain why 9.9 Å
  superconducts, and the mechanism is wrong or incomplete.
* **Internal numerical falsifiers already faced:** the linear deformation potential vanished
  (parity) — reported, model rebuilt on the quadratic vertex; λ_task diverges in the frozen
  regime — reported and replaced by exact matrix elements; the 9.9 Å anchor fails under
  rigid extrapolation — reported, not hidden.

## 8. Honest bottom line

The repo's DFT gives a clean, quantitative two-sided story at the *anhydrous* level: one
critical gallery spacing (c\*_Na = 6.07 Å) where a 2DEG turns on and the Na mode goes
critical, and an electron–phonon model built strictly from that data plus published TB
parameters produces a superconducting window pinned to c\* with Tc of order 10 K — the right
physics shape (dome bounded by "no carriers" and "frozen order") and the right Tc scale.
What the model does *not* do is place the 9.9 Å bilayer hydrate inside that window without
invoking water screening that is not yet computed. The monolayer-hydrate anchor (6.9 Å,
non-SC, model says polaronic/frozen) comes out correctly for γ ≳ 1.2 Å⁻¹. The decisive
missing calculation is E(d) for Na between explicit water layers.

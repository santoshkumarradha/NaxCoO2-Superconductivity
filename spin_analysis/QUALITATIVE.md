# QUALITATIVE.md — experiment-facing predictions of the interlayer-2DEG theory

**Scope.** This is the "even if the numbers are soft" document. The theory under test:
superconductivity in Na_xCoO2·1.3H2O (Tc = 4.5 K) lives in an interlayer (gallery) 2DEG
of alkali sp_z character, created by incomplete Na→CoO2 electron donation when the CoO2
plane-to-plane spacing `c` opens, coupled to an anharmonic Na double-well mode. The 2DEG
mechanics come from Radha & Lambrecht, SciPost Phys. 10, 057 (2021) ["R&L"]: a non-chiral
SSH4 model whose two surface bands (alkali sp_z, CoO2 Wannier) share one electron as

    q_alkali / q_CoO2 = (Γ − 2δ) / (Γ + 2δ)   if Γ > 2δ,  else q_alkali = 0
    Γ = Γ_alkali + Γ_CoO2 = 6·(|t_alkali^xy| + |t_CoO2^xy|)

with 2δ the splitting of the two surface-band centers and Γ_i the band half-widths
(R&L Sec. 3, App. E). Paper values: |t_Li^xy| = 0.6 eV, t_CoO2^xy = 0.09 eV → Γ = 4.14 eV;
Bader charge on Li ≈ 0.4 e ≈ 8×10^14 e/cm² at full coverage.

Each prediction below is stated as **prediction → existing fact → falsifier**. Numbers
marked *ansatz* are assumptions with a range, not measurements; quantitative back-up is
in `stoner_model.py` (run it — all numbers quoted here are printed by it).

---

## (a) Why Li cannot superconduct this way

### Chemical argument
**Prediction.** The mechanism requires the wide-open gallery (bilayer-water spacing,
c ≈ 9.9 Å). Li_xCoO2 does not form the bilayer hydrate (given experimental fact), so the
gallery never opens and Li_xCoO2·yH2O should never superconduct by this mechanism —
regardless of doping.
**Existing fact.** No superconducting Li_xCoO2 hydrate is known; only the Na bilayer
hydrate superconducts, and its monolayer hydrate does not.
**Falsifier.** Any Li-gallery compound with c ≳ 9.9 Å (e.g. organic-co-intercalated) that
still fails to superconduct at comparable carrier density would hurt the theory; a
superconducting Li compound with a *closed* gallery would kill this part outright.

### Model argument: the Li-vs-Na threshold inequality
The carrier condition is q_A > 0 ⟺ Γ(x) > 2δ_A. Diluting the alkali layer to coverage x
stretches the alkali–alkali distance to d(x) = a/√x, and the in-plane hopping decays
exponentially (*ansatz*, motivated by R&L App. C: half coverage in rows → bandwidth ~3×
smaller and the surface band empties; quarter coverage smaller still):

    Γ(x) = 6|t_CoO2^xy| + 6|t_A^xy| · exp[ −a(1/√x − 1)/λ_A ]

so the turn-on coverage is

    x*_A = [ 1 + (λ_A/a) · ln( 6|t_A^xy| / (2δ_A − 6|t_CoO2^xy|) ) ]^(−2)

Inputs (source-tagged; full list in `stoner_model.py`):
- |t_Li^xy| = 0.6 eV, t_CoO2^xy = 0.09 eV — **paper**.
- 2δ_Li = 0.828 eV — **derived**: q_Li(x=1) = 0.4 e (paper Bader) with Γ = 4.14 eV
  gives 2δ = Γ(1−2q).
- |t_Na^xy| = 0.9 eV — **ansatz [0.75, 1.10]**; R&L App. C states only that Na has
  "larger lateral hopping" and a larger electron pocket; the ratio ~1.5 reflects the
  larger 3s/3p orbital extent.
- λ_A = 1/√(2·IP): λ_Li = 0.84 Å, λ_Na = 0.86 Å — **ansatz [±0.12 Å]** (atomic decay
  lengths from ionization potentials 5.39 / 5.14 eV).
- δ_Na = δ_Li — **ansatz [1.00, 1.15]×δ_Li**, see honesty note below.

Result (baseline):

| | 6\|t^xy\| (eV) | 2δ (eV) | x* (turn-on) | Γ(x=0.35) (eV) | q(x=0.35) |
|---|---|---|---|---|---|
| Li | 3.60 | 0.828 | **0.326** | 0.895 | 0.037 |
| Na | 5.40 | 0.828 | **0.280** | 1.094 | 0.122 |

The inequality **x*_Li > x*_Na** holds throughout the ansatz ranges. At the
superconducting composition x ≈ 0.35, Li sits *on top of* its threshold: the baseline
leaves a residual q_Li ≈ 0.04 (within model error of zero), and a Monte-Carlo scan over
the stated ansatz ranges kills Li outright (q_Li = 0) in ≈ 22% of the parameter volume
and leaves a ~2–3× weaker 2DEG otherwise, while Na stays above threshold in ≈ 96% of the
same volume with q_Na ≈ 0.10–0.15 e/Na (n_2DEG ≈ 6×10^13 cm⁻²).

**Honesty notes (required by the derivation, not optional):**
1. The baseline does **not** give exactly q_Li = 0 at x = 0.35 — Li is *marginal*, not
   dead. The robust, falsifiable statements are (i) x*_Li > x*_Na, (ii) q_Na/q_Li ≳ 2–3
   at fixed x, (iii) Li's operating point x ≈ 0.35 lies within one ansatz-σ of its
   threshold while Na's does not.
2. The tasking assumed Li's disadvantage is partly "larger on-site δ". The electron-
   affinity/electronegativity data go the **other way**: Na (χ = 0.93, IP = 5.14 eV) has
   its level slightly *above* Li (χ = 0.98, IP = 5.39 eV), i.e. δ_Na ≳ δ_Li, which works
   *against* Na. Na's advantage is carried entirely by Γ (orbital extent → larger
   prefactor and slower exponential decay under dilution), and it survives a 15% δ
   penalty in the scan. This is also what R&L say qualitatively (App. C).

**Existing fact.** Bulk Li_xCoO2 is not superconducting at any x; Na_xCoO2·1.3H2O
superconducts at x ≈ 0.35.
**Falsifier.** (i) A measured Li-gallery band at E_F at x ≈ 0.35 dilution (STM/ARPES on
dilute Li-terminated surfaces) — R&L already predict the ×10 density collapse at x = 1/2;
(ii) DFT (the runpod Li vs Na job set) showing equal alkali Bader charge for Li and Na at
the same (c, x) — `runpod/postprocess.py` output will show this directly.

---

## (b) The on/off switch: c = 6.9 Å vs c = 9.9 Å

**Prediction.** The gallery band crosses its carrier turn-on between the monolayer-
hydrate spacing (6.9 Å) and the bilayer-hydrate spacing (9.9 Å). Below turn-on there are
no gallery carriers and no double-well: no SC. Above it, carriers + the anharmonic Na
rattler give SC. Water enters *only structurally* — anything that props the gallery open
to ~9.9 Å at x ≈ 0.35 should superconduct near 4.5 K, water or not.
**Existing fact.** Bilayer hydrate (9.9 Å) superconducts, Tc = 4.5 K; monolayer hydrate
(6.9 Å) does not, despite nominally identical x and nearly identical CoO2 sheets. This is
the single strongest qualitative fact in the theory's favor: an electronic-structure-of-
CoO2 explanation must work hard to distinguish 6.9 from 9.9 Å, while a gallery theory
gets it for free.
**Falsifier.** (i) The DFT batch showing alkali Bader charge / gallery-band occupation
essentially equal at c = 6.9 and 9.9 Å; (ii) experimentally, a non-aqueous intercalate
with c ≈ 9.9 Å and x ≈ 0.35 that does *not* superconduct, or any c ≈ 6.9 Å compound that
does; (iii) SC surviving when the gallery is chemically blocked (e.g. large cations
replacing the water network) at fixed c.

---

## (c) The Co-valence offset: +3.4 measured vs +3.65 nominal

**Prediction (as derived — with an honest sign warning).** If Na keeps a fraction q of
its electron in the gallery, a *Co-site-specific* probe (Co L-edge XAS, NQR) should read

    v_Co(site) = 4 − x(1 − q) ≈ 4 − 0.35×0.88 ≈ +3.69   (i.e. *above* nominal 3.65)

while a *charge-counting* probe (redox titration) that digests the whole crystal counts
the gallery electrons as Co reduction and reads lower. The theory therefore predicts a
**probe-dependent valence split** Δv = x·q ≈ 0.04 e/Co between site-specific and
titration valences, present only when the gallery is open (bilayer hydrate), absent in
the anhydrous parent.
**Existing fact.** Measured Co valence ≈ +3.4 vs +3.65 nominal (given). **Honesty note:
incomplete Na donation alone cannot produce this.** Keeping electrons *off* Co can only
push the site valence *up*, and total charge counting can never go below 4 − x = 3.65
without an additional electron reservoir. The offset of ~0.25 e/Co is ~6× larger than
x·q and requires extra donors (the literature's oxonium/H3O+ story). The theory is
*compatible* with this — extra donated electrons land partly in the same gallery band and
help Na clear the carrier threshold — but it does not *explain* the full offset, and we
say so.
**Falsifier.** Site-specific and charge-counting valences agreeing to better than ~0.04
e/Co in the bilayer hydrate (no gallery reservoir), or the valence offset persisting
unchanged in the monolayer hydrate where the theory says the gallery band is empty.

---

## (d) Negative pressure coefficient, dTc/dP < 0

**Prediction.** The c axis (soft, H-bonded gallery) is the most compressible direction;
hydrostatic pressure closes the gallery, moving the system from c = 9.9 Å back toward the
6.9 Å dead state: the Na double-well shallows, the rattler stiffens/harmonizes, q drops —
Tc must fall, and must fall *monotonically to zero* at the pressure where the gallery
band de-occupies (a first-order-like collapse if water is squeezed out). A subtlety the
spin model adds (see SPIN.md): pressure also pushes the system *toward* the Stoner
boundary, so the paramagnon channel briefly strengthens; the observed net negative sign
says the carrier/rattler loss wins — i.e. the pairing is not predominantly
paramagnon-driven.
**Existing fact.** dTc/dP < 0 (given).
**Falsifier.** (i) dTc/dP > 0 in any pressure window; (ii) purely in-plane (a-axis)
strain reproducing the Tc suppression at fixed c — the theory pins Tc to the *gallery*,
so in-plane strain should act only weakly (through Γ's lattice-constant dependence);
(iii) Tc surviving past the pressure where diffraction shows the bilayer→monolayer
gallery collapse.

---

## (e) D2O isotope effect through the water-network stiffness

**Prediction.** In this theory water is not the pairing boson; it is the *container*
whose H-bond network sets the stiffness and asymmetry of the Na double well (the rattler
is Na, mass unchanged by deuteration). D2O substitution changes the cage's librational
and stretch frequencies by ~√2 and its H-bond geometry slightly (stronger effective
H-bonds), which *renormalizes the well parameters* — a second-order effect on the rattler
frequency ω_Na and hence on Tc. Expected signature: a small, possibly *negative-α*
(inverse) isotope shift, |ΔTc/Tc| of at most a few %, i.e. an isotope exponent
α_D2O ≪ 0.5 (*ansatz on magnitude; the sign depends on whether D2O's stiffer network
deepens or shallows the well — undetermined here*).
**Existing fact.** None available to us here — we cite no D2O measurement and predict
blind.
**Falsifier.** α_D2O ≈ 0.5 (water phonons ARE the glue — kills the "water is structural"
reading, though not the gallery-2DEG itself); conversely a strictly null shift at high
precision (ΔTc/Tc < 10⁻³) would say the well is set by geometry alone, which the theory
tolerates but which removes one of its testable handles.

---

## (f) The spin phase diagram (summary; full story in SPIN.md)

**Prediction.** The gallery 2DEG has Stoner factor I_s·N(0) with N(0) ∝ 1/Γ_gallery, and
Γ_gallery grows with gallery opening. Hence a universal ordering of phases along any path
that opens the gallery or adds carriers:

    no 2DEG → **polarized 2DEG** (magnetic order; kills singlet SC; equal-spin f-wave
    triplet allowed) → **paramagnon-boosted paramagnet** (phonon singlet + fluctuation
    channel; this is where SC lives) → plain paramagnetic metal (weak SC at best)

The magnetic phase is a *thin strip hugging the carrier turn-on line* — polarization is
strongest just past turn-on (narrow band, high N(0)) and dies as the band widens. The
bilayer hydrate sits just on the paramagnetic side of the strip (Stoner S ≈ 2–3 in the
baseline model — enhanced but not ordered).
**Existing facts.** (i) Sakurai/Ihara/Takada (Physica C 514, 378 (2015)): a magnetically
ordered phase adjacent to/inside the SC dome vs Co valence — the strip. (ii) The
singlet-vs-triplet NMR contradiction (early powder: triplet-like; later single-crystal
Knight shift: singlet) — samples straddling the Stoner boundary, see SPIN.md. (iii) μSR:
no TRS breaking — consistent with a singlet or a real (unitary, TRS-preserving) f-wave
triplet, the two Mazin–Johannes survivors. (iv) R&L's DFT: the surface 2DEG *is* spin
polarized in the narrow-band (full-coverage surface) limit, LSDA+U- and QSGW-robust.
**Falsifiers.** (i) The runpod DFT batch (`parse_magnetization.py` is ready for it):
if moments do *not* peak just past carrier turn-on (intermediate c ≈ 8.4 Å) and do *not*
decay by c = 9.9 Å, the "polarization dies where SC lives" story fails. (ii) Experiment:
magnetic order found *inside* the optimal-Tc region (not adjacent), or on the wide-band
(low-valence/high-carrier) side of the dome, contradicts the strip topology. (iii) A
Knight-shift drop in *all* directions in a sample that also shows static magnetism would
break the sample-proximity reconciliation.

---

## One-page falsification table

| # | Prediction | Anchoring fact | Kills it |
|---|---|---|---|
| a | Li never clears the carrier threshold at x≈0.35; Na does (x*_Li=0.33 > x*_Na=0.28) | no Li bilayer hydrate; no Li SC | equal Li/Na gallery charge in DFT; dilute-Li surface band at E_F |
| b | Carrier turn-on between c=6.9 and 9.9 Å | bilayer SC / monolayer not | q(6.9)≈q(9.9) in DFT; non-SC 9.9 Å analogue |
| c | Probe-dependent Co valence split Δv = x·q ≈ 0.04 | +3.4 vs +3.65 (needs oxonium too — theory covers only ~1/6 of offset) | XAS/NQR = titration to <0.04 |
| d | dTc/dP < 0, → 0 at gallery collapse | dTc/dP < 0 | dTc/dP > 0; in-plane strain equivalence |
| e | Small anomalous D2O shift, α ≪ 0.5 | (blind prediction) | α ≈ 0.5 |
| f | Magnetic strip hugging turn-on, adjacent to SC dome | Sakurai 2015; NMR contradiction; μSR | moments not peaking at intermediate c; magnetism on wide-band side |

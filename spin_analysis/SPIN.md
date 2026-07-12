# SPIN.md — the spin story of the gallery 2DEG

Companion code: `stoner_model.py` (runnable, numpy/matplotlib) → `spin_phase_diagram.pdf`.
Every assumed parameter is tagged ANSATZ with a range in the script's `PARAMS` block and
printed at runtime; nothing below silently depends on an unlabeled number.

## 1. Where the magnetism comes from

Radha & Lambrecht (SciPost Phys. 10, 057 (2021)) find, in DFT, that the alkali-side 2DEG
on Li(Na)CoO2 layers is **spin polarized**, and attribute it to a **Stoner instability**:
"the cause for spin splitting arises from the avoidance of a high density of states at
the Fermi level" produced by the "rather flat" surface band near Γ (their App. D). The
result survives LSDA+U (U up to 3.4 eV) and QSGW. Their moment decomposition (their
Table 3, LSDA: Co +0.046 μB, Li +0.004, O1 +0.006, smooth/2DEG −0.053 μB) shows the
gallery 2DEG polarized **antiparallel** to the CoO2-side moment — a two-component
magnet, not a single Stoner band.

**Microscopic caveat, stated up front.** In the paper the strongest DOS peak is the flat
*CoO2-side* surface band; the gallery 2DEG polarizes with it (oppositely). In the model
below we follow the simpler effective description — a single gallery band with Stoner
parameter I_s and DOS N(0) ∝ 1/Γ_gallery — because (i) both bands cross E_F together
(they share the same carrier turn-on), (ii) both N(0)'s collapse as the gallery opens
and the bands widen/fill, so the *control-parameter dependence* is the same, and (iii)
the DFT batch (`parse_magnetization.py`) will tell us which site actually carries the
moment as a function of (c, δ). If the batch shows the moment residing ~fully on Co with
no gallery counter-moment, the "2DEG Stoner" language should be replaced by "flat-band
CoO2 Stoner with 2DEG bystander" — the phase-diagram topology below is unchanged, the
pairing implications in §4 are not.

## 2. The model (all in `stoner_model.py`)

Carriers (paper's overlap formula):

    q = (Γ − 2δ)/(2Γ),  Γ = Γ_Co + Γ_A(x, c),  q = 0 if Γ ≤ 2δ
    Γ_Co = 6·t_CoO2^xy = 0.54 eV                       [paper]
    Γ_A(x, c) = Γ0(c) · exp[−a(1/√x − 1)/λ_Na]         [dilution ansatz]
    2δ = 0.828 eV                                       [derived from paper Bader q]

Gallery opening: Γ0(c) linear (ANSATZ), calibrated by exactly two anchors — the x = 0.35
carrier turn-on must fall between the monolayer hydrate (c = 6.9 Å, no SC) and the
bilayer hydrate (c = 9.9 Å, SC), taken at c_on = 8.0 Å (ansatz range 7.2–9.2), and
Γ0(9.9 Å) = 6|t_Na^xy| = 5.4 eV (fully developed orbital). Physical content: opening the
gallery decompresses the Na sp_z orbital and develops the double well, widening the
gallery band and pulling it under E_F.

Stoner criterion on the gallery band. For a 2D triangular-lattice band bottom,
m* = ħ²/(3 t d²), the DOS per spin per alkali site is

    N(0) = ν / Γ_A,   ν = √3/2π = 0.276        [derived, parabolic 2D]

    polarized  ⇔  I_s · N(0) > 1,   I_s = 1.3 eV  [ANSATZ 1.0–1.6 eV]

I_s is the effective exchange of the gallery band (alkali sp_z with Co-d admixture; bare
Co-d Stoner I ≈ 0.9 eV for reference). The script prints the existence condition for the
polarized phase: **the strip exists iff I_s > (2δ − Γ_Co)/ν = 1.04 eV** — i.e. the model
only produces a magnetic phase at all if the effective exchange modestly exceeds 1 eV.
That the paper's DFT *does* find polarization is the empirical reason to take
I_s ≳ 1 eV; a bare-alkali I_s ≈ 0.5 eV would leave the diagram magnetism-free. This is
an honest hinge of the story and exactly what the DFT batch measures.

Stoner enhancement in the paramagnetic phase: S = 1/(1 − I_s N(0)); we call the metal
"paramagnon-boosted" when S > S* = 3 (ANSATZ 2–5).

## 3. The phase diagram (`spin_phase_diagram.pdf`)

Three carrier-bearing regimes plus the empty one, in the (x, c) plane:

1. **No 2DEG** (Γ < 2δ): gallery band empty. Contains the monolayer hydrate
   (x = 0.35, c = 6.9). No carriers → no SC, no gallery magnetism.
2. **Polarized 2DEG** (q > 0, I_s N(0) > 1): a *thin strip hugging the turn-on line*
   (at x = 0.35: c ≈ 8.0–8.5 Å in the baseline). Narrow band ⇒ high N(0) ⇒ Stoner
   order, strongest just past turn-on; in mean field the low-density 2D band polarizes
   completely (half-metallic gallery).
3. **Paramagnon-boosted paramagnet** (q > 0, S > 3): c ≈ 8.5–9.8 Å at x = 0.35.
4. **Plain paramagnetic 2DEG** (S < 3): wide band, weak correlations.

The **bilayer hydrate lands at S ≈ 2.8, q ≈ 0.12** — just on the paramagnetic side of
the strip, carrier-rich, with strong but non-critical ferromagnetic fluctuations. That
placement is the theory's central spin statement: **superconductivity lives one step
beyond the magnetism, and the magnetism lives one step beyond the turn-on.**
(The baseline S ≈ 2.8 vs S* = 3 coincidence is an artifact of round ansatz choices —
the robust statement is the *ordering* of the four regions and the thinness of the
strip, which survive any parameter choice that keeps the strip existing at all.)

## 4. Mapping onto the experimental record

**Mazin–Johannes f-wave survivors.** Symmetry analysis (Mazin & Johannes, Nat. Phys. 1,
91 (2005)) narrowed the pairing to two f-wave *triplet* states. On the paramagnetic side
of a ferromagnetic (q ≈ 0) Stoner boundary, the fluctuation spectrum is FM-like:
paramagnon exchange is attractive in odd-parity equal-spin channels and pair-breaking
for singlets. On the triangular gallery lattice the odd-parity harmonics compatible with
a nodeless-or-symmetry-noded gap are exactly the f-waves. So the theory's region-3
fluctuations *select the same channels Mazin–Johannes kept*: phonon/rattler singlet and
paramagnon f-wave triplet are the two live options, possibly competing within one dome.
Both f-wave states are unitary and TRS-preserving → **consistent with the μSR null**
(no spontaneous fields below Tc).

**The singlet-vs-triplet NMR contradiction.** Early powder NMR suggested triplet
(Knight shift not suppressed); later single-crystal work found a singlet-like Knight
shift drop. In this diagram that is not a contradiction but a *location statement*:
S(x, c) varies steeply near the strip, and (x, c) vary between samples (hydration is
notoriously fragile; c and even x drift with water content and handling). A sample
sitting at S ≫ 1 (closer to the strip) has its spin response dominated by
paramagnons — weak Knight-shift suppression, triplet-favoring — while a sample deeper in
region 3/4 shows the plain singlet drop. Prediction: Knight-shift suppression below Tc
should *anti-correlate* with the normal-state Stoner enhancement (measurable as the
normal-state 1/T1T or χ_spin) across samples. SOC adds a second, smaller contribution to
residual Knight shift (SOC.md).

**Sakurai's magnetic phase (Physica C 514, 378 (2015)).** A magnetically ordered phase
adjacent to/inside the SC dome as a function of Co valence. In the model, moving along
the valence axis at fixed structure moves the gallery filling/width, and the magnetic
strip is *necessarily adjacent* to the SC region because both are pinned to the same
turn-on line — magnetism just past turn-on, SC in the enhanced-paramagnet band beyond
it. Prediction of *side*: the magnetic phase sits on the **low-carrier / high-valence
edge** of the dome (nearer turn-on). Samples inside the dome showing static order are,
in this reading, spatially inhomogeneous — hydration gradients put minority regions
inside the strip; the μSR volume fraction of such order should track water
non-stoichiometry, not Tc.

**R&L's polarized surface 2DEG.** The full-coverage (x = 1) *surface* case is the
narrow-band end-member where DFT already finds the polarization; it anchors I_s ≳ 1 eV.
Note honestly: at x = 1 the alkali band is *wide* (Γ_A = 3.6–5.4 eV), and the paper's
polarization there is driven by the flat CoO2 surface band — see the §1 caveat; the
surface calculation cannot by itself calibrate the gallery strip position.

## 5. What the DFT batch will decide (via `parse_magnetization.py`)

The runpod grid is (element ∈ {Li, Na}) × (c ∈ {5.5, 6.9, 8.4, 9.9 Å}) × (δ ∈ 0–1 Å),
spin-polarized. The model predicts, for the absolute magnetization and the site moments:

1. c = 5.5, 6.9 Å: small/zero gallery moment (band empty); any moment is Co-local.
2. c = 8.4 Å: **maximum** gallery-related moment (just past turn-on in the baseline;
   if turn-on is nearer 9 Å the maximum shifts to between 8.4 and 9.9).
3. c = 9.9 Å: moment *decreasing* again — "polarization dies where SC lives".
4. Na vs Li at fixed (c, δ): Na should show the larger gallery moment onset shift
   (its band turns on earlier in c) — and, per §1, the sign of the smooth/gallery
   moment relative to Co decides the two-component question.

If instead the moment grows monotonically with c, or is element-independent, the spin
story as told here is wrong — report that, don't massage it.

## 6. Parameter table (everything assumed, with ranges)

| Parameter | Value | Range | Status |
|---|---|---|---|
| \|t_Li^xy\| | 0.6 eV | — | paper |
| t_CoO2^xy | 0.09 eV | — | paper |
| \|t_Na^xy\| | 0.9 eV | 0.75–1.10 | ansatz |
| 2δ (Li, full coverage) | 0.828 eV | 0.76–0.90 | derived from paper q=0.4e |
| δ_Na/δ_Li | 1.0 | 1.0–1.15 | ansatz (honest direction: ≥1) |
| λ_Li, λ_Na | 0.84, 0.86 Å | ±0.12 | ansatz (1/√(2·IP)) |
| ν (DOS shape) | 0.276 | — | derived (2D parabolic) |
| I_s | 1.3 eV | 1.0–1.6 | ansatz; strip exists iff > 1.04 |
| S* | 3 | 2–5 | ansatz (labeling threshold) |
| c_on(x=0.35) | 8.0 Å | 7.2–9.2 | ansatz, anchored ∈ (6.9, 9.9) |
| Γ0(c) | linear | — | ansatz (weakest link) |

Weakest links, in order: the Γ0(c) law (shape unconstrained between the two anchors),
I_s (existence of the strip), and the single-band effective treatment (§1 caveat).

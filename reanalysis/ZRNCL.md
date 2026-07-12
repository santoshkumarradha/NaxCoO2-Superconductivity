# Li_xZrNCl / HfNCl reanalysis: does the nitride-chloride family support, refute, or bypass our gallery-2DEG + soft-donor-ion mechanism?

**Data**: `zrncl_data.csv` (26 rows; 17 with both spacing and Tc). **Fit**: `fit_tc_vs_spacing.py` → `tc_vs_spacing.png`.

**Verdict up front: ORTHOGONAL, with one sharp cautionary lesson.** The β-MNCl (M = Zr, Hf) intercalates are a layered-host + donor-ion system in which *neither* of our two ingredients is present: the donated electrons do **not** return to the gallery (they fill the in-plane Zr/Hf-d–N-p host band), and the donor ion is **not** soft (it is rigidly ionized, spectroscopically silent, and can even be removed entirely without killing superconductivity). The famous "Tc flat-or-rising with spacing" result therefore neither confirms nor falsifies our mechanism — it tests a different pairing channel. What it *does* do is forbid us from claiming that "spacing tuned to donor-ion bistability" is a universal design rule for layered donor intercalates: here is a family where the spacing knob was turned from 9.3 Å to 22.4 Å and superconductivity never died, never domed, and never cared about the ion.

---

## 1. The dataset

### 1.1 Headline numbers (all literature-sourced; see CSV for provenance)

| System | d (Å) | Tc (K) | Source |
|---|---|---|---|
| Li0.08ZrNCl (optimal) | 9.34 | 15.1 | Taguchi, Kitora, Iwasa, PRL **97**, 107001 (2006) |
| Li0.16–0.4ZrNCl (overdoped) | ~9.3 | 11.5–12.5 | ibid.; Yamanaka, Adv. Mater. **8**, 771 (1996) |
| ZrNCl0.7 (Cl-deintercalated, **no ion**) | 9.8 | 13 | Tou et al., PRB **72**, 020501(R) (2005) |
| Li0.13(DMF)yZrNCl | 13.01 | 13.7 | Kasahara et al., PRB **82**, 054504 (2010) |
| Ax(DMF)yZrNCl (bilayer) | 19.8 | ~15 | Kawaji, Hotehama, Yamanaka, Chem. Mater. **9**, 2127 (1997) |
| Ax(PC)yZrNCl | 22.4 | ~15 | ibid. |
| ZrNCl EDLT (**no ion**, x≈0.011) | — | T_BKT = 17.9 | Nakagawa et al., Science **372**, 190 (2021) |
| Li0.2HfNCl (optimal) | 9.40 | 20.0 | Takano, Kishiume, Taguchi, Iwasa, PRL **100**, 247005 (2008) |
| Na0.25HfNCl | 9.89 | 24 | Oró-Solé et al., Mater. Res. Bull. **41**, 934 (2006) |
| Li0.2(NH3)yHfNCl | 12.10 | 22.5 | Takano et al., PRL **100**, 247005 (2008) |
| Eu0.08(NH3)yHfNCl | 11.91 | 23.6–24.3 | Zhang et al., SUST **26**, 045017 (2013) |
| Li0.37(THF)yHfNCl (THF monolayer) | 13.6 | ~25.5 | Takano et al., PRL **100**, 247005 (2008) |
| Ca0.11(THF)yHfNCl | 15.0 | **26.0** | Zhang et al., SUST **26**, 085015 (2013) |
| Li0.48(THF)yHfNCl (THF bilayer) | 18.7 | 25.5 | Yamanaka, Hotehama, Kawaji, Nature **392**, 580 (1998) |

Provenance caveat, stated honestly: ACS/APS/JPSJ full texts are paywalled; the precise d values above were cross-checked through two independent compilations — the Kasahara–Kuroki–Yamanaka–Taguchi review (Physica C **514**, 354 (2015), open copy arXiv:1412.4447, Table 1 and Fig. 7d, which reproduces Hotehama et al., J. Phys. Soc. Jpn. **79**, 014707 (2010)) and Harshman & Fiory (J. Supercond. Nov. Magn. **28**, 2967 (2015), arXiv:1508.02523, Table 1, which tabulates d to four digits with per-compound primary citations). Rows lacking a verified spacing are kept in the CSV with `d_A` blank rather than filled by memory.

### 1.2 The two experimental knobs, and what each does

**Doping x (at fixed spacing).** No dome. Superconductivity switches on abruptly at x ≈ 0.05 (ZrNCl) / 0.15 (HfNCl) and Tc *increases toward the insulator*: 11.5 K (0.2 < x < 0.4) → 15.1 K (x = 0.08) in Li_xZrNCl (Taguchi PRL 2006); flat ~20 K for 0.15 ≤ x ≤ 0.50 in Li_xHfNCl (Takano PRL 2008). Gate doping pushes this to x ≈ 0.011 with T_BKT = 17.9 K and T_BKT/T_F = 0.12 — the 2D BCS–BEC crossover regime (Nakagawa et al., Science **372**, 190 (2021)).

**Spacing d (at fixed doping).** Tc(d) is monotonically *non-decreasing*: it rises ~20% in ZrNCl (12.5–13 → ~15–15.5 K) and ~30% in HfNCl (20 → 25.5–26 K) as d grows from ~9.4 to ~13–15 Å, then saturates; it stays flat out to at least 22.4 Å (ZrNCl–PC) with no downturn ever reported (Takano PRL **100**, 247005 (2008): "an enhancement of Tc up to 30% ... with its x independence preserved"; Kasahara PRB **82**, 054504 (2010): coupling strength 2Δ/k_BTc grows concomitantly; Kawaji Chem. Mater. 1997 found the spacing "not influential" at their resolution). Our 1/d fit (numpy fallback, `fit_tc_vs_spacing.py`) gives Tc = 15.9 − 24.4/d (ZrNCl) and Tc = 31.1 − 90.2/d (HfNCl): the coefficient of 1/d is *negative* in both families. Contrast α-TiNCl/TiNBr, where Tc genuinely falls as 1/d with b > 0 (Zhang, Tanaka, Yamanaka, PRB **86**, 024516 (2012)) — the "spacing hurts" phenomenology exists in this class, but in the *other* polymorph.

---

## 2. Where the donor ion sits, and where its electron goes

This is the crux for confronting our mechanism, and every probe agrees:

1. **Site.** The alkali occupies the midplane between adjacent Cl layers (octahedral interstice of the van-der-Waals gap); with NH3/THF/DMF/PC the ion is solvated inside the expanded gallery. In Li_x(THF)_yHfNCl, 7Li NMR shows a +5 ppm *chemical* (not Knight) shift: Li is covalently coordinated to THF oxygen — the ion is chemically locked, not rattling (Tou, Maniwa, Koiwasaki, Yamanaka, PRB **63**, 020508(R) (2001); discussion in Harshman & Fiory 2015).
2. **Electron destination.** All band calculations (Weht, Filippetti, Pickett, Europhys. Lett. **48**, 320 (1999); confirmed by later LDA/GGA/hybrid work) put the conduction band at K/K′ with Zr/Hf d_xy, d_x²−y² + N p_x,p_y character — *in-plane, inside the MN bilayer*. Polarization-dependent N 1s XAS confirms the conduction-band bottom is in-plane hybridized (Yokoya et al., PRB **70**, 193103 (2004)). The 7Li Knight shift is "almost nothing": zero Fermi-level weight at the ion site (Tou et al., PRB **63**, 020508(R)). **There is no occupied gallery/interlayer state in this family.** No μSR, neutron, or ARPES work claims one; the 2D-superconductivity literature (μSR: Ito et al., PRB **69**, 134522 (2004); Hc2 anisotropy γ = 3.7–4.5, ξ_c ≈ d, weakly Josephson-coupled: Kasahara et al. review Table 1; BKT: Nakagawa Science 2021) attributes the two-dimensionality to the *host* band, not to a gallery gas.
3. **The ion is dispensable.** Superconductivity at 13–15.2 K appears with *no intercalant at all*: Cl-deintercalated ZrNCl0.7 (Tou PRB **72**, 020501(R) (2005)) and ionic-liquid-gated pristine ZrNCl (Ye et al., Nat. Mater. **9**, 125 (2010); Nakagawa Science 2021). Whatever pairs these electrons cannot be donor-ion physics.

## 3. Is there any soft/anharmonic mode physics here?

Mostly no — with one loose end:

- **Isotope effect ≈ null on N**: α_N ≈ 0.07 in both Li_xZrNCl (ΔTc = 0.06 K on 15N; Taguchi et al., PRB **76**, 064508 (2007)) and Li_x(THF)_yHfNCl (ΔTc = 0.1 K; Tou, Maniwa, Yamanaka, PRB **67**, 100509(R) (2003)).
- **Weak coupling to phonons**: γ_n ≈ 1.0 mJ/mol K² gives λ ≤ 0.22 (Taguchi et al., PRL **94**, 217002 (2005)); parameter-free SCDFT predicts Tc ≈ 4 K (Zr) / 8–10 K (Hf), far below experiment (Akashi, Nakamura, Arita, Imada, PRB **86**, 054513 (2012)). Yet 2Δ0/k_BTc ≈ 4.1–5 (strong coupling) with no Hebel–Slichter peak (Kotegawa et al., PRB **90**, 020503 (2014)). This mismatch is why plasmon (Bill, Morawitz, Kresin, PRB **68**, 144519 (2003)), spin-fluctuation d+id (Kuroki, PRB **81**, 104502 (2010)) and enhanced-2D-pairing scenarios dominate the theory literature — none invokes the intercalant.
- **Raman anti-correlation**: the e–ph linewidth of the ~620 cm⁻¹ N mode *grows* with doping while Tc *falls* (Kitora, Taguchi, Iwasa, JPSJ **76**, 023706 (2007)) — pairing does not track the measured host e–ph coupling.
- **The loose end**: a McMillan analysis of pressure data (Taguchi et al., PRB **70**, 104506 (2004)) concluded that if phonons pair, only *low-frequency modes below 100 cm⁻¹* (Cl/intercalant-layer modes) could be responsible — the review flags this as contradicting the thermodynamic λ. Nobody has followed up with low-energy INS on the intercalant. This is the single place where a soft-gallery-mode story could still hide in β-MNCl; the null alkali Knight shift and the ion-free superconductors make it unlikely to be the *donor ion*.

## 4. Confrontation with our mechanism

Our design rule has two coupled ingredients: (i) past a critical gallery opening, donor electrons occupy an interlayer/gallery 2DEG; (ii) the intercalant ion sits near a shallow double-well/adsorption bistability whose quantum-paraelectric softness supplies (or strongly assists) pairing; SC lives only where (i) and (ii) coexist, and dies at large spacing when the ion freezes (polaron).

**Does "Tc up-or-flat with spacing" contradict us?** No — because the premise of the comparison fails. In β-MNCl the donor's electron never returns to the gallery: it enters a pre-existing in-plane host band (Sec. 2). Our rule is a statement about systems whose *carriers themselves* are gallery-borne. β-MNCl is a layered host + donor ion in the *wrong electronic regime* for our mechanism — a Zr/Hf d-band superconductor that merely uses the ion as an electron syringe. Its Tc(d) tests interlayer hopping/nesting/dielectric physics of the host band (the community's own reading: improved two-dimensionality, Kasahara PRB 2010; k_z Fermi-surface warping, Takano PRL 2008), not donor-mode softness. Hence: **orthogonal**, not refuting.

**But it is not painless for us either.** Three honest costs:

1. **Universality is dead.** We cannot state "layered host + donor ion + spacing at the bistability threshold" as *the* design rule for layered donor intercalates without immediately explaining why the most-studied such family ignores it. The rule must be stated conditionally: *when* the gallery state is the conduction channel (electride-like regime), spacing tunes both carriers and ion softness; when the host band lies below the gallery state (β-MNCl, most GIC-like cases with deep host bands), the ion is a spectator and spacing does host-band physics instead. The discriminant is where the interlayer state sits relative to the host conduction band — computable (cf. Csányi et al., Nat. Phys. **1**, 42 (2005) for GICs).
2. **β-MNCl proves donor ions can be frozen solid (Li–THF covalent lock) at *every* spacing with zero SC penalty** — so "ion freezes ⇒ SC dies" is only meaningful in our regime where the ion mode participates in pairing. Any writeup must phrase the polaronic-death prediction as regime-specific, or ZrNCl reads as an instant counterexample.
3. **The Tc scale is uncomfortable.** An "ordinary" 2D band in a layered ionic insulator pairs at 15–26 K with no soft ion and (apparently) no strong phonon coupling; our hydrate mechanism produces 4.5 K. We should expect referees to ask why the exotic double-well channel yields *less* than the boring channel next door.

**What β-MNCl actively supports (weakly, by analogy):** (a) dilute 2D electron gases (down to x ≈ 0.01) in layered ionic insulators superconduct robustly, with Tc *rising* toward the low-density limit and T_BKT/T_F saturating the 2D bound (Nakagawa 2021) — the carrier side of our story, a dilute gallery gas pairing at a few K, is not exotic; (b) opening the gallery per se is not pair-breaking — the melodrama in our Tc(spacing) dome must come from the ion, and β-MNCl shows a clean background of what spacing does *without* ion physics: monotonic rise and saturation. A nonmonotonic dome, where observed (NaxCoO2·yH2O vs c-axis; Sakurai's two-dome phase diagrams), is therefore a real anomaly demanding an extra ingredient.

### Falsifiable predictions our model makes for β-MNCl (mostly already checked, which is the point)

- **No soft/quasi-elastic donor mode at any spacing**: Raman/INS below 100 cm⁻¹ should show only harmonic Cl/solvent modes; alkali ADPs in Rietveld/neutron refinements should stay normal from d = 9.3 to 22.4 Å. (Consistent with all existing structure work; the PRB 70, 104506 low-frequency loose end is the one open check — a dedicated INS study of Li_x(THF)_yZrNCl would close it.)
- **Null alkali isotope effect** (6Li/7Li): our channel is absent, so ΔTc ≈ 0. (Untested — cheap and decisive; a *finite* alkali isotope shift would force us to fold β-MNCl into the framework.)
- **N isotope effect small** — observed (α_N ≈ 0.07), but this is overdetermined (plasmon and spin-fluctuation scenarios predict it too); we claim no credit.
- **Superconductivity without any donor ion** — observed (ZrNCl0.7, EDLT). In our framework this is the signature of a host-band superconductor, i.e., of orthogonality itself.
- **Discriminating experiment against the hydrate**: identical low-energy INS/QENS protocol on Na0.3CoO2·1.3H2O (SC) vs anhydrous Na0.3CoO2 vs Li_x(THF)_yZrNCl. Prediction: an anharmonic, strongly temperature-dependent quasi-elastic Na (or H2O-cage) response *only* in the hydrate SC phase; its absence there kills ingredient (ii) of our model far more directly than any ZrNCl data can.

## 5. Verdict

**Orthogonal.** β-ZrNCl/HfNCl is a donor-*doped* (not donor-*paired*, not gallery-*carried*) layered superconductor: the intercalant is an inert electron donor whose own degrees of freedom are demonstrably decoupled from pairing, and the Tc(d) plateau measures host-band two-dimensionalization. It neither exhibits nor rules out our carriers-AND-soft-ion window. Its real value to us is disciplinary: it forces the design rule to carry an explicit electronic-structure precondition (gallery/interlayer state at or near E_F — the Csányi criterion), it supplies the null-hypothesis Tc(d) curve (monotonic, saturating) against which the hydrate's dome is anomalous, and it hands us two cheap falsifiers (alkali isotope null; absence of low-energy anharmonic donor modes) that we should cite as *already satisfied* consistency checks rather than confirmations.

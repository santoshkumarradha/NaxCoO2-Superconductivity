# NaxCoO2·yH2O experimental Tc dome: literature data extraction

**Data**: `experimental_dome.csv` (20 rows, all from full-text tables of three primary papers — no figure-only guessing was needed).

**Bottom line up front**: the experimental "dome" everyone plots is a dome in **Na content x / Co valence**, not in the CoO2–CoO2 spacing. The spacing *does* vary across the superconducting samples, but only within a ~1-2% window (c_hex ≈ 19.4–19.9 Å, i.e. CoO2–CoO2 spacing ≈ 9.7–9.9 Å), and within that narrow window spacing and Tc are **not monotonically related** — the highest-c sample in the Schaak et al. set (x=0.26, c=19.77 Å) has one of the *lowest* Tc's (2.4 K), while the highest-Tc samples (4.7–4.8 K, Milne et al.) sit in the middle of the observed c range. The real control knob, established definitively by Milne et al. (PRL 2004) via direct redox titration, is **Co valence** (equivalently hole doping), with optimal Tc for Co valence 3.24–3.35 and Tc collapsing above 3.35. c-axis expansion is a *correlate* of moving into that optimal valence window (going from anhydrous, c≈11.2 Å, non-SC, to hydrate, c≈19.4-19.9 Å, SC), not an independent tuning axis with its own dome shape.

---

## 1. c-axis → CoO2-CoO2 spacing conversion used

The bilayer-hydrate phase crystallizes in space group **P6₃/mmc** (or the doubled-cell **P6₃/m**, Argyriou et al. 2004) with **two CoO2 sheets per hexagonal c-axis repeat** (standard for this O3-derived stacking, same convention as the parent anhydrous NaxCoO2, c≈11.2 Å ⇒ single-layer spacing ≈5.6 Å). So:

```
d(CoO2-CoO2) = c_hex / 2
```

This gives d ≈ 9.7–9.9 Å for the superconducting hydrate, vs. d ≈ 5.6 Å for non-superconducting anhydrous NaxCoO2 and d ≈ 6 Å for the non-superconducting intermediate (mono-layer) hydrate (c≈12 Å, Foo et al., Solid State Comm. 127, 33 (2003), cited in Milne et al.). This is the formula applied to every `c_hex_A` value in the CSV to produce `coo2_coo2_spacing_A`. (Caveat: Argyriou et al.'s neutron/electron diffraction shows the *true* in-plane cell is doubled, 2a×2a×c — i.e. the a-axis values quoted by Milne/Schaak from simple Rietveld fits to the parent P6₃/mmc cell already correspond to the small (undoubled) hexagonal a; the c-axis, and hence the /2 spacing conversion, is unaffected by the in-plane doubling.)

## 2. What was extracted, and from where

| Source | Type | Rows | What's in the table |
|---|---|---|---|
| Milne, Argyriou, Chemseddine, Aliouane, Veira, Landsgesell, Alber, **PRL 93, 247007 (2004)** (arXiv:cond-mat/0401273, full text obtained) | Primary — Table I | 11 | x (NAA), Co valence (redox titration, most rows), a, c, c/a, Tc (SQUID onset) |
| Schaak, Klimczuk, Foo, Cava, **Nature 424, 527 (2003)** (arXiv:cond-mat/0305450, full text obtained) | Primary — Table 1 | 8 | x (ICP-AES), a, c (Rietveld), Tc (AC susceptibility) — this is the paper that first established the x-dome |
| Argyriou, Milne, Aliouane, Radaelli, Chapon, Chemseddine, Veira, Cox, Mathur, Midgley, arXiv:cond-mat/0403661 (2004) | Primary — Table I | 1 | High-precision NPD Rietveld refinement of one x=0.35, Tc=4.5K sample: a=2.8174 Å, c=19.5456(13) Å |

No values were read off a figure by eye; every number in the CSV comes from a published table with explicit numbers and (where given) uncertainties (uncertainties were dropped from the CSV numeric columns for parseability — see the original PDFs / `notes` field for digit-level precision on request).

### Sources checked but NOT usable for numeric extraction
- **Sakurai, Ihara, Takada, Physica C 514, 378 (2015)** — the requested primary review. Paywalled (ScienceDirect); only the abstract was retrievable, which gives Tc_max = 4.7 K and confirms the line-nodal gap classification but has no table data accessible to this search. Flagged as approximate/abstract-only; not included as CSV rows to avoid inventing precision.
- **Sakurai, Takada, Sasaki, Takayama-Muromachi, J. Phys. Soc. Jpn. 74, 2909 (2005)**, "Phase Diagram of Superconducting NaxCoO2·yH2O" — confirmed to exist and is exactly the "correlation between Tc and c-axis"-type paper the task asked me to verify. It reports that substituting Na+ by H3O+ **at constant Co valence** drives the system from the superconducting phase through a magnetically-ordered phase into a *second, "hidden" superconducting phase* — direct evidence that Co valence (not simple x) is the right abscissa, and that a single-valued Tc(x) or Tc(c) curve is an oversimplification. JPSJ is paywalled; abstract-level claim only, no numeric table extracted.
- **Milne PRL 93, 247007 (2004)** abstract (APS) — 403 Forbidden on direct fetch; full text obtained instead via arXiv, so no loss.
- Mochizuki & Ogata, "CoO2-Layer-Thickness Dependence of Magnetic Properties and Possible Two Different Superconducting States in NaxCoO2·yH2O" (arXiv:cond-mat/0610562, JPSJ 2007) — this is a **theory** paper (multiorbital Hubbard + RPA), not new experimental data; it rationalizes the two-SC-phase picture from Sakurai's JPSJ 2005 result. Included here only as context, not as CSV rows.
- Karppinen, Asako, Motohashi, Yamauchi, "Oxidation State of Cobalt..." Chem. Mater. 16, 1693 (2004) (arXiv:cond-mat/0402484) — gives an independent titration Co valence ≈3.46 for a x≈0.36 sample but no paired Tc/c data in the retrievable excerpt; not included as a CSV row.

## 3. Headline numbers

| Quantity | Range across the 19 structurally-characterized SC samples |
|---|---|
| c_hex (hydrate) | 19.43 – 19.86 Å |
| CoO2–CoO2 spacing (c_hex/2) | 9.72 – 9.93 Å |
| Tc | 2.0 – 4.8 K (0 K just outside this window, in anhydrous and monolayer-hydrate phases) |
| Co valence (titrated, Milne) | 3.24 – 3.45; optimal window 3.24–3.35 |
| Na content x | 0.26 – 0.45 |

Highest Tc in the combined set: **4.8 K**, at x=0.35, Co valence 3.32–3.33, c_hex ≈ 19.66–19.67 Å (Milne et al.). Sakurai's Physica C review abstract separately quotes **Tc = 4.7 K** as the accepted maximum for the compound family.

## 4. The honest caveat on the abscissa (why this matters for overlaying on our Tc(c) panel)

1. **The dome's real x-axis is doping (Co valence / hole count), not spacing.** Both Schaak et al. (2003, plotted vs. Na content x) and Milne et al. (2004, plotted vs. titrated Co valence — the paper this task specifically asked about) draw Tc as a peaked function of a *charge* variable. Milne et al.'s central result is that Schaak's Na-content dome is only approximately right, because Na substoichiometry alone does not fix the charge — oxonium (H3O+) co-intercalation independently dopes holes, so two samples with the same x can have different Co valence and different Tc. The Co-valence dome (their Fig. 3a) is the more physically correct abscissa.
2. **Spacing tracks doping only loosely, and non-monotonically with Tc.** Table I (Milne) shows c/a — and hence c — increasing as Co valence is tuned toward the optimal 3.24–3.35 window, which is the structural fact that motivates "spacing matters" narratives (their explicit claim: "the electronic interlayer coupling between CoO2 sheets becomes more 2D-like" as doping approaches optimal). But within the superconducting samples this is not a clean single-valued Tc(c): e.g. Milne's x=0.29 D2O sample has the *largest* c in the whole set (19.855 Å) with Tc=4.4 K — not the highest Tc — while Schaak's x=0.26 sample has an equally large c (19.77 Å) but Tc=2.4 K, well off the 4.3-4.8 K plateau. c-axis is necessary (SC requires c≳19.4 Å; anhydrous c≈11.2 Å and monolayer-hydrate c≈12 Å are never superconducting above 2K) but not sufficient, and not itself dome-shaped in a well-defined way over the accessible SC range.
3. **The spacing window is narrow by construction.** Because full hydration is essentially an on/off structural transition (the compound either takes up ~1.3-1.4 H2O per formula unit and jumps to c≈19.4-19.9 Å, or it doesn't), there is no continuum of intermediate spacings to scan within the SC phase — unlike, say, a pressure or alloying series that sweeps spacing continuously. The entire experimental Tc variation (0 → 4.8 K) that everyone calls "the dome" happens while c_hex moves by only ≈2% (19.43→19.86 Å out of ~19.6 Å), driven by Na/H3O+/Co-valence chemistry, not by an externally imposed spacing control.
4. **Implication for overlay on our theory's Tc(c) panel.** If our theory predicts Tc as a smooth function of CoO2-CoO2 spacing c over a *wide* range (e.g. spanning the anhydrous 5.6 Å to hydrate 9.7-9.9 Å regime, or beyond), the honest experimental overlay is: (a) two essentially discrete points/clusters — anhydrous (c≈5.6 Å, Tc=0) and hydrate (c≈9.7-9.9 Å, Tc up to 4.7-4.8 K) — rather than a dense curve; and (b) *within* the hydrate cluster, the ~19 data points in the CSV do NOT trace out a clean Tc(c) curve — plotting Tc vs. c_hex directly, without conditioning on Co valence, will look scattered/non-monotonic (see caveat 2), because the real controlling variable (valence/doping) is only partially correlated with c in this dataset. Any Tc(c) panel overlay should either (i) restrict to the narrow c-window and show it as a cloud/error-bar region rather than a fitted curve, or (ii) color/label points by Co valence (or x) to make clear that valence, not c, is what's actually organizing the scatter. Presenting these 19-20 points as if they trace a genuine wide-range Tc(c) dome would misrepresent the source papers.

## 5. Files
- `experimental_dome.csv` — 20 rows: source, sample_id, x_Na, Co_valence_titrated, Co_valence_naive_4minusx (Schaak-style, flagged as superseded by Milne's titration), hydration medium, c_hex_A, coo2_coo2_spacing_A (=c_hex/2), Tc_K, notes, citation.

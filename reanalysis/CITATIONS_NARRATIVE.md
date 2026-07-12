# Citation Research: NaxCoO2·yH2O Mechanism Manuscript

Verified via arXiv/journal pages (July 2026). None of these duplicate the existing bib list.

## 1. Author's prior work — grounding the gallery-band topology in the intro

```bibtex
@article{RadhaLambrecht2021obstructed,
  author  = {Radha, Santosh Kumar and Lambrecht, Walter R. L.},
  title   = {Topological obstructed atomic limit insulators by annihilating {D}irac fermions},
  journal = {Phys. Rev. B},
  volume  = {103},
  pages   = {075435},
  year    = {2021},
  doi     = {10.1103/PhysRevB.103.075435}
}
```
Where it earns its place: intro, one sentence before invoking SciPost(2021) — this is the companion theory paper that establishes the general Zak-phase/obstructed-atomic-limit mechanism (branch cuts in Wannier phase from annihilating Dirac fermions) underlying the SSH4 classification used to justify the Na-derived gallery band as topologically mandated, not incidental.

```bibtex
@article{RadhaLambrecht2021optics,
  author  = {Radha, Santosh Kumar and Lambrecht, Walter R. L. and Cunningham, Brian and
             Gr{\"u}ning, Myrta and Pashov, Dimitar and van Schilfgaarde, Mark},
  title   = {Optical response and band structure of {LiCoO$_2$} including electron-hole interaction effects},
  journal = {Phys. Rev. B},
  volume  = {104},
  pages   = {115120},
  year    = {2021},
  doi     = {10.1103/PhysRevB.104.115120}
}
```
Where it earns its place: methods/intro, to back the QSGW-level electronic-structure treatment of the CoO2/Li(Na)-gallery host material that the SciPost surface-2DEG paper builds on — shows the same group has validated the bulk band structure (including excitonic corrections) that the surface-slab calculation extends.

**Not recommended:** arXiv:2003.00061 is the *preprint* of the already-cited RadhaLambrecht2021 SciPost paper (identical content) — don't cite separately. The 2024 KCoO2 paper (Dernek, Radha, Jackson, Lambrecht, PRB 109, 184442) is a different polymorph (no gallery/2DEG) — skip.

## 2. CRITICAL: D2O isotope-effect experiments — THEY EXIST, ALREADY DONE

**Verdict: found.** Do not present the H→D shift as an untested prediction. Two independent groups measured it in 2003–2004, and both found a **negligible/null Tc shift**:

```bibtex
@article{Jin2003IsotopeD2O,
  author  = {Jin, R. and Sales, B. C. and Khalifah, P. and Mandrus, D.},
  title   = {Observation of bulk superconductivity in {Na$_x$CoO$_2$$\cdot$$y$H$_2$O} and {Na$_x$CoO$_2$$\cdot$$y$D$_2$O} powder and single crystals},
  journal = {Phys. Rev. Lett.},
  volume  = {91},
  pages   = {217001},
  year    = {2003},
  doi     = {10.1103/PhysRevLett.91.217001}
}
```
Exact finding: "The substitution of deuterium for hydrogen has an effect on Tc of less than 0.2 K" — bulk SC with Tc close to 5 K in both H2O and D2O samples (resistivity, susceptibility, specific heat all confirm bulk SC in both). This is a direct H2O↔D2O isotope measurement, not just independent synthesis.

```bibtex
@article{Yokoi2008IsotopeO,
  author  = {Yokoi, Mai and Kobayashi, Yoshiaki and Sato, Masatoshi and Sugai, Shunji},
  title   = {Isotope effect on the superconducting transition temperature of {Na$_x$CoO$_2$$\cdot$$y$H$_2$O}},
  journal = {J. Phys. Soc. Jpn.},
  note    = {arXiv:0803.3254},
  year    = {2008}
}
```
Exact finding (complementary, $^{16}$O$\to^{18}$O not H/D): "we have not found any significant difference of the Tc values between the non-substituted and substituted samples," but the authors caution this "does not necessarily exclude the possibility of the phonon mechanism... because the strong correlation effect suppresses the isotope effect."

**Manuscript-changing consequence:** the claim "deuteration should shift Tc" must be reframed as *already tested and null*. The <0.2 K H/D shift (Jin et al.) is actually *consistent with* a mechanism where water is a dynamical/entropic bath gating Na anharmonicity rather than a mass-coupled phonon participant — the Na double well, not the O–H bond, is the soft pairing mode. Recommend: "the humidity/dynamical fluctuation channel, not O–H(D) vibrational mass, controls Tc — consistent with the near-null isotope shift already reported (Jin et al. 2003; Yokoi et al. 2008), otherwise puzzling for a naive H-bond-phonon picture." Turns a potentially embarrassing null result into supporting evidence.

## 3. Rattling-ion superconductivity precedent

```bibtex
@article{Yonezawa2004KOs2O6,
  author  = {Yonezawa, S. and Muraoka, Y. and Matsushita, Y. and Hiroi, Z.},
  title   = {Superconductivity in a pyrochlore-related oxide {KOs$_2$O$_6$}},
  journal = {J. Phys.: Condens. Matter},
  volume  = {16},
  pages   = {L9--L12},
  year    = {2004},
  doi     = {10.1088/0953-8984/16/3/L01}
}
```
Where it earns its place: discussion, precedent paragraph — first clean report (Tc = 9.6 K) of a rattling alkali-metal ion (K in an oversized Os6O cage) driving strong-coupling s-wave SC; direct structural/conceptual analogue to a loosely-bound Na ion softening in the CoO2 gallery. (Note: the companion JPSJ 73, 1651 (2004) "Unprecedented Superconductivity..." paper was later *retracted* by the authors, JPSJ 74, 3399 (2005) — use this J. Phys.: Condens. Matter report instead, which was not retracted.)

```bibtex
@article{Hiroi2012RattlingReview,
  author  = {Hiroi, Zenji and Yamaura, Jun-ichi and Hattori, Kazumasa},
  title   = {Rattling good superconductor: {$\beta$}-pyrochlore oxides {$A$Os$_2$O$_6$}},
  journal = {J. Phys. Soc. Jpn.},
  volume  = {81},
  pages   = {011012},
  year    = {2012},
  doi     = {10.1143/JPSJ.81.011012}
}
```
Where it earns its place: same paragraph, as the review citation summarizing how rattling anharmonicity and electron-rattler coupling strength (not rattler mass alone) set Tc across K/Rb/Cs — directly supports framing anharmonicity, not simple isotope mass, as the tuning knob, which parallels the H2O/D2O argument in item 2.

## 4. Anharmonic-phonon / quantum-paraelectric pairing theory beyond Edge2015

```bibtex
@article{vanderMarel2019SrTiO3,
  author  = {van der Marel, D. and Barantani, F. and Rischau, C. W.},
  title   = {Possible mechanism for superconductivity in doped {SrTiO$_3$}},
  journal = {Phys. Rev. Research},
  volume  = {1},
  pages   = {013003},
  year    = {2019},
  doi     = {10.1103/PhysRevResearch.1.013003}
}
```
Where it earns its place: theory/discussion section — proposes pairing via exchange of *two* transverse-optical (soft ferroelectric) phonons, quantitatively matched to measured phonon spectral weight; gives a concrete, testable template for how a soft polar/anharmonic mode (directly analogous to the Na double well) can mediate pairing even where a single-phonon linear-coupling channel is weak/vanishes at low density.

```bibtex
@article{Yu2022QuantumParaelectric,
  author  = {Yu, Yue and Hwang, Harold Y. and Raghu, S. and Chung, Suk Bum},
  title   = {Theory of superconductivity in doped quantum paraelectrics},
  journal = {npj Quantum Mater.},
  volume  = {7},
  pages   = {63},
  year    = {2022},
  doi     = {10.1038/s41535-022-00466-2}
}
```
Where it earns its place: same section — shows single soft-TO-phonon pairing reproduces the observed BCS-ratio-preserving dome shape in doped SrTiO3, giving a second, complementary microscopic mechanism to cite alongside van der Marel for "why a soft anharmonic mode near a structural instability can pair electrons in a narrow density/parameter window" — maps onto the manuscript's "narrow spacing window" dome claim.

(Foundational Müller & Burkard 1979 reference for "quantum paraelectric" itself is listed under item 7 below.)

## 5. Water dynamics in confinement, beyond Jalarvo2008

Note: Jalarvo2008 (already in bib) already covers water dynamics *in the hydrate itself* (Jalarvo, Bordallo, Aliouane, Adams, Pieper, Argyriou, J. Phys. Chem. B 112, 703 (2008)) — so the natural complement is a same-technique study in a related layered intercalate, showing the phenomenon (cation-gated, sub-bulk-rate confined water reorientation on ps timescales) generalizes:

```bibtex
@article{Bordallo2008Clay,
  author  = {Bordallo, H. N. and Aldridge, L. P. and Churchman, G. J. and Gates, W. P. and
             Telling, M. T. F. and Kiefer, K. and Fouquet, P. and Seydel, T. and Kimber, S. A. J.},
  title   = {Quasi-elastic neutron scattering studies on clay interlayer-space highlighting
             the effect of the cation in confined water dynamics},
  journal = {J. Phys. Chem. C},
  volume  = {112},
  number  = {36},
  pages   = {13982--13991},
  year    = {2008},
  doi     = {10.1021/jp712305c}
}
```
Where it earns its place: discussion of the water-dynamics mechanism — same lead author (Bordallo) as Jalarvo2008, using the identical QENS methodology on clay interlayers, showing the interlayer cation (not the host lattice) sets the confined-water relaxation timescale — directly supports the "water's dynamical fluctuations, gated by the alkali ion, not water's static presence" argument by showing this is a general phenomenon in cation-intercalated layered galleries, not an idiosyncrasy of the cobaltate.

## 6. NaxCoO2 hydrate phase diagram vs c-axis / Co valence

```bibtex
@article{Sakurai2004Correlation,
  author  = {Sakurai, H. and Takada, K. and Sasaki, T. and Izumi, F. and
             Dilanian, R. A. and Takayama-Muromachi, E.},
  title   = {Correlation between {$T_c$} and lattice parameters of novel superconducting
             {Na$_x$CoO$_2$$\cdot y$H$_2$O}},
  journal = {J. Phys. Soc. Jpn.},
  volume  = {73},
  number  = {9},
  pages   = {2590--2596},
  year    = {2004},
  doi     = {10.1143/JPSJ.73.2590}
}
```
Where it earns its place: this is the single best experimental anchor for the "phonon dome lives in a narrow spacing window" claim — five independently synthesized batches with Tc ranging 3.2–4.6 K, shown to correlate directly with lattice-parameter (not just nominal composition) variation, i.e. Tc tracks the CoO2/Na spacing itself rather than x alone. Pair it with the already-cited Foo2003/Schaak2003 (Na-content dome, x ≈ 0.26–0.35, optimal Co valence ≈ +3.70) to make the two-axis (x-dome × c-axis/hydration-state-dome) argument explicit.

## 7. Additional strengthening references (judgment call, max 3)

```bibtex
@article{OhtomoHwang2004,
  author  = {Ohtomo, A. and Hwang, H. Y.},
  title   = {A high-mobility electron gas at the {LaAlO$_3$/SrTiO$_3$} heterointerface},
  journal = {Nature},
  volume  = {427},
  pages   = {423--426},
  year    = {2004},
  doi     = {10.1038/nature02308}
}
```
Why: the canonical, most-cited precedent for "a polarity/termination discontinuity at an oxide interface creates a 2DEG" — useful in the intro as the general phenomenon class (polar-catastrophe-driven 2DEG) that the manuscript's "opening the gallery creates a 2DEG" claim is a variant of, giving referees an immediately recognizable anchor before the more specialized SSH4/topological argument.

```bibtex
@article{Muller1979QuantumParaelectric,
  author  = {M{\"u}ller, K. A. and Burkard, H.},
  title   = {{SrTiO$_3$}: An intrinsic quantum paraelectric below 4 {K}},
  journal = {Phys. Rev. B},
  volume  = {19},
  pages   = {3593--3602},
  year    = {1979},
  doi     = {10.1103/PhysRevB.19.3593}
}
```
Why: defines "quantum paraelectric" at its source; strengthens item 4's framing and is a one-line, unimpeachable citation for the general concept of a soft, anharmonic, temperature/quantum-fluctuation-suppressed lattice mode — the same conceptual category being claimed for the Na double well.

No third addition recommended — the reference set above (KOs2O6 rattling precedent + SrTiO3 quantum-paraelectric pairing theory + polar-2DEG precedent) already covers the three distinct physical analogies the manuscript leans on; adding more risks diluting rather than strengthening the narrative.

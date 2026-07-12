# Prior computational and experimental literature on water in Na_xCoO2·yH2O

Annotated bibliography compiled to honestly source the "water section" of our DFT+model
paper (claims: (a) interlayer 2DEG in the gallery past a critical CoO2-plane spacing,
(b) Na sits in a soft anharmonic double-well there, (c) water's essential role is
**dynamical** — the solvation cage prevents Na chemisorption onto the CoO2 sheet and
keeps the Na mode soft, rather than water acting as a passive spacer/electrostatic
dielectric). All entries below were checked against arXiv abstract pages, journal
metadata, or (for the Lynn et al. paper) the primary PDF text itself. Where I could not
verify exact wording from a primary source, I say so explicitly rather than paraphrase
a secondary summary as a quote.

Ordering is chronological by first (arXiv) submission date.

---

## 1. The 2003 "rush" — first theoretical attempts to explain the hydrate superconductivity

**G. Baskaran**, "Electronic Model for CoO2 Layer Based Systems: Chiral Resonating
Valence Bond Metal and Superconductivity," *Phys. Rev. Lett.* **91**, 097003 (2003).
arXiv:cond-mat/0303649 (submitted ~24 Mar 2003).
Models the neutral CoO2 layer as a triangular-lattice S=1/2 Mott insulator and
Na_xCoO2·yH2O as an electron-doped Mott insulator (effective t-J model). Proposes a
chiral-RVB metal with d-wave pairing at moderate doping and possible weak
ferromagnetism/p-wave pairing at higher doping. **Water is not modeled at all** —
hydration enters only implicitly as "the doping mechanism that produces the metal";
no water atoms, no lattice dynamics. Generally regarded as the first theoretical paper
responding to the discovery. *Relation to our claim: orthogonal* (purely electronic/
magnetic RVB mechanism; water absent as a physical species).

**B. Kumar and B. S. Shastry**, "Superconductivity in CoO2 Layers and the Resonating
Valence Bond Mean Field Theory of the Triangular Lattice t-J Model," *Phys. Rev. B*
**68**, 104508 (2003). arXiv:cond-mat/0304210 (submitted ~9 Apr 2003).
RVB mean-field theory of the triangular t-J model; finds a time-reversal-symmetry-
breaking complex order parameter with a fully gapped quasiparticle spectrum away from
half filling, and compares an LDA Fermi surface to an effective tight-binding model.
No water atoms in the calculation; water's role is not addressed beyond providing the
observed doping/2D geometry. *Relation to our claim: orthogonal.*

**Q.-H. Wang, D.-H. Lee, and P. A. Lee**, "Doped t-J Model on a Triangular Lattice:
Possible Application to Na_xCoO2·yH2O and Na_{1-x}TiO2," *Phys. Rev. B* **69**, 092504
(2004). arXiv:cond-mat/0304377 (submitted ~16 Apr 2003).
Slave-boson mean-field treatment of the triangular t-J model; finds a time-reversal-
breaking d_{x²-y²}+id_{xy} superconducting ground state. Purely a doped-Mott-insulator
model; water enters only as the historical reason the material is 2D/doped, no
explicit water physics. *Relation to our claim: orthogonal.*

**A. Tanaka and X. Hu**, "Possible Spin Triplet Superconductivity in Na_xCoO2·yH2O,"
*Phys. Rev. Lett.* **91**, 257006 (2003). arXiv:cond-mat/0304409 (submitted ~17 Apr
2003).
Proposes spin-triplet pairing from antiferromagnetic fluctuations on the hexagonal
Fermi surface, drawing an analogy to Sr2RuO4. Model-only (no DFT, no water atoms).
*Relation to our claim: orthogonal.*

---

## 2. First neutron-scattering structure/dynamics study — and the first explicit
##    "water dynamics may matter for pairing" speculation

**J. W. Lynn, Q. Huang, C. M. Brown, V. L. Miller, M. L. Foo, R. E. Schaak, C. Y.
Jones, E. A. Mackey, and R. J. Cava**, "Structure and Dynamics of Superconducting
Na_xCoO2 Hydrate and Its Unhydrated Analog," *Phys. Rev. B* **68**, 214516 (2003).
arXiv:cond-mat/0307263 (submitted 11 Jul 2003, revised 24 Oct 2003).

Neutron diffraction + inelastic neutron scattering on Na_0.3CoO2·1.4(H/D)2O and the
unhydrated parent. Determines the crystal structure (c-axis expands 11.16 Å → 19.5 Å;
water forms an ice-like bilayer between the Na and CoO2 sheets; Na sites shift relative
to the parent structure to accommodate the water layer) and measures the phonon
density of states, finding the lattice dynamical scattering of the superconductor is
**dominated by hydrogen modes** with librational and bending character "quite similar
to ice." Establishes the strong inverse correlation between CoO2-layer thickness
(O–Co–O bond angle) and T_c that our own spacing-vs-T_c analysis builds on.

Verified verbatim quote from the paper's Introduction (checked against the primary PDF,
not a secondary summary):

> "the recent discovery of superconductivity in hydrated Na_xCoO2 — the first
> superconductor where the presence of water is critical to superconductivity — has
> been of particular interest with regard to the superconducting cuprates."

> "On the other hand, the traditional electron-phonon interaction may be establishing
> conventional s-wave pairing, with the possibility that the anharmonic motion of the
> hydrogen and oxygen ions might be playing a role in enhancing the superconducting
> properties, in a manner similar to MgB2."

This is, as far as I can determine, **the earliest explicit statement in the
literature that anharmonic motion of the water/hydrogen sublattice — not just its
static structural presence — might be mechanistically important for the
superconductivity**, drawn directly by analogy to phonon-mediated (MgB2-like)
pairing. It is a *speculative aside in the introduction*, not a calculation: the
paper's own inelastic data are used only to characterize the phonon density of states
and to show the H-dominated modes resemble ice, not to compute any electron-phonon
coupling or to identify a soft anharmonic Na mode. It does **not** discuss Na
chemisorption, a double-well Na potential, or a 2DEG in the gallery.
*Relation to our claim: anticipates* (the general idea that water dynamics — not
water-as-spacer — could matter for the mechanism) *but does not develop it*; the
proposed channel (H/O phonon-mediated conventional s-wave pairing) is a different
physical mechanism from our claim (water dynamically stabilizing a soft anharmonic Na
mode/gallery 2DEG by suppressing Na–CoO2 chemisorption). We must cite this paper
prominently and be precise that we are citing a one-sentence speculative aside, not a
worked-out theory.

---

## 3. First DFT calculations that actually include explicit H2O atoms

**C. A. Marianetti, G. Kotliar, and G. Ceder**, "Role of Hybridization in NaxCoO2 and
the Effect of Hydration," *Phys. Rev. Lett.* **92**, 196405 (2004).
arXiv:cond-mat/0312514 (submitted 19 Dec 2003 — **this appears to be the earliest
arXiv-dated first-principles calculation that explicitly includes H2O atoms in the
structure**, predating Johannes & Singh by about six weeks).

DFT-LDA on Na_{1/3}CoO2 and Na_{1/3}CoO2(H2O)_{4/3}, plus a Hubbard+DMFT model for the
rehybridization physics. Finds that Na doping does not simply add electrons to the Co
t2g manifold — the added electron is "dressed" by O 2p/Co e_g rehybridization — and
separately compares LDA band structures of the hydrated and unhydrated compounds,
concluding explicitly that **"hydration does cause the electronic structure to become
more two-dimensional."** Water's role here is structural/electrostatic — increasing
interlayer spacing and thereby 2D character — with no discussion of water dynamics,
Na double-well potentials, or gallery states. *Relation to our claim: partially
anticipates* the "critical CoO2-plane spacing/2D character" part of claim (a) (water
as the agent that drives 2D-ization by expansion), but is silent on claims (b) and (c)
and treats water purely statically. Should be cited as an early, and possibly the
earliest, DFT paper to explicitly connect hydration-driven spacing to electronic
dimensionality — a precursor to, not an anticipation of, the 2DEG/soft-Na-mode/dynamical-water claims.

**M. D. Johannes and D. J. Singh**, "Comparison of the Electronic Structures of
Hydrated and Unhydrated Na_xCoO_2: The Effect of H_2O," *Phys. Rev. B* **70**, 014507
(2004). arXiv:cond-mat/0401646 (submitted 30 Jan 2004).

DFT electronic-structure calculation for Na_{1/3}CoO2·1.33H2O compared directly to
the parent Na_{0.3}CoO2. Verified abstract (via arXiv):

> "We find that the intercalation of water into the parent compound has little effect
> on the Fermi surface outside of the predictable effects expansion, in particular
> increased two-dimensionality. This implies an intimate connection between the
> electronic properties of the hydrated and unhydrated phases."

This is the paper named as the priority target in the task brief. Its conclusion is
essentially the null result for any *electronic* role of water beyond expansion: no
gallery/interlayer states, no new bands at E_F attributable to water orbitals, and by
implication no support for water acting as an electronic (rather than structural)
actor. Note this is the *opposite* framing from Marianetti/Kotliar/Ceder's more
positive-sounding "hydration does cause... more two-dimensionality" — both groups
essentially agree on the physics (water only matters via expansion/2D-ization) but
Johannes–Singh phrase it as "little effect... outside of predictable expansion
effects" while M/K/C emphasize the 2D consequence more. *Relation to our claim:
contradicts* the idea that water plays any special electronic role at the Fermi
level beyond passive spacing — directly supports our framing that water is NOT doing
interesting electronic-structure work as a static spacer, which is exactly why we
need to argue its role is dynamical instead. This is a load-bearing citation: it is
the strongest existing evidence that the *static* DFT picture of water is a dead end,
motivating the dynamical picture.

---

## 4. Later DFT / model studies of hydration (2004–2006) — mostly doping-only or
##    static-water refinements

**K.-W. Lee, J. Kuneš, and W. E. Pickett**, "Charge Disproportionation and Spin
Ordering Tendencies in Na_xCoO2," *Phys. Rev. B* **70**, 045104 (2004).
arXiv:cond-mat/0403018 (submitted ~1 Mar 2004).
First-principles band structure of **unhydrated** Na_xCoO2 (x = 0.33, 0.61) plus a
Heisenberg-model analysis of charge/spin ordering and Wannier-function hopping
integrals. **Does not include water atoms** — this is a Na-ordering/charge-
disproportionation study of the anhydrous compound, not a hydrate calculation, despite
being frequently cited alongside the hydrate literature. (Pickett's group has several
related papers on Na ordering and charge disproportionation in Na_xCoO2, e.g. a later
Na_0.5CoO2 charge/spin-ordering paper, arXiv:cond-mat/0510555; none of those found
include explicit H2O either.) *Relation to our claim: orthogonal* — useful only for
the anhydrous Na-ordering baseline, not for water physics. **Must not be cited as a
hydrate-with-water DFT study** — verify this distinction carefully since several
secondary sources blur it.

**R. Arita**, "Electronic Structure of Sodium Cobalt Oxide: Comparing Mono- and
Bilayer-Hydrate," *Phys. Rev. B* **71**, 132503 (2005). arXiv:cond-mat/0502256
(submitted 10 Feb 2005).
Full structural optimization + DFT of both the monolayer-hydrate (MLH) and
bilayer-hydrate (BLH) phases (the BLH is the actual superconductor; MLH is not).
Central finding: the a1g band has negligible c-axis dispersion in BLH but ~0.1 eV
dispersion in MLH, i.e., **BLH is much more 2D than MLH**, offered as the likely reason
superconductivity is absent in the MLH phase. This is effectively a "critical spacing"
result — the extra water layer in BLH vs MLH is what pushes the system into the more
2D, superconducting regime. Does not discuss gallery states, an interlayer 2DEG, or Na
double-well/soft modes, and treats water statically (structure only, no dynamics).
*Relation to our claim: anticipates* the "critical spacing" part of claim (a) — this is
probably the closest prior DFT result to our critical-spacing claim — but is silent on
(b) and (c).

**Y. Tanaka, Y. Yanase, and M. Ogata**, "Superconductivity in Na_xCoO2·yH2O due to
Charge Fluctuation," *J. Phys. Soc. Jpn.* **73**, 319 (2004). arXiv:cond-mat/0311266
(submitted ~Nov 2003). Single-band extended Hubbard model (RPA) on the triangular
lattice; finds f-wave triplet pairing near a charge-density-wave instability. Model
only, no water atoms, water enters only as "why the layers are 2D."
*Relation to our claim: orthogonal.*

**Y. Yanase, M. Mochizuki, and M. Ogata**, "Multi-Orbital Analysis on the
Superconductivity in Na_xCoO2·yH2O," *J. Phys. Soc. Jpn.* **74**, 430 (2005).
arXiv:cond-mat/0407563 (submitted Jul 2004). Multi-orbital Hubbard model with
tight-binding parameters fit to DFT band structure; Hund's-coupling-driven
ferromagnetic fluctuations favor triplet pairing on the hole-pocket (a1g) band. This
and the related Mochizuki–Ogata series (ferromagnetic/triplet-pairing instabilities
controlled by trigonal distortion, CoO2-layer-thickness dependence of magnetic
properties, etc., through ~2007) are **model Hamiltonians parameterized from DFT band
structures of the doped material**, not direct DFT+water calculations — water enters
only implicitly through the assumed CoO2-layer geometry/doping level, never as
explicit atoms or dynamics. I could not find a single canonical "Mochizuki & Ogata
review" article in JPSJ matching the brief's description; their contributions are
better characterized as a multi-paper program (2004–2007) rather than one review.
*Relation to our claim: orthogonal* (electronic/magnetic mechanism on a doping-only
model).

**R. J. Xiao, H. X. Yang, and J. Q. Li**, "Influence of Water Intercalation on the
Electronic Structure of the Hydrated Na0.3CoO2·yH2O Using a Local Spin Density
Approximation," *Phys. Rev. B* **73**, 092517 (2006). arXiv:cond-mat/0602629
(submitted 27 Feb 2006).
LSDA calculation comparing parent Na_0.3CoO2 and hydrated Na_0.3CoO2·1.3H2O; reports
hydration-induced changes in valence charge density and O 2p–Co a1g orbital
hybridization (decreased a1g occupation, rehybridization), validated against EELS
data. This is in some tension with Johannes–Singh's "little effect... outside of
expansion" conclusion — Xiao/Yang/Li emphasize real orbital-hybridization changes near
E_F, not just a passive 2D-ization. Static DFT only; no dynamics, no gallery states,
no Na double-well discussion. *Relation to our claim: orthogonal-to-mildly relevant*
— shows the "water has no electronic effect" conclusion of Johannes–Singh is not
universally reproduced, so we should not over-claim consensus on that point; still,
no group in this cluster proposes a dynamical role for water.

---

## 5. Phonon/electron-phonon proposals specific to Na_xCoO2·yH2O (lattice-based
##    pairing mechanisms)

**K. Yada and H. Kontani**, "Electron-Phonon Mechanism for Superconductivity in
Na_{0.35}CoO2: Valence-Band Suhl-Kondo Effect Driven by Shear Phonons,"
arXiv:cond-mat/0512440 (submitted Dec 2005; published as *J. Phys. Soc. Jpn.* series
of papers), and the fuller **K. Yada and H. Kontani**, "s-wave superconductivity due
to Suhl-Kondo mechanism in Na_xCoO2·yH2O: Effect of Coulomb interaction and trigonal
distortion," *Phys. Rev. B* **77**, 184521 (2008), arXiv:0801.3495.
Computes electron-phonon coupling for CoO2-layer **shear and breathing phonons**
(inter-orbital hopping modulation from trigonal CoO6 distortion) and finds Suhl-Kondo-
enhanced s-wave pairing that survives strong Coulomb repulsion, by analogy to MgB2.
This is the clearest fully worked-out *phonon-mediated* pairing proposal specific to
this material family that I found. However — important caveat — **these are CoO2-layer
lattice phonons (Co–O bond distortions), not water/hydrogen modes**; water's ice-like
librational/bending modes (identified by Lynn et al.) are a separate, higher-energy
part of the spectrum and are not the phonons Yada–Kontani couple to electrons. So while
this is "phonon pairing proposed for the hydrate specifically," it is *not* a
water-dynamics proposal — it is orthogonal to our water-dynamics claim even though it
answers the "any phonon-based pairing proposal" sub-question. *Relation to our claim:
orthogonal* (right material, right general idea — anharmonic lattice dynamics
matter — but the wrong sublattice; good to cite as the extant "phonon pairing in this
material" precedent to contrast against, since it targets CoO2 phonons, not the
water/Na gallery dynamics we invoke).

---

## 6. NMR / QENS / µSR studies of Na and water dynamics in the hydrate (experimental)

**N. Jalarvo, H. N. Bordallo, N. Aliouane, M. A. Adams, J. Pieper, and D. N.
Argyriou**, "Dynamics of Water in Na_xCoO2·yH2O," *J. Phys. Chem. B* **112**(3),
703–709 (2008).
Incoherent inelastic + quasielastic neutron scattering (QENS) on Na_0.7CoO2 and
Na_0.28CoO2·1.3H2O to directly probe the hydrogen-bond network dynamics. Identifies
**two distinct proton dynamical processes**: a fast, localized process attributed to
water strongly bound within the sodium-cobalt-oxyhydrate cage, and a slower process
attributed to more collective motion; above 310 K the fast water follows a random
jump-diffusion model with D ≈ 0.9×10⁻⁹ m²/s (well below bulk-water diffusion), and at
room temperature the sodium ions were found to have **no measurable influence** on the
fast-water rotational dynamics. This is a direct, quantitative experimental
characterization of the "water is dynamic, not static" picture and is the single
most important experimental water-dynamics reference for our claim. It does not
propose any link to superconductivity or to a Na double-well/2DEG — it is a pure
dynamics-characterization paper. *Relation to our claim: anticipates* (establishes
empirically that water in the gallery is dynamically active on relevant timescales,
which is a necessary experimental precondition for our "water = lubricant" claim) but
does **not** connect this to the Na potential, chemisorption, or superconductivity —
we add that connection.

Complementary/context: J. W. Lynn et al. (Section 2 above) separately report, in
follow-up high-resolution backscattering work referenced on the NCNR summary pages,
that "most of the water is static on an energy scale of 1 meV" — i.e., the ice-like
bulk of the water network is largely rigid on sub-meV (long-time) scales even though
Jalarvo et al. resolve faster (meV–tens-of-meV) local dynamics. These two results are
not contradictory (different energy/time windows) but should both be cited so we do
not overstate "water is dynamic" without the important caveat that the bulk hydrogen-
bonded ice-like network is largely static; the dynamically relevant species for our
argument would need to be the more mobile component (or the Na motion itself, which is
the object of the NMR studies below), not the bulk ice lattice.

**H. Ohta, Y. Itoh, C. Michioka, and K. Yoshimura**, "23Na NMR study of
non-superconducting double-layer hydrate Na_xCoO2·yH2O," *Physica C: Superconductivity*
445–448, 69–72 (2006). arXiv:cond-mat/0605487 (submitted 19 May 2006).
23Na NMR on hydrated vs. dehydrated samples. Finds the hyperfine field and electric
field gradient at the Na site are **markedly reduced by hydration** — i.e., water
electronically screens/shields the Na ion — while the spin-lattice relaxation
behavior otherwise resembles the non-hydrated compound, arguing the SC-relevant
physics is magnetic-fluctuation-driven rather than water-dynamics-driven in their
framing. *Relation to our claim: partially anticipates* the "water screens/shields Na"
piece of claim (b)/(c) (consistent with water suppressing a strong Na–CoO2 Coulomb/
chemisorption interaction) but the authors' own interpretive framing points away from
water dynamics and toward magnetic quantum-criticality — should be cited carefully,
noting we are using their *data* (Na hyperfine screening by water) but not endorsing
their magnetic-QCP framing as the whole story.

**Y. Itoh, H. Ohta, C. Michioka, M. Kato, and K. Yoshimura**, "59Co, 23Na, and 1H NMR
Studies of Double-Layer Hydrated Superconductors Na_xCoO2·yH2O," in *Advances in Solid
State Physics* Vol. 47, pp. 329–341 (Springer, 2008); also arXiv:0803.0366 (submitted
4 Mar 2008). Proceedings review of the group's NMR program; high-resolution 1H NMR
identifies the presence of H3O⁺ oxonium ions (not pure H2O) in the gallery, and 23Na
T1⁻¹ measurements track local field fluctuations correlated with T_c. Useful as a
compact summary/review of the NMR evidence for gallery composition (Na⁺/H3O⁺/H2O
coexistence) but again frames the physics magnetically rather than dynamically.
*Relation to our claim: orthogonal-to-supporting* (useful structural/compositional
context — H3O⁺ presence matters for our "solvation cage" picture — but not itself a
dynamical-water argument).

**W. Higemoto et al.**, "Muon Spin Relaxation Measurements in Na_xCoO2·yH2O,"
arXiv:cond-mat/0311427 (2003), and related µSR papers (e.g. Physica C 2005,
ScienceDirect S0921452605012949) find zero-field µSR relaxation rate temperature-
independent down to 2 K (no static magnetism in the SC phase; weak magnetism seen only
in an off-stoichiometric, non-SC, water-excess sample). These probe local magnetic
(not structural/vibrational) dynamics and provide evidence against broken time-reversal
symmetry in the SC state. *Relation to our claim: orthogonal* (magnetic probe, not
water-dynamics probe), included for completeness since µSR was explicitly named as a
technique to check.

---

## Citation obligations — the papers we MUST cite in the water section

1. **Lynn et al., PRB 68, 214516 (2003)** [cond-mat/0307263] — because it is the
   earliest paper to explicitly float the idea that anharmonic water/hydrogen motion
   (not static water) could matter for the pairing mechanism, by direct analogy to
   MgB2. We must cite it and be precise that we are extending a one-sentence aside into
   a worked-out mechanism, not rediscovering something already established.

2. **Johannes & Singh, PRB 70, 014507 (2004)** [cond-mat/0401646] — the canonical
   "water has little electronic effect beyond expansion" result. This is a load-bearing
   negative result for us: it is the strongest existing DFT evidence that a *static*
   picture of water fails to explain anything interesting, which is precisely the
   motivation for our dynamical reframing. Must be cited accurately (not overstated —
   they do find increased 2D character, they just attribute it to trivial expansion).

3. **Marianetti, Kotliar & Ceder, PRL 92, 196405 (2004)** [cond-mat/0312514] — the
   earliest identified water-inclusive DFT calculation (Dec 2003), and the first to
   explicitly state hydration drives 2D-ization of the electronic structure. Must be
   cited both for priority (it predates Johannes–Singh) and because it is the closest
   existing anticipation of our "critical spacing" claim (a).

4. **Arita, PRB 71, 132503 (2005)** [cond-mat/0502256] — direct DFT comparison of
   monolayer- vs. bilayer-hydrate showing the extra water layer is what pushes the a1g
   band into the 2D (dispersionless along c) regime associated with superconductivity.
   This is the closest prior computational statement of a "critical spacing" effect
   and must be cited and clearly distinguished from our 2DEG/soft-Na-mode extension of
   it.

5. **Jalarvo et al., J. Phys. Chem. B 112, 703 (2008)** — the only paper found that
   directly and quantitatively measures water dynamics (QENS) in the gallery,
   establishing the empirical fact that gallery water is dynamically active (fast
   jump-diffusion component) rather than a rigid spacer. This is the strongest
   experimental anchor for the "water = lubricant, not spacer" claim and must be cited
   as the empirical basis, paired with the caveat (from the Lynn et al. follow-up
   backscattering work) that the bulk ice-like network is largely static on sub-meV
   scales — so we should be precise about which component of the water dynamics we are
   invoking.

(Honorable mention, recommended but not mandatory: Xiao, Yang & Li, PRB 73, 092517
(2006), to show the "water has no electronic effect" conclusion is contested in the
static-DFT literature, and Ohta et al., Physica C 445-448, 69 (2006), for the NMR
evidence that water electronically screens the Na site, which supports the
chemisorption-suppression mechanism.)

---

## Notes on searches performed / limitations

- Search coverage was thorough for 2003–2010 for the specific author names and topics
  named in the task brief (Johannes & Singh; K.-W. Lee & Pickett; Mochizuki & Ogata;
  Arita/Kuroki; Lechermann/Piefke; Baskaran; Kumar & Shastry; Tanaka & Hu; Wang/Lee/Lee;
  Ogata's group; Yada & Kontani), plus NMR/QENS/µSR water-dynamics literature.
- I did **not** find a single canonical "Mochizuki & Ogata review" article matching the
  brief's description of a JPSJ review; their program is a multi-paper series
  (2004–2007+) on multi-orbital Hubbard-model triplet pairing, none of which include
  explicit water atoms in a DFT sense.
- I did **not** find a Lechermann/Piefke paper specifically treating the *hydrate* with
  explicit water in a DFT+DMFT calculation; Lechermann's group's DFT+DMFT work located
  (Boehnke & Lechermann; Lechermann/Boehnke/Grieger/Piefke) treats anhydrous Na_xCoO2
  spectral functions and thermopower, not the hydrate. This should be treated as a gap
  rather than a citation — if such a paper exists it was not surfaced by these
  searches, and I flag this rather than fabricate a reference.
- Some search-tool outputs (via the WebSearch summarizer) paraphrase content and could
  contain small inaccuracies; wherever a claim in this document is presented as a
  "verified quote," it was cross-checked against an arXiv abstract page, journal
  metadata, or (for the Lynn et al. paper) the actual PDF text extracted via a direct
  file read — not taken solely from a search-engine summary. Abstracts for Johannes &
  Singh, Marianetti/Kotliar/Ceder, Xiao/Yang/Li, and Arita were fetched from their
  arXiv abstract pages directly. All arXiv IDs, journal volumes, and page numbers above
  should still be spot-checked against the publisher/arXiv record before being placed
  in the final manuscript's bibliography, standard practice for any citation.

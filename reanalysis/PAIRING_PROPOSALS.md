# Proposed pairing states for Na_xCoO2·yH2O and their experimental status

Compiled from the verified literature research in this project (see also PRIOR_DFT_WATER.md).
Source for the "critical assessment": Mazin & Johannes, Nature Physics 1, 91 (2005) — eliminated
all symmetry-allowed states except two f-wave candidates, using the constraints below.

## Experimental constraints (the gauntlet)
- Knight shift (single-crystal 59Co NMR): decreases below Tc both field orientations → singlet
  (arXiv:cond-mat/0306264; contradicted by early powder NMR claiming no drop → triplet claims).
- 1/T1: no Hebel-Slichter coherence peak → against conventional isotropic s-wave.
- Specific heat: T² term (line nodes) in some samples (Yang et al., cond-mat/0308031); other
  samples fit two-gap / extended-s; explicit sample-age dependence (cond-mat/0503690).
- μSR (Higemoto et al., PRB 70, 134508 (2004)): no TRS breaking, no static magnetism → rules out
  chiral TRSB states (d+id, chiral p) in their simplest form.
- SC survives nonmagnetic disorder to a degree unusual for sign-changing states (sample dependence
  muddies this).

## Proposals (chronological)
| State / mechanism | Proposers | Status vs constraints |
|---|---|---|
| RVB singlet (t-J, d+id family) | Baskaran PRL 91, 097003 (2003); Kumar & Shastry; Wang, Lee & Lee | d+id breaks TRS → tension with μSR; singlet OK with Knight shift |
| Spin-triplet p/f-wave (ferromagnetic fluctuations, e_g' pockets) | Tanaka & Hu; Kuroki, Tanaka, Arita | relies on e_g' pockets ARPES never sees; triplet contradicts single-crystal Knight shift |
| Chiral d+id (multiorbital fRG) | Kiesel, Platt, Hanke, Thomale PRL 110 (2013) | TRSB → μSR tension; post-dates and contradicts Mazin-Johannes elimination |
| f-wave (two candidates surviving elimination) | Mazin & Johannes, Nat. Phys. 1, 91 (2005) | consistent with no-HS-peak + nodes; triplet f needs Knight-shift reconciliation |
| Conventional s-wave el-ph (H/O anharmonic phonons, MgB2-like aside) | Lynn et al. PRB 68, 214516 (2003) (one-sentence speculation) | no HS peak + nodal specific heat argue against simple s |
| Magnetic-fluctuation pairing near magnetic phase | Sakurai, Ihara, Takada, Physica C 514, 378 (2015) | phase diagram support (magnetic phase adjacent); mechanism not microscopically derived |

## Where our mechanism sits
- Gallery-2DEG + soft Na cage mode, phonon-mediated: natural ground state is singlet s±-like on the
  2DEG band (anisotropic two-gap with the a1g band), TRS-preserving → passes μSR and single-crystal
  Knight shift.
- No Hebel-Slichter peak: strong-coupling/anharmonic broadening of the soft mode suppresses the
  coherence peak (standard for soft anharmonic pairing modes) — state plainly, not oversold.
- "Line nodes" in specific heat: small gap on the 2DEG band easily suppressed by gallery disorder
  (dehydration!) mimics nodes; explains the sample-age dependence of gap structure directly.
- Near the Stoner boundary of the polarized 2DEG strip (see spin_analysis/), an equal-spin f-wave
  channel exists — connecting to the Mazin-Johannes survivors; sample proximity to the boundary
  explains the historical singlet/triplet contradiction.
- Falsifier: observation of TRSB, or a Knight shift that stays flat in a fresh single crystal,
  would exclude our singlet channel.

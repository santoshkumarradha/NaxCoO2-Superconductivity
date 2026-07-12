# STORY_NOTES.md — how Feynman would structure this paper

Physics to be told: (i) DFT shows that opening the gallery between CoO2 layers past a
critical spacing *simultaneously* creates an interlayer 2D electron gas and a soft Na
double-well mode; (ii) gas = carriers, soft mode = glue, superconductivity only in the
window between "no gas" and "Na frozen"; (iii) the one experimental mismatch (bilayer
hydrate at 9.9 Å lies outside the static window) becomes the prediction that water's job
is to re-soften the Na well — dynamical lubricant, not spacer; (iv) the picture explains
why Li never superconducts, why the monolayer hydrate fails, and the magnetic phase at the
carrier turn-on edge.

Style basis: `FEYNMAN_STYLE.md` (register: written Feynman, 1948, not lecture Feynman).

## 1. The first paragraph

Feynman opens with the *fact that wants explaining*, not with the field ("It is a curious
historical fact that..." [RMP48]). Here the curious fact is a gift: **this superconductor
only works wet.** Paragraph one states that fact flatly, states the claim of the paper in
one sentence, and — the 1948 move — immediately deflates: the calculations are ordinary
DFT, nothing exotic; what is new is only the point of view (the action is *between* the
layers, not in them). By the end of paragraph one the reader knows the whole mechanism in
words: gallery opens → gas condenses and sodium starts to rattle → gas carries, rattle
glues → superconductivity lives in a window.

## 2. Where the mismatch goes: early, and honestly

Feynman puts the defect in the main text with its own heading ("11. Inadequacies of the
Formulation" [RMP48]; the thesis conclusion is a list of failures). But here the mismatch
is better than a defect — it is the *hinge of the argument* — so it goes even earlier than
Feynman's inadequacies sections: **flagged in the introduction, resolved in its own section
in the middle of the paper.** The 9.9 Å bilayer hydrate sits outside the static window; a
lesser paper hides that in a supplement. This paper says, in the introduction: "One
measured number does not fit the static picture, and it is the most important number in
the paper." Resolving it is what *forces* the reinterpretation of water as a dynamical
lubricant — the paper's real discovery. Feynman's Nobel lecture is one long demonstration
that the wrong idea, honestly tracked, is the story ("...to describe how you had the wrong
idea first" is exactly what journals don't allow; do it anyway).

## 3. Section plan (Feynman logic: picture → calculation → verdict → difficulty → checks)

1. **Introduction.** The wet fact; the claim; the deflation; the flagged mismatch; the
   roadmap in verbs ("We first show... We then ask... Finally we check the picture against
   the three facts any mechanism must explain.").
2. **What happens when the gallery opens.** The static DFT result, told as an event: pull
   the layers apart, and at a critical spacing two things happen *at once* — the oxygen
   charges bifurcate (an electron gas condenses on the interior surfaces) and the sodium
   potential splits into a double well. One event, two faces. Concrete picture first
   (rattling ion, gas condensing on an internal surface), band structures and Bader
   charges after. Equations introduced in words ("the sodium sits in a potential which is,
   to an excellent approximation, a pair of wells separated by a barrier of a few meV...").
3. **The pairing window.** Carriers from the gas, glue from the soft mode; T_c estimated
   from the electron-phonon calculation. The window stated as two walls: close the gallery
   and the gas evaporates (no carriers); open it too far, or cool the sodium into one
   well, and the glue sets (no pairing). *The money plot lives here*: T_c (or λ, the
   coupling) versus layer spacing, a dome with both walls marked, experimental compounds
   pinned on it. Every later section points back at this one figure — Feynman's "sort of
   bird's-eye view" [RMP48 §13] in plot form.
4. **The difficulty with the bilayer hydrate.** Its own section, Feynman-plain heading.
   The 9.9 Å point does not sit where the static dome says it should. One cause, named:
   the calculation holds the water still, and the water does not hold still. The water's
   hydrogen-bond network re-softens the sodium well at spacings where the bare double well
   would have frozen — the water is not a spacer, it is a lubricant. This converts the
   mismatch into the sharpest prediction of the paper (isotope/deuteration shift, pressure
   dependence, water-dynamics correlation with T_c).
5. **Three checks.** The picture must explain, and does: (a) *Li never superconducts* —
   lithium holds its charge and its well never softens in the accessible range;
   (b) *the monolayer hydrate fails* — gallery open enough for neither gas nor soft mode
   in the right regime (one wall of the dome); (c) *the magnetic phase at the carrier
   turn-on edge* — right at the wall where the gas first condenses, the carriers are few
   and the gas orders instead of pairing. Each check is one paragraph: fact, then why the
   picture requires it.
6. **Inadequacies and predictions.** In the main text, cheerfully, thesis-conclusion
   style: what the static calculation cannot do; which numbers are soft (the barrier
   height, the anharmonic frequency); "The final test of any physical theory lies, of
   course, in experiment" — then the ranked list of experiments, most lethal first.

## 4. Candidate titles (his register: plain nominal phrase, no colon)

Feynman's actual titles are flatly descriptive: "Space-Time Approach to Non-Relativistic
Quantum Mechanics", "Space-Time Approach to Quantum Electrodynamics", "Simulating Physics
with Computers", "The Principle of Least Action in Quantum Mechanics".

- **The Role of Water in the Superconductivity of Sodium Cobaltate**
- **Interlayer Approach to Superconductivity in Hydrated Sodium Cobaltate** (deliberate
  echo of "Space-Time Approach to...")
- **An Electron Gas and a Rattling Ion Between the Layers of NaxCoO2**
- **Why Sodium Cobaltate Superconducts Only When Wet** (bolder; Feynman-adjacent rather
  than Feynman — he asked questions like this in talks, not in Phys. Rev. titles)
- **Superconductivity from a Soft Sodium Well in Hydrated NaxCoO2**

Avoid: anything with "novel", "emergent", "unconventional", a colon, or "insights".

## 5. Opening paragraph, drafted in the extracted style

> It is a curious fact that sodium cobaltate becomes a superconductor only when it is wet.
> Dry Na_xCoO2 is an ordinary metal down to the lowest temperatures measured; slip two
> layers of water between its CoO2 sheets and it superconducts near 4.5 K; and twenty
> years of looking for the reason inside the CoO2 sheet itself have not produced an
> accepted mechanism. This paper takes the water seriously as a mechanical object and asks
> what it does. We shall show, by ordinary density-functional calculations, that pulling
> the layers apart past a critical separation does two things at once: the charge on the
> oxygen layers bifurcates, so that an electron gas condenses in the gallery between the
> sheets, and the sodium ion, which had been sitting at the bottom of one stiff well,
> finds itself with two shallow wells and begins to rattle between them. The gas supplies
> the carriers and the rattle supplies the glue, and superconductivity can live only in
> the window between the layers too close, where there is no gas, and the sodium frozen,
> where there is no glue. There is nothing exotic in the calculations, and the ingredients
> — a two-dimensional electron gas, an anharmonic mode — are old ones; what is new is only
> where we look for them, which is between the layers rather than in them. One measured
> number does not fit this static picture: the superconducting bilayer hydrate holds its
> layers 9.9 Å apart, wider than the window we calculate, and we shall meet this
> difficulty early, in Section IV, because resolving it is what teaches us the actual role
> of the water. The water is not a spacer. It is a lubricant, and that statement can be
> tested.

(Word-level debts, all to retrieved text: "It is a curious fact" [RMP48 opening]; "We
shall show" [RMP48 roadmap voice]; the deflation move "There are, therefore, no
fundamentally new results... a pleasure in recognizing old things from a new point of
view" [RMP48]; the flagged-difficulty move from "11. Inadequacies of the Formulation"
[RMP48] and the thesis Conclusion; the two short verdict sentences at the end follow the
rhythm of "In classical mechanics it is always true. In quantum mechanics it is often
false." [RMP48].)

## 6. Running reminders while drafting

- The dome figure is the paper. Every section should end by placing itself on that figure.
- The mismatch is an asset; never apologize for it, spend it.
- Enthusiasm goes to nature (the coincidence of gas and soft mode is *nature's* trick, "one
  event seen two ways"), never to the method or the authors.
- One hedge per claim, at the uncertain clause exactly.
- If a paragraph could open an AI survey of the field, delete it.

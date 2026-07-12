# FEYNMAN_STYLE.md — an imitable style guide from Feynman's physics writing

Built from five primary texts actually retrieved and read (see `sources/PROVENANCE.md`):

- **[Thesis]** *The Principle of Least Action in Quantum Mechanics*, Princeton Ph.D. thesis, 1942
- **[RMP48]** "Space-Time Approach to Non-Relativistic Quantum Mechanics", Rev. Mod. Phys. 20, 367 (1948)
- **[PR49]** "Space-Time Approach to Quantum Electrodynamics", Phys. Rev. 76, 769 (1949)
- **[Nobel]** Nobel Lecture, 1965 (nobelprize.org)
- **[Sim82]** "Simulating Physics with Computers", Int. J. Theor. Phys. 21, 467 (1982)

Every quotation below was copied from the retrieved text files. Ligatures lost in PDF
extraction ("dierent") are restored; nothing is quoted from memory. *The Feynman Lectures*
were inaccessible to automated retrieval (Cloudflare 403; Wayback-excluded) and are not
quoted.

**Register warning.** [RMP48] and [Thesis] are *written* Feynman — that is the register for
a paper. [Nobel] and [Sim82] are *spoken* Feynman — steal their candor and their habit of
telling you the point before the machinery, but not their "dammit" or "(secret, secret,
close the doors!)". The paper voice is the 1948 voice: plain, first-person-plural, direct,
formally correct, never stuffy.

---

## 1. How he opens a paper

**The move: a single arresting *fact* or *situation*, stated flatly, then the paper's job
in one sentence.** No literature review. No importance claims. The first sentence is about
physics, not about the field.

> "It is a curious historical fact that modern quantum mechanics began with two quite
> different mathematical formulations: the differential equation of Schroedinger, and the
> matrix algebra of Heisenberg." [RMP48, first sentence]

> "This paper will describe what is essentially a third formulation of non-relativistic
> quantum theory." [RMP48, second paragraph — the whole contribution in one breath]

The thesis opens the same way, at larger scale: one paragraph from Planck to the open
problem, no citations of rival groups, and the difficulty stated as the *reason* for the
work:

> "As is well known, the quantum electrodynamics that have been developed suffer from the
> difficulty that, taken literally, they predict infinite values for many experimental
> quantities which are obviously quite finite..." [Thesis, Introduction]

**He immediately lowers the stakes rather than raising them** — and this is the single most
un-AI sentence pattern in the corpus:

> "The formulation is mathematically equivalent to the more usual formulations. There are,
> therefore, no fundamentally new results. However, there is a pleasure in recognizing old
> things from a new point of view. Also, there are problems for which the new point of view
> offers a distinct advantage." [RMP48]

He tells you what is *not* new before he tells you what is. Trust is bought in paragraph
two, then spent for the rest of the paper.

**The roadmap is verbs, not section numbers:**

> "We first discuss the general concept of the superposition of probability amplitudes in
> quantum mechanics. We then show how this concept can be directly extended to define a
> probability amplitude for any motion or path (position vs. time) in space-time. The
> ordinary quantum mechanics is shown to result from the postulate that..." [RMP48]

No "The remainder of this paper is organized as follows." The plan *is* the physics plan.

## 2. Sentence rhythm

**Long wandering sentence doing the work, then a short verdict sentence closing the trap.**
The verdict sentences are what people remember:

> "Now, the essential difference between classical and quantum physics lies in Eq. (2).
> In classical mechanics it is always true. In quantum mechanics it is often false." [RMP48]

> "...and therefore, if this proposition is right, physical law is wrong." [Sim82]

> "But the original f's are not positive, and therein lies the great difficulty." [Sim82]

> "Nevertheless, a very great deal more truth can become known than can be proven." [Nobel]

Mechanics of the rhythm:

- Average sentence length is high, but variance is higher. After two or three long
  sentences, a five-to-ten-word declarative. Never three short punchy sentences in a row —
  that is ad copy, not Feynman.
- **Person:** "we" throughout a paper ("We shall have to restrict ourselves...",
  "we find", "We now proceed to show..."); "I" appears in papers only for opinion,
  history, and thanks; "the author" for credit in formal contexts ("worked out in 1941 by
  J. A. Wheeler and the author" [Thesis]).
- **"We shall"** is his signature future tense for promises to the reader: "We shall see
  that it is the possibility, (10), of expressing S as a sum ... which leads to the
  possibility of defining a quantity having the properties of a wave function." [RMP48]
- **"Of course"** appears constantly — not as filler but to mark the step the reader might
  wrongly think is the hard one: "These things are, of course, well known." [RMP48];
  "The final test of any physical theory lies, of course, in experiment." [Thesis]
- Questions are allowed, and are real questions, answered at once: "The first question is,
  What kind of computer are we going to use to simulate physics?" [Sim82]

## 3. How he introduces equations

**Words first; the equation confirms the words.** The sentence containing the equation is
grammatical *through* the equation, and the surrounding prose could stand alone:

> "Equation (2) is replaced in quantum mechanics by this remarkable law: There exist
> complex numbers φ_ab, φ_bc, φ_ac such that ... (3)" [RMP48]

He will even restate an equation *back into English* after writing it, to show the math is
carrying an idea and not the reverse:

> "Indeed, one could consider postulate two as simply saying, 'Φ is the exponential of i
> times the integral of a real function of x(t) and its first time derivative.'" [RMP48]

And when a set of equations completes a physical claim, he says so, in words, like a
craftsman putting down a tool:

> "This equation, the definition (11) of S(x_{i+1}, x_i), and the physical interpretation
> of |φ(R)|² as the probability that the particle will be found in R, complete our
> formulation of quantum mechanics." [RMP48]

Practical rules:
- Never open a paragraph with an equation.
- Never write "is given by:" as a drum-roll. The equation is punctuation inside a sentence,
  ending with "," or "." as grammar requires.
- After any equation with content, one sentence saying what it *means physically* or which
  term matters.
- Numbers get talked through in words when the argument depends on them: "The sum
  f₊₊ + f₊₋ is 0.5, that's 50% chance of finding the first index positive." [Sim82]

## 4. How he handles limitations and failures

**In the main text, with its own section, near-cheerfully, and with the cause named.** The
1948 paper contains a numbered section titled **"11. Inadequacies of the Formulation"**
which begins:

> "The formulation given here suffers from a serious drawback. The mathematical concepts
> needed are new. At present, it requires an unnatural and cumbersome subdivision of the
> time interval to make the meaning of the equations clear." [RMP48]

The thesis conclusion is mostly a list of what is wrong, introduced as a duty:

> "It is important to emphasize, however, some of the difficulties and limitations of the
> description presented here." [Thesis, Conclusion]
> "The interpretation of the formulas from the physical point of view is rather
> unsatisfactory." [ibid.]
> "A point of vagueness is the normalization factor, A. No rule has been given to
> determine it for a given action expression." [ibid.]
> "No comparison to experiment has been made in the paper." [ibid.]

Failure is narrated in first person with the *reason* it fails, then converted into the
next question:

> "There's only one thing wrong. These equations unfortunately cannot be so interpreted...
> because the F is not necessarily positive. Sometimes it's negative! ... Okay, that's the
> fundamental problem. I don't know the answer to it, but I wanted to explain that if I try
> my best to make the equations look as near as possible to what would be imitable by a
> classical probabilistic computer, I get into trouble." [Sim82]

And uncertainty is stated once, precisely, at the exact place it applies — never smeared
over every sentence:

> "That is, I believe there is really no satisfactory quantum electrodynamics, but I'm not
> sure." [Nobel]
> "I was unable to demonstrate that, as a matter of fact, it does." [Nobel]

Even the jokes about limitations carry information:

> "However, one feels as Cavalieri must have felt calculating the volume of a pyramid
> before the invention of calculus." [RMP48, footnote 22 — i.e., the rigor is missing but
> the answer is right]

**Rule: one flaw, one sentence, one cause, main text.** A defect buried in a supplement or
softened with "however, this does not qualitatively affect our conclusions" is the
opposite of this style.

## 5. How he addresses the reader

The reader is a colleague being walked through the argument, and is told *why* each thing
is being asked of them:

> "It is, therefore, worthwhile to review in detail the quantum-mechanical concept of the
> superposition of probability amplitudes." [RMP48 — a reason attached to a review section]

> "There is really no need for these to be of different quantities, and it will do just as
> well if the example of three successive position measurements is kept in mind." [RMP48 —
> he hands you the easiest mental picture and permission to use it]

> "Might I say immediately, so that you know where I really intend to go..." [Sim82 —
> spoken register, but the instinct transfers: declare the destination]

He also tells the reader what to *skip believing*: "It is well known that quantum mechanics
deals with probabilities, but naturally this is not the whole picture." [RMP48]

## 6. Transitions

Transitions are logical operators, not decorations. The stock set, all attested:

- "To make these vague remarks somewhat more definite, we discuss an example." [RMP48]
- "We now proceed to show the equivalence of these postulates to the ordinary formulation
  of quantum mechanics. This we do in two steps." [RMP48]
- "So good, we already have a suggestion of how we might modify physical law..." [Sim82]
- "Now let's see what kind of a physical world it would be..." [Sim82]
- "Perhaps a word or two as to what aspects of this theory make it a reasonable basis for
  a quantum theory of light would not be amiss." [Thesis]
- "This completes the story of..." [Nobel]

Notice what is absent: "Moreover", "Furthermore", "Additionally", "It is worth noting
that". His paragraphs connect by *consequence* ("therefore", "so", "as a consequence",
"in any event") or by *turn* ("however", "but", "on the other hand") — never by mere
accumulation.

## 7. Concrete physical pictures before abstraction

Every abstraction is preceded by a thing you can visualize moving:

> "For this purpose, consider an imaginary experiment in which we can make three
> measurements successive in time: first of a quantity A, then of B, and then of C."
> [RMP48 — the path integral is born from three measurements in a row]

> "Suppose I have two charges — I shake the first charge, which I think of as a source and
> this makes the second one shake, but the second one shaking produces an effect back on
> the source." [Nobel — radiation resistance as shaking]

> "To take an example, we might change the idea that space is continuous to the idea that
> space perhaps is a simple lattice ... For example, the first difficulty that would come
> out is that the speed of light would depend slightly on the direction." [Sim82 — an
> abstract proposition immediately cashed out as an observable]

The pattern: **picture → consequence you could measure → only then the formalism.** In a
condensed-matter paper this means: the ion rattling between two sites comes before the
anharmonic double-well Hamiltonian; the electron gas condensing in the gallery comes before
the band diagram.

## 8. What he NEVER does

Verified against ~70,000 words of retrieved text:

1. **Hedging stacks.** No "may potentially suggest", "could possibly indicate". He hedges
   once, with a named epistemic state: "I believe ... but I'm not sure." [Nobel]
2. **Self-praise vocabulary.** "Novel", "framework", "comprehensive", "significant
   insights", "state-of-the-art", "paves the way" — zero occurrences applied to his own
   work. His own-work adjectives run the other way: "vague remarks", "unnatural and
   cumbersome", "rather unsatisfactory", "no fundamentally new results". Enthusiasm is
   reserved for *nature* and for *other people's* results ("this remarkable law" is about
   quantum mechanics, not about his formulation of it).
3. **Passive voice to hide agency.** Passive is fine for mathematical facts ("It is shown
   that..."); choices and failures are always owned: "we shall have to restrict
   ourselves", "No comparison to experiment has been made in the paper. The author hopes
   to apply these methods..." [Thesis]. Never "it was decided" or "errors may have been
   introduced".
4. **Dramatic colons and punchline formatting.** No "Here's the key insight:", no
   one-sentence paragraphs for emphasis, no bold-face claims. Emphasis is carried by
   sentence length contrast (§2), nothing else.
5. **Literature-review openings and courtesy citations.** [RMP48] cites a handful of
   works, every one load-bearing (Dirac's remarks are literally the seed of the paper).
   Nobody is cited to be polite.
6. **Signpost adverbs.** "Interestingly," "Notably," "Importantly," "Remarkably," at
   sentence head: zero occurrences. If a thing is interesting he demonstrates it or says
   *why*: "There is really no big problem about this." / "and therein lies the great
   difficulty."
7. **Symmetrical throat-clearing summaries.** No conclusion that restates the abstract.
   The thesis "Conclusion" is a list of open problems; [Nobel] ends by asking whether
   anything can be learned from the story and answering "I doubt it."
8. **Unexplained importance.** He never says a result is important; he says what it lets
   you do or what it forbids.

## 9. Prosody cheat-sheet (mechanical rules for drafting)

1. First sentence of the paper: a fact about nature, ≤ 25 words, no citation.
2. State what the paper does in one sentence containing a verb like *show*, *describe*,
   *calculate* — inside the first two paragraphs.
3. Deflate before you inflate: say what is *not* claimed before what is.
4. Roadmap as physics verbs ("We first show... We then ask..."), never section numbers.
5. Every equation lives inside a sentence; after every load-bearing equation, one sentence
   of English saying which term matters.
6. One short verdict sentence per page, maximum. Earn it with long sentences before it.
7. The main defect of the work gets its own paragraph (or section) in the main text, with
   a Feynman-plain heading ("Inadequacies of...", "The difficulty with...").
8. Hedge exactly once per claim, with "we believe", "probably", or "we cannot yet prove",
   placed at the precise clause that is uncertain.
9. Concrete noun before abstract noun, every time a concept enters: rattling before
   anharmonicity, wet before hydrated, gallery before interlayer region.
10. Kill on sight: novel, framework, comprehensive, insights, notably, interestingly,
    moreover, furthermore, "plays a crucial role", "sheds light on", "paves the way",
    "taken together", "it is worth noting".

## 10. Translation table: AI-paper sentence → Feynman-voiced rewrite

| # | Typical AI/boilerplate sentence | Feynman-voiced rewrite |
|---|---|---|
| 1 | In this work, we present a novel framework for understanding superconductivity in hydrated sodium cobaltate. | This paper describes a mechanism for the superconductivity of hydrated sodium cobaltate. |
| 2 | Superconductivity in Na_xCoO2·yH2O has attracted significant attention due to its unique structural properties. | Sodium cobaltate superconducts only when it is wet. That is a strange fact, and it wants a mechanical explanation. |
| 3 | It is worth noting that these results may potentially suggest a phonon-mediated pairing scenario. | The results suggest — we cannot yet prove — that the pairing is phonon-mediated. |
| 4 | Leveraging state-of-the-art first-principles calculations, we systematically investigate the interlayer region. | We calculate, by ordinary density-functional methods, what happens in the gallery when the layers are pulled apart. |
| 5 | The remainder of this paper is organized as follows. Section II describes our methods... | We first show that opening the gallery creates the electron gas and the soft sodium mode together. We then ask what kills the superconductivity on either side of the window. |
| 6 | These findings provide crucial insights into the underlying pairing mechanism. | This tells us where the pairing comes from. |
| 7 | However, a comprehensive understanding of the hydration dependence remains elusive. | We do not understand what the water does. That is the question of this paper. |
| 8 | The critical role of hydration in enabling superconductivity cannot be overstated. | Without the water there is no superconductivity; the question is what the water actually does. |
| 9 | Interestingly, the isostructural Li compound exhibits no superconducting transition. | Lithium cobaltate never superconducts, and any picture of the pairing had better explain why. Ours does, and for a simple reason. |
| 10 | Our results are in reasonable agreement with experiment, although some discrepancies remain. | One measured number disagrees with our picture, and it is important. We take it up at once. |
| 11 | This discrepancy may be attributed to various factors, including dynamical effects not captured at the static level. | The discrepancy has, we believe, a single cause: the calculation holds the water still, and the water does not hold still. |
| 12 | Future work will explore the broader implications of these findings for related layered systems. | The obvious next experiment is to squeeze the hydrate and watch T_c fall when the sodium freezes. |
| 13 | It was found that the Na ion occupies a symmetric double-well potential beyond the critical spacing. | Past the critical spacing we find that the sodium ion no longer sits in one well; it has two, and it rattles between them. |
| 14 | This novel mechanism paves the way for the discovery of new superconducting materials. | If the picture is right, any layered oxide whose gallery can be opened this far and kept soft should do the same thing. |
| 15 | Notably, the emergence of the 2DEG coincides exactly with the softening of the Na mode — a striking correlation. | The electron gas and the soft mode appear at the same critical spacing. They must: they are one event seen two ways. |

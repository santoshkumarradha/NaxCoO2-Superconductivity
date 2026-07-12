# Provenance of source texts

All files retrieved 2026-07-12. Every quote in `../FEYNMAN_STYLE.md` and `../STORY_NOTES.md`
was checked against these local text files; nothing is quoted from memory.

## Retrieved

| File | Work | Retrieved from | Notes |
|---|---|---|---|
| `feynman_1942_thesis.txt` | R. P. Feynman, *The Principle of Least Action in Quantum Mechanics*, Ph.D. thesis, Princeton (1942) | `files.untiredwithloving.org/thesis.pdf` (PDF of the World Scientific 2005 reprint, ed. L. M. Brown) | Thesis text only was extracted (~24,500 words); the rest of the book (editor's preface, 1948 reprint, Dirac 1933 paper) and the book PDF were deleted. The hosting of the book PDF is third-party, not the publisher; the official edition is World Scientific ISBN 981-256-366-0 (paywalled). Kept locally for private style analysis; quote only briefly. |
| `feynman_1948_rmp.pdf`, `feynman_1948_spacetime_rmp.txt` | "Space-Time Approach to Non-Relativistic Quantum Mechanics", Rev. Mod. Phys. **20**, 367 (1948) | `jontalle.web.engr.illinois.edu/uploads/Feynman.48.pdf` (retyped course copy, Univ. of Illinois) | © APS; open copies circulate on university course pages. Official: `link.aps.org/doi/10.1103/RevModPhys.20.367`. |
| `feynman_1949_qed.pdf`, `feynman_1949_spacetime_qed.txt` | "Space-Time Approach to Quantum Electrodynamics", Phys. Rev. **76**, 769 (1949) | `alternatievewiskunde.nl/QED/eng.pdf` | Bonus source (found while locating the 1948 paper). © APS. |
| `feynman_1982.pdf`, `feynman_1982_simulating.txt` | "Simulating Physics with Computers", Int. J. Theor. Phys. **21**, 467 (1982) | `s2.smu.edu/~mitch/class/5395/papers/feynman-quantum-1981.pdf` (SMU course copy) | OCR of the journal scan; contains OCR artifacts ("qua~atum", "y o u"). Quotes were cleaned of obvious OCR noise only. |
| `feynman_nobel_1965.html`, `feynman_nobel_1965.txt` | Nobel Lecture, "The Development of the Space-Time View of Quantum Electrodynamics" (11 Dec 1965) | `nobelprize.org/prizes/physics/1965/feynman/lecture/` (official, freely readable) | © The Nobel Foundation 1965. |

## Not retrieved

- **The Feynman Lectures on Physics** (`feynmanlectures.caltech.edu`): free to read in a
  browser, but Caltech's Cloudflare configuration returns 403 to automated fetching, and the
  domain is **excluded from the Wayback Machine**. Not quoted anywhere in the style guide,
  per the honesty rule. If Lectures quotes are wanted later, copy them manually from a browser.
- The **original 1942 typescript** of the thesis (Princeton archives) was not accessed; the
  text used is the World Scientific 2005 typesetting.

## Extraction caveats

- `pdftotext` dropped "ff/fi/ffi" ligatures in the 1948 and 1942 files ("dierent" =
  "different", "denition" = "definition"). Quotes used in the style guide restore the
  ligatures; check against the PDF if in doubt.
- The 1948 file is a retyped (LaTeX) copy, not the APS scan; equation numbering matches
  the original but hyphenation/line breaks do not.

#!/usr/bin/env python
"""Apply WRL (2026-07) revisions to main.tex and supplement.tex.

Substantive edits are wrapped in `changes`-package markup (\replaced/\added/
\deleted) so the review PDF shows exactly what moved; flipping the preamble
line to \\usepackage[final]{changes} yields the clean version in one step.
Copyedits (US spelling, N(0)->N(E_F) in the SM) are applied silently.
Every replacement is asserted to match exactly once.
"""
import re
import sys

MAIN = "main.tex"
SUP = "supplement.tex"
BIB = "refs.bib"


def apply(path, edits, spelling=False):
    t = open(path).read()
    if spelling:
        t = t.replace("centred", "centered").replace("centre", "center")
    for old, new, tag in edits:
        n = t.count(old)
        if n != 1:
            print(f"!! [{tag}] match count {n} in {path}:\n   {old[:90]}...")
            sys.exit(1)
        t = t.replace(old, new)
    open(path, "w").write(t)
    print(f"{path}: {len(edits)} edits applied")


# ============================================================ main.tex ======
PREAMBLE = (
    "\\usepackage{xcolor}",
    "\\usepackage{xcolor}\n\n"
    "% --- revision markup (WRL comments, Jul 2026): the `changes` package ----\n"
    "% Compile as-is for the annotated review version; for the final clean\n"
    "% version replace the line below with \\usepackage[final]{changes} -- a\n"
    "% one-line toggle that hides all markup.\n"
    "\\usepackage{changes}",
    "preamble",
)

AUTHOR = (
    """\\author{Walter R. L. Lambrecht}
\\affiliation{Department of Physics, Case Western Reserve University, Cleveland, Ohio 44106-7079, USA}

""",
    "% WRL declined co-authorship (private communication, 14 Jul 2026):\n"
    "% \"you should be a single author on this\" -- thanked in the acknowledgments.\n\n",
    "author",
)

ACK = (
    "% TODO: acknowledgments / funding.",
    "\\added{S.K.R. thanks W.~R.~L.~Lambrecht for extensive discussions, for a\n"
    "critical reading of the manuscript, and for several of the clarifications\n"
    "reported here.}",
    "ack",
)

DIRAC = (
    """that the topology of an obstructed atomic limit, seeded by an annihilating Dirac
cone~\\cite{RadhaLambrecht2021,RadhaLambrecht2021obstructed}, requires the half-open gallery at the
surface to host a spin-split two-dimensional electron gas residing entirely outside the CoO$_2$
sheet, its density centered on the surface alkali ions~\\cite{RadhaLambrecht2021optics}.""",
    "\\replaced{that an obstructed atomic limit hosts, at the half-open gallery\n"
    "of the terminated surface, a chiral SSH$_4$ band pair whose splitting by\n"
    "the ionic asymmetry of the surface---the term that breaks the chiral\n"
    "(positive--negative-energy) symmetry of the SSH$_4$ chain---produces a\n"
    "metallic, spin-split two-dimensional electron gas dispersing parallel to\n"
    "the layers and residing entirely outside the CoO$_2$ sheet, its density\n"
    "centered on the surface alkali\n"
    "ions~\\cite{RadhaLambrecht2021,RadhaLambrecht2021obstructed,RadhaLambrecht2021optics}}"
    "{that the topology of an obstructed atomic limit, seeded by an annihilating Dirac\n"
    "cone~\\cite{RadhaLambrecht2021,RadhaLambrecht2021obstructed}, requires the half-open gallery at the\n"
    "surface to host a spin-split two-dimensional electron gas residing entirely outside the CoO$_2$\n"
    "sheet, its density centered on the surface alkali ions~\\cite{RadhaLambrecht2021optics}.}",
    "dirac",
)

BARREL = (
    "matches the Fermi-surface\nbarrel computed independently at 9.9~\\AA{} (SM~S4, S5)",
    "\\replaced{matches the Fermi-surface\nbarrel---a cylinder warped outward about $k_z=0$, which is its computed\n$k_z$ dispersion (SM~S4)---obtained independently at 9.9~\\AA{} (SM~S4, S5)}"
    "{matches the Fermi-surface\nbarrel computed independently at 9.9~\\AA{} (SM~S4, S5)}",
    "barrel",
)

GATING = (
    "and the well stiffness responds directly to the oxygen-height channel that carries the\ntransfer, $d\\alpha/dz_{\\rm O}\\approx-0.52$~eV/\\AA$^2$ per~\\AA.",
    "\\replaced{and, in an independent structural\ncontrol at fixed charge, the well stiffness responds to the oxygen height\n$z_{\\rm O}$, the CoO$_2$ sheet thickness through which the transfer acts,\nsoftening linearly as the sheet thickens, $d\\alpha/dz_{\\rm O}\\approx-0.52$~eV/\\AA$^2$\nper~\\AA}"
    "{and the well stiffness responds directly to the oxygen-height channel that carries the\ntransfer, $d\\alpha/dz_{\\rm O}\\approx-0.52$~eV/\\AA$^2$ per~\\AA.}",
    "gating",
)

FIG2PARA = (
    """Figure~\\ref{fig:money} collects $\\alpha(c)$, $N(0)$, and the superconducting dome constructed from them on a
common axis with a broken scale spanning all four measured spacings.
Inserting $\\lambda(c)$, the boson energy, and $N(0)$ into the Allen--Dynes formula (coupling model
below, $\\mu^\\ast=0.10$; carrier gate and polaron cut, SM~S1) yields a dome pinned just above
$c^\\ast$,""",
    "\\replaced{The pairing boson of the mechanism is the even overtone of the\n"
    "quantized sodium well and the coupling the quadratic electron--phonon\n"
    "vertex of Eq.~\\eqref{eq:lambda} below. Figure~\\ref{fig:money} collects\n"
    "$\\alpha(c)$, $N(E_F)$, and the superconducting dome constructed from them\n"
    "on a common axis with a broken scale spanning all four measured spacings.\n"
    "Inserting $\\lambda(c)$, the boson energy, and $N(E_F)$ into the\n"
    "Allen--Dynes formula (coupling model below, $\\mu^\\ast=0.10$; carrier gate\n"
    "and polaron cut, SM~S1) yields a dome pinned just above $c^\\ast$,}"
    "{Figure~\\ref{fig:money} collects $\\alpha(c)$, $N(0)$, and the superconducting dome constructed from them on a\n"
    "common axis with a broken scale spanning all four measured spacings.\n"
    "Inserting $\\lambda(c)$, the boson energy, and $N(0)$ into the Allen--Dynes formula (coupling model\n"
    "below, $\\mu^\\ast=0.10$; carrier gate and polaron cut, SM~S1) yields a dome pinned just above\n"
    "$c^\\ast$,}",
    "fig2para",
)

FIG3CAP = (
    "$\\Gamma$--M--K--$\\Gamma$ at (a) 5.5, (b) 6.9, and (c) 9.9~\\AA, referenced to $E_F$ (dashed);",
    "$\\Gamma$--M--K--$\\Gamma$ at (a) 5.5, (b) 6.9, and (c) 9.9~\\AA, referenced to $E_F$ (dashed)\\added{, placed by the code at\n  mid-gap of the insulating 5.5~\\AA{} cell; the chemical potential of a doped\n  semiconductor would shift within the gap with temperature, doping, and\n  band-edge masses};",
    "fig3cap",
)

EG = (
    """The absence of the $e_g'$ pockets,
on which several in-sheet pairing proposals were built, is an experimental fact already at the
parent level.""",
    "\\replaced{The absence of the $e_g'$ pockets,\n"
    "on which several in-sheet pairing proposals were built, is an experimental\n"
    "fact already at the parent level, traced to the Na charge transfer and the\n"
    "surface termination rather than to\n"
    "correlations~\\cite{Pillay2008PRB,Pillay2008PRL}.}"
    "{The absence of the $e_g'$ pockets,\n"
    "on which several in-sheet pairing proposals were built, is an experimental fact already at the\n"
    "parent level.}",
    "eg",
)

INTERCALANT = (
    "it is the interlayer gas predicted at the surface~\\cite{RadhaLambrecht2021}, now\nrealized in the bulk gallery.",
    "it is the interlayer gas predicted at the surface~\\cite{RadhaLambrecht2021}, now\nrealized in the bulk gallery\\added{; that an intercalant-derived band descends toward $E_F$ as the\nlayers separate was anticipated for the surface case~\\cite{Pillay2008PRL}}.",
    "intercalant",
)

FORTO = (
    "cobalt for 91\\% sodium",
    "\\replaced{cobalt to 91\\% sodium}{cobalt for 91\\% sodium}",
    "forto",
)

SSH4 = (
    """The two bands sharing the gallery electron, alkali
$sp_z$ and CoO$_2$, form the SSH-type pair of our earlier surface model~\\cite{RadhaLambrecht2021},
split by the same Bader charge that fills the gallery band in Fig.~\\ref{fig:money}(b).""",
    "\\replaced{The two orbitals sharing the gallery electron, the alkali\n"
    "$sp_z$ and its CoO$_2$ partner, are the two sublattices of the chiral\n"
    "SSH$_4$ surface model of Ref.~\\cite{RadhaLambrecht2021}; the\n"
    "Na$\\to$CoO$_2$ charge transfer is the ionic term of that model, breaking\n"
    "its chiral (positive--negative-energy) symmetry, so the same Bader charge\n"
    "that fills the gallery band in Fig.~\\ref{fig:money}(b) splits the\n"
    "topological band pair}"
    "{The two bands sharing the gallery electron, alkali\n"
    "$sp_z$ and CoO$_2$, form the SSH-type pair of our earlier surface model~\\cite{RadhaLambrecht2021},\n"
    "split by the same Bader charge that fills the gallery band in Fig.~\\ref{fig:money}(b).}",
    "ssh4",
)

MIDPLANE = (
    """The sodium displacement
$\\delta$ is odd under reflection through the gallery midplane, so the linear electron--phonon
coupling vanishes identically at the symmetric point (verified at several $k_z$).""",
    "\\replaced{The sodium displacement\n"
    "$\\delta$ is measured from the gallery midplane and is odd under reflection\n"
    "through it; because the ion oscillates between two minima symmetric about\n"
    "that plane, its probability density remains centered there, and the linear\n"
    "electron--phonon coupling vanishes identically at the symmetric point\n"
    "(verified at several $k_z$)}"
    "{The sodium displacement\n"
    "$\\delta$ is odd under reflection through the gallery midplane, so the linear electron--phonon\n"
    "coupling vanishes identically at the symmetric point (verified at several $k_z$).}",
    "midplane",
)

GDEF = (
    "with $N(0)$ from the calculation, not fitted, and $g$ the quadratic vertex (SM~S3).",
    "\\replaced{with $N(E_F)$ from the calculation, not fitted, and $g$ the\n"
    "quadratic electron--phonon vertex that couples the gallery band to this\n"
    "mode (SM~S3).}"
    "{with $N(0)$ from the calculation, not fitted, and $g$ the quadratic vertex (SM~S3).}",
    "gdef",
)

ROADMAP = (
    "The position of the bilayer hydrate is the central difficulty, which we address next.",
    "\\replaced{The static model developed so far already predicts\n"
    "superconductivity, but at the wrong spacing and the wrong scale,\n"
    "$T_c\\approx14$~K near 6.4~\\AA{} against the observed 4.5~K at 9.9~\\AA;\n"
    "why the only superconducting member of the family lies beyond the static\n"
    "window is the central difficulty, which we address next.}"
    "{The position of the bilayer hydrate is the central difficulty, which we address next.}",
    "roadmap",
)

INERT = (
    """This is the
inert-water null result of Refs.~\\cite{JohannesSingh2004,Arita2005} reappearing at the level of
the lattice dynamics, and it indicates that the water must be treated dynamically.""",
    "\\replaced{That static water leaves the electronic structure unchanged was\n"
    "established in Refs.~\\cite{JohannesSingh2004,Arita2005}; the present null\n"
    "result is distinct and stronger: treated statically, the water pushes the\n"
    "only superconducting member of the family deeper into the polaronic\n"
    "regime. The water must therefore be treated dynamically.}"
    "{This is the\n"
    "inert-water null result of Refs.~\\cite{JohannesSingh2004,Arita2005} reappearing at the level of\n"
    "the lattice dynamics, and it indicates that the water must be treated dynamically.}",
    "inert",
)

TIMESCALE = (
    """the sodium mode has a
period near 200~fs, while the network reorganizes on the picosecond scale measured by quasielastic
neutron scattering~\\cite{Jalarvo2008}.""",
    "\\replaced{the sodium mode is a local oscillation in its shallow well,\n"
    "with a period near 200~fs, whereas the water does not relax molecule by\n"
    "molecule: its motion is the cooperative reorganization of a hydrogen-bonded\n"
    "network, measured by quasielastic neutron scattering on the picosecond\n"
    "scale~\\cite{Jalarvo2008}. The heavier ion is therefore, counter to\n"
    "intuition, the faster degree of freedom: it moves alone, while the network\n"
    "must move together.}"
    "{the sodium mode has a\n"
    "period near 200~fs, while the network reorganizes on the picosecond scale measured by quasielastic\n"
    "neutron scattering~\\cite{Jalarvo2008}.}",
    "timescale",
)

MDOPEN = (
    "and the decisive open test is a finite-temperature Born--Oppenheimer or path-integral\nmolecular-dynamics sampling of the coupled sodium--water system (SM~S7).",
    "\\replaced{and the decisive open test---proposed\n"
    "here, not yet performed---is a finite-temperature Born--Oppenheimer or\n"
    "path-integral molecular-dynamics sampling of the coupled sodium--water\n"
    "system (SM~S7).}"
    "{and the decisive open test is a finite-temperature Born--Oppenheimer or path-integral\nmolecular-dynamics sampling of the coupled sodium--water system (SM~S7).}",
    "mdopen",
)

STONER = (
    """Immediately above the carrier onset the gallery band is narrow enough to be
Stoner-polarized, producing a thin magnetic strip adjacent to the dome (SM~S8, where the
computed regimes are assembled into a phase diagram); this accounts""",
    "\\replaced{Immediately above the carrier onset the gallery band is at its\n"
    "narrowest and its $N(E_F)$ at its largest, so the Stoner product\n"
    "$N(E_F)\\,I$ exceeds threshold in a thin strip adjacent to the dome;\n"
    "further out the band widens, the product drops, and the strip closes\n"
    "(SM~S8, where the computed regimes are assembled into a phase diagram).\n"
    "This accounts}"
    "{Immediately above the carrier onset the gallery band is narrow enough to be\n"
    "Stoner-polarized, producing a thin magnetic strip adjacent to the dome (SM~S8, where the\n"
    "computed regimes are assembled into a phase diagram); this accounts}",
    "stoner",
)

NACOSEO = (
    """Na$_2$CoSe$_2$O~\\cite{Cheng2024} superconducts at 6.3~K through its Co--Se
network; we find its gallery band unoccupied, so it demonstrates the distinction the picture
requires: a cobaltate can superconduct without the gallery, but then not through it (SM~S12).""",
    "\\replaced{As a negative control on the mechanism we computed\n"
    "Na$_2$CoSe$_2$O~\\cite{Cheng2024}, which superconducts at 6.3~K through its\n"
    "Co--Se network; its gallery band is unoccupied, so it demonstrates the\n"
    "distinction the picture requires: a cobaltate can superconduct without the\n"
    "gallery, but then not through it (SM~S12).}"
    "{Na$_2$CoSe$_2$O~\\cite{Cheng2024} superconducts at 6.3~K through its Co--Se\n"
    "network; we find its gallery band unoccupied, so it demonstrates the distinction the picture\n"
    "requires: a cobaltate can superconduct without the gallery, but then not through it (SM~S12).}",
    "nacoseo",
)

SELFEN = (
    """together with a
20--30~meV self-energy feature from the confined sodium mode""",
    "\\replaced{together with a\n"
    "20--30~meV self-energy kink---a mass-enhancement feature imprinted on the\n"
    "gallery band by the confined sodium mode, the photoemission signature of\n"
    "the coupling in Eq.~\\eqref{eq:lambda}}"
    "{together with a\n"
    "20--30~meV self-energy feature from the confined sodium mode}",
    "selfen",
)

X13 = (
    """(full-coverage cell; the physical $x=1/3$ cell shifts the
crossing by only $\\approx0.2$~\\AA, SM~S2).""",
    """(full-coverage cell; the physical $x=1/3$ cell shifts the
crossing by only $\\approx0.2$~\\AA, SM~S2)\\added{. At that composition the identification
is direct: in the ordered $x=1/3$ cell the sodium-projected band at $\\Gamma$ is
empty at 8.4~\\AA{} and occupied at 9.9~\\AA{} (92\\% filling, 1.5~eV of
in-plane broadening), and the charge the sodium retains in the Bader count
grows accordingly (SM~S5)}""",
    "x13",
)

main_edits = [PREAMBLE, AUTHOR, ACK, DIRAC, BARREL, GATING, FIG2PARA, FIG3CAP,
              EG, INTERCALANT, FORTO, SSH4, MIDPLANE, GDEF, ROADMAP, INERT,
              TIMESCALE, MDOPEN, STONER, NACOSEO, SELFEN, X13]

# ======================================================== supplement.tex ====
SUP_FIG = (
    """in Fig.~\\ref{fig:s3fold} should be read as folded.""",
    """in Fig.~\\ref{fig:s3fold} should be read as folded.

\\added{Figure~\\ref{fig:galleryx13} quantifies the onset at this composition.
The sodium-projected gallery band, 40--55\\% Na in character at $\\Gamma$, lies
$+1.4$~eV above $E_F$ in the closed gallery and descends through it as the
spacing opens [panel~(a)]; at $\\Gamma$ it crosses between 8.4 and 9.9~\\AA{}
[panel~(b)], where it is 92\\% occupied. Its in-plane dispersion width, 2.9~eV
at 5.5~\\AA{} and 1.5~eV at 9.9~\\AA{} [panel~(a)], is the broadening the
surface-model picture requires of the split SSH$_4$ band pair. The Bader count
tracks the same event [panel~(c)]: with the gallery closed the sodium donates
$0.85$--$0.94\\,e$ to the CoO$_2$ sheets, but only $\\approx0.65\\,e$ once the
pocket is open, the difference being the 2DEG charge residing in the gallery
basin. There is no finite-$\\delta$ threshold at 9.9~\\AA: the pocket is
present at $\\delta=0$ and is partially emptied as the ion moves off-center
(occupancy 0.92 at $\\delta=0$ to 0.74 at $\\delta=0.75$~\\AA), a negative
feedback that assists in bounding the double-well minimum [panel~(d)].}

% ---- FIG S: gallery-band onset, Bader, and wells at x = 1/3 ----
\\begin{figure}[t]
  \\centering
  \\figorbox{figS_gallery_x13}{figS_gallery_x13}{\\columnwidth}
  \\caption{\\added{The gallery-band onset at the physical composition, $x=1/3$.
  (a) Na-projected gallery band along $\\Gamma$--M--K--$\\Gamma$ (spin up,
  marker size proportional to Na weight) at 5.5 and 9.9~\\AA, referenced to
  $E_F$; the band drops $\\approx5$~eV as the gallery opens and crosses $E_F$
  only at 9.9~\\AA{} (widths $W$ in the panel). (b) The same band at $\\Gamma$
  against spacing for $\\delta=0$ and $0.75$~\\AA; the pocket opens between 8.4
  and 9.9~\\AA{} (occupancy 0.92 at $\\delta=0$). (c) Bader charge donated by
  the sodium to the CoO$_2$ sheets against displacement (9 valence $e$ less
  the Na-basin count, semicore-inclusive PAW). (d) $E(\\delta)$ wells at the
  expanded spacings; the 5.5~\\AA{} well lies far above the frame
  ($+2.56$~eV at $\\delta=0.75$~\\AA). All quantities re-parsed from the set-B
  scan outputs; no new calculations.}}
  \\label{fig:galleryx13}
\\end{figure}""",
    "supfig",
)

SUP_04 = (
    """\\subsection{The 0.4 eV configurational spread and the sixty-step artifact}

The single relaxed configuration is only the deepest point of a landscape over which the four-water
hydrogen-bond network fluctuates. A relaxation of the same $\\delta=0.30$~\\AA{} cell,
stopped at sixty ionic steps rather than carried to full convergence, comes to rest 394~meV above
the fully converged minimum, a different, shallower local minimum of the cage at the identical
sodium position (the open star in the Letter's Fig.~4(a)). The hydrogen-bond network therefore
presents the sodium with a family of potentials rather than a single one, spanning at least 0.4~eV
between local minima at fixed $\\delta$; the sodium potential is set by the instantaneous cage
configuration, of which the deep adiabatic surface is only the lower envelope. (This also explains
why an earlier interpretation of the sixty-step run as a near-flat 5~meV well was an artifact of
incomplete relaxation rather than a physical softening; the corrected, fully converged comparison
is the one reported throughout.)""",
    """\\subsection{The 0.4 eV configurational spread}

The single relaxed configuration is only the deepest point of a landscape over which the four-water
hydrogen-bond network fluctuates. \\replaced{A relaxation of the same $\\delta=0.30$~\\AA{} cell,
stopped at sixty ionic steps rather than carried to full convergence, comes to rest 394~meV above
the fully converged minimum. We do not claim this partially relaxed point as a
distinct local minimum---an unconverged run stopping high is expected---but as
an indication of the scale of the landscape; the direct demonstration of
multiple converged minima is the relaxation ensemble of the next subsection.}{A relaxation of the same $\\delta=0.30$~\\AA{} cell,
stopped at sixty ionic steps rather than carried to full convergence, comes to rest 394~meV above
the fully converged minimum, a different, shallower local minimum of the cage at the identical
sodium position (the open star in the Letter's Fig.~4(a)).} The hydrogen-bond network therefore
presents the sodium with a family of potentials rather than a single one, spanning at least 0.4~eV
between local minima at fixed $\\delta$; the sodium potential is set by the instantaneous cage
configuration, of which the deep adiabatic surface is only the lower envelope.
\\deleted{(This also explains
why an earlier interpretation of the sixty-step run as a near-flat 5~meV well was an artifact of
incomplete relaxation rather than a physical softening; the corrected, fully converged comparison
is the one reported throughout.)}

\\subsection{\\added{The relaxation ensemble}}
\\label{sec:wensemble}

\\added{[To be filled on completion of the running set-W calculation: sixteen
fully converged BFGS relaxations ($n_{\\rm step}=250$, force threshold
$10^{-3}$~Ry/bohr) of the same $\\delta=0.30$~\\AA{} hydrate cell, each started
from a different seeded-random orientation of the four water molecules
(water-oxygen sites fixed at the Sec.~\\ref{sec:water} ansatz), plus four
controls at $\\delta=0$. Distinct converged energies across the starts
demonstrate a rough multi-minimum landscape directly, replacing the
sixty-step indication above; identical energies would falsify it. The
converged energy spread and the final Na height of each basin will be
reported here.]}""",
    "sup04",
)

SUP_STONER_S = (
    "(b) Cuts at $x=0.35$: the carrier fraction $q$ per Na and the Stoner enhancement $S$ against $c$.",
    "(b) Cuts at $x=0.35$: the carrier fraction $q$ per Na and \\replaced{the Stoner enhancement\n  $S=N(E_F)\\,I$, with the Stoner $I$ taken from the surface model and the\n  dashed line $S^\\ast=3$ an ansatz threshold above which we call the gas\n  polarized (sensitivity range 2--5, source-tagged in the generating\n  script),}{the Stoner enhancement $S$} against $c$.",
    "supstoners",
)

SUP_STONER_G2D = (
    "in the plane of sodium content $x$ and CoO$_2$ spacing $c$, the gallery band is empty\n  below the carrier turn-on line,",
    "in the plane of sodium content $x$ and CoO$_2$ spacing $c$, the gallery band is empty\n  below the carrier turn-on line \\added{(the model condition $\\Gamma=2\\delta$: the\n  Na$\\to$CoO$_2$ level splitting of the SSH$_4$-type pair reached by the\n  transfer)},",
    "supstonerg2d",
)

SUP_DFTM = (
    "(c) The DFT cell-integrated moment $\\langle\\lvert m\\rvert\\rangle$ for Na and Li $1\\times1$,\n  switching on at the same spacing the well transition does;",
    "(c) The DFT cell-integrated moment $\\langle\\lvert m\\rvert\\rangle$ for Na and Li $1\\times1$\n  \\added{(the ``DFT $\\lvert m\\rvert$ onset'' marker of the panel)},\n  switching on at the same spacing the well transition does;",
    "supdftm",
)

SUP_WIDEN = (
    "strongest just past\nit and weakening as the band widens,",
    "strongest just past\nit and weakening as the band widens \\added{(with increasing spacing a larger\nfraction of the gallery band drops below $E_F$, so the occupied bandwidth\ngrows and the Stoner product falls back below threshold)},",
    "supwiden",
)

SUP_LSDA = (
    """One quantity is excluded as
unreliable: the moment of the $\\sqrt3\\times\\sqrt3$ cell, over-polarized in the local-spin-density
approximation and omitted from the diagram; the well energies on which the mechanism rests shift by
less than 6\\% under the spin treatment (Sec.~\\ref{sec:correlations}), so the exclusion does not
affect the structural conclusions.""",
    """\\replaced{One quantity is excluded from the diagram on stated grounds: the
moment of the $\\sqrt3\\times\\sqrt3$ cell. The local-spin-density approximation
systematically over-stabilizes magnetism in the narrow-band limit---it
underestimates the exchange narrowing of local moments, moving the Stoner
boundary too far into the paramagnetic side---and the dilute $x=1/3$ cell,
with the narrowest gallery band of the set, is where that error concentrates;
its moment is therefore omitted rather than plotted. The omission is
conservative and immaterial to the structure: the well energies on which the
mechanism rests shift by less than 6\\% under the spin treatment
(Sec.~\\ref{sec:correlations}).}
{One quantity is excluded as
unreliable: the moment of the $\\sqrt3\\times\\sqrt3$ cell, over-polarized in the local-spin-density
approximation and omitted from the diagram; the well energies on which the mechanism rests shift by
less than 6\\% under the spin treatment (Sec.~\\ref{sec:correlations}), so the exclusion does not
affect the structural conclusions.}""",
    "suplsda",
)

SUP_HS = (
    "There is no\nHebel--Slichter peak because the pairing boson is a soft, strongly anharmonic mode, and the",
    "There is \\replaced{no\nHebel--Slichter peak (the NMR $1/T_1$ coherence peak just below $T_c$,\nthe classic signature of a singlet phonon-mediated gap)}{a\nHebel--Slichter peak} because the pairing boson is a soft, strongly anharmonic mode, and the",
    "suphs",
)

sup_edits = [SUP_FIG, SUP_04, SUP_STONER_S, SUP_STONER_G2D, SUP_DFTM,
             SUP_WIDEN, SUP_LSDA, SUP_HS]

# ================================================================ refs.bib ==
bib_add = """
@article{Pillay2008PRB,
  author = {Pillay, D. and Johannes, M. D. and Mazin, I. I. and Andersen, O. K.},
  title = {Origin of $a_{1g}$ and $e_g'$ orderings in Na$_x$CoO$_2$},
  journal = {Phys. Rev. B},
  volume = {78},
  pages = {012501},
  year = {2008}
}

@article{Pillay2008PRL,
  author = {Pillay, D. and Johannes, M. D. and Mazin, I. I.},
  title = {Electronic Structure of the Na$_x$CoO$_2$ Surface},
  journal = {Phys. Rev. Lett.},
  volume = {101},
  pages = {246808},
  year = {2008}
}
"""

if __name__ == "__main__":
    apply(MAIN, main_edits, spelling=True)
    apply(SUP, sup_edits, spelling=True)
    # silent notation copyedit in the SM: N(0) -> N(E_F) (marked in the Letter)
    t = open(SUP).read()
    t = t.replace("$N(0)$", "$N(E_F)$")
    open(SUP, "w").write(t)
    # SM preamble: add the changes package if absent
    t = open(SUP).read()
    if "usepackage{changes}" not in t:
        for anchor in ("\\usepackage{xcolor}", "\\usepackage{booktabs}",
                       "\\usepackage{amssymb}"):
            if anchor in t:
                t = t.replace(anchor, anchor + "\n\\usepackage{changes} "
                              "% revision markup; [final] option hides all",
                              1)
                break
        open(SUP, "w").write(t)
        print("supplement.tex: changes package added")
    b = open(BIB).read()
    if "Pillay2008PRB" not in b:
        open(BIB, "a").write(bib_add)
        print("refs.bib: Pillay entries added")
    print("done")

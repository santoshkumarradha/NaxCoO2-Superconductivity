#!/usr/bin/env python
"""
Figure 3 - Theory vs experiment scoreboard.
A paper-style summary matrix: the six experimental facts about Na_xCoO2.yH2O
(rows) against the verdict of the soft-mode / charge-bifurcation mechanism
(explained = green check, prediction = amber). Glyphs are drawn as vector
paths (no glyph-font dependency). References cited in the footer.
"""
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

import _style as S

S.use_house_style()

# (fact, anchor tag, verdict, verdict-word, one-line reason [hard-wrapped])
ROWS = [
    ("Anhydrous (5.5 Å) is not SC", "Foo 2003", "ok", "explained",
     r"$\alpha>0$: single stiff well ($\hbar\omega\!\approx\!30$ meV),"
     "\n" r"$N(0)\!\approx\!0$ — no 2DEG, nothing to pair."),
    ("Monolayer hydrate (6.9 Å) is not SC", "Foo 2003", "ok", "explained",
     r"Past the window: $\lambda\!\approx\!20\!\gg\!2$ — the alkali"
     "\n" r"self-traps (polaronic), so SC is killed."),
    (r"Bilayer hydrate (9.9 Å) is SC, $T_c\!=\!4.5$ K", "Takada 2003",
     "pred", "prediction",
     r"Bare geometry is polaronic here; SC needs water to"
     "\n" r"soften the Na well back into the window (untested)."),
    ("No superconducting Li analogue", "", "ok", "explained",
     r"$c^\ast_{\mathrm{Li}}\!\leq\!5.5$ Å and lighter mass $\Rightarrow$ stiffer mode;"
     "\n" r"Li never enters the soft-mode window."),
    ("Magnetic phase borders the dome", "expt.", "ok", "explained",
     r"Local moment turns on exactly at the $\alpha<0$ well"
     "\n" r"transition (Fig. 4), coincident with the dome."),
    ("Superconducting $T_c$ is only a few K", "Takada 2003", "ok",
     "explained",
     r"Allen$-$Dynes with DFT $\omega_{\mathrm{eff}},\lambda$ gives peak"
     "\n" r"$T_c\!\approx\!14$ K (band 11$-$23) — the right scale."),
]

STAT = {"ok": S.C_GOOD, "pred": S.C_WARN}

fig, ax = plt.subplots(figsize=(6.9, 4.5))
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.axis("off")

# column layout (axes fraction)
X_FACT, X_VERD, X_REAS = 0.015, 0.395, 0.545
top, bot = 0.90, 0.09
n = len(ROWS)
rh = (top - bot) / n

# ---- header ----
ax.text(X_FACT, 0.955, "EXPERIMENTAL FACT", fontsize=8.5, fontweight="bold",
        color=S.C_INK, va="center")
ax.text(X_VERD, 0.955, "OUR VERDICT", fontsize=8.5,
        fontweight="bold", color=S.C_INK, va="center", ha="left")
ax.text(X_REAS, 0.955, "WHY IT FOLLOWS", fontsize=8.5, fontweight="bold",
        color=S.C_INK, va="center")
ax.plot([X_FACT, 0.985], [0.925, 0.925], color=S.C_INK, lw=1.0)


for i, (fact, tag, stat, word, reason) in enumerate(ROWS):
    yc = top - (i + 0.5) * rh
    if i % 2 == 1:
        ax.add_patch(Rectangle((X_FACT - 0.005, yc - rh / 2), 0.995,
                               rh, facecolor="#f4f3ef", edgecolor="none",
                               zorder=0))
    col = STAT[stat]
    # fact
    ax.text(X_FACT, yc + 0.010, fact, fontsize=7.8, color=S.C_INK,
            va="center", ha="left", fontweight="bold")
    if tag:
        ax.text(X_FACT, yc - 0.028, tag, fontsize=6.2, color=S.C_MUT,
                va="center", ha="left", style="italic")
    # verdict chip: round marker (display-space) + drawn white glyph
    cx, cy = X_VERD + 0.045, yc
    ax.plot(cx, cy, "o", ms=15, color=col, mew=0, zorder=3)
    if stat == "ok":
        ax.plot([cx - 0.008, cx - 0.002, cx + 0.009],
                [cy + 0.002, cy - 0.008, cy + 0.011],
                color="white", lw=1.6, solid_capstyle="round",
                solid_joinstyle="round", zorder=4)
    else:
        ax.plot([cx, cx], [cy + 0.011, cy - 0.002], color="white", lw=1.8,
                solid_capstyle="round", zorder=4)
        ax.plot(cx, cy - 0.010, "o", ms=2.4, color="white", mew=0, zorder=4)
    ax.text(cx + 0.026, cy, word, fontsize=7.2, color=col, va="center",
            ha="left", fontweight="bold")
    # reason
    ax.text(X_REAS, yc, reason, fontsize=7.2, color=S.C_SEC, va="center",
            ha="left", linespacing=1.45)
    # row separator
    ax.plot([X_FACT, 0.985], [yc - rh / 2, yc - rh / 2], color=S.C_GRID,
            lw=0.6, zorder=1)

# footer references
ax.text(X_FACT, 0.035,
        r"Foo $et\ al.$, Solid State Commun. $\mathbf{127}$, 33 (2003)   $\cdot$   "
        r"Takada $et\ al.$, Nature $\mathbf{422}$, 53 (2003)",
        fontsize=6.3, color=S.C_MUT, va="center")

fig.subplots_adjust(left=0.015, right=0.985, top=0.99, bottom=0.01)
S.save(fig, "fig3")
print("fig3 written")

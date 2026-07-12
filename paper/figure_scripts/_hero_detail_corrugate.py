#!/usr/bin/env python
"""Deterministic corrugation of the 2DEG band in fig0_hero_detail.png.

Physics target (Radha & Lambrecht 2021; legacy/monolayer.cube planar average):
the gallery electron gas is maximal ON the Na plane and AROUND each ion
(91% Na character), corrugated, thinning between ions.  So: densify/widen the
dotted band into a halo hugging the Na+ rim, taper the band thinner/fainter
away from the ion.  Midline stays at the Na plane.  No text or geometry is
touched: all ops are darken-only dot stamps gated by a light-pixel rule, plus
a white-multiply taper confined to band rows away from the ion and away from
the label boxes.

Measured geometry (full-res 6336x2688):
  Na center (3172.5, 1344.5), sphere radius 188 (+ ring ~10)
  dot rows y = 1212.2, 1277.6, 1344.3, 1410.7, 1477.5   (spacing 66.3)
  x pitch 189.5; row phases (x mod pitch): 64, 132, 94.5, 132, 64
  dot: core departure-from-wash (71, 58, 34), radius ~24 core / ~29 with AA
  band rows ~1100-1590; wash bg ~(231,237,247) at band core
  "2DEG" label bbox x 5907-6335, y 1285-1401; left c-axis arrow/label x < 200
"""
import numpy as np
from PIL import Image

SRC = "../figures/fig0_hero_detail.png"
OUT = "../figures/fig0_hero_detail.png"          # in-place (backup kept in archive)

img = np.asarray(Image.open(SRC).convert("RGB")).astype(np.float64)
H, W, _ = img.shape
assert (H, W) == (2688, 6336), (H, W)

CX, CY = 3172.5, 1344.5          # Na center
PITCH, ROWDY = 189.5, 66.3
ROWS = [1212.2, 1277.6, 1344.3, 1410.7, 1477.5]
PHASES = [64.0, 132.0, 94.5, 132.0, 64.0]

yy = np.arange(H)[:, None]
xx = np.arange(W)[None, :]

def smoothstep(e0, e1, v):
    t = np.clip((v - e0) / (e1 - e0), 0.0, 1.0)
    return t * t * (3 - 2 * t)

# ---------------------------------------------------------------- 1) taper --
# White-multiply the band away from the ion: fainter overall, fringes fade
# more than the midline (=> band reads thinner).  Zero effect for |dx|<550.
band = (yy >= 1085) & (yy <= 1605)
dx = xx - CX
dyc = yy - CY
T = 0.52 * smoothstep(550.0, 2400.0, np.abs(dx))            # lateral taper
V = 0.62 + 0.38 * smoothstep(55.0, 195.0, np.abs(dyc))      # thin the fringes
alpha = np.where(band, T * V, 0.0)

# guard boxes: never touch the "2DEG" label or the left c-axis arrow/label
guard = np.ones((H, W))
guard[:, :210] = 0.0                                   # left arrow + c = 9.9 A
guard[1240:1450, 5860:] = 0.0                          # 2DEG label (padded)
# feather guards so no seam cuts through dots
from scipy.ndimage import gaussian_filter
guard = gaussian_filter(guard, 25)
alpha *= guard
# only lighten light-ish pixels (dots/wash); leave any dark ink alone
lightness = img.min(axis=2)
alpha *= smoothstep(150.0, 200.0, lightness)

img = img * (1 - alpha[..., None]) + 255.0 * alpha[..., None]

# --------------------------------------------------------- 2) halo stamps --
# elliptical "isosurface" metric around the ion
d_ell = np.sqrt(dx**2 + (2.0 * dyc) ** 2)
r_circ = np.sqrt(dx**2 + dyc**2)

DOT_A = np.array([71.0, 58.0, 34.0])     # departure amplitude per channel
R_IN, R_OUT = 20.0, 29.0                 # dot profile radii

def stamp_dots(centers, weights):
    """Darken-only stamp of synthetic dots (phase-matched clones)."""
    global img
    dep = np.zeros((H, W))
    for (px, py), w in zip(centers, weights):
        if w <= 0.003:
            continue
        x0, x1 = int(px - 32), int(px + 33)
        y0, y1 = int(py - 32), int(py + 33)
        if x0 < 0 or y0 < 0 or x1 > W or y1 > H:
            continue
        gy = np.arange(y0, y1)[:, None] - py
        gx = np.arange(x0, x1)[None, :] - px
        r = np.sqrt(gx**2 + gy**2)
        prof = 1.0 - smoothstep(R_IN, R_OUT, r)
        dep[y0:y1, x0:x1] = np.maximum(dep[y0:y1, x0:x1], w * prof)
    # light-pixel gate: never stamp over sphere/ring/arrow/dashes/waters/text
    gate = smoothstep(185.0, 215.0, lightness)
    dep *= gate
    dep[r_circ < 214] = 0.0              # hard exclusion: Na sphere + ring
    img = np.clip(img - dep[..., None] * DOT_A[None, None, :], 0, 255)

def halo_w(px, py, full, zero):
    d = np.sqrt((px - CX) ** 2 + (2.0 * (py - CY)) ** 2)
    return float(1.0 - smoothstep(full, zero, d))

# (a) interleaved half-pitch dots on the five existing rows -> denser near ion
cent, wts = [], []
for ry, ph in zip(ROWS, PHASES):
    for m in range(-40, 40):
        px = ph + (m + 0.5) * PITCH
        if abs(px - CX) > 1000:
            continue
        w = halo_w(px, ry, 260.0, 900.0)
        cent.append((px, ry)); wts.append(w)
stamp_dots(cent, wts)

# (b) two extra rows just outside the band -> halo bulges around the rim
cent, wts = [], []
for ry, ph in [(ROWS[0] - ROWDY, 132.0), (ROWS[-1] + ROWDY, 132.0)]:
    for m in range(-40, 40):
        for half in (0.0, 0.5):          # full + interleaved on bulge rows
            px = ph + (m + half) * PITCH
            if abs(px - CX) > 750:
                continue
            w = halo_w(px, ry, 210.0, 640.0)
            cent.append((px, ry)); wts.append(w)
stamp_dots(cent, wts)

# (c) darken the existing lattice dots near the ion (denser look, not bigger)
cent, wts = [], []
for ry, ph in zip(ROWS, PHASES):
    for m in range(-40, 40):
        px = ph + m * PITCH
        if abs(px - CX) > 1000:
            continue
        w = 0.5 * halo_w(px, ry, 320.0, 850.0)
        cent.append((px, ry)); wts.append(w)
stamp_dots(cent, wts)

# (d) gentle wash glow hugging the rim (soft isosurface, not smoke)
glow = np.exp(-np.clip(d_ell - 200.0, 0, None) ** 2 / (2 * 230.0**2))
glow *= smoothstep(165.0, 210.0, lightness)
glow[r_circ < 206] = 0.0
WASH_A = np.array([26.0, 19.0, 8.0])
img = np.clip(img - glow[..., None] * WASH_A[None, None, :], 0, 255)

Image.fromarray(img.astype(np.uint8)).save(OUT)
print("wrote", OUT)

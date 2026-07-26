#!/usr/bin/env python
"""Deterministic Gaussian modulation of the dotted 2DEG fill in
fig0_hero_boost.png (row-1 triptych), scenes ii and iii only.

The uniform dotted gallery fill is reshaped with a soft vertical Gaussian
centred on the Na plane (gallery midplane): dots fade toward the CoO2 sheets
and are subtly intensified around the amber ions.

The gallery dots are small dark slate specks whose color matches the
octahedra, so they are separated MORPHOLOGICALLY: dark structures that
survive a binary opening (octahedra, thick strokes) are protected; the tiny
specks that vanish under opening are the dots and may be modulated.  Labels,
the ion-pair arrow and the water cage are additionally covered by feathered
guard boxes; amber/red pixels are excluded by a color rule.

Measured geometry (full-res 6336x2688):
  scene ii: dotted x 2130-4077, y ~917-1435; midplane y0=1170; ions at
            (3137,1052) solid / (3137,1285) ghost; "2DEG" x3805-4037 y1135-1217
  scene iii: dotted x 4228-6184, y ~687-1618; midplane y0=1165; Na at
            (5215,1172); water cage bbox x~4830-5600 y~790-1560;
            "2DEG" x5886-6130 y1130-1221
"""
import numpy as np
from PIL import Image
from scipy.ndimage import gaussian_filter, binary_opening, binary_dilation

SRC = "../archive/ai_mockups/hero/boost_pre_corrugation.png"   # clean backup
OUT = "../figures/fig0_hero_boost.png"

img = np.asarray(Image.open(SRC).convert("RGB")).astype(np.float64)
H, W, _ = img.shape
assert (H, W) == (2688, 6336), (H, W)
r, g, b = img[:, :, 0], img[:, :, 1], img[:, :, 2]

yy = np.arange(H)[:, None].astype(float)
xx = np.arange(W)[None, :].astype(float)

def smoothstep(e0, e1, v):
    t = np.clip((v - e0) / (e1 - e0), 0.0, 1.0)
    return t * t * (3 - 2 * t)

def disk(rad):
    d = np.arange(-rad, rad + 1)
    return (d[:, None] ** 2 + d[None, :] ** 2) <= rad * rad

lightness = img.min(axis=2)
darkish = lightness < 205

# large dark structures (octahedra bodies + thick outlines) survive opening;
# the tiny gallery dots (~10 px specks) do not
large = binary_opening(darkish, structure=disk(7))
protected = binary_dilation(large, structure=disk(11))
protect_soft = gaussian_filter(protected.astype(float), 4)

# feathered guard boxes: labels, ion-pair arrow, water cage
guard = np.ones((H, W))
guard[1105:1245, 3775:4065] = 0.0        # scene ii "2DEG"
guard[1100:1250, 5855:6160] = 0.0        # scene iii "2DEG"
guard[1080:1265, 3075:3205] = 0.0        # scene ii double-headed arrow
guard_boost = gaussian_filter(guard.copy(), 18)   # boost: no cage guard
guard[780:1570, 4820:5610] = 0.0         # scene iii H2O cage + Na (whitening)
guard = gaussian_filter(guard, 18)

# color rule: never touch amber/red (ions, ghost, waters); dots/white pass
colorok = smoothstep(-30.0, -8.0, b - r)

# vertical Gaussian fade profile per scene
alpha = np.zeros((H, W))
A = 0.80
for x0, x1, y0lo, y0hi, ymid, sig in (
        (2115, 4095, 820, 1500, 1170.0, 120.0),   # scene ii
        (4212, 6200, 620, 1700, 1165.0, 195.0)):  # scene iii
    dy = yy - ymid
    prof = A * (1.0 - np.exp(-dy**2 / (2.0 * sig**2)))
    reg = np.zeros((H, W))
    reg[y0lo:y0hi, x0:x1] = 1.0
    reg = gaussian_filter(reg, 12)
    alpha = np.maximum(alpha, prof * reg)

alpha *= (1.0 - protect_soft) * guard * colorok
img = img * (1 - alpha[..., None]) + 255.0 * alpha[..., None]

# ---- subtle intensification of dots around the amber ions -----------------
boost = np.zeros((H, W))
for cx, cy, sigr, amp in (
        (3137.0, 1052.0, 150.0, 0.30),   # scene ii solid ion
        (3137.0, 1285.0, 150.0, 0.18),   # scene ii ghost (lower well) ion
        (5215.0, 1172.0, 185.0, 0.30)):  # scene iii Na
    rad2 = (xx - cx) ** 2 + (yy - cy) ** 2
    boost = np.maximum(boost, amp * np.exp(-rad2 / (2.0 * sigr**2)))
boost *= (1.0 - protect_soft) * guard_boost * colorok * (lightness < 246)
dep = 255.0 - img
img = np.clip(img - boost[..., None] * dep, 0, 255)

Image.fromarray(img.astype(np.uint8)).save(OUT)
print("wrote", OUT)

"""E(delta, c) scan for NaCoO2: Na displacement between CoO2 layers vs layer spacing.

Produces the raw data for three curves vs c:
  alpha(c)  — quadratic coefficient of E = alpha*d^2 + beta*d^4
  n_Na(c)   — Bader/band charge on Na (2DEG turn-on)
  dE/dd     — deformation potential input for lambda(c)

Usage: python na_double_well.py [--smoke]
"""
import argparse, json, os
import numpy as np
from ase import Atoms
from gpaw import GPAW, PW, FermiDirac

p = argparse.ArgumentParser()
p.add_argument("--smoke", action="store_true", help="tiny cutoff/kpts sanity run")
p.add_argument("--element", default="Na", choices=["Na", "Li"])
args = p.parse_args()

A = 2.888 if args.element == "Na" else 2.82   # in-plane lattice constant (Ang)
Z_O = 0.96                                     # O height above/below Co plane (Ang)
ecut = 300 if args.smoke else 600
kpts = (4, 4, 2) if args.smoke else (12, 12, 4)

# O3-like stack, 1 formula unit: CoO2 layer + alkali in the gallery.
# c = CoO2-plane to CoO2-plane distance (one gallery per cell, periodic in z).
c_values = [5.2, 5.46, 5.75, 6.1, 6.5, 6.94, 7.5] if not args.smoke else [5.46, 6.5]
d_values = np.linspace(0.0, 1.0, 9) if not args.smoke else np.array([0.0, 0.5])

os.makedirs("data", exist_ok=True)
results = []
for c in c_values:
    for d in d_values:
        cell = [[A, 0, 0], [-A / 2, A * np.sqrt(3) / 2, 0], [0, 0, c]]
        # Co at origin, O above/below, alkali mid-gallery displaced by d along z
        atoms = Atoms(
            f"CoOO{args.element}",
            cell=cell, pbc=True,
            scaled_positions=[
                (0, 0, 0),
                (1 / 3, 2 / 3, Z_O / c),
                (2 / 3, 1 / 3, -Z_O / c % 1),
                (1 / 3, 2 / 3, (c / 2 + d) / c),
            ],
        )
        calc = GPAW(mode=PW(ecut), xc="PBE", kpts=kpts,
                    occupations=FermiDirac(0.05), spinpol=True,
                    txt=f"data/{args.element}_c{c:.2f}_d{d:.2f}.txt")
        atoms.calc = calc
        e = atoms.get_potential_energy()
        calc.write(f"data/{args.element}_c{c:.2f}_d{d:.2f}.gpw")
        results.append({"c": c, "d": float(d), "E": float(e)})
        print(f"c={c:.2f} d={d:.2f} E={e:.4f} eV", flush=True)
        json.dump(results, open(f"data/{args.element}_scan.json", "w"), indent=1)
print("done")

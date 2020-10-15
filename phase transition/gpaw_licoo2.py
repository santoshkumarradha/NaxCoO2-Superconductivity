import pymatgen.io.ase as pase
from ase import Atoms
from gpaw import GPAW, PW
from tqdm import tqdm

###----load structures
structures=pickle.load(open("min_structures_licoo2", "rb"))
total_energies=[]
for i in tqdm(structures):
    fname=str(np.round(i.structure.lattice.c,3))
    bulk=pase.AseAtomsAdaptor().get_atoms(structures[0].structure)
    k = 4
	calc = GPAW(mode=PW(300),       # cutoff
	            kpts=(k, k, k),     # k-points
	            txt=name + '.txt')  # output file
	bulk.set_calculator(calc)
	energy = bulk.get_potential_energy()
	total_energies.append(energy)
	print('Energy for {} = '.format(fname), energy, 'eV')
pickle.dump(total_energies, open("energies", "wb"))
print("done")
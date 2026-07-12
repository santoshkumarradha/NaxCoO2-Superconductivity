import pymatgen.io.ase as pase
from ase import Atoms
from tqdm import tqdm
import pickle
import time
import numpy as np
from pathlib import Path


@contextmanager
def cd(newdir):
    prevdir = os.getcwd()
    os.chdir(os.path.expanduser(newdir))
    try:
        yield
    finally:
        os.chdir(prevdir)



###----load structures
structures=pickle.load(open("min_structures_licoo2", "rb"))
total_energies=[]

###---calculate
for i in tqdm(structures):
	start_time = time.time()
	fname = "./" + str(np.round(i.structure.lattice.c,3)) + "_caxis/" + str(d)
	Path(fname).mkdir(parents=True, exist_ok=True)
    print("running c={} ".format(np.round(i.structure.lattice.c,3)))
    with cd(fname):
        calculator = lmf(nkabc=[5, 5, 5], ctrl="temp", p=num_processors)
        pot_energy = calculator.get_potential_energy(atoms)
        #lmf().clean()

    #---time calculations
    elapsed_time = time.time() - start_time
    total_time = time.strftime("%H:%M:%S", time.gmtime(elapsed_time))
    program_time_ellapsed_unformated = time.time() - program_start_time
    program_time_ellapsed = time.strftime(
        "%H:%M:%S", time.gmtime(program_time_ellapsed_unformated))
    current_time=time.strftime(
        "%H:%M:%S", time.gmtime(time.time()))
    print("completed in {} total time {} current timr {}\n =================\n".format(
                total_time, program_time_ellapsed,current_time))
import sys, os
sys.path.append(os.getcwd())
from lmf import *
from pymatgen.io.cif import CifParser, CifWriter
from pymatgen import Structure, Lattice
from matplotlib import pyplot as plt
from pymatgen.io.ase import AseAtomsAdaptor as p2ase
import ase
from IPython.core.display import Image
from contextlib import contextmanager
from pathlib import Path
import time
import json
import os

#----for logging into file------------
import subprocess, os, sys
tee = subprocess.Popen(["tee", "log.txt"], stdin=subprocess.PIPE)
os.dup2(tee.stdin.fileno(), sys.stdout.fileno())
os.dup2(tee.stdin.fileno(), sys.stderr.fileno())


def get_structure(d=10, a=0):
    species = []
    coords = []
    for i in structure:
        if i.species_string == "Li":
            species.append(i.species)
            coords.append(i.coords + [0, 0, +.5 * d + a])
        else:
            species.append(i.species)
            coords.append(i.coords)
    a = structure.lattice.matrix.copy()
    a[2][2] = a[2][2] - d
    struc = Structure(Lattice(a), species, coords, coords_are_cartesian=True)
    #     CifWriter(struc).write_file(fname+'licoo2_'+str(struc.lattice.c)[:4]+'.cif')
    #     CifWriter(struc).write_file(fname+'licoo2_big.cif')
    #     os.system("open "+fname+'licoo2_big.cif')
    return struc


@contextmanager
def cd(newdir):
    prevdir = os.getcwd()
    os.chdir(os.path.expanduser(newdir))
    try:
        yield
    finally:
        os.chdir(prevdir)


def append_record(record, fname):
    with open(fname, 'a') as f:
        json.dump(record, f)
        f.write(os.linesep)


program_start_time = time.time()

#----------------variables--------
num_processors = 12
nc = 1
nd = 10
c_axis = np.linspace(.5, 1.2, nc)
d_move = np.linspace(0, .18, nd)
structure_data_fname = "structures_done.txt"
#--------load previous calculation and see if its already done.------------------
try:
    with open(structure_data_fname) as f:
        my_list = [json.loads(line) for line in f]
    done_structures = [Structure.from_dict(i) for i in my_list]
except:
    done_structures = []
    print("no old calculation found")

#-----------load Main structure---------------------------------------------------
parser = CifParser("structures_files/licoo2.cif")
structure = parser.get_structures()[0]

#-------------Main calculation---------------------------------------------------
for c in c_axis:
    for d in d_move:
        struc = get_structure(c, d)
        dictonary = struc.as_dict()
        atoms = p2ase().get_atoms(struc)
        if struc not in done_structures:
            start_time = time.time()
            fname = "./" + str(c) + "_caxis/" + str(d)
            Path(fname).mkdir(parents=True, exist_ok=True)
            print("running c={} and d= {}".format(c, d))
            with cd(fname):
                calculator = lmf(nkabc=[4, 4, 4], ctrl="temp", p=num_processors)
                pot_energy = calculator.get_potential_energy(atoms)
                dictonary["energy"] = pot_energy
                #lmf().clean()

            #---time calculations
            elapsed_time = time.time() - start_time
            total_time = time.strftime("%H:%M:%S", time.gmtime(elapsed_time))
            program_time_ellapsed_unformated = time.time() - program_start_time
            program_time_ellapsed = time.strftime(
                "%H:%M:%S", time.gmtime(program_time_ellapsed_unformated))

            #---load the structure and energy to completed dictonray
            append_record(dictonary, structure_data_fname)
            print("completed in {} total time {}\n =================\n".format(
                total_time, program_time_ellapsed))
        else:
            print("!!!!!!!!!!!!!!!!!")
            print("calculation for c={} and d= {} alread done \n".format(c, d))

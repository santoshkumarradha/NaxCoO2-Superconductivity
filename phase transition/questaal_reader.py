from dotmap import DotMap
from pymatgen import Structure
import numpy as np
import os
ry2ev = 13.605662285137


class reader:

	def __init__(self, output_file="output"):
		self.fname = output_file
		self.data = get_data(self.fname)
		self.iterations = make_iterations(self.data)
		# self.natoms=get_nbas(self.data)
		# self.nbas=self.natoms

		#---finding and setting up structure
		self.species = get_species(self.data)
		self.atoms = self.species
		self.structure = make_structure(self.data)

		#---setting energy data
		self.niter = len(self.iterations)
		set_iteration_energy(self.data, self.iterations)
		get_charges(self.data, self.iterations)

		#---setting energy data
		self.energy = self.iterations[-1].ehf
		self.ehf = self.iterations[-1].ehf
		self.ehk = self.iterations[-1].ehk
		#---setting band data in eV
		set_band_data(self.data, self.iterations)
		self.gap = self.iterations[-1].gap
		self.valance_band_max = self.iterations[-1].valance_band_max
		self.conduction_band_min = self.iterations[-1].conduction_band_min

	def get_variables(self):
		lst = {
			"data":
				"raw data string",
			"Iterations":
				"Iteration object, further contains other data about iterations",
			"structure":
				"calculations structure (returns in pymatgen format",
			"energy":
				"total energy in (ev) of the final iteration (ehf)",
			"ehf":
				"final energy ehf",
			"ehk":
				"final energy ehk",
			"atoms":
				"species name in lmf codes",
			"gap":
				"final band gap from given k mesh (eV)",
			"valance_band_max":
				"valance band max energy in code(eV)",
			"conduction_band_min":
				"conduction band min in code (eV)"
		}
		for i in lst:
			print(i, "-", lst[i])


def get_data(fname="output"):
	try:
		with open(fname, 'r') as f:
			data = f.read()
		return data
	except IOError:
		print("file not found")


def get_lines(data, key, num_lines=0, return_index=False):
	'''
	get the num_lines along with line matching "key" in data text
	if return_index returns the index
	'''
	index = [i for i, s in enumerate(data.splitlines()) if key in s]
	values = [data.splitlines()[i:i + num_lines + 1] for i in index]
	index_vals = [list(range(i, i + num_lines + 1)) for i in index]
	values = [list(filter(lambda name: name.strip(), i)) for i in values]
	# if return_index:
	#   return values,index_vals
	# else:
	return values


def make_structure(data):
	screen_charge,charge_dict=get_final_charge(data)

	return Structure(lattice=get_lattice(data),
					 species=get_z(data, get_nbas(data)),
					 coords=get_atomicpos(data, frac=True),
					 coords_are_cartesian=False,
					 charge=screen_charge,
					 site_properties=charge_dict)


def get_nbas(data):
	''' 
	Extract number of atoms
	'''
	try:
		nbas = int(
			get_lines(data, "nbas")[0][0].partition('nbas = ')[2].split()[0])
	except ValueError:
		print("unable to fing number of atoms")
	return nbas


def get_species(data):
	'''
	return species type for all atoms present
	'''
	natoms = get_nbas(data)
	return [
		''.join(i[0].split()[-1].split(':')[::-1]).lower()
		for i in get_lines(data, "   species  ")[:natoms]
	]


def make_iterations(data):
	'''
	make an iteration object to hold information
	'''
	iterations = []
	for i in range(get_iter(data)):
		dummy = DotMap()
		dummy.niter = i
		iterations.append(dummy)
	return iterations


def get_iter(data):
	'''
	returns total num of iterations
	'''
	return int(
		get_lines(data, "iteration",
				  return_index=False)[-1][0].split("iteration ")[-1].split()[0])


def get_z(data, natoms):
	'''get the Z's of atoms'''
	z = []
	for j in range(natoms):
		z.append(
			float(
				get_lines(data,
						  " site  {}  z=".format(j + 1))[0][0].split()[3]))
	z = np.array(z).astype(np.int)
	return z


def get_atomicpos(data, frac=True):
	'''
	data: string with output
	frac: bool returns frac coordinates 
	'''
	natoms = get_nbas(data)
	if frac:
		start = 5
		frac_pos = [
			i.split()[start:start + 3] for i in get_lines(
				data, "pos (Cartesian coordinates)", num_lines=natoms)[0][1:]
		]
		frac_pos = np.array(frac_pos).astype(np.float)
		return frac_pos
	else:
		start = 3
		cart_pos = [
			i.split()[start:start + 3] for i in get_lines(
				data, "pos (Cartesian coordinates)", num_lines=natoms)[0][1:]
		]
		cart_pos = np.array(cart_pos).astype(np.float) / 1.8897259885789
		return cart_pos


def get_lattice(data):
	'''gets the lattice vector'''
	lat = get_lines(data, "Plat", num_lines=5)
	plat = np.array([i.split()[:3] for i in lat[0][1:4]]).astype(np.float)
	alat = float(lat[0][-1].split("alat = ")[1].split()[0])
	lattice = alat * plat / 1.8897259885789
	return lattice


def get_energy(data):
	'''gets the energy at each iteration'''
	ehf = [i[0].split()[2].split("=")[-1] for i in get_lines(data, " nit=")]
	ehk = [i[0].split()[3].split("=")[-1] for i in get_lines(data, " nit=")]
	ehf = np.array(ehf).astype(np.float) * ry2ev
	ehk = np.array(ehk).astype(np.float) * ry2ev
	return ehf, ehk


def set_iteration_energy(data, iterations):
	'''set energy to iteration object'''
	ehf, ehk = get_energy(data)
	for j in range(len(iterations)):
		iterations[j].ehf = ehf[j]
		iterations[j].ehk = ehf[j]
		iterations[j].energy = ehf[j]


def get_charges(data, iterations=None):
	'''
	get charge data for each iterations
	'''
	if iterations is None:
		iterations = make_iterations(data)
	for j in range(len(iterations)):
		iter_i = iterations[j].niter
		key = "charges:       old"
		natoms = get_nbas(data)
		iter_data_txt = get_lines(data, key, natoms + 1,
								  return_index=False)[iter_i]
		natoms = get_nbas(data)
		species = get_species(data)
		iteration_charge_data = DotMap()
		charge = DotMap()
		niter = get_iter(data)
		for i in range(natoms + 1):
			line_data = iter_data_txt[1:][i].split()
			if i == 0:
				name = "smooth"
				tmp = 0
			else:
				name = species[i - 1]
				tmp = 1
			iteration_charge_data[name].old_charge = float(line_data[tmp + 1])
			iteration_charge_data[name].new_charge = float(line_data[tmp + 2])
			iteration_charge_data[name].screened_charge = float(line_data[tmp +
																		  3])
			iteration_charge_data[name].rms_charge = float(line_data[tmp + 4])
			iteration_charge_data[name].diff_charge = float(line_data[tmp + 5])
			iteration_charge_data[name].charge = float(line_data[tmp + 2])
		iterations[j].charge = iteration_charge_data

	if iterations is None:
		return iterations

def get_final_charge(data):
	natoms=get_nbas(data)
	key = "charges:       old"
	natoms = get_nbas(data)
	charge_data= get_lines(data, key, natoms + 1,
							  return_index=False)[-1]
	screen_charge=float(charge_data[1:][0].split()[2])
	final_charges=[float(i.split()[3]) for i in charge_data[1:][1:]]
	charge_dict={"charge":final_charges}
	return screen_charge,charge_dict

def get_band_data(data):
	'''get band gap data'''
	vals = [[i[0].split()[2], i[0].split()[5], i[0].split()[11]]
			for i in get_lines(data, "gap")]
	vals = np.array(vals).astype(float)
	valance_band_max = vals.T[0]
	conduction_band_min = vals.T[1]
	gap = vals.T[2]
	return valance_band_max, conduction_band_min, gap


def set_band_data(data, iterations):
	'''set to iterations data'''
	valance_band_max, conduction_band_min, gap = get_band_data(data)
	for i in range(len(iterations)):
		iterations[i].valance_band_max = valance_band_max[i] * ry2ev
		iterations[i].conduction_band_min = conduction_band_min[i] * ry2ev
		iterations[i].gap = gap[i]

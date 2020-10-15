from dotmap import DotMap
import numpy as np


class questaal_data():
	def __init__(self,output_file="output"):
		self.fname=output_file
		self.data=self.get_data()
		self.iterations=self.make_iterations(self.data)
		# self.get_charges(self.data,self.iterations)


	def get_data(self):
		with open(self.fname, 'r') as f:
			data = f.read()
		return data

	def get_lines(data,key,num_lines=0,return_index=False):
		'''
		get the num_lines along with line matching "key" in data text
		if return_index returns the index
		'''
		index=[i for i, s in enumerate(data.splitlines()) if key in s]
		values=[data.splitlines()[i:i+num_lines+1] for i in index]
		index_vals=[list(range(i,i+num_lines+1)) for i in index]
		values=[list(filter(lambda name: name.strip(), i)) for i in values]
		if return_index: 
			return values,index_vals
		else: 
			return values

	def get_nbas(data):
		''' 
		Extract number of atoms
		'''
		try:
			nbas=int(get_lines(data,"nbas")[0][0].partition('nbas = ')[2].split()[0])
		except ValueError:
			print("unable to fing number of atoms")
		return nbas
	def get_iter(data):
		'''
		returns total num of iterations
		'''
		return int(get_lines(data,"iteration",return_index=False)[-1][0].split("iteration ")[-1].split()[0])
	def get_species(data):
		'''
		return species type for all atoms present
		'''
		natoms=get_nbas(data)
		return [''.join(i[0].split()[-1].split(':')[::-1]).lower() for i in get_lines(data,"   species  ")[:natoms]]


	def make_iterations(data):
		'''
		make an iteration object to hold information
		'''
		iterations=[]
		for i in range(get_iter(data)):
			dummy=DotMap()
			dummy.niter=i
			iterations.append(dummy)
		return iterations

	def get_charges(data,iterations=None):
		'''
		get charge data for each iterations
		'''
		if iterations is None:
			iterations=make_iterations(data)
		for j in range(len(iterations)):
			iter_i=iterations[j].niter
			key="charges:       old"
			natoms=get_nbas(data)
			iter_data_txt=get_lines(data,key,natoms+1,return_index=False)[iter_i]
			natoms=get_nbas(data)
			species=get_species(data)
			iteration_charge_data=DotMap()
			charge=DotMap()
			niter=get_iter(data)
			for i in range(natoms+1):
				line_data=iter_data_txt[1:][i].split()
				if i==0:
					name="smooth"
					tmp=0
				else:name=species[i-1];tmp=1
				iteration_charge_data[name].old_charge=float(line_data[tmp+1])
				iteration_charge_data[name].new_charge=float(line_data[tmp+2])
				iteration_charge_data[name].screened_charge=float(line_data[tmp+3])
				iteration_charge_data[name].rms_charge=float(line_data[tmp+4])
				iteration_charge_data[name].diff_charge=float(line_data[tmp+5])
				iteration_charge_data[name].charge=float(line_data[tmp+2])
			iterations[j].charge=iteration_charge_data
		if iterations is None:
			return iterations

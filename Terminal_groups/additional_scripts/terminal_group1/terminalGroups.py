''' this program computes the terminal groups for a range of structures.
	groups include OH, methyl group, vinyl group, phenyl group '''

import copy
import numpy
import os

#-----------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------#
# utility functions 

def Distance(dx, dy, dz):
	return (dx*dx + dy*dy + dz*dz)**0.5

def AverageMax(data, N):
	# data must be a numpy array
	data = numpy.sort(data)
	return numpy.mean(data[-N:])

def AverageMin(data, N):
	# data must be a numpy array
	data = numpy.sort(data)
	return numpy.mean(data[:N])

def ConvertAtomType(atomtype):
	convertTypes = {1: 'Si', 2: 'O', 3: 'C'}
	return convertTypes[atomtype]

#-----------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------#

class TerminalGroups:

	def __init__(self, inputfile, moltype):
		# data attributes of the network 
		self.Natoms = 0 # number of total atoms
		self.NSi = 0 # number of Si atoms 
		self.NC = 0 # number of C atoms
		self.NO = 0 # number of O atoms
		self.Lx = 0 # cell length in x dimension
		self.Ly = 0 # cell length in y dimension
		self.Lz = 0	# cell length in z dimension

		# analyze the network and compute the terminal groups 
		self.data = self.LoadData(inputfile)
		print 'Data loaded.'

		print 'Computing Si-O bonds...'
		AM_SiO = self.ComputeAMSiO(self.data)

		print 'Computing Si-C and C-C bonds...'
		AM_SiC = self.ComputeAMSiC(self.data)

		print 'Finding precursors...'
		precursors, histogram = self.GetPrecursors(AM_SiC)

		print 'Finding terminal groups...'
		self.terminalOH = self.ComputeTerminalOH(AM_SiO)
		self.terminalGroups = self.ComputeTerminalGroups(AM_SiC, precursors, moltype) 

		# get the number of bonds 
		self.Nbonds = self.ComputeNbonds(AM_SiO, AM_SiC)
		self.Nbonds_term_min,\
		self.Nbonds_term_max = self.ComputeNbondsTerminalGroups(self.terminalGroups, moltype)

		# rename the inputfile
		inputfile = self.renameFile(inputfile)

		print 'Writing terminal groups to a file...'
		self.WriteTerminalOH(self.data, inputfile)
		self.WriteTerminalGroups(self.data, inputfile)
		self.WriteTerminalGroupStats(precursors, histogram, inputfile)


	# set getr methods to retrieve the data from the class
	def getTerminalGroups(self):
		# send out the data 
		return (self.terminalOH, self.terminalGroups)

	def getCellDimensions(self):
		return (self.Lx, self.Ly, self.Lz)

	def getAtomData(self):
		return self.data

	def getNbonds(self):
		return (self.Nbonds, self.Nbonds_term_min, self.Nbonds_term_max)

	#------------------------------------------------------------------------#
	def renameFile(self, inputfile):
		# set a custom directory for the inputfile 
		base_dir = '%s/' %inputfile[:-4]
		if not os.path.exists(base_dir):
			os.makedirs(base_dir)

		# rename the inputfile to store the results in a common directory
		dir_name = '%sterminal_groups/' %base_dir
		if not os.path.exists(dir_name):
			os.makedirs(dir_name)
		inputfile = '{}{}'.format(dir_name, inputfile)
		return inputfile

	#------------------------------------------------------------------------#
	def LoadData(self,inputfile):
		""" reads LAMMPS .xyz file, sorts data, and computes cell dimensions """ 
		f = open(inputfile)
		self.Natoms = int(f.readline())

		f.readline()

		data = []
		while True:
			fields = f.readline().strip().split()
			if fields:
				atomtype = int(fields[0])
				xcoord = float(fields[1])
				ycoord = float(fields[2])
				zcoord = float(fields[3])

				# determine type of atom
				if atomtype == 1:
					self.NSi += 1
				elif atomtype == 2:
					self.NO += 1
				elif atomtype == 3:
					self.NC += 1
				else: 
					raise RuntimeError, "Incorrect atom type."

				# populate the data array
				data.append([atomtype,xcoord,ycoord,zcoord])

			else:
				break
		f.close()

		data = numpy.array(data)
		b = numpy.zeros((self.Natoms,4))

		# sort the atoms by Si, O, C
		# decreases computation time knowing which types of atoms will bond
		Siidx = 0
		Oidx = self.NSi
		Cidx = self.NSi + self.NO
		for i in range(self.Natoms):
			if data[i][0] == 1:
				b[Siidx] = data[i]
				Siidx += 1
			elif data[i][0] == 2:
				b[Oidx] = data[i]
				Oidx += 1
			elif data[i][0] == 3:
				b[Cidx] = data[i]
				Cidx += 1

		data = b

		# find the dimensions of the simulation cell
		# average the smallest 50 points and largest 50 points
		x_coords = data[:,1]
		y_coords = data[:,2]
		z_coords = data[:,3]

		self.Lx = AverageMax(x_coords, 50) - AverageMin(x_coords, 50)
		self.Ly = AverageMax(y_coords, 50) - AverageMin(y_coords, 50)
		self.Lz = AverageMax(z_coords, 50) - AverageMin(z_coords, 50)

		return data

	#------------------------------------------------------------------------#
	# methods to compute adjacency matricies for given sets of atoms
	def AMSiCtoData(self, idx):
		# mapping between AMSiO idx and data idx
		# order of data array is NSi, NO, NC 
		if idx < self.NSi:
			return idx
		elif idx >= self.NSi:
			# idx already contains NSi so shift by NO to get to C IDs
			return idx + self.NO

	def AMSiOtoData(self, idx):
		# mapping between AMSiO idx and data idx
		# order of data array is NSi, NO, NC 
		# shift idx by NSi to O IDs
		return idx + self.NSi


	def ComputeAMSiO(self, data):
		""" Computes all Si-O bonds and populates an adjacency matrix""" 
		# adjacency matrix for Si-O bonding 
		AM_SiO = numpy.zeros((self.NSi,self.NO))

		# note that the matrix is not symmetric
		symmetric = False

		# set cutoff to be 2.3 A
		cutoff = 2.3

		for i in range(self.NSi):
			for j in range(self.NSi, self.NSi+self.NO):

				dx = abs(data[j][1] - data[i][1])
				dy = abs(data[j][2] - data[i][2])
				dz = abs(data[j][3] - data[i][3])

				d = Distance(dx,dy,dz)
				if d < cutoff:
					AM_SiO[i][j-self.NSi] = d
					continue
				
				# consider PBCs
				if dx > self.Lx - cutoff:
					d = Distance(dx + self.Lx, dy, dz)
					if d < cutoff:
						AM_SiO[i][j-self.NSi] = d
						continue

					d = Distance(dx - self.Lx, dy, dz)
					if d < cutoff:
						AM_SiO[i][j-self.NSi] = d
						continue

				if dy > self.Ly - cutoff:
					d = Distance(dx ,dy + self.Ly, dz)
					if d < cutoff:
						AM_SiO[i][j-self.NSi] = d
						continue

					d = Distance(dx, dy - self.Ly, dz)
					if d < cutoff:
						AM_SiO[i][j-self.NSi] = d
						continue

				if dz > self.Lz - cutoff:
					d = Distance(dx, dy, dz + self.Lz)
					if d < cutoff:
						AM_SiO[i][j-self.NSi] = d
						continue

					d = Distance(dx, dy, dz - self.Lz)
					if d < cutoff:
						AM_SiO[i][j-self.NSi] = d
		return AM_SiO


	def ComputeAMSiC(self, data):
		""" Computes all Si-C and C-C bonds and populates an adjacency matrix""" 
		# symmetric array is arranged by (Si C) x (Si C) 
		AM_SiC = numpy.zeros((self.NSi+self.NC, self.NSi+self.NC))
		symmetric = True

		# set Si-C cutoff at 2.3 and C-C cutoff at 2.3
		cutoff = 2.3

		for i in range(self.NSi+self.NC):
			for j in range(i+1, self.NSi+self.NC):

				# get mapping between AMSiC id and data array id 
				idx_i = self.AMSiCtoData(i)
				idx_j = self.AMSiCtoData(j)

				# no Si-Si bonds so check atom type and skip if necessary
				# 1 = Si, 2 = O, 3 = C
				if (data[idx_i][0]==1) and (data[idx_j][0]==1):
					continue

				dx = abs(data[idx_j][1] - data[idx_i][1])
				dy = abs(data[idx_j][2] - data[idx_i][2])
				dz = abs(data[idx_j][3] - data[idx_i][3])

				d = Distance(dx,dy,dz)
				if d < cutoff:
					AM_SiC[i][j] = d
					AM_SiC[j][i] = d
					continue

				# consider PBCs
				if dx > self.Lx - cutoff:
					d = Distance(dx + self.Lx, dy, dz)
					if d < cutoff:
						AM_SiC[i][j] = d
						AM_SiC[j][i] = d
						continue

					d = Distance(dx - self.Lx, dy, dz)
					if d < cutoff:
						AM_SiC[i][j] = d
						AM_SiC[j][i] = d
						continue

				if dy > self.Ly - cutoff:
					d = Distance(dx ,dy + self.Ly, dz)
					if d < cutoff:
						AM_SiC[i][j] = d
						AM_SiC[j][i] = d
						continue

					d = Distance(dx, dy - self.Ly, dz)
					if d < cutoff:
						AM_SiC[i][j] = d
						AM_SiC[j][i] = d
						continue

				if dz > self.Lz - cutoff:
					d = Distance(dx, dy, dz + self.Lz)
					if d < cutoff:
						AM_SiC[i][j] = d
						AM_SiC[j][i] = d
						continue

					d = Distance(dx, dy, dz - self.Lz)
					if d < cutoff:
						AM_SiC[i][j] = d
						AM_SiC[j][i] = d
		return AM_SiC


	def ComputeNbondsTerminalGroups(self, terminalGroups, moltype):
		# get the number of terminal groups
		Nterm_groups = len(terminalGroups)

		# get the size of each terminal group
		Natoms_term_group = len(terminalGroups[0])

		# each precursor has a terminal group with a specific number of bonds
		# we use the number of atoms that are bonded to correct for over counting
		# of steric interactions within a cluster. 
		# set the minimum and maximum interactions of the terminal group

		# minimum interactions is n-1
		Nints_term_group_min = Natoms_term_group - 1
		Nbonds_term_min = Nterm_groups*Nints_term_group_min

		# Nbonds_term_min = {	'EtOCSMethyl' : 0*Nterm_groups,
		# 					'EtOCSVinyl'  : 1*Nterm_groups,
		# 					'EtOCSPhenyl' : 5*Nterm_groups}

		# maximum interactions is n*(n-1)/2
		Nints_term_group_max = Natoms_term_group*(Natoms_term_group - 1)/2.0
		Nbonds_term_max = Nterm_groups*Nints_term_group_max

		# Nbonds_term_max = {	'EtOCSMethyl' : 0*Nterm_groups,
		# 					'EtOCSVinyl'  : 1*Nterm_groups,
		# 					'EtOCSPhenyl' : 15*Nterm_groups}

		return (Nbonds_term_min, Nbonds_term_max)


	def ComputeNbonds(self, AM_SiO, AM_SiC):
		# get the total number of bonds in the system
		Nbonds = 0

		# note: nonzero() gives the indicies in a flattened version
		# of the input array -> index 0 is i-coords, index 1 is j-coords
		N_SiO_bonds = len(numpy.nonzero(AM_SiO)[0])

		# AM_SiC is a symmetric matrix so each bond is repeated
		# i.e divide by 2 
		N_SiC_CC_bonds = len(numpy.nonzero(AM_SiC)[0])/2.0

		Nbonds = N_SiO_bonds + N_SiC_CC_bonds
		return Nbonds

	# def deltaCoordinatePBC(self, data, i, j, coord):
	# 	dx, dy, dz = self.deltaCoordinate(data, i, j)

	# 	deltaCoords = { 'x':  (dx + self.Lx, dy, dz),
	# 					'x-': (dx - self.Lx, dy, dz),
	# 					'y':  (dx, dy + self.Ly, dz),
	# 					'y-': (dx, dy - self.Ly, dz),
	# 					'z':  (dx, dy, dz + self.Lz),
	# 					'z-': (dx, dy, dz - self.Lz)}

	# 	return deltaCoords[coord]


	# def deltaCoordinate(self, data, i, j):
	# 	dx = abs(data[j][1] - data[i][1])
	# 	dy = abs(data[j][2] - data[i][2])
	# 	dz = abs(data[j][3] - data[i][3])	
	# 	return (dx, dy, dz)	

	# def updateAM(self, d, i, j, AM, symmetric):
	# 	AM[i][j] = d
	# 	if symmetric:
	# 		AM[j][i] = d
	# 	return AM

	# def nearBoundary(self, dx, dy, dz, cutoff, coord):
	# 	cellLengths = { 'x' : self.Lx,
	# 					'y' : self.Ly,
	# 					'z' : self.Lz}
	# 	deltacoord = {  'x' : dx,
	# 					'y' : dy,
	# 					'z' : dz}

	# 	if deltacoord[coord] > cellLengths[coord] - cutoff:
	# 		return True
	# 	return False


	# def BondPBC(self, cutoff, data, i, j, coord):
	# 	# try adding L to the change in coordinate (e.g. dx)
	# 	dx, dy, dz = self.deltaCoordinatePBC(data, i, j, coord)
	# 	d = Distance(dx, dy, dz)
	# 	if d < cutoff:
	# 		return d

	# 	# try subtracting L from the change in coordinate (e.g. dx)
	# 	coord = coord + '-'
	# 	dx, dy, dz = self.deltaCoordinatePBC(data, i, j, coord)
	# 	d = Distance(dx, dy, dz)
	# 	if d < cutoff:
	# 		return d

	# 	# if no bonds found, return 0  
	# 	return 0


	# def addBond(self, i, j, idx_i, idx_j, cutoff, AM, data, symmetric):
	# 	# if the atoms are within a given cutoff distance, then they are bonded 
	# 	dx, dy, dz = self.deltaCoordinate(data, idx_i, idx_j)
	# 	d = Distance(dx, dy, dz)

	# 	if d < cutoff:
	# 		return self.updateAM(d, i, j, AM, symmetric)

	# 	# consider PBC if bond not previously found
	# 	for coord in ['x', 'y', 'z']: 
	# 		if self.nearBoundary(dx, dy, dz, cutoff, coord):
	# 			d = self.BondPBC(cutoff, data, idx_i, idx_j, coord)
	# 			if d:
	# 				return self.updateAM(d, i, j, AM, symmetric)

	# 	# if no bonds found, return the original AM
	# 	return AM


	# def ComputeAMSiO(self, data):
	# 	""" Computes all Si-O bonds and populates an adjacency matrix""" 
	# 	# adjacency matrix for Si-O bonding 
	# 	AM_SiO = numpy.zeros((self.NSi,self.NO))

	# 	# note that the matrix is not symmetric
	# 	symmetric = False

	# 	# set cutoff to be 2.3 A
	# 	cutoff = 2.3

	# 	for i in range(self.NSi):
	# 		for j in range(self.NO):
	# 			# get mapping between AMSiO id and data array id 
	# 			idx_i = i
	# 			idx_j = self.AMSiOtoData(j)

	# 			# update the AM if there is a bond between the atoms
	# 			AM_SiO = self.addBond(i, j, idx_i, idx_j, cutoff, AM_SiO, data, symmetric)

	# 	return AM_SiO

	# def ComputeAMSiC(self, data):
	# 	""" Computes all Si-C and C-C bonds and populates an adjacency matrix""" 
	# 	# symmetric array is arranged by (Si C) x (Si C) 
	# 	AM_SiC = numpy.zeros((self.NSi+self.NC, self.NSi+self.NC))
	# 	symmetric = True

	# 	# set Si-C cutoff at 2.3 and C-C cutoff at 2.3
	# 	cutoff = 2.3

	# 	for i in range(self.NSi+self.NC):
	# 		for j in range(i+1, self.NSi+self.NC):

	# 			# get mapping between AMSiC id and data array id 
	# 			idx_i = self.AMSiCtoData(i)
	# 			idx_j = self.AMSiCtoData(j)

	# 			# no Si-Si bonds so check atom type and skip these interactions
	# 			# 1 = Si, 2 = O, 3 = C
	# 			if (data[idx_i][0]==1) and (data[idx_j][0]==1):
	# 				continue

	# 			# update the AM if there is a bond between the atoms
	# 			AM_SiC = self.addBond(i, j, idx_i, idx_j, cutoff, AM_SiC, data, symmetric)

	# 	return AM_SiC

	#------------------------------------------------------------------------#
	# get the atoms involved in each precursor 

	def addAtoms(self, bondedTo, precursor):
		precursor_updated = copy.deepcopy(precursor)
		for atom in bondedTo:
			# python set will not repeat atom id
			# so no need to check for unique atoms
			precursor_updated.add(atom)
		return precursor_updated


	def generateHistogram(self, precursors):
		precursor_sizes = []
		for precursor in precursors:
			precursor_sizes.append(len(precursor))

		precursor_sizes = numpy.array(precursor_sizes)
		max_size = max(precursor_sizes)

		bins = numpy.linspace(0,max_size,max_size+1)
		histogram = numpy.histogram(precursor_sizes,bins=bins,normed=False)
		return histogram		

	def GetPrecursors(self, AM_SiC):
		""" Compute the precursors given the AM """ 
		# store list of atom IDs in each precursors
		precursors = []

		for i in range(len(AM_SiC)):
			# find all atoms in the precurosr 
			precursor = set([])

			# seed the precursor with the list of atoms that current atom is bonded to
			bondedTo = numpy.nonzero(AM_SiC[i])[0]
			precursor = self.addAtoms(bondedTo, precursor)

			while True:
				precursor_updated = set([])
				for atom in precursor:
					bondedTo = numpy.nonzero(AM_SiC[atom])[0]
					new_atoms = self.addAtoms(bondedTo, precursor)
					precursor_updated = precursor_updated.union(new_atoms)

				# if no additional atoms found, break
				if precursor_updated == precursor:
					break
				else:
					precursor = precursor_updated

			# update precursor list
			precursors.append(precursor)

		# clean list of precursors to contain only uniuqe precursors
		unique_precursors = []
		for precursor in precursors:
			if len(precursor) == 0:
				continue 
			if precursor not in unique_precursors:
				# note: no need to convert IDs to data since only interested in Si atoms
				# which have the same IDs
				unique_precursors.append(precursor)

		# generate histogram of number of atoms per precursor
		histogram = self.generateHistogram(unique_precursors)

		return (unique_precursors, histogram)

	#------------------------------------------------------------------------#
	# compute the terminal groups 

	def ComputeTerminalOH(self, AM_SiO):
		""" Compute the terminal OH groups given the AM""" 
		terminalOH = []

		for i in range(self.NO):
			Nbonds = len(numpy.nonzero(AM_SiO[:,i])[0])

			# print numpy.nonzero(AM_SiO[:,i])[0]
			# print Nbonds

			if Nbonds == 1: # terminal O
				# convert the index to the data idx 
				idx = self.AMSiOtoData(i)
				terminalOH.append(idx)
		return terminalOH


	def ComputeTerminalMethyl(self, AM_SiC, precursors):
		# strategy: all the C atoms are in a bonded precursor; 
		# if the C atom is only bonded to one other atoms, 
		# then it's a terminal methyl group (H atoms modeld implicitly)
		terminalMethyl = []

		# organsilicate with terminal methyl has C with one bond (H modeled implicitly)
		for i in range(self.NSi, self.NSi+self.NC):
			Nbonds = len(numpy.nonzero(AM_SiC[:,i])[0])

			# print numpy.nonzero(AM_SiC[:,i])[0]
			# print Nbonds

			if Nbonds == 1:
				# convert the index to the data idx
				idx = self.AMSiCtoData(i)
				# store a list of the id to be consistent with other terminal groups
				terminalMethyl.append([idx])
		return terminalMethyl


	def ComputeTerminalVinyl(self, AM_SiC, precursors):
		# strategy: the precursors that we currently use with a vinyl group only 
		# have 1 C=C double bond (from the vinyl group). Look for the C=C bond
		# in each precursors and store the atom IDs
		terminalVinyl = []

		cutoff = 1.40 # double bonded C=C cutoff

		for precursor in precursors:
			group = set([])
			for idx in precursor:
				bonds = numpy.nonzero(AM_SiC[:,idx])[0]
				for idx_2 in bonds:
					d = AM_SiC[idx][idx_2]
					if d < cutoff:
						idx_1 = self.AMSiCtoData(idx)
						idx_2 = self.AMSiCtoData(idx_2)

						# ensure only unique groups are saved by 
						# using a python set (only unique values stored)
						group.add(idx_1)
						group.add(idx_2)

			if group:
				terminalVinyl.append(group)

		return terminalVinyl


	def ComputeTerminalPhenyl(self, AM_SiC, precursors):
		# strategy: the precursors that we currently use with phenyl groups
		# only have aromatic carbon in the phenyl group. Look for thc Car-Car
		# bonds in the precursor and store the atom IDs
		terminalPhenyl = []

		cutoff = 1.43 # double bonded C=C cutoff

		for precursor in precursors:
			aromatic_group = set([])
			for idx in precursor:
				bonds = numpy.nonzero(AM_SiC[:,idx])[0]
				for idx_2 in bonds:
					d = AM_SiC[idx][idx_2]
					if d < cutoff:
						idx_1 = self.AMSiCtoData(idx)
						idx_2 = self.AMSiCtoData(idx_2)

						# note: python sets only store unique values, 
						# so no need to check if ID is repeated
						aromatic_group.add(idx_1)
						aromatic_group.add(idx_2)

			if aromatic_group:
				terminalPhenyl.append(aromatic_group)

		return terminalPhenyl

	def ComputeTerminalGroups(self, AM_SiC, precursors, moltype):
		terminalGroups = {  'EtOCSMethyl' : self.ComputeTerminalMethyl,
							'EtOCSVinyl'  : self.ComputeTerminalVinyl,
							'EtOCSPhenyl' : self.ComputeTerminalPhenyl}

		return terminalGroups[moltype](AM_SiC, precursors)


	#------------------------------------------------------------------------#
	# write out terminal groups in VMD format 

	def WriteTerminalOH(self, data, inputfile):
		filename = '%s_terminalOH_VMD.xyz' %inputfile[:-4]
		f = open(filename, 'w')
		f.write('%d\n' %len(self.terminalOH))
		f.write('Terminal OH atoms\n')
		for atomID in self.terminalOH:
			atomtype, x, y, z = data[atomID]
			atomtype = ConvertAtomType(atomtype)
			f.write('%s  %.4f  %.4f  %.4f\n' %(atomtype, x, y, z))
		f.close()

	def WriteTerminalGroups(self, data, inputfile):
		Natoms = 0
		for group in self.terminalGroups:
			for atom in group:
				Natoms += 1

		filename = '%s_terminalGroups_VMD.xyz' %inputfile[:-4]
		f = open(filename, 'w')
		f.write('%d\n' %Natoms)
		f.write('Terminal Groups atoms\n')	
		for group in self.terminalGroups:
			for atomID in group:
				atomtype, x, y, z = data[atomID]
				atomtype = ConvertAtomType(atomtype)
				f.write('%s  %.4f  %.4f  %.4f\n' %(atomtype, x, y, z))		
		f.close()		

	# get the terminal groups stats and write them out 
	def WriteTerminalGroupStats(self, precursors, histogram, inputfile):
		filename = '%s_terminalGroup_stats.txt' %inputfile[:-4]
		f = open(filename, 'w')
		f.write('Terminal Groups Stats\n\n')
		f.write('N terminal OH = %d\n' %len(self.terminalOH))

		f.write('\nN precursors = %d\n' %len(precursors))
		f.write('Precursor histogram:\n')
		hist = histogram[0]
		bins = histogram[1][1:]
		for i in range(len(hist)):
			f.write('%d, %d\n' %(bins[i], hist[i]))

		f.write('\nN terminal groups (methyl, vinyl, phenyl) = %d\n' %len(self.terminalGroups))
		f.write('Terminal groups histogram:\n')
		term_groups_histogram = self.generateHistogram(self.terminalGroups)
		hist = term_groups_histogram[0]
		bins = term_groups_histogram[1][1:]
		for i in range(len(hist)):
			f.write('%d, %d\n' %(bins[i], hist[i]))
		
		f.close()


#-----------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------#

# import glob
# import sys
# import time


# if len(sys.argv) < 2:
# 	print 'Usage:'
# 	print '  python %s <precursor type>' %sys.argv[0]
# 	exit()

# moltype = sys.argv[1]

# t0 = time.time()
		
# # get all the relevant files and process each network
# inputfiles = glob.glob('{}_*.xyz'.format(moltype))
# for inputfile in inputfiles:

# 	print '\nWorking with %s...' %inputfile
# 	TerminalGroups(inputfile, moltype)

# print 'Analyzed network in %.4f seconds.' %(time.time()-t0)

''' This class determines the min-cut of a material given a LAMMPS .xyz file '''

import copy
import numpy
import scipy.sparse
import time

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

class MinCut:

	def __init__(self,inputfile,molType):

		if molType not in set(['EtOCS', 'EtOCSMe', '135Benz','SiO2']):
			raise RuntimeError, 'Invalid precursor type.'
		else:
			self.molType = molType # precursor type: Et-OCS, Et-OCS(Me), 135-Benzene, SiO2

		self.Natoms = 0
		self.NSi = 0
		self.NC = 0
		self.NO = 0

		self.Lx = 0 
		self.Ly = 0
		self.Lz = 0
		self.y_avg = 0 

		self.E_SiC_bond = 5000001 # weight of Si-C bonds
		self.E_SiO_bond = 6000100 # weight of Si-O bonds
		self.E_CC_bond = 7000001 # weight of C-C bonds
		self.E_bd = 100000000 # weight of network bonds and atoms in source and sink

		self.AM = [] # adjacency matrix
		self.Clist = [] # list of C atoms
		self.Olist = [] # list of O atoms
		self.mollist = [] # list of molecules

		self.network = [] # idx of network bonds for calculating fracture energy
		self.source = [] # idx of atoms in upper side clamp
		self.sink = [] # idx of atoms in lower side clamp

		self.Nnetwork = 0 # number of atoms in network
		self.Nsource = 0 # number of atoms in source
		self.Nsink = 0 # number of atoms in sink

		delta = 15 # min-cut height

		# self.O_stats = numpy.zeros(6) # freeO, non-bridging O, bridging O 
		# self.Si_stats = numpy.zeros(6) # T0, T1, T2, T3
		# self.Si_statsT = numpy.zeros(6) # only considers bridging O 

		# self.p = 0 # Si-X-Si network connectivity
		# self.q = 0 # condensation degree
		# self.mSi = 0 # mean Si network connectivity
		# self.rho = 0 # density 

		# self.data = []
		# self.CM = [] # Si-O bonding matrix
		# self.CM_T = [] # only considers bridging O

		# analysis of network 
		self.LoadData(inputfile)
		print 'Data loaded.'

		print 'Computing the adjacency matrix...'
		self.ComputeAM()

		# print 'Analyzing network...'
		# self.AnalyzeConnectivity()
		# self.ComputeDensity() # must compute density after connecitivty
		# self.WriteSolution(inputfile)

	#------------------------------------------------------------------------------#
	''' utility methods ''' 

	def Distance(self,dx,dy,dz):
		return (dx*dx + dy*dy + dz*dz)**0.5

	#------------------------------------------------------------------------------#

	def LoadData(self,inputfile):
		""" reads LAMMPS file, sorts data, and computes cell dimensions """ 
		f = open(inputfile)
		self.Natoms = int(f.readline())

		f.readline()

		self.data = numpy.zeros((self.Natoms,4))

		for i in range(self.Natoms):
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
				self.data[i] = [atomtype,xcoord,ycoord,zcoord]
		f.close()

		# sort the rows by atom type: Si (1), O (2), then C (3)
		self.data = numpy.array(sorted(self.data, key=lambda a_entry: a_entry[0]))

		# find the dimensions of the simulation cell 
		self.Lx = max(self.data[:,1]) - min(self.data[:,1])
		self.Lz = max(self.data[:,3]) - min(self.data[:,3])

		y_min = min(self.data[:,2])
		y_max = max(self.data[:,2])
		self.Ly = self.y_max - self.y_min
		self.y_avg = 0.5*(y_max + y_min)
		
	#------------------------------------------------------------------------------#

	def isBond(self,i,j):
		dx = abs(self.data[j][1] - self.data[i][1])
		dy = abs(self.data[j][2] - self.data[i][2])
		dz = abs(self.data[j][3] - self.data[i][3])

		if self.Distance(dx,dy,dz) < 2.0:
			return True
		
		# consider PBCs in x and z directions; not in y-dir due to source and sink
		if dx > self.Lx - 2.0:
			if self.Distance(dx + self.Lx, dy, dz) < 2.0:
				return True

			if self.Distance(dx - self.Lx, dy, dz) < 2.0:
				return True

		if dz > self.Lz - 2.0:
			if self.Distance(dx, dy, dz + self.Lz) < 2.0:
				return True

			if self.Distance(dx, dy, dz - self.Lz) < 2.0:
				return True

		return False


	def updateList(self,i,j,bond_energy,row,col,val):
		row.append(i)
		col.append(j)
		val.append(bond_energy)

		row.append(j)
		col.append(i)
		val.append(bond_energy)

		return (row,col,val)


	def ComputeAM(self):
		'''
		AM is a symmetric sparse block matrix consisting of all bonds
		
		      NSi    NO    NC 
		NSi |  0  |     |     |
		NO  |     |  0  |  0  | 
		NC  |     |  0  |     |

		only Si-O, Si-C and C-C bonds are possible
		'''

		# an AM array is too large for memory so use COO sparse matrix format
		row = []
		col = []
		val = []

		# first compute the Si-O bonds
		for i in range(self.NSi):
			Nbonds = 0
			for j in range(self.NSi,self.NSi+self.NO):

				if self.isBond(i,j):
					row, col, val = self.updateList(i,j,self.E_SiO_bond,row,col,val)
					Nbonds += 1

				# can have a maximum of 3 Si-O bonds
				if Nbonds == 3:
					break

		# remove the non-bonded and non-bridging (terminal) O atoms
		# col_idx must appear at least 2 times for briding O
		row_bridge = []
		col_bridge = []
		val_bridge = []

		Nbridging = 0
		for i in range(len(col)):
			Nbonds = col.count(col[i])

			if Nbonds >= 2:
				row.append(row[i])
				col.append(col[i])
				val.append(val[i])

				Nbridging += 1

		row = row_bridge
		col = col_bridge
		val = val_bridge

		# compute the Si-C bonds
		for i in range(self.NSi):
			for j in range(self.NSi+self.NO,self.Natoms):

				if self.isBond(i,j):
					row, col, val = self.updateList(i,j,self.E_SiC_bond,row,col,val)
					break # only 1 possible Si-C bond

		# compute the C-C bonds
		for i in range(self.NSi+self.NO,self.Natoms):
			for j in range(i+1,self.Natoms):

				if self.isBond(i,j):
					row, col, val = self.updateList(i,j,self.E_CC_bond,row,col,val)
					break # only 1 possible C-C bond

		# convert AM to sparse coo format
		row = numpy.array(row)
		col = numpy.array(col)
		val = numpy.array(val)
		self.AM = scipy.sparse.coo_matrix((val,(row,col)),shape=(self.Natoms,self.Natoms))


	#------------------------------------------------------------------------------#

	def formFlowNetwork(self,delta):
		for i in range(self.Natoms):
			ycoord = self.data[i][2]

			# the source is all atoms above yavg + delta/2 
			if ycoord > self.y_avg + delta*0.5:
				self.source.append(i)

			# the sink is all atoms below yavg - delta/2 
			elif ycoord < self.y_avg - delta*0.5:
				self.sink.append(i)

			# all other atoms are within the network
			else:
				self.network.append(i)

		self.Nsource = len(self.source)
		self.Nsink = len(self.sink)
		self.Nnetwork = len(self.network)







"""This class determines the network connectivity of the material given a LAMMPS .xyz file"""

import copy
import numpy

from params import ValidMolTypes, CondensationDegree, Connectivity, MeanSiConnectivity,\
					MassPrecursor, Nprecursors, Mass

#-----------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------#
class Network:

	def __init__(self, inputfile, PrecursorParams, q_relative):

		# validate the precursor type
		molType = PrecursorParams['molType']
		if molType not in ValidMolTypes:
			raise RuntimeError, 'Invalid precursor type.'
		else:
			self.molType = molType 

		self.Natoms = 0
		self.NSi = 0
		self.NC = 0
		self.NO = 0

		self.Lx = 0
		self.Ly = 0
		self.Lz = 0

		self.O_stats = numpy.zeros(6) # freeO, non-bridging O, bridging O 
		self.Si_stats = numpy.zeros(6) # Si atom coordination with O (CN1, CN2, CN3)
		self.Si_statsT = numpy.zeros(6) # only considers bridging O (T0, T1, T2, T3) 
		self.Si_statsT_rel = numpy.zeros(6) # removes influence of bridging O within a precursor

		self.p = 0 # Si-X-Si network connectivity
		self.q = 0 # condensation degree
		self.mSi = 0 # mean Si network connectivity
		self.p_rel = 0 # relative network connectivity
		self.q_rel = 0 # relative condensation degree which accounts for Si-O-Si bonds in the precursors
		self.mSi_rel = 0 # relative mean Si network connectivity
		self.rho = 0 # density 

		self.data = []
		self.CM = [] # Si-O bonding matrix
		self.CM_T = [] # Si-O bonding that only considers bridging O

		self.AMSiC = [] # adjacency matrix for Si-C and C-C bonding 
		self.precursors = [] # list of sets with atomIDs in each precursor
		self.precursor_histogram = [] # histogram of the precursor size
		self.Nprec_w_bridging = 0 # number of precursors with bridging Si-O-Si

		# if a dual precursor model is being used, set the appropriate parameters
		if PrecursorParams['molType1']:
			self.mixed = True
			self.molType1 = PrecursorParams['molType1']
			self.molType2 = PrecursorParams['molType2']
			self.Nprec1 = PrecursorParams['Nprecursors1']
			self.Nprec2 = PrecursorParams['Nprecursors2']
		else:
			self.mixed = False

		# analysis of network 
		self.LoadData(inputfile)
		print 'Data loaded.'

		print 'Computing Si-O bonds...'
		self.ComputeCM()

		if q_relative:
			print 'Computing Si-C and C-C bonds...'
			self.ComputeAMSiC()
			print 'Finding precursors...'
			self.GetPrecursors()

		print 'Analyzing network...'
		self.AnalyzeConnectivity()
		self.ComputeDensity() # must compute density after connecitivty
		self.WriteSolution(inputfile)

	def LoadData(self,inputfile):
		""" reads LAMMPS file, sorts data, and computes cell dimensions """ 
		f = open(inputfile)
		self.Natoms = int(f.readline())

		f.readline()

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
				self.data.append([atomtype,xcoord,ycoord,zcoord])

			else:
				break
		f.close()

		self.data = numpy.array(self.data)
		b = numpy.zeros((self.Natoms,4))

		# sort the atoms by Si, O, C
		# decreases computation time knowing which types of atoms will bond
		Siidx = 0
		Oidx = self.NSi
		Cidx = self.NSi + self.NO
		for i in range(self.Natoms):
			if self.data[i][0] == 1:
				b[Siidx] = self.data[i]
				Siidx += 1
			elif self.data[i][0] == 2:
				b[Oidx] = self.data[i]
				Oidx += 1
			elif self.data[i][0] == 3:
				b[Cidx] = self.data[i]
				Cidx += 1

		self.data = b

		# find the dimensions of the simulation cell
		# average the lowest 50 points and largest 50 points
		x_coords = self.data[:,1]
		y_coords = self.data[:,2]
		z_coords = self.data[:,3]

		self.Lx = self.AverageMax(x_coords, 50) - self.AverageMin(x_coords, 50)
		self.Ly = self.AverageMax(y_coords, 50) - self.AverageMin(y_coords, 50)
		self.Lz = self.AverageMax(z_coords, 50) - self.AverageMin(z_coords, 50)


	def Distance(self,dx,dy,dz):
		return (dx*dx + dy*dy + dz*dz)**0.5

	def AverageMax(self, data, N):
		data = numpy.sort(data)
		return numpy.mean(data[-N:])

	def AverageMin(self, data, N):
		data = numpy.sort(data)
		return numpy.mean(data[:N])

	def ComputeCM(self):
		""" Computes all Si-O bonds and populates a CM matrix""" 
		# adjacency matrix for Si-O bonding 
		self.CM = numpy.zeros((self.NSi,self.NO))

		# set cutoff to be 2.3 A
		cutoff = 2.3

		for i in range(self.NSi):
			for j in range(self.NSi, self.NSi+self.NO):

				dx = abs(self.data[j][1] - self.data[i][1])
				dy = abs(self.data[j][2] - self.data[i][2])
				dz = abs(self.data[j][3] - self.data[i][3])

				if self.Distance(dx,dy,dz) < cutoff:
					self.CM[i][j-self.NSi] = 1
				
				else:
					# consider PBCs
					if dx > self.Lx - cutoff:
						if self.Distance(dx + self.Lx, dy, dz) < cutoff:
							self.CM[i][j-self.NSi] = 1

						if self.Distance(dx - self.Lx, dy, dz) < cutoff:
							self.CM[i][j-self.NSi] = 1

					if dy > self.Ly - cutoff:
						if self.Distance(dx ,dy + self.Ly, dz) < cutoff:
							self.CM[i][j-self.NSi] = 1

						if self.Distance(dx, dy - self.Ly, dz) < cutoff:
							self.CM[i][j-self.NSi] = 1

					if dz > self.Lz - cutoff:
						if self.Distance(dx, dy, dz + self.Lz) < cutoff:
							self.CM[i][j-self.NSi] = 1

						if self.Distance(dx, dy, dz - self.Lz) < cutoff:
							self.CM[i][j-self.NSi] = 1

		# adjacency matrix with only bridging O 
		self.CM_T = copy.deepcopy(self.CM)


	def AMSiCIdxToDataIdx(self, idx):
		# mapping between AMSiC and data
		# order of data array is NSi, NO, NC 
		if idx < self.NSi:
			return idx
		elif idx >= self.NSi:
			# idx already contains NSi so shift by NO to get to C IDs
			return idx + self.NO

	def ComputeAMSiC(self):
		# compute adjacency matrix for Si-C and C-C bonding to determine the 
		# atoms associated with each precursor.
		# symmetric array is arranged by (Si C) x (Si C) 
		self.AMSiC = numpy.zeros((self.NSi+self.NC, self.NSi+self.NC))

		# set Si-C cutoff at 2.3 and C-C cutoff at 2.3
		cutoff_Si_C = 2.3
		cutoff_C_C = 2.3

		for i in range(self.NSi+self.NC):
			for j in range(i+1, self.NSi+self.NC):

				# get mapping between AMSiC id and data array id 
				idx_i = self.AMSiCIdxToDataIdx(i)
				idx_j = self.AMSiCIdxToDataIdx(j)

				# no Si-Si bonds so check atom type and skip if necessary
				# 1 = Si, 2 = O, 3 = C
				if (self.data[idx_i][0]==1) and (self.data[idx_j][0]==1):
					continue

				dx = abs(self.data[idx_j][1] - self.data[idx_i][1])
				dy = abs(self.data[idx_j][2] - self.data[idx_i][2])
				dz = abs(self.data[idx_j][3] - self.data[idx_i][3])

				if (self.data[idx_i][0]==3) and (self.data[idx_j][0]==3):
					cutoff = cutoff_C_C
				else:
					cutoff = cutoff_Si_C

				if self.Distance(dx,dy,dz) < cutoff:
					self.AMSiC[i][j] = 1
					self.AMSiC[j][i] = 1

				else:
					# consider PBCs
					# if dx > 0.25*self.Lx:
					if dx > self.Lx - cutoff:
						if self.Distance(dx + self.Lx, dy, dz) < cutoff:
							self.AMSiC[i][j] = 1
							self.AMSiC[j][i] = 1

						if self.Distance(dx - self.Lx, dy, dz) < cutoff:
							self.AMSiC[i][j] = 1
							self.AMSiC[j][i] = 1

					# if dy > 0.25*self.Ly:
					if dy > self.Ly - cutoff:
						if self.Distance(dx ,dy + self.Ly, dz) < cutoff:
							self.AMSiC[i][j] = 1
							self.AMSiC[j][i] = 1

						if self.Distance(dx, dy - self.Ly, dz) < cutoff:
							self.AMSiC[i][j] = 1
							self.AMSiC[j][i] = 1

					# if dz > 0.25*self.Lz:
					if dz > self.Lz - cutoff:
						if self.Distance(dx, dy, dz + self.Lz) < cutoff:
							self.AMSiC[i][j] = 1
							self.AMSiC[j][i] = 1

						if self.Distance(dx, dy, dz - self.Lz) < cutoff:
							self.AMSiC[i][j] = 1
							self.AMSiC[j][i] = 1


	def addAtoms(self, bondedTo, precursor):
		precursor_updated = copy.deepcopy(precursor)
		for atom in bondedTo:
			# python set will not repeat atom id
			# so no need to check for unique atoms
			precursor_updated.add(atom)
		return precursor_updated


	def GetPrecursors(self):
		# store list of atom IDs in each precursors
		precursors = []

		for i in range(len(self.AMSiC)):
			# find all atoms in the precurosr 
			precursor = set([])

			# seed the precursor with the list of atoms that current atom is bonded to
			bondedTo = numpy.nonzero(self.AMSiC[i])[0]
			precursor = self.addAtoms(bondedTo, precursor)

			while True:
				precursor_updated = set([])
				for atom in precursor:
					bondedTo = numpy.nonzero(self.AMSiC[atom])[0]
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

		self.precursors = unique_precursors

		# generate histogram of number of atoms per precursor
		precursor_sizes = []
		for precursor in unique_precursors:
			precursor_sizes.append(len(precursor))

		precursor_sizes = numpy.array(precursor_sizes)
		max_size = max(precursor_sizes)

		bins = numpy.linspace(0,max_size,max_size+1)
		self.precursor_histogram = numpy.histogram(precursor_sizes,bins=bins,normed=False)


	def AnalyzeConnectivity(self):
		""" Determines the connectivity of the network """ 

		# determine O bonding statistics
		for i in range(self.NO):
			NO = numpy.sum(self.CM[:,i])
			if  NO == 0:
				# free O
				self.O_stats[0] += 1
			elif NO == 1:
				# non-bridging O; set CM column to O since connectivity depends on bridigng O
				self.O_stats[1] += 1
				# the adjacnecy matrix with only bridging oxygen
				self.CM_T[:,i] = 0
			elif NO == 2:
				# bridging O
				self.O_stats[2] += 1
			elif NO == 3:
				self.O_stats[3] += 1
			elif NO == 4:
				self.O_stats[4] += 1
			elif NO == 5:
				self.O_stats[5] += 1

		# print self.O_stats

		# determine Si bonding statistics
		# T0 -> Si has all terminal groups or no Si-O bonds
		# T1 -> Si has 1 bridging O
		# T2 -> Si has 2 bridging O 
		# T3 -> Si has 3 bridging O 
		for i in range(self.NSi):
			NSi = numpy.sum(self.CM[i,:])
			
			# Si atom coordination
			if NSi == 0:
				self.Si_stats[0] += 1
			elif NSi == 1:
				self.Si_stats[1] += 1
			elif NSi == 2:
				self.Si_stats[2] += 1
			elif NSi == 3:
				self.Si_stats[3] += 1
			elif NSi == 4:
				self.Si_stats[4] += 1
			elif NSi == 5:
				self.Si_stats[5] += 1

			# Si atom configuration with Si-O-Si bridging (T groups)
			NSi_T = numpy.sum(self.CM_T[i,:])
			if NSi_T == 0:
				self.Si_statsT[0] += 1
			elif NSi_T == 1:
				self.Si_statsT[1] += 1
			elif NSi_T == 2:
				self.Si_statsT[2] += 1
			elif NSi_T == 3:
				self.Si_statsT[3] += 1
			elif NSi_T == 4:
				self.Si_statsT[4] += 1
			elif NSi_T == 5:
				self.Si_statsT[5] += 1

		# normalize the number of Si bond types 
		self.Si_stats = self.Si_stats/numpy.sum(self.Si_stats)
		self.Si_statsT = self.Si_statsT/numpy.sum(self.Si_statsT)

		# print self.Si_stats
		# print self.Si_statsT

		# compute the condensation degree
		NT0 = self.Si_statsT[0]
		NT1 = self.Si_statsT[1]
		NT2 = self.Si_statsT[2]
		NT3 = self.Si_statsT[3]
		NT4 = self.Si_statsT[4]

		self.q = CondensationDegree[self.molType](NT0, NT1, NT2, NT3, NT4)	

		# compute the relative condensation degree
		if len(self.precursors):
			self.RelativeCondensationDegree()

		# compute the Si-X-Si network connectivity
		if self.mixed:
			self.p = Connectivity[self.molType](self.q, self.Nprec1, self.Nprec2,\
												self.molType1, self.molType2)
			self.p_rel = Connectivity[self.molType](self.q_rel, self.Nprec1, self.Nprec2,\
													self.molType1, self.molType2)
		else:
			self.p = Connectivity[self.molType](self.q)
			self.p_rel = Connectivity[self.molType](self.q_rel)

		# comptue the mean Si network connectivity
		if self.mixed:
			self.mSi = MeanSiConnectivity[self.molType](self.p, self.Nprec1, self.Nprec2,\
														self.molType1, self.molType2)
			self.mSi_rel = MeanSiConnectivity[self.molType](self.p_rel, self.Nprec1, self.Nprec2,\
															self.molType1, self.molType2)
		else:
			self.mSi = MeanSiConnectivity[self.molType](self.p)
			self.mSi_rel = MeanSiConnectivity[self.molType](self.p_rel)

	def RelativeCondensationDegree(self):
		# adjust CM_T to remove bridging O 
		for i in range(self.NO):
			bondedTo = numpy.nonzero(self.CM_T[:,i])[0]
			if self.BridgingOinPrecursor(bondedTo):
				self.Nprec_w_bridging += 1
				self.CM_T[:,i] = 0

		# Si atom configurations with relevant Si-O-Si bridging (exludes Si-O-Si within the precursor)
		for i in range(self.NSi):
			NSi = numpy.sum(self.CM_T[i,:])
			if NSi == 0:
				self.Si_statsT_rel[0] += 1
			elif NSi == 1:
				self.Si_statsT_rel[1] += 1
			elif NSi == 2:
				self.Si_statsT_rel[2] += 1
			elif NSi == 3:
				self.Si_statsT_rel[3] += 1
			elif NSi == 4:
				self.Si_statsT_rel[4] += 1
			elif NSi == 5:
				self.Si_statsT_rel[5] += 1

		# normalize the number of Si bond types
		self.Si_statsT_rel = self.Si_statsT_rel/numpy.sum(self.Si_statsT_rel)

		# compute the relative condensation degree
		NT0 = self.Si_statsT_rel[0]
		NT1 = self.Si_statsT_rel[1]
		NT2 = self.Si_statsT_rel[2]
		NT3 = self.Si_statsT_rel[3]
		NT4 = self.Si_statsT_rel[4]

		self.q_rel = CondensationDegree[self.molType](NT0, NT1, NT2, NT3, NT4)	


	def BridgingOinPrecursor(self, bondedTo):
		# track the number of Si-O bonds in the precursor 
		# if the Nbonds > 1 then the precursor forms a ring 
		# via a Si-O-Si bond and that oxygen does not contribute 
		# to the network connectivity
		NSiO_bonds = 0
		for precursor in self.precursors:
			for atomID in bondedTo:
				if atomID in precursor:
					NSiO_bonds += 1

			if NSiO_bonds > 1:
				return True
			else:
				NSiO_bonds = 0
		return False


	def ComputeDensity(self):
		""" computes the density of the material """ 
	
		V = (self.Lx*1e-10)*(self.Ly*1e-10)*(self.Lz*1e-10)
		An = 6.02e23
		cc = 1e6

		NbridgingO = self.O_stats[2]
		NbridgingO_2 = self.O_stats[3] # over-coordinated O
		NterminalO = self.O_stats[1]
		NSiT2 = self.Si_stats[2]

		# print 'NbridgingO = %d' %NbridgingO
		# print 'Nnon_bridgingO = %d' %Nnon_bridgingO
		# print 'frac of T2 Si = %.4f' %NSiT2

		# compute the density depending on the precursor type
		if self.mixed:
			Nprecs = Nprecursors[self.molType](self.NSi, self.Nprec1, self.Nprec2,\
												self.molType1, self.molType2)
			Mprecs = MassPrecursor[self.molType](self.Nprec1, self.Nprec2,\
												self.molType1, self.molType2)
		else:
			Nprecs = Nprecursors[self.molType](self.NSi)
			Mprecs = MassPrecursor[self.molType]

		self.rho = (Nprecs*Mprecs + NbridgingO*Mass['O'] + NbridgingO_2*Mass['O'] +\
					NterminalO*Mass['OH'])/(V*An*cc)

		# self.rho = (Nprecursors[self.molType]*MassPrecursor[self.molType] + NbridgingO*Mass['O'] +\
		# 			+ NbridgingO_2*Mass['O'] + NterminalO*Mass['OH'])/(V*An*cc)

		# # compute the density depending on the precursor type
		# if self.molType == 'EtOCS' or 'OCSEt':
		# 	self.rho = (self.NSi*(29+14) + NbridgingO*16 + NterminalO*17 + NSiT2*self.NSi*17)/(V*An*cc)
		# elif self.molType == 'EtOCSMe':
		# 	self.rho = (self.NSi*(29+7) + NbridgingO*16 + NterminalO*17 + NSiT2*self.NSi*17)/(V*An*cc)
		# elif self.molType == '135Benz':
		# 	self.rho = (self.NSi*(29+12+13) + NbridgingO*16 + NterminalO*17 + NSiT2*self.NSi*17)/(V*An*cc)
		# elif self.molType == 'SiO2':
		# 	self.rho = (self.NSi*29 + NbridgingO*16 + NterminalO*17 + NSiT2*self.NSi*17)/(V*An*cc)


	def WriteSolution(self,inputfile):
		""" write out the network analysis """ 
		outfile = '%s_results.txt' %(inputfile[:-4])
		f = open(outfile, 'w')
		f.write('Precursor structure: %s \n\n' %self.molType)
		f.write('Condensation, q = %.4f \n' %self.q)
		f.write('Connectivity, p = %.4f \n' %self.p)
		f.write('Mean silicon network connectivity, mSi = %.4f \n' %self.mSi)
		f.write('Density, rho = %.4f \n\n' %self.rho)

		if len(self.precursors):
			f.write('Relative structure values considering Si-O-Si bridging in the precursors:\n')
			f.write('Relative condensation, q_rel = %.4f\n' %self.q_rel)
			f.write('Relative connectivity, p_rel = %.4f\n' %self.p_rel)
			f.write('Relative mean Si network connectivity, mSi_rel = %.4f\n\n' %self.mSi_rel)

			f.write('Precursor statistics:\n')
			f.write('Number of precursors = %d\n' %len(self.precursors))
			f.write('Fraction of precursors with Si-O-Si bridging = %.4f\n\n' %(float(self.Nprec_w_bridging)/len(self.precursors)))

		f.write('Si bonding statstics:\n')
		f.write('fraction N0 = %.4f \n' %self.Si_stats[0])
		f.write('fraction N1 = %.4f \n' %self.Si_stats[1])
		f.write('fraction N2 = %.4f \n' %self.Si_stats[2])
		f.write('fraction N3 = %.4f \n' %self.Si_stats[3])
		f.write('fraction N4 = %.4f \n' %self.Si_stats[4])
		f.write('fraction N5 = %.4f \n\n' %self.Si_stats[5])
		f.write('fraction T0 = %.4f \n' %self.Si_statsT[0])
		f.write('fraction T1 = %.4f \n' %self.Si_statsT[1])
		f.write('fraction T2 = %.4f \n' %self.Si_statsT[2])
		f.write('fraction T3 = %.4f \n' %self.Si_statsT[3])
		f.write('fraction T4 = %.4f \n' %self.Si_statsT[4])
		f.write('fraction T5 = %.4f \n\n' %self.Si_statsT[5])

		f.write('O bonding statistics:\n')
		f.write('free O = %d \n' %self.O_stats[0])
		f.write('non-bridging O = %d \n' %self.O_stats[1])
		f.write('bridging O = %d \n' %self.O_stats[2])
		f.write('3 coord O = %d \n' %self.O_stats[3])
		f.write('4 coord O = %d \n' %self.O_stats[4])
		f.write('5 coord O = %d \n' %self.O_stats[5])
		f.close()

		# write the precursor histogram to a separate file
		if len(self.precursors):
			histfile = '%s_histogram.txt' %inputfile[:-4]
			fhist = open(histfile, 'w')

			hist = self.precursor_histogram[0]
			# shift the bin centers to the upper bound since only integers used
			bins = self.precursor_histogram[1][1:]

			for i in range(len(hist)):
				fhist.write('%d, %d\n' %(hist[i], bins[i]))
			fhist.close()










import copy
import numpy
import os

''' finds the amount of stiffening in the matrix due to interacting terminal groups 
assumes that all bonds are weighted with equal stiffness '''


class Stiffening:

	def __init__(self,inputfiles,clustercutoff,cutoff_OC):
		self.clustercutoff = clustercutoff
		self.cutoff_OC = cutoff_OC

		self.Natoms = 0
		self.NSi = 0
		self.NC = 0
		self.NO = 0

		self.Lx = 0
		self.Ly = 0
		self.Lz = 0

		self.data = []
		self.CM = [] # Si-O bonding matrix
		self.CM_T = [] 

		self.Si_statsT = numpy.zeros(6) # only considers bridging O 

		self.p = 0 # Si-X-Si network connectivity

		self.terminalO = []
		self.all_clusters = []

		self.bins = []
		self.counts = []

		self.fname1 = 'stiffening_min_cutoff_{}.txt'.format(clustercutoff)
		self.fname2 = 'stiffening_max_cutoff_{}.txt'.format(clustercutoff)

		self.ProcessData(inputfiles)

	def ProcessData(self,inputfiles):
		f1 = open(self.fname1, 'a')
		f2 = open(self.fname2, 'a')

		for i in range(len(inputfiles)):
			# re-initialize the data attributes
			self.Natoms = 0
			self.NSi = 0
			self.NC = 0
			self.NO = 0

			self.Lx = 0
			self.Ly = 0
			self.Lz = 0

			self.data = []
			self.CM = [] # Si-O bonding matrix
			self.CM_T = [] 

			self.Si_statsT = numpy.zeros(6) # only considers bridging O 

			self.p = 0 # Si-X-Si network connectivity

			self.terminalO = []
			self.all_clusters = []

			self.bins = []
			self.counts = []


			# analyze the inputfile
			inputfile = inputfiles[i]

			self.LoadData(inputfile)
			print '\nData loaded for %s.' %inputfile

			print 'Computing Si-O bonds...'
			self.ComputeCM()

			print 'Computing the network connectivity...'
			self.ComputeConnectivity()

			# if there is no histogram files, then compute the terminal O clusters
			print 'Finding terminal O...'
			self.GetTerminalO()

			print 'Finding terminal O clusters...'
			self.ComputeClusters()

			print 'Computing the cluster stats...'
			self.ClusterStats()

			print 'Computing the total number of bonds...'
			Nbonds = self.GetNumberBonds()

			print 'Computing the number of steric bridges ...'
			Nbridges = self.GetStericBridges()

			print 'Computing the degree of stiffening from terminal O bridges...'
			degreeStiffening = self.ComputeStiffening(Nbonds, Nbridges)
			degreeStiffening_min = degreeStiffening[0]
			degreeStiffening_max = degreeStiffening[1]

			print 'Writing out results...\n'
			self.WriteResults(degreeStiffening_min, inputfile, f1)
			self.WriteResults(degreeStiffening_max, inputfile, f2)

		f1.close()
		f2.close()


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
		self.Lx = max(self.data[:,1] - min(self.data[:,1]))
		self.Ly = max(self.data[:,2] - min(self.data[:,2]))
		self.Lz = max(self.data[:,3] - min(self.data[:,3]))

		# switch numpy array back to list
		self.data = self.data.tolist()



	def Distance(self,dx,dy,dz):
		return (dx*dx + dy*dy + dz*dz)**0.5

	def Volume(self,r):
		return (4.0/3)*numpy.pi*r**3


	def ComputeCM(self):
		""" Computes all Si-O bonds and populates a CM matrix and an ID matrix""" 

		self.CM = numpy.zeros((self.NSi,self.NO))
	
		for i in range(self.NSi):
			for j in range(self.NSi,self.NSi+self.NO-1):

				dx = abs(self.data[j][1] - self.data[i][1])
				dy = abs(self.data[j][2] - self.data[i][2])
				dz = abs(self.data[j][3] - self.data[i][3])

				if self.Distance(dx,dy,dz) < 2.3:
					self.CM[i][j-self.NSi] = 1
				
				# consider PBCs
				if dx > self.Lx - 2.3:
					if self.Distance(dx + self.Lx, dy, dz) < 2.3:
						self.CM[i][j-self.NSi] = 1

					if self.Distance(dx - self.Lx, dy, dz) < 2.3:
						self.CM[i][j-self.NSi] = 1

				if dy > self.Ly - 2.3:
					if self.Distance(dx ,dy + self.Ly, dz) < 2.3:
						self.CM[i][j-self.NSi] = 1

					if self.Distance(dx, dy - self.Ly, dz) < 2.3:
						self.CM[i][j-self.NSi] = 1

				if dz > self.Lz - 2.3:
					if self.Distance(dx, dy, dz + self.Lz) < 2.3:
						self.CM[i][j-self.NSi] = 1

					if self.Distance(dx, dy, dz - self.Lz) < 2.3:
						self.CM[i][j-self.NSi] = 1


	def ComputeConnectivity(self):
		""" Determines the connectivity of the network """ 

		self.CM_T = copy.deepcopy(self.CM)

		# determine the bridging O
		for i in range(self.NO):
			NO = numpy.sum(self.CM[:,i])
			if NO == 1:
				# non-bridging O; set CM column to O since connectivity depends on bridigng O
				self.CM_T[:,i] = 0

		# determine Si bonding statistics
		# T0 -> Si has all terminal groups or no Si-O bonds
		# T1 -> Si has 1 bridging O
		# T2 -> Si has 2 bridging O 
		# T3 -> Si has 3 bridging O 
		for i in range(self.NSi):
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
		self.Si_statsT = self.Si_statsT/numpy.sum(self.Si_statsT)

		# compute the condensation degree
		NT0 = self.Si_statsT[0]
		NT1 = self.Si_statsT[1]
		NT2 = self.Si_statsT[2]
		NT3 = self.Si_statsT[3]
		NT4 = self.Si_statsT[4]
		q = (1.0/3)*(0*NT0 + 1*NT1 + 2*NT2 + 3*NT3)

		# compute the Si-X-Si network connectivity
		self.p = (1.0/4)*(1.0 + 3*q)



	def GetTerminalO(self):
		""" Compute free and terminal O and write out""" 

		for i in range(self.NO):
			NO = numpy.sum(self.CM[:,i])
			if  NO == 1: # terminal O
				# terminal O atom from data array
				idx = i + self.NSi
				self.terminalO.append(self.data[idx])



	def region(self,orig_cluster,atom_guage,cutoff):
		new_cluster = []
		for atom1 in orig_cluster:
			for atom2 in atom_guage:
				dx = abs(atom1[1]-atom2[1])
				dy = abs(atom1[2]-atom2[2])
				dz = abs(atom1[3]-atom2[3])

				if self.Distance(dx,dy,dz) <= cutoff:
					new_cluster.append(atom2)
					atom_guage.remove(atom2)

				#account for PBCs
				if dx > self.Lx-cutoff:
					if self.Distance(dx+self.Lx, dy, dz) <= cutoff:
						new_cluster.append(atom2)
						atom_guage.remove(atom2)

					if self.Distance(dx-self.Lx, dy, dz) <= cutoff:
						new_cluster.append(atom2)
						atom_guage.remove(atom2)

				if dy > self.Ly-cutoff:
					if self.Distance(dx, dy+self.Ly, dz) <= cutoff:
						new_cluster.append(atom2)
						atom_guage.remove(atom2)

					if self.Distance(dx, dy-self.Ly, dz) <= cutoff:
						new_cluster.append(atom2)
						atom_guage.remove(atom2)

				if dz > self.Lz-cutoff:
					if self.Distance(dx, dy, dz+self.Lz) <= cutoff:
						new_cluster.append(atom2)
						atom_guage.remove(atom2)

					if self.Distance(dx, dy, dz-self.Lz) <= cutoff:
						new_cluster.append(atom2)
						atom_guage.remove(atom2)

		new_cluster = new_cluster+orig_cluster
		return new_cluster,atom_guage


	def ComputeClusters(self):
		# get the atom data from the terminalO indicies
		atom_guage = copy.deepcopy(self.terminalO)

		# all_clusters = []

		if self.p < 0.8:
			cutoff = 3.1
		else:
			cutoff = self.clustercutoff

		# cutoff = self.clustercutoff

		while len(atom_guage)>0:

			cluster = [atom_guage[0]]
			atom_guage.remove(atom_guage[0])

			while True:
				size_orig_cluster = len(cluster)
				cluster, atom_guage = self.region(cluster,atom_guage,cutoff)
				size_new_cluster = len(cluster)

				if size_new_cluster==size_orig_cluster:
					self.all_clusters.append(cluster)
					break


	def ClusterStats(self):
			#get the cluster sizes
			cluster_sizes = numpy.arange(len(self.all_clusters))
			for i in range(len(self.all_clusters)):
				cluster_sizes[i] = len(self.all_clusters[i])

			# generate the bins for collecting the data frequency
			# generate the normed and counts histograms for the data set
			# note that len(hist) = len(data)-1
			self.bins = numpy.linspace(0,200,201)

			histdata = numpy.histogram(cluster_sizes,bins=self.bins,normed=False)
			self.counts = histdata[0]



	def ComputeSteric_OC(self):
		''' compute additional steric interactions of terminal O and C atoms in silane '''

		Nbridges = 0

		for cluster in self.all_clusters:
			if len(cluster) == 1:
				x1 = cluster[0][1]
				y1 = cluster[0][2]
				z1 = cluster[0][3]

				# compare each terminal O atom to the C atoms
				for i in range(self.NSi + self.NO, self.Natoms):
					x2 = self.data[i][1]
					y2 = self.data[i][2]
					z2 = self.data[i][3]

					dx = abs(x1-x2)
					dy = abs(y1-y2)
					dz = abs(z1-z2)

					# find if the C atom is within the radial cutoff
					d = self.Distance(dx,dy,dz)

					if d < self.cutoff_OC:
						Nbridges += 1

		# print 'Terminal O - C bridges = %d' %Nbridges
		return Nbridges 



	def ComputeSteric_OO(self):
		''' finds the bounds for the number of terminal O bridges;
		the bounds for each cluster are (Natoms-1) <= Nbridges <= (Natoms*(Natoms-1)/2);
		these are effectively linear ineractions and all interacting'''

		stericbonds_min = 0
		stericbonds_max = 0

		for i in range(len(self.counts)):
			# only clusters of 2 or more terminal O can form steric bridges
			Natoms_in_cluster = self.bins[i]
			frequency = self.counts[i]

			stericbonds_min += frequency*(Natoms_in_cluster-1)
			stericbonds_max += frequency*(Natoms_in_cluster*(Natoms_in_cluster-1)/2)

		return (stericbonds_min, stericbonds_max)


	def GetStericBridges(self):
		''' the connectivity model has three regimes; connectivity just above 
		the percolation threshold is a very different matrix, wherein the terminal O
		consume the majority of the free volume; at intermediate connectivity values,
		we have the porosity model, where we compute minimal and maximal terminal O interactions;
		at high connectivity values, the isolated terminal O are in a dense enviornment and can 
		interact with other atoms sterically, not just terminal O'''

		# compute the terminal O-O interactions
		bonds_min, bonds_max = self.ComputeSteric_OO()

		# compute additional interactions from terminal O and C atoms
		bonds_additional = self.ComputeSteric_OC()
		
		bonds_min += bonds_additional
		bonds_max += bonds_additional

		# pc = 0.6 #percolation threshold 

		# beta_max = 0.04706 + 3*(self.p-pc)
		# beta_min = 0.516 + 5.46*(self.p - pc) - 11.72*(self.p - pc)**2

		# return (bonds_min, bonds_max, beta_min, beta_max)

		return (bonds_min, bonds_max)



	def GetNumberBonds(self):
		''' number of Si-O bonds formed plus the number of bonds from silane precursors ''' 
		Nsilane = self.NC/2
		Nbonds_silane = 3*Nsilane

		Nbonds_SiO = numpy.sum(self.CM)

		return Nbonds_silane + Nbonds_SiO


	def ComputeStiffening(self, Nbonds, Nbridges):
		''' the degree of stiffening is related to the increase in the number of bonds;
		in our case, the terminal O are forming steric bridges and we weight a steric interaction
		and bond to have the same stiffness '''

		Nbridges_min = Nbridges[0]
		Nbridges_max = Nbridges[1]

		# degreeStiffening_min = float(Nbridges_min + Nbonds)/Nbonds
		# degreeStiffening_max = float(Nbridges_max + Nbonds)/Nbonds

		if Nbridges_max > 0.5*Nbonds:
			# assume a composite type modulus, where we weight by bond fraction
			beta_max = Nbridges_max/(Nbridges_max + Nbonds)
			beta_min = Nbridges_min/(Nbridges_min + Nbonds)

			# degreeStiffening_min = float(beta_min*(0.33*Nbridges_min) + Nbonds)/Nbonds
			# degreeStiffening_max = float(beta_max*(0.33*Nbridges_max) + Nbonds)/Nbonds

			degreeStiffening_min = float(0.5*(Nbridges_min) + Nbonds)/Nbonds
			degreeStiffening_max = float(0.5*(Nbridges_max) + Nbonds)/Nbonds

		else:
			degreeStiffening_min = float(Nbridges_min + Nbonds)/Nbonds
			degreeStiffening_max = float(Nbridges_max + Nbonds)/Nbonds

		return (degreeStiffening_min, degreeStiffening_max) 

	def WriteResults(self, degreeStiffening, inputfile, f):		
		if degreeStiffening and inputfile:
			trial = inputfile[-6:-4]
			f.write('%s,%.4f\n' %(trial,degreeStiffening))

	def WriteLine(self, string, f):
		f.write('%s,'%string)
		














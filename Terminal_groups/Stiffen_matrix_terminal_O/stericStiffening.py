import copy
import numpy
import os

''' finds the amount of stiffening in the matrix due to interacting terminal groups 
assumes that all bonds are weighted with equal stiffness '''


class Stiffening:

	def __init__(self,inputfiles,histfiles,clustercutoff,cutoff_OC):
		fname1 = 'steric_stiffening1_min_cutoff_{}.txt'.format(clustercutoff)
		fname2 = 'steric_stiffening1_max_cutoff_{}.txt'.format(clustercutoff)
		fname3 = 'steric_stiffening2_min_cutoff_{}.txt'.format(clustercutoff)
		fname4 = 'steric_stiffening2_max_cutoff_{}.txt'.format(clustercutoff)
		self.ProcessData(inputfiles,histfiles,clustercutoff,cutoff_OC,fname1,fname2,fname3,fname4)

	def ProcessData(self,inputfiles,histfiles,clustercutoff,cutoff_OC,fname1,fname2,fname3,fname4):
		f1 = open(fname1, 'a')
		f2 = open(fname2, 'a')
		f3 = open(fname3, 'a')
		f4 = open(fname4, 'a')

		f1.write('Porosity model\n')
		f2.write('Porosity model\n')
		f3.write('Connectivity model\n')
		f4.write('Connectivity model\n')

		for i in range(len(inputfiles)):
			inputfile = inputfiles[i]

			# histfile = histfiles[i]

			# bins, counts = self.LoadDataHist(histfile)
			# print '\n\nData loaded for %s.' %histfile

			data, specs = self.LoadData(inputfile)
			self.Lx = specs[4]
			self.Ly = specs[5]
			self.Lz = specs[6]
			print '\nData loaded for %s.' %inputfile

			print 'Computing Si-O bonds...'
			CM = self.ComputeCM(data,specs)

			# if there is no histogram files, then compute the terminal O clusters
			print 'Finding terminal O...'
			terminalO = self.GetTerminalO(CM,data,specs)

			print 'Finding terminal O clusters...'
			all_clusters = self.ComputeClusters(terminalO,clustercutoff)

			print 'Computing the cluster stats...'
			histdata, clusterstats = self.ClusterStats(all_clusters)
			bins = histdata[0]
			counts = histdata[1]

			print 'Computing the total number of bonds...'
			NC = specs[3]
			Nbonds = self.GetNumberBonds(NC,CM)

			print 'Computing the number of steric bridges (porosity model)...'
			Nbridges1 = self.GetStericBridgesPorosity(bins,counts)

			print 'Computing the number of steric bridges (connectivity model)...'
			NbridgesI, NbridgesII, NbridgesIII = self.GetStericBridgesConnectivity(bins,counts,data,CM,specs,all_clusters,cutoff_OC,terminalO)

			print 'Computing the degree of stiffening from terminal O bridges...'
			degreeStiffening1 = self.ComputeStiffening(Nbonds, Nbridges1, 1)
			degreeStiffening_min1 = degreeStiffening1[0]
			degreeStiffening_max1 = degreeStiffening1[1]

			# (1-Vterm/Vtot)*((Nbonds + Nsteric)/Nbonds)
			degreeStiffeningI = self.ComputeStiffening(Nbonds, NbridgesI, NbridgesI[2])
			degreeStiffening_minI = degreeStiffeningI[0]
			degreeStiffening_maxI = degreeStiffeningI[1]

			# (Nbonds + Nsteric)/Nbonds
			degreeStiffeningII = self.ComputeStiffening(Nbonds, NbridgesII,  NbridgesII[2])
			degreeStiffening_minII = degreeStiffeningII[0]
			degreeStiffening_maxII = degreeStiffeningII[1]	

			# (Nbonds + Nsteric + Nsteric_additional)/Nbonds
			degreeStiffeningIII = self.ComputeStiffening(Nbonds, NbridgesIII,  NbridgesIII[2])
			degreeStiffening_minIII = degreeStiffeningIII[0]
			degreeStiffening_maxIII = degreeStiffeningIII[1]			

			print 'Writing out results...\n'
			self.WriteResults(degreeStiffening_min1, inputfile, f1)
			self.WriteResults(degreeStiffening_max1, inputfile, f2)

			self.WriteLine('modelI',f3)
			self.WriteLine('modelI',f4)
			self.WriteResults(degreeStiffening_minI, inputfile, f3)
			self.WriteResults(degreeStiffening_maxI, inputfile, f4)

			self.WriteLine('modelII',f3)
			self.WriteLine('modelII',f4)
			self.WriteResults(degreeStiffening_minII, inputfile, f3)
			self.WriteResults(degreeStiffening_maxII, inputfile, f4)

			self.WriteLine('modelIII',f3)
			self.WriteLine('modelIII',f4)
			self.WriteResults(degreeStiffening_minIII, inputfile, f3)
			self.WriteResults(degreeStiffening_maxIII, inputfile, f4)

		f1.close()
		f2.close()
		f3.close()
		f4.close()


	def LoadDataHist(self,inputfile):
			""" reads histogram file """ 
			bins = []
			counts = []

			f = open(inputfile)

			while True:
				fields = f.readline().strip().split(',')
				if fields[0]:
					bin = float(fields[0])
					count = float(fields[1])

					bins.append(bin)
					counts.append(count)
				else:
					break
			f.close()

			return bins, counts

	def LoadData(self,inputfile):
		""" reads LAMMPS file, sorts data, and computes cell dimensions """ 
		data = []

		NSi = 0
		NC = 0
		NO = 0

		f = open(inputfile)
		Natoms = int(f.readline())

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
					NSi += 1
				elif atomtype == 2:
					NO += 1
				elif atomtype == 3:
					NC += 1
				else: 
					raise RuntimeError, "Incorrect atom type."

				# populate the data array
				data.append([atomtype,xcoord,ycoord,zcoord])

			else:
				break
		f.close()

		data = numpy.array(data)
		b = numpy.zeros((Natoms,4))

		# sort the atoms by Si, O, C
		Siidx = 0
		Oidx = NSi
		Cidx = NSi + NO
		for i in range(Natoms):
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
		Lx = max(data[:,1] - min(data[:,1]))
		Ly = max(data[:,2] - min(data[:,2]))
		Lz = max(data[:,3] - min(data[:,3]))

		# switch numpy array back to list
		data = data.tolist()

		return data, (Natoms, NSi, NO, NC, Lx, Ly, Lz)


	def Distance(self,dx,dy,dz):
		return (dx*dx + dy*dy + dz*dz)**0.5

	def Volume(self,r):
		return (4.0/3)*numpy.pi*r**3


	def ComputeCM(self,data,specs):
		""" Computes all Si-O bonds and populates a CM matrix and an ID matrix""" 
		Natoms, NSi, NO, NC, Lx, Ly, Lz = specs

		CM = numpy.zeros((NSi,NO))
	
		for i in range(NSi):
			for j in range(NSi,NSi+NO-1):

				dx = abs(data[j][1] - data[i][1])
				dy = abs(data[j][2] - data[i][2])
				dz = abs(data[j][3] - data[i][3])

				if self.Distance(dx,dy,dz) < 2.3:
					CM[i][j-NSi] = 1
				
				# consider PBCs
				if dx > Lx - 2.3:
					if self.Distance(dx + Lx, dy, dz) < 2.3:
						CM[i][j-NSi] = 1

					if self.Distance(dx - Lx, dy, dz) < 2.3:
						CM[i][j-NSi] = 1

				if dy > Ly - 2.3:
					if self.Distance(dx ,dy + Ly, dz) < 2.3:
						CM[i][j-NSi] = 1

					if self.Distance(dx, dy - Ly, dz) < 2.3:
						CM[i][j-NSi] = 1

				if dz > Lz - 2.3:
					if self.Distance(dx, dy, dz + Lz) < 2.3:
						CM[i][j-NSi] = 1

					if self.Distance(dx, dy, dz - Lz) < 2.3:
						CM[i][j-NSi] = 1
		return CM


	def ComputeConnectivity(self,CM,specs):
		""" Determines the connectivity of the network """ 
		NSi = specs[1]
		NO = specs[2]

		CM_T = copy.deepcopy(CM)

		Si_statsT = numpy.zeros(6) # T0, T1, T2, T3, T4, T5

		# determine the bridging O
		for i in range(NO):
			NO = numpy.sum(CM[:,i])
			if NO == 1:
				# non-bridging O; set CM column to O since connectivity depends on bridigng O
				CM_T[:,i] = 0

		# determine Si bonding statistics
		# T0 -> Si has all terminal groups or no Si-O bonds
		# T1 -> Si has 1 bridging O
		# T2 -> Si has 2 bridging O 
		# T3 -> Si has 3 bridging O 
		for i in range(NSi):
			NSi_T = numpy.sum(CM_T[i,:])
			if NSi_T == 0:
				Si_statsT[0] += 1
			elif NSi_T == 1:
				Si_statsT[1] += 1
			elif NSi_T == 2:
				Si_statsT[2] += 1
			elif NSi_T == 3:
				Si_statsT[3] += 1
			elif NSi_T == 4:
				Si_statsT[4] += 1
			elif NSi_T == 5:
				Si_statsT[5] += 1


		# normalize the number of Si bond types 
		Si_statsT = Si_statsT/numpy.sum(Si_statsT)

		# compute the condensation degree
		NT0 = Si_statsT[0]
		NT1 = Si_statsT[1]
		NT2 = Si_statsT[2]
		NT3 = Si_statsT[3]
		NT4 = Si_statsT[4]
		q = (1.0/3)*(0*NT0 + 1*NT1 + 2*NT2 + 3*NT3)

		# compute the Si-X-Si network connectivity
		p = (1.0/4)*(1.0 + 3*q)
		print 'Connectivity = %.4f' %p
		return p


	def GetTerminalO(self,CM,data,specs):
		""" Compute free and terminal O and write out""" 
		NSi = specs[1]
		NO = specs[2]

		terminalO = []

		for i in range(NO):
			NO = numpy.sum(CM[:,i])
			if  NO == 1: # terminal O
				# terminal O atom from data array
				idx = i + NSi
				terminalO.append(data[idx])

		return terminalO


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


	def ComputeClusters(self,atoms,cutoff):
		# get the atom data from the terminalO indicies
		atom_guage = copy.deepcopy(atoms)

		all_clusters = []

		while len(atom_guage)>0:

			cluster = [atom_guage[0]]
			atom_guage.remove(atom_guage[0])

			while True:
				size_orig_cluster = len(cluster)
				cluster, atom_guage = self.region(cluster,atom_guage,cutoff)
				size_new_cluster = len(cluster)

				if size_new_cluster==size_orig_cluster:
					all_clusters.append(cluster)
					break

		return all_clusters


	def ClusterStats(self,all_clusters):
			#get the cluster sizes
			cluster_sizes = numpy.arange(len(all_clusters))
			for i in range(len(all_clusters)):
				cluster_sizes[i] = len(all_clusters[i])

			# get average with and without isolated terminal groups 
			lg_cluster = max(cluster_sizes)
			avg_cluster = numpy.mean(cluster_sizes)
			std_cluster = numpy.std(cluster_sizes)

			steric_clusters = []
			for size in cluster_sizes:
				if size != 1:
					steric_clusters.append(size)

			steric_clusters = numpy.array(steric_clusters)

			avg_steric_cluster = numpy.mean(steric_clusters)
			std_steric_cluster = numpy.std(steric_clusters)

			# generate the bins for collecting the data frequency
			# generate the normed and counts histograms for the data set
			# note that len(hist) = len(data)-1
			bins = numpy.linspace(0,200,201)

			histdata = numpy.histogram(cluster_sizes,bins=bins,normed=False)
			hist = histdata[0]

			histdata_norm = numpy.histogram(cluster_sizes,bins=bins,normed=True)
			hist_norm = histdata_norm[0]

			# get the fraction of isolated (non-interacting) terminal O
			for i in range(len(hist_norm)):
				if bins[i] == 1:
					isolated_terminalO = hist_norm[i]

			return ((bins, hist, hist_norm), (lg_cluster, avg_cluster, std_cluster, avg_steric_cluster, std_steric_cluster, isolated_terminalO))


	def ComputeVolumeFraction(self,specs,terminalO):
		Natoms, NSi, NO, NC, Lx, Ly, Lz = specs
		NtermO = len(terminalO)

	
		gamma = 1 - NtermO/Natoms

		Vsim = Lx*Ly*Lz
		Vatoms = NO*self.Volume(1.52) + NSi*self.Volume(2.1) + NC*self.Volume(1.7)
		gamma = abs(Vsim - Vatoms)/Vsim

		return gamma


	def ComputeStericInteractions(self,data,specs,all_clusters,cutoff):
		''' compute additional steric interactions of terminal O and C atoms in silane '''
		Natoms = specs[0]
		NSi = specs[1]
		NO = specs[2]

		Nbridges = 0

		for cluster in all_clusters:
			if len(cluster) == 1:
				x1 = cluster[0][1]
				y1 = cluster[0][2]
				z1 = cluster[0][3]

				# compare each terminal O atom to the C atoms
				for i in range(NSi + NO, Natoms):
					x2 = data[i][1]
					y2 = data[i][2]
					z2 = data[i][3]

					dx = abs(x1-x2)
					dy = abs(y1-y2)
					dz = abs(z1-z2)

					# find if the C atom is within the radial cutoff
					d = self.Distance(dx,dy,dz)

					if d < cutoff:
						Nbridges += 1

		# print 'Terminal O - C bridges = %d' %Nbridges
		return Nbridges 



	def GetStericBridgesPorosity(self,bins,counts):
		''' finds the bounds for the number of terminal O bridges;
		the bounds for each cluster are (Natoms-1) <= Nbridges <= (Natoms*(Natoms-1)/2);
		these are effectively linear ineractions and all interacting'''

		stericbonds_min = 0
		stericbonds_max = 0

		for i in range(len(counts)):
			# only clusters of 2 or more terminal O can form steric bridges
			Natoms_in_cluster = bins[i]
			frequency = counts[i]

			stericbonds_min += frequency*(Natoms_in_cluster-1)
			stericbonds_max += frequency*(Natoms_in_cluster*(Natoms_in_cluster-1)/2)

		return (stericbonds_min, stericbonds_max)


	def GetStericBridgesConnectivity(self,bins,counts,data,CM,specs,all_clusters,cutoff,terminalO):
		''' the connectivity model has three regimes; connectivity just above 
		the percolation threshold is a very different matrix, wherein the terminal O
		consume the majority of the free volume; at intermediate connectivity values,
		we have the porosity model, where we compute minimal and maximal terminal O interactions;
		at high connectivity values, the isolated terminal O are in a dense enviornment and can 
		interact with other atoms sterically, not just terminal O'''

		# compute the terminal O-O interactions
		bonds_min, bonds_max = self.GetStericBridgesPorosity(bins,counts)


		# # model I: larger volume fraction of terminal O reduces strength of interaction
		# # (1-Vterm/Vtot)*((Nbonds + Nsteric)/Nbonds)
		# gamma = self.ComputeVolumeFraction(specs,terminalO)
		# NbridgesI = (bonds_min, bonds_max, 1-gamma)


		# model I: terminal O-O interactions 
		# (Nbonds + Nsteric)/Nbonds
		NbridgesI = (bonds_min, bonds_max, 1)

		# model II: terminal O-O interactions in addition to terminal O sterically interacting with C atoms
		# (Nbonds + Nsteric + Nsteric_additional)/Nbonds
		bonds_additional = self.ComputeStericInteractions(data,specs,all_clusters,cutoff)
		
		bonds_min2 = bonds_additional + bonds_min
		bonds_max2 = bonds_additional + bonds_max
		NbridgesII = (bonds_min2, bonds_max2, 1)

		# model III: model II with consideration of the percolation threshold; 
		# close to the percolation threshold, the stiffening does not have a large effect
		p = self.ComputeConnectivity(CM,specs)
		pc = 0.6 #percolation threshold 

		beta = 3*(p-pc)
		NbridgesIII = (bonds_min, bonds_max, beta)

		return NbridgesI, NbridgesII, NbridgesIII




	def GetNumberBonds(self,NC,CM):
		''' number of Si-O bonds formed plus the number of bonds from silane precursors ''' 
		Nsilane = NC/2
		Nbonds_silane = 3*Nsilane

		Nbonds_SiO = numpy.sum(CM)

		return Nbonds_silane + Nbonds_SiO


	def ComputeStiffening(self, Nbonds, Nbridges, beta):
		''' the degree of stiffening is related to the increase in the number of bonds;
		in our case, the terminal O are forming steric bridges and we weight a steric interaction
		and bond to have the same stiffness '''

		Nbridges_min = Nbridges[0]
		Nbridges_max = Nbridges[1]

		degreeStiffening_min = beta * float(Nbridges_min + Nbonds)/Nbonds
		degreeStiffening_max = beta * float(Nbridges_max + Nbonds)/Nbonds

		return (degreeStiffening_min, degreeStiffening_max) 

	def WriteResults(self, degreeStiffening, inputfile, f):		
		if degreeStiffening and inputfile:
			trial = inputfile[-6:-4]
			f.write('%s,%.4f\n' %(trial,degreeStiffening))

	def WriteLine(self, string, f):
		f.write('%s,'%string)
		














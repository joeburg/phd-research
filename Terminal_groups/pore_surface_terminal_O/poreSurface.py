import copy
import numpy
import os

''' finds the terminal O atoms situated at pore surfaces '''

class PoreSurface:

	def __init__(self,inputfiles,porefiles,clustercutoff,porecutoff):
		self.Lx = 0
		self.Ly = 0
		self.Lz = 0

		self.ProcessData(inputfiles,porefiles,clustercutoff,porecutoff)

	def ProcessData(self,inputfiles,porefiles,clustercutoff,porecutoff):
		for i in range(len(inputfiles)):
			inputfile = inputfiles[i]
			porefile = porefiles[i]

			# analysis of network 
			data, specs = self.LoadData(inputfile)
			self.Lx = specs[4]
			self.Ly = specs[5]
			self.Lz = specs[6]
			print '\nData loaded for %s.' %inputfile

			data_wo_pp, specs_pp = self.LoadData(porefile)
			print 'Data loaded for %s.' %porefile

			print 'Computing Si-O bonds...'
			CM = self.ComputeCM(data,specs)

			print 'Finding terminal O...'
			terminalO = self.GetTerminalO(CM,data,specs)

			print 'Finding porogen atoms...'
			porogenAtoms = self.GetPorogenAtoms(data, specs, data_wo_pp, specs_pp)

			print 'Finding terminal O on pore surface...'
			poresurface = self.ComputeOonSurface(terminalO, porogenAtoms, porecutoff)

			print 'Finding terminal O clusters on pore surface...'
			all_clusters = self.ComputeClusters(poresurface,clustercutoff)

			print 'Writing out results...'
			self.AnalyzeResults(poresurface,terminalO,all_clusters,porefile,clustercutoff,porecutoff)


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

		return data, (Natoms, NSi, NO, NC, Lx, Ly, Lz)


	def Distance(self,dx,dy,dz):
		return (dx*dx + dy*dy + dz*dz)**0.5


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


	def GetPorogenAtoms(self, data, specs, data_wo_pp, specs_pp):
		''' Gets porogen atoms by comparing the two files 
		(one with porogen and one without) '''

		Natoms, NSi, NO, NC, Lx, Ly, Lz = specs
		Natoms_pp, NSi_pp, NO_pp, NC_pp, Lx_pp, Ly_pp, Lz_pp = specs_pp

		Cidx = NSi + NO
		Cidx_pp = NSi_pp + NO_pp

		# the number of porogen atoms is the difference between the 
		# filled matrix and porous matrix
		porogenAtoms = numpy.zeros((NC - NC_pp,4))
		
		n = 0
		for i in range(Cidx, Cidx+NC):
			x1 = data[i][1]
			y1 = data[i][2]
			z1 = data[i][3]

			N = 0
			for j in range(Cidx_pp, Cidx_pp + NC_pp):
				x2 = data_wo_pp[j][1]
				y2 = data_wo_pp[j][2]
				z2 = data_wo_pp[j][3]

				if (x1==x2) and (y1==y2) and (z1==z2):
					break
				else:
					N += 1

			# if there are no matches, then it's a porogen atom
			if N == NC_pp:
				porogenAtoms[n] = data[i]
				n += 1

		return porogenAtoms


	def ComputeOonSurface(self, terminalO, porogenAtoms, cutoff):
		''' find terminal O that are within a cutoff distance of porogen atoms; 
		we assume this to mean they are on the pore surface '''

		poresurface = []

		for Oatom in terminalO:
			for Catom in porogenAtoms:

				dx = abs(Oatom[1] - Catom[1])
				dy = abs(Oatom[2] - Catom[2])
				dz = abs(Oatom[3] - Catom[3])

				if self.Distance(dx, dy, dz) < cutoff:
					# switch numpy array back to list for clustering algo
					poresurface.append(Oatom.tolist())
					break

		return poresurface


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


	def WriteHistData(self,bins,hist,fname):
		f = open(fname,'w')
		for i in range(len(hist)):
			if hist[i]:
				f.write('%.2f,%.2f\n'%(bins[i],hist[i]))
		f.close()

	def WriteClusterStats(self,stats,all_clusters,fname,cutoff):
		lg_cluster, avg_cluster, std_cluster,avg_steric_cluster, std_steric_cluster, isolated_terminalO = stats

		f = open(fname,'w')
		f.write('Radial cutoff = %.4f\n\n' %cutoff)
		f.write('Number of clusters = %d\n' % len(all_clusters))
		f.write('Largest cluster = %d\n' % lg_cluster)
		f.write('Fraction of isolated clusters = %.4f\n' % isolated_terminalO)
		f.write('Average cluster = %.4f\n' % avg_cluster)
		f.write('Standard deviation = %.4f\n\n' % std_cluster)
		f.write('Average steric cluster = %.4f\n' % avg_steric_cluster)
		f.write('Standard deviation steric clusters = %.4f' % std_steric_cluster)
		f.close()


	def WriteAtomsVMD(self,poresurface,fname):
		f = open(fname, 'w')
		f.write('%d\n' %len(poresurface))
		f.write('Atoms\n')
		for atom in poresurface:
			atomtype = atom[0] 
			if atomtype == 2:
				atomtype = 'O'
			elif (atomtype == 1) or (atomtype == 3):
				raise RuntimeError, "Must be a terminal O atom!"

			x = atom[1]
			y = atom[2]
			z = atom[3]

			f.write('%s  %.4f  %.4f  %.4f\n' %(atomtype, x, y, z))
		f.close()


	def WritePoreSurfaceStats(self,poresurface,terminalO,fname,cutoff):
		NatomsOnsurface = len(poresurface)
		NterminalO = len(terminalO)

		f = open(fname, 'w')
		f.write('Radial cutoff for pore surface = %.4f\n\n' %cutoff)
		f.write('Number terminal O = %d\n' %NterminalO)
		f.write('Number of terminal O on pore surface = %d\n' %NatomsOnsurface)
		f.write('Fraction terminal O on pore surface = %.6f\n' %(float(NatomsOnsurface)/NterminalO))
		f.close()


	def AnalyzeResults(self,poresurface,terminalO,all_clusters,inputfile,clustercutoff,porecutoff):
		# write out the fraction of terminal O on pore surfaces
		fname0 = 'pore_surface_cutoff_{}_{}.txt'.format(porecutoff, inputfile[:-4])

		self.WritePoreSurfaceStats(poresurface,terminalO,fname0,porecutoff)

		# get the cluster stats 
		histdata, clusterstats = self.ClusterStats(all_clusters)
		bins = histdata[0]
		hist = histdata[1]
		hist_norm = histdata[2]

		# write out the histogram results
		fname1 = 'clusters_hist_cutoff_{}_{}.csv'.format(clustercutoff, inputfile[:-4])
		fname2 = 'clusters_histnorm_cutoff_{}_{}.csv'.format(clustercutoff, inputfile[:-4])

		self.WriteHistData(bins,hist,fname1)
		self.WriteHistData(bins,hist_norm,fname2)

		# write out the cluster stats
		fname3 = 'clusters_stats_cutoff_{}_{}.txt'.format(clustercutoff, inputfile[:-4])

		self.WriteClusterStats(clusterstats, all_clusters, fname3, clustercutoff)

		# write out the atoms on the pore surfaces to VMD files in a visualization directory
		# make the top directory if it dosent exit
		vis_dir = 'visualizations/'
		if not os.path.exists(vis_dir):
			os.makedirs(vis_dir)

		fname4 = '{}{}_pore_surfaceO_VMD.xyz'.format(vis_dir,inputfile[:-4])

		self.WriteAtomsVMD(poresurface,fname4)













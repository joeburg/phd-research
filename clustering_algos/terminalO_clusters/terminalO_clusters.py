import copy
import numpy
import os


class Clusters:

	def __init__(self,inputfiles,cutoff):
		self.ProcessData(inputfiles,cutoff)

	def ProcessData(self,inputfiles,cutoff):
		for inputfile in inputfiles:
			self.Natoms = 0
			self.NSi = 0
			self.NC = 0
			self.NO = 0

			self.Lx = 0
			self.Ly = 0
			self.Lz = 0

			self.data = []
			self.CM = [] # Si-O bonding matrix

			# analysis of network 
			self.LoadData(inputfile)
			print '\nData loaded for %s.' %inputfile

			print 'Computing Si-O bonds...'
			self.ComputeCM()

			print 'Finding terminal O...'
			terminalO = self.GetTerminalO()

			print 'Finding clusters...'
			all_clusters = self.ComputeClusters(terminalO,cutoff)

			print 'Writing cluster statistics...'
			self.WriteClusterStats(all_clusters,inputfile,cutoff)

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

		self.CM_T = copy.deepcopy(self.CM)

	def GetTerminalO(self):
		""" Compute free and terminal O and write out""" 
		terminalO = []

		for i in range(self.NO):
			NO = numpy.sum(self.CM[:,i])
			if  NO == 1: # terminal O
				# index corresponding to data[] array
				idx = i + self.NSi
				terminalO.append(idx)

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

	def ComputeClusters(self,terminalO,cutoff):
		# get the atom data from the terminalO indicies
		atom_guage = []
		for i in range(len(terminalO)):
			atom_guage.append(self.data[terminalO[i]])

		all_clusters = []
		# t=0
		while len(atom_guage)>0:
			# t += 1
			# print t

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

	def WriteAtomsVMD(self,cluster,f):
		for atom in cluster:
			atomtype = atom[0] 
			if atomtype == 2:
				atomtype = 'O'
			elif (atomtype == 1) or (atomtype == 3):
				raise RuntimeError, "Must be a terminal O atom!"

			x = atom[1]
			y = atom[2]
			z = atom[3]

			f.write('%s  %.4f  %.4f  %.4f\n' %(atomtype, x, y, z))


	def WriteClusterStats(self,all_clusters,inputfile,cutoff):
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
		# generate the histogram for the data set
		# note that len(hist) = len(data)-1
		bins = numpy.linspace(0,200,201)

		histdata = numpy.histogram(cluster_sizes,bins=bins,normed=False)
		hist = histdata[0]

		# write out the histogram (only write out counts > 0)
		fname1 = 'clusters_hist_cutoff_{}_{}.csv'.format(cutoff, inputfile[:-4])
		f1 = open(fname1,'w')
		for i in range(len(hist)):
			if hist[i]:
				f1.write('%.2f,%.2f\n'%(bins[i],hist[i]))
		f1.close()

		# get normalized histogram and writeout as well
		histdata_norm = numpy.histogram(cluster_sizes,bins=bins,normed=True)
		hist_norm = histdata_norm[0]

		fname1 = 'clusters_histnorm_cutoff_{}_{}.csv'.format(cutoff, inputfile[:-4])
		f1 = open(fname1,'w')
		for i in range(len(hist_norm)):
			if hist_norm[i]:
				f1.write('%.2f,%.8f\n'%(bins[i],hist_norm[i]))
		f1.close()

		# get the fraction of isolated terminal groups and other sizes
		cluster_gte_3 = 0
		cluster_gte_4 = 0
		cluster_gte_5 = 0
		cluster_gte_6 = 0
		for i in range(len(hist_norm)):
			if bins[i] == 1:
				isolated_terminalO = hist_norm[i]

			if bins[i] >= 3:
				if bins[i] >= 4:
					if bins[i] >= 5:
						if bins[i] >= 6:
							cluster_gte_6 += hist_norm[i]
							cluster_gte_5 += hist_norm[i]
							cluster_gte_4 += hist_norm[i]
							cluster_gte_3 += hist_norm[i]
						else:	
							cluster_gte_5 += hist_norm[i]
							cluster_gte_4 += hist_norm[i]
							cluster_gte_3 += hist_norm[i]
					else:
						cluster_gte_4 += hist_norm[i]
						cluster_gte_3 += hist_norm[i]
				else:
					cluster_gte_3 += hist_norm[i]

		fname2 = 'clusters_stats_cutoff_{}_{}.txt'.format(cutoff, inputfile[:-4])
		f2 = open(fname2,'w')
		f2.write('Radial cutoff = %.4f\n\n' %cutoff)
		f2.write('Number of clusters = %d\n' % len(all_clusters))
		f2.write('Largest cluster = %d\n' % lg_cluster)
		f2.write('Fraction of isolated clusters = %.4f\n' % isolated_terminalO)
		f2.write('Fraction clusters of at least size 3 = %.6f\n' % cluster_gte_3)
		f2.write('Fraction clusters of at least size 4 = %.6f\n' % cluster_gte_4)
		f2.write('Fraction clusters of at least size 5 = %.6f\n' % cluster_gte_5)
		f2.write('Fraction clusters of at least size 6 = %.6f\n\n' % cluster_gte_6)
		f2.write('Average cluster = %.4f\n' % avg_cluster)
		f2.write('Standard deviation = %.4f\n\n' % std_cluster)
		f2.write('Average steric cluster = %.4f\n' % avg_steric_cluster)
		f2.write('Standard deviation steric clusters = %.4f' % std_steric_cluster)
		f2.close()

		# write out cluster VMD files into a visualization directory
		# make the top directory if it dosent exit
		vis_dir = 'visualizations/'
		if not os.path.exists(vis_dir):
			os.makedirs(vis_dir)

		fname3 = '{}cutoff{}_{}_clusters_VMD.xyz'.format(vis_dir,cutoff,inputfile[:-4])
		f3 = open(fname3, 'w')

		N = 0
		f3.write('%d\n' %sum(cluster_sizes))
		f3.write('Atoms\n')
		for cluster in all_clusters:
			fname4 = '{}cutoff{}_{}_clusters_VMD_{}.xyz'.format(vis_dir,cutoff,inputfile[:-4],N)

			# only write separate files for large clusters
			if len(cluster) > 10:
				N += 1
				f4 = open(fname4, 'w')
				f4.write('%d\n' % len(cluster))
				f4.write('Atoms\n')
				self.WriteAtomsVMD(cluster,f4)
				f4.close()

			# write all atoms to single file
			self.WriteAtomsVMD(cluster,f3)
		f3.close()










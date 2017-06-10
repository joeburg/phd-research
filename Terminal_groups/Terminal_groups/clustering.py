import copy
import numpy
import os

#-----------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------#
# utility functions 

def Distance(dx, dy, dz):
	return (dx*dx + dy*dy + dz*dz)**0.5


#-----------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------#

class Clusters:

	def __init__(self, moltype, data, dimensions, terminalOH, terminalGroups, inputfile, cutoff):
		# rename the file 
		inputfile = self.renameFile(inputfile)

		print '\nFinding terminal group (from precursor) clusters...'
		clusters = self.ComputeClusters(data, dimensions, terminalGroups, cutoff)
		print 'Writing cluster statistics...'
		self.histogram_groups = self.WriteClusterStats(moltype, clusters, inputfile, cutoff, 'Prec_Groups')

		print '\nFinding terminal OH clusters...'
		clusters = self.ComputeClusters(data, dimensions, terminalOH, cutoff)
		print 'Writing cluster statistics...'
		self.histogram_OH = self.WriteClusterStats(moltype, clusters, inputfile, cutoff, 'OH_Groups')
		
		print '\nFinding all clusters...'
		terminalGroups_mapped = self.MapDataStructure(terminalGroups)
		terminalOH_mapped = self.MapDataStructure(terminalOH)
		allterminalGroups = terminalGroups_mapped + terminalOH_mapped
		clusters = self.ComputeClusters(data, dimensions, allterminalGroups, cutoff)
		print 'Writing cluster statistics...'
		self.histogram_all_groups = self.WriteClusterStats(moltype, clusters, inputfile, cutoff, 'All_Groups')

	# getr methods to retrive data attributes of the class
	def getClusterHistograms(self):
		return (self.histogram_OH, self.histogram_groups,\
				self.histogram_all_groups)

	#------------------------------------------------------------------------#
	def renameFile(self, inputfile):
		# set a custom directory for the inputfile 
		base_dir = '%s/' %inputfile[:-4]
		if not os.path.exists(base_dir):
			os.makedirs(base_dir)

		# rename the inputfile to store the cluster results in a common directory
		dir_name = '%sclusters/' %base_dir
		if not os.path.exists(dir_name):
			os.makedirs(dir_name)
		inputfile = '{}{}'.format(dir_name, inputfile)
		return inputfile

	#------------------------------------------------------------------------#
	# clustering of terminal groups

	def MapDataStructure(self, terminalGroups):
		# check to see if the terminal groups contain a collection of atoms
		group = terminalGroups[0]
		data_structures = set([type(set([])), type(list()), type(tuple())])

		if type(group) in data_structures:
			# for clustering purposes, map the list of sets/lists/tuples to a 
			# tuple of the indicies
			terminalGroups_mapped = set([])
			for group in terminalGroups:
				for atom in group:
					terminalGroups_mapped.add(atom)
		else:
			terminalGroups_mapped = terminalGroups

		# convert to tuple so it can be indexed
		return tuple(terminalGroups_mapped)

	def AveragePosition(self, group, data):
		# return the average position of a group of atoms
		atom_type = None
		coords = numpy.array([0.0, 0.0, 0.0])
		for atom_idx in group:
			atom = data[atom_idx]
			if not atom_type:
				atom_type = atom[0]
			coords[0] += atom[1]
			coords[1] += atom[2]
			coords[2] += atom[3]
		coords = coords/len(group)
		return [atom_type, coords[0], coords[1], coords[2]]


	def GetAtomGuageData(self, data, terminalGroups):
		# convert the data array to a list so items can be removed 
		data = data.tolist()

		# get the atom positions based on the terminalGroups indicies
		atom_guage = []

		# if the group contains more than 1 atom it will be a data structure 
		data_structures = set([type(set([])), type(list()), type(tuple())])
		for group in terminalGroups: break # gets an element from the set/list/tuple
		if type(group) in data_structures:
			if len(group) > 1:
				# when the terminal groups has more than 1 atom, 
				# get the average postion of the terminal group atoms
				for group in terminalGroups:
					atom_guage.append(self.AveragePosition(group, data))

			else:
				# there's only 1 atom in the terminal group 
				for atom in terminalGroups:
					atom_guage.append(data[atom[0]])

			return atom_guage
			
		# populate the atom guage with the terminal groups
		for i in range(len(terminalGroups)):
			atom_guage.append(data[terminalGroups[i]])

		return atom_guage


	def region(self, dimensions, orig_cluster, atom_guage, cutoff):
		Lx, Ly, Lz = dimensions
		new_cluster = []
		for atom1 in orig_cluster:
			for atom2 in atom_guage:
				dx = abs(atom1[1]-atom2[1])
				dy = abs(atom1[2]-atom2[2])
				dz = abs(atom1[3]-atom2[3])

				if Distance(dx,dy,dz) <= cutoff:
					new_cluster.append(atom2)
					atom_guage.remove(atom2)
					continue

				#account for PBCs
				if dx > Lx-cutoff:
					if Distance(dx+Lx, dy, dz) <= cutoff:
						new_cluster.append(atom2)
						atom_guage.remove(atom2)
						continue

					if Distance(dx-Lx, dy, dz) <= cutoff:
						new_cluster.append(atom2)
						atom_guage.remove(atom2)
						continue

				if dy > Ly-cutoff:
					if Distance(dx, dy+Ly, dz) <= cutoff:
						new_cluster.append(atom2)
						atom_guage.remove(atom2)
						continue

					if Distance(dx, dy-Ly, dz) <= cutoff:
						new_cluster.append(atom2)
						atom_guage.remove(atom2)
						continue

				if dz > Lz-cutoff:
					if Distance(dx, dy, dz+Lz) <= cutoff:
						new_cluster.append(atom2)
						atom_guage.remove(atom2)
						continue

					if Distance(dx, dy, dz-Lz) <= cutoff:
						new_cluster.append(atom2)
						atom_guage.remove(atom2)
						continue

		new_cluster = new_cluster+orig_cluster
		return new_cluster,atom_guage


	def ComputeClusters(self, data, dimensions, terminalGroups, cutoff):
		# # ensure the data structure is in the right format 
		# terminalGroups = self.MapDataStructure(terminalGroups)

		# # convert the data array to a list so items can be removed 
		# data = data.tolist()

		# # get the atom data from the terminalGroups indicies
		# atom_guage = []
		# for i in range(len(terminalGroups)):
		# 	atom_guage.append(data[terminalGroups[i]])

		# get the atom data from the terminalGroups indicies
		atom_guage = self.GetAtomGuageData(data, terminalGroups)

		all_clusters = []
		while len(atom_guage)>0:

			cluster = [atom_guage[0]]
			atom_guage.remove(atom_guage[0])

			while True:
				size_orig_cluster = len(cluster)
				cluster, atom_guage = self.region(dimensions, cluster, atom_guage, cutoff)
				size_new_cluster = len(cluster)

				if size_new_cluster==size_orig_cluster:
					all_clusters.append(cluster)
					break

		return all_clusters

	#------------------------------------------------------------------------#
	# write out results 

	def WriteAtomsVMD(self,cluster,f):
		atomType = {1: 'Si', 2: 'O', 3: 'C'}
		for atom in cluster:
			atomtype, x, y, z = atom
			atomtype = atomType[atomtype]
			f.write('%s  %.4f  %.4f  %.4f\n' %(atomtype, x, y, z))

	def WriteHistogram(self, filename, histogram):
		# writes out non-zero values of a histogram
		hist = histogram[0]
		bins = histogram[1][:-1]
		f = open(filename,'w')
		for i in range(len(hist)):
			if hist[i]:
				f.write('%.2f,%.8f\n'%(bins[i],hist[i]))
		f.close()	

	def WriteClusterStats(self, moltype, all_clusters, inputfile, cutoff, clustertype):
		#get the cluster sizes
		cluster_sizes = numpy.zeros(len(all_clusters))
		for i in range(len(all_clusters)):
			cluster_sizes[i] = len(all_clusters[i])

		# get average with and without isolated terminal groups 
		lg_cluster = max(cluster_sizes)
		avg_cluster = numpy.mean(cluster_sizes)
		std_cluster = numpy.std(cluster_sizes)

		# clusters of at least size 2 experience steric effects
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
		bins = numpy.linspace(0,1000,1001)
		fname = '{}_{}_cluster_hist_cutoff_{}.txt'.format(inputfile[:-4], clustertype, cutoff)
		histogram = numpy.histogram(cluster_sizes,bins=bins,normed=False)
		self.WriteHistogram(fname, histogram)

		# get normalized histogram and writeout as well
		fname = '{}_{}_cluster_histnorm_cutoff_{}.txt'.format(inputfile[:-4], clustertype, cutoff)
		histogram_norm = numpy.histogram(cluster_sizes,bins=bins,normed=True)
		self.WriteHistogram(fname, histogram_norm)

		# get isolated groups
		isolated_terminalO = histogram_norm[0][1]

		# get the fraction of isolated terminal groups and other sizes
		fname = '{}_{}_cluster_stats_cutoff_{}.txt'.format(inputfile[:-4], clustertype, cutoff)
		f = open(fname,'w')
		f.write('Radial cutoff = %.4f\n\n' %cutoff)
		f.write('Number of clusters = %d\n' % len(all_clusters))
		f.write('Isolated clusters frac = %.4f\n' % isolated_terminalO)
		f.write('Largest cluster = %d\n' % lg_cluster)
		f.write('Average cluster, std = %.4f, %.4f\n' % (avg_cluster, std_cluster))
		f.write('Average steric cluster, std = %.4f, %.4f\n' % (avg_steric_cluster, std_steric_cluster))
		f.close()

		# write out cluster VMD files into a visualization directory
		# # make the top directory if it dosent exit
		# vis_dir = 'clusters_VMD/'
		# if not os.path.exists(vis_dir):
		# 	os.makedirs(vis_dir)

		# set cluster size cutoffs 
		cluster_cutoff = {	'EtOCSMethyl' : 6,
							'EtOCSVinyl'  : 16,
							'EtOCSPhenyl' : 100}

		# fname = '{}{}_{}_clusters_cutoff_{}_VMD.xyz'.format(vis_dir, inputfile[:-4], clustertype, cutoff)
		fname = '{}_{}_clusters_cutoff_{}_VMD.xyz'.format(inputfile[:-4], clustertype, cutoff)
		f = open(fname, 'w')
		N = 0
		f.write('%d\n' %sum(cluster_sizes))
		f.write('Atoms\n')
		for cluster in all_clusters:
			fname2 = '{}_{}.xyz'.format(fname[:-4], N)
			# only write separate files for large clusters
			if len(cluster) > cluster_cutoff[moltype]:
				N += 1
				f2 = open(fname2, 'w')
				f2.write('%d\n' % len(cluster))
				f2.write('Atoms\n')
				self.WriteAtomsVMD(cluster, f2)
				f2.close()

			# write all atoms to single file
			self.WriteAtomsVMD(cluster, f)
		f.close()

		return histogram

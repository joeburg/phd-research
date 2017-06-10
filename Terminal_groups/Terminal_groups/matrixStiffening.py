import os

#-----------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------#

class Stiffening: 
	def __init__(self, inputfile, Nbonds, Nbonds_term_min, Nbonds_term_max, histogram, beta, clustertype):

		# reanme the files to store in a common directory
		inputfile = self.renameFile(inputfile)

		print '\nComputing the number of steric interactions...'
		gamma = self.ComputeStericInteractions(histogram, Nbonds_term_min, Nbonds_term_max)

		print 'Computing the stiffening coefficients from terminal group interactions...'
		Gamma = self.ComputeStiffening(Nbonds, gamma, beta)

		print 'Writing out results...'
		self.WriteStiffeningResults(inputfile, gamma, Nbonds, Gamma, beta, clustertype)

	#------------------------------------------------------------------------#
	def renameFile(self, inputfile):
		# set a custom directory for the inputfile 
		base_dir = '%s/' %inputfile[:-4]
		if not os.path.exists(base_dir):
			os.makedirs(base_dir)

		# rename the inputfile to store the cluster results in a common directory
		dir_name = '%ssteric_interactions/' %base_dir
		if not os.path.exists(dir_name):
			os.makedirs(dir_name)
		inputfile = '{}{}'.format(dir_name, inputfile)
		return inputfile

	#------------------------------------------------------------------------#

	def ComputeStericInteractions(self, histogram, Nbonds_term_min, Nbonds_term_max):
		''' finds the bounds for the number of steric interactions;
		the bounds for each cluster are (Natoms-1) <= Nbridges <= (Natoms*(Natoms-1)/2);
		these are effectively linear ineractions and all interacting'''
		counts, bins = histogram

		# the minimum and maximum number of interactions per cluster 
		gamma_min = 0
		gamma_max = 0

		for i in range(len(counts)):
			# only clusters of 2 or more terminal O can form steric bridges
			Natoms_cluster = bins[i]
			Nclusters = counts[i]

			gamma_min += Nclusters*(Natoms_cluster-1)
			gamma_max += Nclusters*(Natoms_cluster*(Natoms_cluster-1)/2.0)

		# remove the number of interactions from bonded atoms in the terminal 
		# groups so they are not counted as steric interactions 
		gamma_min -= Nbonds_term_min
		gamma_max -= Nbonds_term_max

		return (gamma_min, gamma_max)


	def ComputeStiffening(self, Nbonds, gamma, beta):
		''' the degree of stiffening is related to the number of steric interactions;
		in our case, the terminal groups form sterically interact and we weight a steric interaction
		and bond to have the same stiffness '''

		# get the bounds on the number of steric interactions for all the clusters 
		gamma_min, gamma_max = gamma

		# compute the stiffening coefficients: Gamma_min, Gamma_max
		Gamma_min = beta * float(gamma_min + Nbonds)/Nbonds
		Gamma_max = beta * float(gamma_max + Nbonds)/Nbonds

		return (Gamma_min, Gamma_max)

	#------------------------------------------------------------------------#

	def WriteStiffeningResults(self, inputfile, gamma, Nbonds, Gamma, beta, clustertype):
		filename = '{}_stiffening_results_{}.txt'.format(inputfile[:-4], clustertype)
		f = open(filename, 'w')
		f.write('Cluster Type: %s\n' %clustertype)
		f.write('gamma min = %.4f\n' %gamma[0])
		f.write('gamma max = %.4f\n' %gamma[1])
		f.write('Nbonds = %d\n' %Nbonds)
		f.write('\nFor Beta = %.4f:\n' %beta)
		f.write('Gamma min = %.4f\n' %Gamma[0])
		f.write('Gamma max = %.4f\n' %Gamma[1])
		f.close()
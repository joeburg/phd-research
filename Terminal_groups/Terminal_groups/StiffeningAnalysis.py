''' The previous code does not correcly compute the stiffening coefficients 
This program uses the clustering data to re-compute the stiffening coefficients '''

import glob
import sys
import time

#-----------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------#

def LoadClusterHistogram(inputfile):
	f = open(inputfile)

	data = []
	while True:
		fields = f.readline().strip().split(',')
		if len(fields)>1: 
			nAtomsInCluster = float(fields[0])
			nClusters = float(fields[1])
			data.append((nAtomsInCluster, nClusters))
		else:
			break
	return data


def NIntsBetweenTerminalGroupsMax(nGroups):
	return nGroups*(nGroups-1)*0.5

def NIntsBetweenTerminalGroupsMin(nGroups):
	return nGroups - 1

def NTerminalGroupsInCluster(nAtomsInCluster, moltype):
	nAtomsPerGroup = {'EtOCSMethyl': 1.0, 'EtOCSVinyl': 2.0, 'EtOCSPhenyl': 6.0}
	return int(nAtomsInCluster/nAtomsPerGroup[moltype])

def ComputeStiffening(data, moltype):
	# the min and max number of interactions between pairs of terminal groups 
	nAtomIntsPerPairOfGroupsMin = {'EtOCSMethyl': 1, 'EtOCSVinyl': 1, 'EtOCSPhenyl': 4}
	nAtomIntsPerPairOfGroupsMax = {'EtOCSMethyl': 1, 'EtOCSVinyl': 4, 'EtOCSPhenyl': 36} 

	nStericInteractionsMin = 0 # gamma_min
	nStericInteractionsMax = 0 # gamma_max
	for cluster in data:
		nAtomsInCluster, nClusters = cluster

		nTerminalGroups = NTerminalGroupsInCluster(nAtomsInCluster, moltype)
		
		nGroupIntsMin = NIntsBetweenTerminalGroupsMin(nTerminalGroups)
		nGroupIntsMax = NIntsBetweenTerminalGroupsMax(nTerminalGroups)

		nStericInteractionsMin += nGroupIntsMin * nAtomIntsPerPairOfGroupsMin[moltype] * nClusters
		nStericInteractionsMax += nGroupIntsMax * nAtomIntsPerPairOfGroupsMax[moltype] * nClusters

	return (nStericInteractionsMin, nStericInteractionsMax)

def ComputeStiffeningOH(data):
	nStericInteractionsMin = 0 # gamma_min
	nStericInteractionsMax = 0 # gamma_max
	for cluster in data:
		nAtomsInCluster, nClusters = cluster

		nStericInteractionsMin += (nAtomsInCluster-1)*nClusters
		nStericInteractionsMax += (nAtomsInCluster*(nAtomsInCluster-1)*0.5)*nClusters

	return (nStericInteractionsMin, nStericInteractionsMax)


def ComputeStiffeningCoeffs(data):
	nStericInteractionsMin = 0 # gamma_min
	nStericInteractionsMax = 0 # gamma_max
	for cluster in data:
		nAtomsInCluster, nClusters = cluster

		nStericInteractionsMin += (nAtomsInCluster-1)*nClusters
		nStericInteractionsMax += (nAtomsInCluster*(nAtomsInCluster-1)*0.5)*nClusters

	return (nStericInteractionsMin, nStericInteractionsMax)

#-----------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------#

if len(sys.argv) < 2:
	print 'Usage:'
	print '  python %s <precursor type> [OH - False]' %sys.argv[0]
	exit()

moltype = sys.argv[1]

if len(sys.argv) > 2: 
	OHGroups = True
else:
	OHGroups = False

t0 = time.time()
		
# get all the relevant files and process each network
inputfiles = glob.glob('{}_*.txt'.format(moltype))

# write all the results to the same file
f = open('steric_interactions.txt', 'w')
f.write('Filename : gamma_min, gamma_max\n')

for inputfile in inputfiles:

	print 'Working with %s...' %inputfile
	
	data = LoadClusterHistogram(inputfile)
	gamma_min, gamma_max = ComputeStiffeningCoeffs(data)

	# if OHGroups:
	# 	gamma_min, gamma_max = ComputeStiffeningOH(data)
	# else:
	# 	gamma_min, gamma_max = ComputeStiffening(data, moltype)

	f.write('%s : %.4f, %.4f\n' %(inputfile, gamma_min, gamma_max))

f.close()

print 'Analyzed network in %.4f seconds.' %(time.time()-t0) 
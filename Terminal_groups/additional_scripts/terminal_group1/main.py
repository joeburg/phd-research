import glob
import sys
import time

from clustering import Clusters
from matrixStiffening import Stiffening
from poreSurface import PoreSurface
from terminalGroups import TerminalGroups

#-----------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------#



if len(sys.argv) < 3:
	print 'Usage:'
	print '  python %s <precursor type> <cluster cutoff> [pore surface = False]' %sys.argv[0]
	exit()
elif len(sys.argv) == 4:
	pore_surface = True
else:
	pore_surface = False

moltype = sys.argv[1]
cutoff = float(sys.argv[2])


# get all the relevant files and process each network
inputfiles = glob.glob('{}_*.xyz'.format(moltype))

for inputfile in inputfiles:
	t0 = time.time()
	print '\nWorking with %s...' %inputfile

	# compute the terminal groups 
	terminalgroups = TerminalGroups(inputfile, moltype)
	terminalOH, terminalGroups = terminalgroups.getTerminalGroups()
	dimensions = terminalgroups.getCellDimensions()
	atomData = terminalgroups.getAtomData()
	Nbonds, Nbonds_term_min, Nbonds_term_max = terminalgroups.getNbonds()

	# compute the terminal groups on the pore surface 
	if pore_surface:
		PoreSurface(inputfile, atomData, terminalOH, terminalGroups, dimensions)

	else:
		# compute the clustering of the terminal groups 
		clusters = Clusters(moltype, atomData, dimensions, terminalOH,\
							terminalGroups, inputfile, cutoff)
		histogram_OH, histogram_groups,\
					histogram_all_groups = clusters.getClusterHistograms()

		# compute the stiffening coefficients
		beta = 1  
		Stiffening(inputfile, Nbonds, 0, 0, histogram_OH, beta, 'OH_Groups')
		Stiffening(inputfile, Nbonds, Nbonds_term_min, Nbonds_term_max, histogram_groups, beta, 'Prec_Groups')
		Stiffening(inputfile, Nbonds, Nbonds_term_min, Nbonds_term_max, histogram_all_groups, beta, 'All_Groups')


	print '\nAnalyzed network in %.4f seconds.\n' %(time.time()-t0)
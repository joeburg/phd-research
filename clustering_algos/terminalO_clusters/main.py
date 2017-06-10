import terminalO_clusters
import glob
import sys
import time

##---------------------------------------------------------------------------##

if len(sys.argv) < 2:
	print 'Usage:'
	print '  python %s <radial cutoff>' %sys.argv[0]
	exit()

try:
	t0 = time.time()
	filenames = glob.glob('OCSEtpore_*.xyz')
	cutoff = float(sys.argv[1])

	terminalO_clusters.Clusters(filenames,cutoff)

	# computeDistTerminalO.TerminalODist(inputfile)
	print 'Found terminal O clusters in %.4f seconds.' %(time.time()-t0)
except Exception, e:
	print 'ERROR: %s' % e 
	exit()
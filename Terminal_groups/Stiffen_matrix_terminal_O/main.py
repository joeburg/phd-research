import stericStiffening
import glob
import sys
import time

##---------------------------------------------------------------------------##

if len(sys.argv) < 3:
	print 'Usage:'
	print '  python %s <radial cluster cutoff> <radial C-O cutoff>' %sys.argv[0]
	exit()

try:
	t0 = time.time()

	clustercutoff = float(sys.argv[1])
	cutoff_CO = float(sys.argv[2])
	
	inputfiles = glob.glob('OCSEt_140000_*.xyz')
	histfiles = glob.glob('clusters_hist_cutoff_{}_OCSEt_140000_*.csv'.format(sys.argv[1]))


	stericStiffening.Stiffening(inputfiles,histfiles,clustercutoff,cutoff_CO)

	print 'Found the degree of stiffening due to terminal O in %.4f seconds.' %(time.time()-t0)
except Exception, e:
	print 'ERROR: %s' % e 
	exit()
import poreSurface
import glob
import sys
import time

##---------------------------------------------------------------------------##

if len(sys.argv) < 3:
	print 'Usage:'
	print '  python %s <radial cluster cutoff> <radial pore cutoff>' %sys.argv[0]
	exit()

try:
	t0 = time.time()
	
	porefilenames = glob.glob('OCSEtpore_*.xyz')
	filenames = glob.glob('OCSEt_*.xyz')
	clustercutoff = float(sys.argv[1])
	porecutoff = float(sys.argv[2])

	poreSurface.PoreSurface(filenames,porefilenames,clustercutoff,porecutoff)

	print 'Found terminal O on pore surfaces in %.4f seconds.' %(time.time()-t0)
except Exception, e:
	print 'ERROR: %s' % e 
	exit()
import glob
import network
import sys
import time

if len(sys.argv) < 2:
	print 'Usage:'
	print '  python %s <precursor type (EtOCS, EtOCSMethyl, 135Benz, etc.) [relative q = False]>' %sys.argv[0]
	exit()

molType = sys.argv[1]

if len(sys.argv) == 3:
	q_relative = True
else:
	q_relative = False

inputfiles = glob.glob('{}_*.xyz'.format(molType))

for inputfile in inputfiles:
	try:
		t0 = time.time()
		print '\nWorking with %s...' %inputfile
		network.Network(inputfile, molType, q_relative)
		print 'Analyzed network in %.4f seconds.' %(time.time()-t0)
	except Exception, e:
		print 'ERROR: %s' % e 
		exit()
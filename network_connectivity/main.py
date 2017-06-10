import glob
import network
import sys
import time

from utils import getPrecursorParams
#-----------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------#


if len(sys.argv) < 2:
	print 'Usage:'
	print '  python %s <precursor type> [relative q = False] [Nprecursors1] [Nprecursors2]' %sys.argv[0]
	exit()

molType = sys.argv[1]

# check if the relative condensation degree will be computed
if len(sys.argv) == 3:
	q_relative = sys.argv[2]
	if q_relative in ['True', 'true', 'TRUE']:
		q_relative = True
	else:
		q_relative = False
else:
	q_relative = False

# analyze the molType
if len(sys.argv) > 3:
	PrecursorParams = getPrecursorParams(molType, int(sys.argv[3]), int(sys.argv[4]))
else:
	PrecursorParams = getPrecursorParams(molType)

# get all the relevant files and process each network
inputfiles = glob.glob('{}_*.xyz'.format(molType))
for inputfile in inputfiles:
	# t0 = time.time()
	# print '\nWorking with %s...' %inputfile
	# network.Network(inputfile, PrecursorParams, q_relative)
	# print 'Analyzed network in %.4f seconds.' %(time.time()-t0)

	try:
		t0 = time.time()
		print '\nWorking with %s...' %inputfile
		network.Network(inputfile, PrecursorParams, q_relative)
		print 'Analyzed network in %.4f seconds.' %(time.time()-t0)
	except Exception, e:
		print 'ERROR: %s' % e 
		exit()
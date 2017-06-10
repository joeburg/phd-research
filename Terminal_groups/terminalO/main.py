import terminalO
import sys
import time

if len(sys.argv) < 2:
	print 'Usage:'
	print '  python %s <inputfile>' %sys.argv[0]
	exit()

inputfile = sys.argv[1]

try:
	t0 = time.time()
	terminalO.TerminalO(inputfile)
	print 'Found terminal O in %.4f seconds.' %(time.time()-t0)
except Exception, e:
	print 'ERROR: %s' % e 
	exit()
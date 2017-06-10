import AtomDistance
import glob
import matplotlib.pyplot as plt
import pylab
import sys
import time


if len(sys.argv) < 1:
	print 'Usage:'
	print '  python %s' %sys.argv[0]
	exit()


try:
	t0 = time.time()

	filenames = glob.glob('OCSEtpore_*.xyz')
	colors = ['g','b','r','#228b22','k','c','m','y','#6a5acd']

	AtomDistance.AtomDistance(filenames,colors)	
	
	print 'Found terminal O in %.4f seconds.' %(time.time()-t0)
except Exception, e:
	print 'ERROR: %s' % e 
	exit()
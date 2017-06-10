# remove bins with not counts from the histogram data

if len(sys.argv) < 2:
	print 'Usage:'
	print '  python %s <radial cutoff>' %sys.argv[0]
	exit()

try:
	t0 = time.time()
	filenames = glob.glob('OCSEtpore_*.xyz')

	
	

	# computeDistTerminalO.TerminalODist(inputfile)
	print 'Removed bins with 0 counts in %.4f seconds.' %(time.time()-t0)
except Exception, e:
	print 'ERROR: %s' % e 
	exit()
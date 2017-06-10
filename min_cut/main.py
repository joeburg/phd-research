import min_cut
import sys
import time

if len(sys.argv) < 3:
	print 'Usage:'
	print '  python %s <inputfile> <min-cut height (in Angstroms)>' %sys.argv[0]
	# print '  python %s <inputfile> <precursor type (EtOCS, EtOCSMe, 135Benz, SiO2)>' %sys.argv[0]
	exit()

inputfile = sys.argv[1]
delta = float(sys.argv[2])


t0 = time.time()
min_cut.MinCut(inputfile,delta)
print 'Analyzed network in %.4f seconds.' %(time.time()-t0)

# try:
# 	t0 = time.time()
# 	min_cut.MinCut(inputfile,molType)
# 	print 'Analyzed network in %.4f seconds.' %(time.time()-t0)
# except Exception, e:
# 	print 'ERROR: %s' % e 
# 	exit()
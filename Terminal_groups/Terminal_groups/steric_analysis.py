import glob
import sys
import time

#-----------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------#



#-----------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------#

if len(sys.argv) < 2:
	print 'Usage:'
	print '  python %s <precursor type> [OH - False]' %sys.argv[0]
	exit()

moltype = sys.argv[1]

if len(sys.argv) > 2: 
	OHGroups = True
else:
	OHGroups = False

t0 = time.time()
		
# get all the relevant files and process each network
inputfiles = glob.glob('{}_*.txt'.format(moltype))

# write all the results to the same file
f = open('steric_interactions.txt', 'w')
f.write('Filename : gamma_min, gamma_max\n')

for inputfile in inputfiles:

	print 'Working with %s...' %inputfile
	
	data = LoadClusterHistogram(inputfile)

	if OHGroups:
		gamma_min, gamma_max = ComputeStiffeningOH(data)
	else:
		gamma_min, gamma_max = ComputeStiffening(data, moltype)

	f.write('%s : %.4f, %.4f\n' %(inputfile, gamma_min, gamma_max))

f.close()

print 'Analyzed network in %.4f seconds.' %(time.time()-t0) 
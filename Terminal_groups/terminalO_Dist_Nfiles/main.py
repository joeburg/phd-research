import computeDistTerminalO
import glob
import matplotlib.pyplot as plt
import pylab
import sys
import time

# if len(sys.argv) < 2:
# 	print 'Usage:'
# 	print '  python %s <inputfile>' %sys.argv[0]
# 	exit()
# 
# inputfile = sys.argv[1]


if len(sys.argv) < 1:
	print 'Usage:'
	print '  python %s' %sys.argv[0]
	exit()


try:
	t0 = time.time()

	filenames = glob.glob('OCSEtpore_*.xyz')
	colors = ['g','b','r','#228b22','k','c','m','y','#6a5acd']

	# set font type and axes thickness
	arialfont = {'fontname':'Arial'}
	plt.rcParams['axes.linewidth'] = 2
	plt.rcParams.update({'fontname':'Arial','font.size': 16})

	plt.figure(1)
	for i in range(len(filenames)):
		file = filenames[i]
		color = colors[i]
		computeDistTerminalO.TerminalODist(file,color)

	plt.xlabel('Terminal O-O Distance, $\delta$ ($\AA$)',size=20,**arialfont)
	plt.ylabel('Frequency',size=20,**arialfont)
	plt.xlim((0,120))
	plt.axes().set_aspect(1./plt.axes().get_data_ratio())
	plt.savefig('terminalODist.pdf',bbox_inches='tight')
	plt.close(1)		

	# computeDistTerminalO.TerminalODist(inputfile)
	print 'Found terminal O in %.4f seconds.' %(time.time()-t0)
except Exception, e:
	print 'ERROR: %s' % e 
	exit()
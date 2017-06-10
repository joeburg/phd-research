import glob
import sys
import time

from analyzeGeneralStructure import analyzeStructure
from generateDataFiles import createDataFile
from generateInputFiles import createInputFiles
from utils import getSimulationData

#-------------------------------------------------------------------------------------------------#
if len(sys.argv) < 2:
	print 'Usage:'
	print '  python %s <Usage Info: True/False>' %sys.argv[0]
	exit()

t0 = time.time()

usageInfo = sys.argv[1]

if usageInfo == 'True':
	print '''\nUse the 'simulation_starter.yml' file in 'inputfiles/' to 
	setup the simulation.

	IMPORTANT: Place all files in the 'inputfiles/' directory.
	IMPORTANT: Do not change the name of the .yml file! 
	'''
	exit()

# get the simulation data 
SimulationData = getSimulationData()

# analyze the precursor structure 
Structure, Potentials = analyzeStructure(SimulationData)

# create the data file 
createDataFile(SimulationData, Structure, Potentials)

# create the input file 
createInputFiles(SimulationData, Structure, Potentials)
	
print '\nCreated simulation files in %.4f seconds.\n' %(time.time()-t0)
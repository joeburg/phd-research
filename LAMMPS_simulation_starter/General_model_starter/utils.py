import glob
import os
import yaml

#-------------------------------------------------------------------------------------------------#

def getSimulationData():
	input_dir = 'inputfiles/'
	simulationFile = '%ssimulation_starter.yml' %input_dir
	with open(simulationFile) as f:
		SimulationData = yaml.load(f)

	# adjust the path names with the input directory
	SimulationData['atomfile'] = '%s%s' %(input_dir, SimulationData['atomfile'])
	SimulationData['atomTypesFile'] = '%s%s' %(input_dir, SimulationData['atomTypesFile'])
	SimulationData['bondTypesFile'] = '%s%s' %(input_dir, SimulationData['bondTypesFile'])
	SimulationData['angleTypesFile'] = '%s%s' %(input_dir, SimulationData['angleTypesFile'])
	SimulationData['dihedralTypesFile'] = '%s%s' %(input_dir, SimulationData['dihedralTypesFile'])
	SimulationData['pairPotentialFile'] = '%s%s' %(input_dir, SimulationData['pairPotentialFile'])

	# create the output directory 
	sim_dir = SimulationData['sim_dir']
	SimulationData['sim_dir'] = MakeDirectory(sim_dir)

	return SimulationData

# def GetFileTypes(filename):
# 	atomTypesFile = glob.glob('inputfiles/*_atomTypes.txt')[0]
# 	# atomTypesFile = '%s_atomTypes.txt' %filename[:-4]
# 	bondTypesFile = '%s_bondTypes.txt' %filename[:-4]
# 	angleTypesFile = '%s_angleTypes.txt' %filename[:-4]
# 	dihedralTypesFile = '%s_dihedralTypes.txt' %filename[:-4]

# 	return (atomTypesFile, bondTypesFile,\
# 			angleTypesFile, dihedralTypesFile)

# def MakeResultsDirectory():
# 	# make nice path slug if spaces are given
# 	dir_name = 'structure_data/'
# 	# make the results directory if it dosent exit
# 	if not os.path.exists(dir_name):
# 		os.makedirs(dir_name)
# 	return dir_name

def MakeDirectory(sim_name):
	# make nice path slug if spaces are given
	slug_list = sim_name.strip().split()
	slug = ''
	for item in slug_list:
		slug += '{}_'.format(item)
	# replace the last _ with /
	slug = slug[:-1]+'/'

	# make the top directory if it dosent exit
	if not os.path.exists(slug):
		os.makedirs(slug)

	return slug
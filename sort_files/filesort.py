import glob
import os

#------------------------------------------------------------------------------#
def getSimulationMetrics(logfiles):
	conditions = []
	trials = []
	for logfile in logfiles:
		# format of a logfile name is log_1.2.lammps
		# get the condition index (i.e. '1' in the example name)
		condition = int(logfile[logfile.find('_')+1 : logfile.find('.')])
		conditions.append(condition)

		# get the trial number of the condition (i.e. '2' in the example)
		trial_per_condition = int(logfile[logfile.find('.')+1 : -7])
		trials.append(trial_per_condition)

	Nconditions = max(conditions)
	Ntrials = max(trials)
	return (Nconditions, Ntrials)


def sortRelevantFiles(fileList, trial, file_dir):
	# add directory to trial directory
	if not os.path.exists(file_dir):
		os.makedirs(file_dir)

	relevant_files = []
	for dfile in fileList:
		if trial in dfile:
			relevant_files.append(dfile)

	for rfile in relevant_files:
		if os.path.exists(rfile):
			os.rename(rfile, file_dir+rfile)


def sortSimulationFiles(simulation_files):
	sim_dir = 'simulation_files/'

	for simfile in simulation_files:
		if not os.path.exists(sim_dir):
			os.makedirs(sim_dir)

		# rename the file
		if os.path.exists(simfile):
			os.rename(simfile, '{}{}'.format(sim_dir, simfile))


def createModulsDir(path):
	modulus_dir = path + 'modulus/'
	if not os.path.exists(modulus_dir):
		os.makedirs(modulus_dir)

def createConnectivityDir(path):
	network_conn_dir = path + 'network_connectivity/'
	if not os.path.exists(network_conn_dir):
		os.makedirs(network_conn_dir)


#------------------------------------------------------------------------------#

moltype = raw_input('Give the precursor name: ')

# get lists of the data files
datafiles = glob.glob('data.{}_*_*'.format(moltype))
inputfile = 'in.{}'.format(moltype)
hostfile = 'hostfile'
logfile = 'log.lammps'
executable = 'lmp_rhd'
potential = 'OCS.sw'
sim_desc = 'simulation_description.txt'
simulation_files = [inputfile, hostfile, logfile, executable, potential, sim_desc] + datafiles

atomPostionfiles = [f for f in glob.glob('{}_*_*.xyz'.format(moltype)) if 'matrix' not in f]
logfiles = glob.glob('log_*.lammps')
rdffiles = glob.glob('{}_rdf_*.txt'.format(moltype))
matrixfiles = glob.glob('{}_matrix_*_*.xyz'.format(moltype))
porogenfiles = glob.glob('porogen_*_*.xyz')

# initialize the trial directory
trial_dir = 'trial_{}/'

# sort simulation files into the respective directory
sortSimulationFiles(simulation_files)

# sort each trial into the respective directory
Nconditions, Ntrials = getSimulationMetrics(logfiles)

for i in range(Nconditions):
	top_trial_dir = trial_dir.format(i+1)
	if not os.path.exists(top_trial_dir):
		os.makedirs(top_trial_dir)

	for j in range(Ntrials):
		trial = '{}.{}'.format(i+1,j+1)
		sub_trial_dir = trial_dir.format(trial)
		path = '{}{}'.format(top_trial_dir, sub_trial_dir)

		# sort atom data files 
		sortRelevantFiles(atomPostionfiles, trial, path+'atom_data/')

		# sort log files
		sortRelevantFiles(logfiles, trial, path)

		# sort rdf files
		sortRelevantFiles(rdffiles, trial, path+'rdf_data/')

		# sort porogen files
		sortRelevantFiles(matrixfiles+porogenfiles, trial, path+'porogen_data/')

		# create connectivity and modulus dir 
		createConnectivityDir(path)
		createModulsDir(path)



print 'Files have been sorted!'


























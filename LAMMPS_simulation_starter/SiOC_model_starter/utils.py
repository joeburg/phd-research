import os
import shutil
import yaml
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#

def CopyFileToNewFolder(filename, new_dir):
	src_dir= os.curdir
	dst_dir= os.path.join(src_dir, new_dir)
	src_file = os.path.join(src_dir, filename)
	shutil.copy(src_file, dst_dir)
        

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


def GetMolType(molID):
	saved_precursors = {1: 'EtOCS', 
						2: 'EtOCSMethyl', 
						3: 'EtOCSVinyl',
						4: 'EtOCSPhenyl', 
						5: '135Benzene', 
						6: '13Benzene',
						7: '14Benzene', 
						8: 'MSSQ', 
						9: 'SiO2',
						10: 'TVSR1',
						11: 'TVSR2',
						12: 'TVSR3',
						13: 'TVSR4',
						14: 'TVSR6'}
	return saved_precursors[molID]


def GetSimulationParams(sim_type):
	Nconditions = 1
	Ntrials_per_cond = 1
	O_Si_ratio = None 
	porosity = None
	O_Si_ratio_start = None
	O_Si_ratio_end = None
	porosity_start = None
	porosity_end = None 

	# general data file takes no params (sim_type == 1)

	# single simulation without porosity 
	if sim_type == 2:
		O_Si_ratio = input('Give the O:Si ratio: ')
		porosity = 0

	# single simulation with porosity 
	elif sim_type == 3:
		O_Si_ratio = input('Give the O:Si ratio: ')
		porosity = input('Give the volume percent porosity [e.g. 30]: ')

	# network connectivity scan
	elif sim_type == 4:
		Ntrials_per_cond = input('How many trials per condition? ')
		O_Si_ratio_start = input('Give the starting O:Si ratio (Note: ascending range): ')
		O_Si_ratio_end = input('Give the ending O:Si ratio: ')
		Nconditions = input('How many O:Si ratio steps? ')
		porosity = input("Give the volume percent porosity [e.g. 30]: ")

		if O_Si_ratio_start > O_Si_ratio_end:
			raise RuntimeError, "The O:Si ratio range must be ascending!"

	# porosity scan 
	elif sim_type == 5:
		Ntrials_per_cond = input('How many trials per condition? ')
		O_Si_ratio = input('Give the O:Si ratio: ')
		porosity_start = input('Give the starting porosity (Note: ascending range) [e.g. 0]: ')
		porosity_end = input('Give the ending porosity [e.g. 50]: ')
		Nconditions = input('How many porosity steps? ')

		if porosity_start > porosity_end:
			raise RuntimeError, "The porosity range must be ascending!"

	# one connectivity, multiple trials
	elif sim_type == 6: 
		O_Si_ratio = input('Give the O:Si ratio: ')
		Ntrials_per_cond = input('How many trials? ')
		porosity = input("Give the volume percent porosity [e.g. 30]: ")

	return {'Nconditions' : Nconditions,
			'Ntrials_per_cond' : Ntrials_per_cond, 
			'O_Si_ratio' : O_Si_ratio,
			'porosity' : porosity,
			'O_Si_ratio_start' : O_Si_ratio_start,
			'O_Si_ratio_end' : O_Si_ratio_end,
			'porosity_start' : porosity_start,
			'porosity_end' : porosity_end}


def GetSimulationDescriptions(sim_type):
	sim_desc = {1: 'generate general data file',
				2: 'single simulation (no porosity)',
				3: 'single simulation (with porosity)',
				4: 'network connectivity scan',
				5: 'porosity scan',
				6: 'one connectivity, multiple trials'}
	return sim_desc[sim_type]


def ChooseDynamics():
	return input('The dynamics of your precursor are not saved!\n'+\
					'Choose dynamics based on a similar precursor structure\n'+\
					'(Note: this will get you started with an inputfile but'+\
					'you will need to test the dyanmics and add to the database - \n'+\
					'sim_dynamics.yml):\n\n'+\
					'\t1) EtOCS, MSSQ, SiO2, TVSR2, TVSR4\n'+\
					'\t2) EtOCSMethyl, EtOCSVinyl, TVSR6\n'+\
					'\t3) EtOCSPhenyl\n'+\
					'\t4) 135Benzene, 13Benzene, 14Benzene\n\n'+\
					'\tchoice: ')


def LoadSimulationDynamics(sim_type, moltype1, moltype2=False):
	filename = 'database/sim_dynamics.yml'
	with open(filename) as f:
		SimDynamics = yaml.load(f)

	dynamics = {'EtOCS' : SimDynamics['sim_dynam_1'],
				'MSSQ' : SimDynamics['sim_dynam_1'],
				'SiO2' : SimDynamics['sim_dynam_1'],
				'TVSR2' : SimDynamics['sim_dynam_1'],
				'TVSR4' : SimDynamics['sim_dynam_1'],
				'EtOCSMethyl' : SimDynamics['sim_dynam_2'],
				'EtOCSVinyl' : SimDynamics['sim_dynam_2'],
				'TVSR6' : SimDynamics['sim_dynam_2'],
				'EtOCSPhenyl_no_porosity' : SimDynamics['sim_dynam_3'],
				'EtOCSPhenyl_porosity' : SimDynamics['sim_dynam_5'],
				'135Benzene' : SimDynamics['sim_dynam_4'],
				'13Benzene' : SimDynamics['sim_dynam_4'],
				'14Benzene' : SimDynamics['sim_dynam_4'],
				1 : SimDynamics['sim_dynam_1'],
				2 : SimDynamics['sim_dynam_2'],
				3 : SimDynamics['sim_dynam_3'],
				4 : SimDynamics['sim_dynam_4'],
				5 : SimDynamics['sim_dynam_5'],}

	# if dual precursors, then prompt the user to choose the anneal schedule
	if moltype2:
		sim_dynam = ChooseDynamics()
		return dynamics[sim_dynam]

	# the EtOCSPhenyl simulations may need different anneal dynamics
	# either sim_dynam_3 (no porosity) or sim_dynam_5 (with porosity)
	if moltype1 == 'EtOCSPhenyl':
		print moltype1
		# sim type 5 is with porosity 
		if sim_type == 5:
			moltype1 = 'EtOCSPhenyl_porosity'
		else:
			moltype1 = 'EtOCSPhenyl_no_porosity'


	# if the dynamics are not found, prompt the user to choose
	try:
		return dynamics[moltype1]
	except:
		sim_dynam = ChooseDynamics()
		return dynamics[sim_dynam]

def getSimulationDynamicsType(sim_dynamics):
	filename = 'database/sim_dynamics.yml'
	with open(filename) as f:
		SimDynamics = yaml.load(f)

	for key in SimDynamics:
		if sim_dynamics == SimDynamics[key]:
			return key


def GetDihedralFile(moltype):
	saved_dihedral_files = {'TVSR1': 'TVSR1_dihedrals_all.txt',
							'TVSR2': 'TVSR2_dihedrals_all.txt',
							'TVSR3': 'TVSR3_dihedrals_all.txt',
							'TVSR4': 'TVSR4_dihedrals_all.txt',
							'TVSR6': 'TVSR6_dihedrals_all.txt'}
	return saved_dihedral_files[moltype]






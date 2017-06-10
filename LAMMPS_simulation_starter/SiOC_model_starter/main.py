import sys
import time

from createSimulation import CreateSimulation, CreateMixedSimulation

from utils import GetMolType, MakeDirectory, GetSimulationParams

#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
# simualation setup

print '\n\nWelcome to the LAMMPS starter program!\n\n'

sim_name = raw_input('\nName the simulation folder: ')
sim_dir = MakeDirectory(sim_name)

mixed_precursors = False

precursor_type = input('\nNew precursor or saved precursor:\n\n'+\
								'\t1) saved\n'+\
								'\t2) new\n\n'+\
								'\tchoice: ')

if precursor_type == 2:
	filename = raw_input("\nGive the new precursor's filename [format = molecule_name.xyz]: ")
	moltype = filename[:-4]

else:
	molID = input('\nType the number of the saved precursor:\n\n'+\
						'\t1) EtOCS\n'+\
						'\t2) EtOCSMethyl\n'+\
						'\t3) EtOCSVinyl\n'+\
						'\t4) EtOCSPhenyl\n'+\
						'\t5) 135Benzene\n'+\
						'\t6) 13Benzene\n'+\
						'\t7) 14Benzene\n'+\
						'\t8) MSSQ\n'+\
						'\t9) SiO2\n'+\
						'\t10) TVSR1\n'+\
						'\t11) TVSR2\n'+\
						'\t12) TVSR3\n'+\
						'\t13) TVSR4\n'+\
						'\t14) TVSR6\n'+\
						'\t15) mixed precursors\n\n'+\
						'\tchoice: ')

	if molID == 15:
		mixed_precursors = True

		molID1 = input('\nType the number of the first precursor:\n\n'+\
						'\t1) EtOCS\n'+\
						'\t2) EtOCSMethyl\n'+\
						'\t3) EtOCSVinyl\n'+\
						'\t4) EtOCSPhenyl\n'+\
						'\t5) 135Benzene\n'+\
						'\t6) 13Benzene\n'+\
						'\t7) 14Benzene\n'+\
						'\t8) MSSQ\n'+\
						'\t9) SiO2\n'+\
						'\t10) TVSR1\n'+\
						'\t11) TVSR2\n'+\
						'\t12) TVSR3\n'+\
						'\t13) TVSR4\n'+\
						'\t14) TVSR6\n'+\
						'\tchoice 1: ')

		molID2 = input('\nType the number of the second precursor:\n\n'+\
						'\t1) EtOCS\n'+\
						'\t2) EtOCSMethyl\n'+\
						'\t3) EtOCSVinyl\n'+\
						'\t4) EtOCSPhenyl\n'+\
						'\t5) 135Benzene\n'+\
						'\t6) 13Benzene\n'+\
						'\t7) 14Benzene\n'+\
						'\t8) MSSQ\n'+\
						'\t9) SiO2\n'+\
						'\t10) TVSR1\n'+\
						'\t11) TVSR2\n'+\
						'\t12) TVSR3\n'+\
						'\t13) TVSR4\n'+\
						'\t14) TVSR6\n'+\
						'\tchoice 2: ')
		moltype1 = GetMolType(molID1)
		moltype2 = GetMolType(molID2)
		filename1 = 'database/{}.xyz'.format(moltype1)
		filename2 = 'database/{}.xyz'.format(moltype2)

	else:
		moltype = GetMolType(molID)
		filename = 'database/{}.xyz'.format(moltype)

dihedral = input('\nSaved or custom dihedrals?\n\n'+\
							'\t1) auto-generate\n'+\
							'\t2) saved to database\n'+\
							'\t3) custom\n\n'+\
							'\tchoice: ')
if dihedral == 1:
	fname_dihedrals = False
elif dihedral == 2:
	fname_dihedrals = 'database/'
elif dihedral == 3:
	fname_dihedrals = raw_input("\nGive the dihedral filename (must be in lammps format): ")
else:
	raise RuntimeError, "Error: choice '%d' is not an option. Please select a valid option." %dihedral

if mixed_precursors:
	Nprecursors1 = input('\nNumber of %s precursors: ' %moltype1)
	Nprecursors2 = input('\nNumber of %s precursors: ' %moltype2)
else:
	Nprecursors = input('\nNumber of %s precursors: '%moltype)

sim_type = input('\nType the number of the simulation type:\n\n'+\
						'\t1) generate general data file\n'+\
						'\t2) single simulation (no porosity)\n'+\
						'\t3) single simulation (with porosity)\n'+\
						'\t4) network connectivity scan\n'+\
						'\t5) porosity scan\n'+\
						'\t6) one connectivity, multiple trials\n\n'+\
						'\tchoice: ')
sim_params = GetSimulationParams(sim_type)


#-------------------------------------------------------------------------------------------------#
# generate the files 

t0 = time.time()

if mixed_precursors:
	CreateMixedSimulation(sim_dir, filename1, moltype1, filename2,\
							moltype2, fname_dihedrals, Nprecursors1,\
							Nprecursors2, sim_type, sim_params)
else:
	CreateSimulation(sim_dir, filename, moltype, fname_dihedrals,\
						Nprecursors, sim_type, sim_params)

print 'Analyzed the precursor and generated the simulation files in %.4f seconds.' %(time.time()-t0)
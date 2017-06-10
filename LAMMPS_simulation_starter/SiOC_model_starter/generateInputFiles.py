from utils import LoadSimulationDynamics, getSimulationDynamicsType
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#

def WriteSimDescription(sim_dir, moltype, Nprecursors, sim_description, O_Si_ratios,\
						porosities, Nconditions, Ntrials_per_cond, sim_type):
	filename = '%ssimulation_description.txt' %sim_dir
	f = open(filename, 'w')
	f.write('Precursor: %s\n' %moltype)
	f.write('Number of precursors: %d\n' %Nprecursors)
	f.write('Simulation type: %s\n' %sim_description)
	f.write('Simulation dymamics: %s' %sim_type)
	for i in range(Nconditions):
		for j in range(Ntrials_per_cond):
			trial = '%s.%s' %(i+1, j+1)
			O_Si_ratio = O_Si_ratios[i]
			porosity = porosities[i]
			f.write('\nTrial %s\n' %trial)
			f.write('\tO:Si ratio: %.4f\n' %O_Si_ratio)
			f.write('\tPorosity: %.4f%%\n' %porosity)
	f.close()

def WriteMixedSimDescription(sim_dir, moltype1, Nprecursors1, moltype2, Nprecursors2,\
								sim_description, O_Si_ratios, porosities,\
								Nconditions, Ntrials_per_cond, sim_type):
	filename = '%ssimulation_description.txt' %sim_dir
	f = open(filename, 'w')
	f.write('Precursor 1: %s\n' %moltype1)
	f.write('Number of %s precursors: %d\n' %(moltype1, Nprecursors1))
	f.write('Precursor 2: %s\n' %moltype2)
	f.write('Number of %s precursors: %d\n' %(moltype2, Nprecursors2))
	f.write('Simulation type: %s\n' %sim_description)
	f.write('Simulation dymamics: %s' %sim_type)
	for i in range(Nconditions):
		for j in range(Ntrials_per_cond):
			trial = '%s.%s' %(i+1, j+1)
			O_Si_ratio = O_Si_ratios[i]
			porosity = porosities[i]
			f.write('\nTrial %s\n' %trial)
			f.write('\tO:Si ratio: %.4f\n' %O_Si_ratio)
			f.write('\tPorosity: %.4f%%\n' %porosity)
	f.close()
#-------------------------------------------------------------------------------------------------#

def WriteSimVars(f, inputfile_params, moltype1, moltype2=False):
	f.write('\n# trial parameter\n')
	f.write('variable a index ')
	for params in inputfile_params:
		f.write('%s ' %params['trial'])

	# write precursor for double (mixed) and single precursor simulations
	if moltype2:
		f.write('\n# %s precursor start moleculeID parameter\n' %moltype1)
		f.write('variable p index ')
		for params in inputfile_params:
			f.write('%d ' %params['precursor1_start'])

		f.write('\n# %s precursor end moleculeID parameter\n' %moltype1)
		f.write('variable q index ')
		for params in inputfile_params:
			f.write('%d ' %params['precursor1_end'])		

		f.write('\n# %s precursor start moleculeID parameter\n' %moltype2)
		f.write('variable x index ')
		for params in inputfile_params:
			f.write('%d ' %params['precursor2_start'])

		f.write('\n# %s precursor end moleculeID parameter\n' %moltype2)
		f.write('variable y index ')
		for params in inputfile_params:
			f.write('%d ' %params['precursor2_end'])	

	else:
		f.write('\n# %s precursor start moleculeID parameter\n' %moltype1)
		f.write('variable p index ')
		for params in inputfile_params:
			f.write('%d ' %params['precursor_start'])

		f.write('\n# %s precursor end moleculeID parameter\n' %moltype1)
		f.write('variable q index ')
		for params in inputfile_params:
			f.write('%d ' %params['precursor_end'])	

	f.write('\n# free O start atomID parameter\n')
	f.write('variable r index ')
	for params in inputfile_params:
		f.write('%d ' %params['freeO_start'])

	f.write('\n# free O end atomID parameter\n')
	f.write('variable s index ')
	for params in inputfile_params:
		f.write('%d ' %params['freeO_end'])

	f.write('\n# porogen start moleculeID parameter\n')
	f.write('variable t index ')
	for params in inputfile_params:
		f.write('%d ' %params['porogen_start'])

	f.write('\n# porogen end moleculeID parameter\n')
	f.write('variable v index ')
	for params in inputfile_params:
		f.write('%d ' %params['porogen_end'])

	f.write('\n\nlog log_$a.lammps\n\n')


def WriteSimSetup(f, inputfile_params):
	f.write('#########################################################################\n')
	f.write('dimension 3\n')
	f.write('atom_style full\n\n')

	f.write('# bonded potentials\n')
	f.write('bond_style harmonic\n')
	f.write('angle_style harmonic\n')
	f.write('dihedral_style hybrid opls class2\n\n')

	f.write('boundary p p p\n')
	f.write('units metal\n')
	f.write('neighbor 2 bin\n')
	f.write('neigh_modify check yes\n\n')

	f.write('timestep 0.001\n\n')
	f.write('read_data %s$a\n\n' %inputfile_params[0]['datafile_base'])


def WriteSimGroups(f, moltype1, moltype2=False):
	f.write('#########################################################################\n')

	if moltype2:
		f.write('# number of %s precursors (referenced by moleculeID)\n' %moltype1)
		f.write('group precursor1 molecule <> $p $q\n\n')

		f.write('# number of %s precursors (referenced by moleculeID)\n' %moltype2)
		f.write('group precursor2 molecule <> $x $y\n\n')

		f.write('group precursors union precursor1 precursor2\n\n')

	else:
		f.write('# number of %s precursors (referenced by moleculeID)\n' %moltype1)
		f.write('group precursors molecule <> $p $q\n\n')

	f.write('# number of free O atoms (referenced by atomID)\n')
	f.write('group freeO id <> $r $s\n\n')

	f.write('# number of porogen molecules (referenced by moleculeID)\n')
	f.write('group porogen molecule <> $t $v\n\n')

	f.write('group active union precursors freeO porogen\n')
	f.write('group activematrix union precursors freeO\n')
	f.write('group inactive subtract all active\n')
	f.write('delete_atoms group inactive\n\n')


def WriteSoftPotential(f, sim_dynamics):
	f.write('#########################################################################\n')
	f.write('#soft potential to push the atoms apart\n')
	f.write('pair_style soft 3\n')
	f.write('pair_coeff	* * 0 500 3\n\n')

	f.write('thermo_style custom step temp tave press pave lx ly lz vol epair etotal\n\n')

	nvt = sim_dynamics['soft_potential']
	f.write('fix 1\t\tactive nvt %d %d %.4f drag %.2f\n' %(nvt['Tstart'], nvt['Tstop'],\
															nvt['Tdamp'], nvt['drag']))
	f.write('thermo\t\t1000\n')
	f.write('run\t\t%d\n\n' %nvt['run'])

	f.write('thermo\t\t1\n')
	f.write('minimize\t\t1e-4 10000 10000\n\n')


def WriteSimulatedAnnealing(f, sim_dynamics, moltype1, moltype2=False):
	f.write('#########################################################################\n')
	f.write('# simulated annealing to create structure\n')
	f.write('pair_style sw\n')
	f.write('pair_coeff * * OCS.sw Si O C\n\n')

	f.write('thermo\t\t100\n\n')

	npt = sim_dynamics['sim_anneal_1']
	f.write('unfix\t\t1\n')
	f.write('fix\t\t1 active npt %d %d %.4f xyz %d %d %.4f drag %.2f\n' %(npt['Tstart'], npt['Tstop'],\
																			npt['Tdamp'], npt['Pstart'],\
																			npt['Pstop'], npt['Pdamp'],\
																			npt['drag']))
	f.write('run\t\t%d\n\n' %npt['run'])

	npt = sim_dynamics['sim_anneal_2']
	f.write('unfix\t\t1\n')
	f.write('fix\t\t1 active npt %d %d %.4f xyz %d %d %.4f drag %.2f\n' %(npt['Tstart'], npt['Tstop'],\
																			npt['Tdamp'], npt['Pstart'],\
																			npt['Pstop'], npt['Pdamp'],\
																			npt['drag']))
	f.write('run\t\t%d\n\n' %npt['run'])

	if moltype2:
		moltype = '%s_%s' %(moltype1, moltype2)
	else:
		moltype = moltype1

	f.write('unfix\t\t1\n')
	f.write('fix\t\t1 all npt 300 300 0.01 xyz 1 1 0.1 drag 2\n')
	f.write('fix\t\t2 active rdf 2000 %s_rdf_$a.txt 100 1 1 1 2 2 1 1 3 3 1 2 3 3 2 2 2 3 3\n' %moltype)
	f.write('dump\t\t1 porogen xyz 10000 porogen_*_$a.xyz\n')
	f.write('dump\t\t2 active xyz 10000 %s_*_$a.xyz\n' %moltype)
	f.write('dump\t\t3 activematrix xyz 10000 %s_matrix_*_$a.xyz\n' %moltype)
	f.write('run\t\t10000\n')
	f.write('undump\t\t1\n')
	f.write('undump\t\t2\n')
	f.write('undump\t\t3\n\n')


def WriteRemovePorogen(f, moltype1, moltype2=False):
	f.write('#########################################################################\n')
	f.write('# remove the porogen molecules from the simulation and equilibrate the structure\n\n')
	f.write('delete_atoms group porogen\n\n')

	if moltype2:
		moltype = '%s_%s' %(moltype1, moltype2)
	else:
		moltype = moltype1

	f.write('dump\t\t1 active xyz 5000 %s_*_$a.xyz\n\n' %moltype)	

	f.write('unfix\t\t1\n')
	f.write('fix\t\t1 all npt 300 300 0.01 xyz 1 1 0.1 drag 2\n')
	f.write('run\t\t20000\n\n')	


def WriteBulkModulus(f, sim_dynamics):
	f.write('#########################################################################\n')
	f.write('#apply setpwise pressure to measure modulus\n\n')

	npt = sim_dynamics['bulk_modulus']
	Plist = npt['Plist']
	for press in Plist:
		f.write('unfix\t\t1\n')
		f.write('fix\t\t1 all npt %d %d %.4f xyz %d %d %.4f drag %.2f\n'%(npt['Tstart'], npt['Tstop'],\
																			npt['Tdamp'], press,\
																			press, npt['Pdamp'],\
																			npt['drag']))
		f.write('run\t\t%d\n\n' %npt['run'])


def WriteSimJump(f, filename_base, mixed=False):
	f.write('#########################################################################\n')
	f.write('clear\n')
	f.write('next a\n')
	f.write('next p\n')
	f.write('next q\n')

	if mixed:
		f.write('next x\n')
		f.write('next y\n')

	f.write('next r\n')
	f.write('next s\n')
	f.write('next t\n')
	f.write('next v\n')
	f.write('jump %s' %filename_base)

#-------------------------------------------------------------------------------------------------#

def WriteInputFile(sim_dir, moltype, sim_type, sim_description, inputfile_params):
	# get the simulation dynamics based on the precursor
	sim_dynamics = LoadSimulationDynamics(sim_type, moltype)

	filename = '%sin.%s' %(sim_dir, moltype)
	filename_base = 'in.%s' %moltype
	f = open(filename, 'w')

	f.write('# %s -- %s -- with bulk modulus simulation\n' %(moltype, sim_description))
	
	WriteSimVars(f, inputfile_params, moltype)
	WriteSimSetup(f, inputfile_params)
	WriteSimGroups(f, moltype)
	WriteSoftPotential(f, sim_dynamics)
	WriteSimulatedAnnealing(f, sim_dynamics, moltype)
	WriteRemovePorogen(f, moltype)	
	WriteBulkModulus(f, sim_dynamics)
	WriteSimJump(f, filename_base)
	f.close()


def WriteMixedInputFile(sim_dir, moltype1, moltype2, sim_type, sim_description, inputfile_params):
	# get the simulation dynamics based on the precursor
	sim_dynamics = LoadSimulationDynamics(sim_type, moltype1, moltype2)

	filename = '%sin.%s_%s' %(sim_dir, moltype1, moltype2)
	filename_base = 'in.%s_%s' %(moltype1, moltype2)
	f = open(filename, 'w')

	f.write('# %s and %s -- %s -- with bulk modulus simulation\n' %(moltype1, moltype2, sim_description))
	
	WriteSimVars(f, inputfile_params, moltype1, moltype2)
	WriteSimSetup(f, inputfile_params)
	WriteSimGroups(f, moltype1, moltype2)
	WriteSoftPotential(f, sim_dynamics)
	WriteSimulatedAnnealing(f, sim_dynamics, moltype1, moltype2)
	WriteRemovePorogen(f, moltype1, moltype2)
	WriteBulkModulus(f, sim_dynamics)
	WriteSimJump(f, filename_base, mixed=True)
	f.close()


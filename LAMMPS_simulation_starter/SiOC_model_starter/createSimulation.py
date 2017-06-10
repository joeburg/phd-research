import math
import numpy
import sys
import time

from generateDataFiles import LoadData, LoadPorogen, AnalyzePrecursor,\
								ReplicateStructure, ReplicateMixedStructure,\
								WriteDataFile, WriteMixedDataFile, sortBondedPotentials

from generateInputFiles import WriteInputFile, WriteMixedInputFile,\
								WriteSimDescription, WriteMixedSimDescription

from utils import GetSimulationDescriptions, CopyFileToNewFolder
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#

def getOSiRatioList(sim_type, sim_params):
	O_Si_ratios = []

	if sim_params['O_Si_ratio']:
		# Note cannot have O:Si = 0
		O_Si_ratios.append(sim_params['O_Si_ratio'])

		if sim_params['Ntrials_per_cond'] > 1:
			O_Si_ratios = O_Si_ratios*sim_params['Ntrials_per_cond']

	elif sim_params['O_Si_ratio_start']:
		O_Si_start = sim_params['O_Si_ratio_start']
		O_Si_end = sim_params['O_Si_ratio_end']
		Nconds = sim_params['Nconditions']

		step = abs(O_Si_end - O_Si_start)/float(Nconds-1)

		for i in range(Nconds):
			O_Si_ratios.append(O_Si_start + i*step)

	elif sim_type == 1:
		# datafile only
		# no ratios given, so assume O:Si = 3
		O_Si_ratios.append(3)

	return O_Si_ratios


def getPorosityList(sim_type, sim_params):
	# given in terms of volume % porosity
	porosities = []

	if sim_params['porosity']:
		porosities.append(sim_params['porosity'])

	# note: many times the starting porosity wil be 0
	elif sim_params['porosity_end']:
		porosity_start = sim_params['porosity_start']
		porosity_end = sim_params['porosity_end']
		Nconds = sim_params['Nconditions']

		step = abs(porosity_end - porosity_start)/float(Nconds-1)

		for i in range(Nconds):
			porosities.append(porosity_start + i*step)

	elif sim_params['porosity'] == 0:
		porosities.append(0)

	elif sim_type == 1:
		# datafile only
		# no porosity given, so generate enough porogen to have 99% porosity
		porosities.append(99)

	return porosities


def NormalizeParamsList(O_Si_ratios, porosities):
	if len(O_Si_ratios)>1 and len(porosities)==1:
		porosities = porosities*len(O_Si_ratios)

	elif len(O_Si_ratios)==1 and len(porosities)>1:
		O_Si_ratios = O_Si_ratios*len(porosities)

	elif len(O_Si_ratios)>1 and len(porosities)>1:
		# this simulation is for a porosity scan with multiple
		# trials per condition 
		# the O:Si ratio should be less than the porosity
		if len(O_Si_ratios) < len(porosities):
			O_Si_ratios = [O_Si_ratios[0]]*len(porosities)
		else:
			porosities = [porosities[0]]*len(O_Si_ratios)

		# O_Si_ratios_norm = O_Si_ratios*len(porosities)

		# porosities_norm = []
		# for i in range(len(porosities)):
		# 	trial = []
		# 	trial.append(porosities[i])
		# 	trial = trial*len(O_Si_ratios)
		# 	porosities_norm += trial

		# O_Si_ratios = O_Si_ratios_norm
		# porosities = porosities_norm

	return (O_Si_ratios, porosities)


def getPrecursorAtoms(atom_metrics, Nprecursors):
	NSi = atom_metrics['NSi_per_precursor']*Nprecursors
	NO = atom_metrics['NO_per_precursor']*Nprecursors
	NC = atom_metrics['NC_per_precursor']*Nprecursors
	return (NSi, NO, NC)

def getNoxygen(O_Si_ratio, NSi):
	return int(O_Si_ratio*NSi)

def VolSphere(r):
	return 4.0/3*math.pi*r**3

def getNporogen(NSi, NO, NC, porosity):
	# porosity given in terms of % (convert to decimal)
	porosity = 0.01*porosity

	VSi = VolSphere(2.1)
	VC = VolSphere(1.7)
	VO = VolSphere(1.52)

	NC_porogen = 40.0 # extend if porogen molecule changes

	Nporogen = (1.0/(VC*NC_porogen))*(porosity/(1.0-porosity))*\
				(VSi*NSi + VC*NC + VO*NO)
	return int(Nporogen)


#-------------------------------------------------------------------------------------------------#
def CreateSimulation(sim_dir, precursor_filename, moltype, fname_dihedrals, Nprecursors, sim_type, sim_params):

	O_Si_ratios = getOSiRatioList(sim_type, sim_params)
	porosities = getPorosityList(sim_type, sim_params)

	# ensure the O:Si and porosities have the same length
	O_Si_ratios, porosities = NormalizeParamsList(O_Si_ratios, porosities)

	# load the new file or saved files
	data, atom_metrics = LoadData(precursor_filename)
	NSi_prec, NO_prec, NC_prec = getPrecursorAtoms(atom_metrics, Nprecursors)

	# analyze the precursor structure
	# precursor_structure = AnalyzePrecursor(moltype, data, fname_dihedrals)
	PrecursorStructure, BondedPotentials = AnalyzePrecursor(moltype, data, fname_dihedrals)

	# load the porogen template (Note: the bond coeffs may change)
	# porogen_templates, bond_coeffs = LoadPorogen(PrecursorStructure['bond_coeffs'])
	# PrecursorStructure['bond_coeffs'] = bond_coeffs
	PorogenTemplates, BondedPotentials = LoadPorogen(BondedPotentials)

	# sort the bonded potentials 
	BondedPotentials = sortBondedPotentials(BondedPotentials)

	# get total number of simulations to run
	Nconditions = sim_params['Nconditions']
	Ntrials_per_cond = sim_params['Ntrials_per_cond']

	inputfile_params = []
	# genereate and write the LAMMPS data file for each simulation
	for i in range(Nconditions):
		for j in range(Ntrials_per_cond):
			trial = '%d.%d' %(i+1, j+1)

			NO_free = getNoxygen(O_Si_ratios[i], NSi_prec)
			NO = NO_prec + NO_free
			Nporogen = getNporogen(NSi_prec, NO, NC_prec, porosities[i])

			# replicate the basic structures (precursor, free O, porogen)
			atoms, bonds, angles,\
			dihedrals, cell_length = ReplicateStructure(PrecursorStructure,\
														PorogenTemplates,\
														Nprecursors,\
														NO, Nporogen) 	

			# write out the LAMMPS data file
			datafile_base = WriteDataFile(sim_dir, Nprecursors, moltype, trial,\
											BondedPotentials, atoms, bonds,\
											angles, dihedrals, cell_length)

			# parameters from the data file that will be used to generate the input file
			# data file order: precursors, free O and then porogen molecules 
			# NOTE: for simulations with 0% porosity, using a start molecleID of Nprec+1 and an 
			# end moleculeID of Nprec (since Nporogen=0), leads to 0 porogen molecules in LAMMPS,
			# so current logic works! 
			params = {'precursor_start': 1,
						'precursor_end': Nprecursors,
						'freeO_start': Nprecursors*PrecursorStructure['metrics']['Natoms']+1,
						'freeO_end': Nprecursors*PrecursorStructure['metrics']['Natoms']+NO_free,
						'porogen_start': Nprecursors+1,
						'porogen_end': Nprecursors+Nporogen,
						'trial': trial,
						'datafile_base': datafile_base}

			inputfile_params.append(params)

	# write the LAMMPS input file
	if sim_type != 1: # dont write an input file when just generating a general data file 
		sim_description = GetSimulationDescriptions(sim_type)
		WriteInputFile(sim_dir, moltype, sim_type, sim_description, inputfile_params)

		# write the simulation description file
		WriteSimDescription(sim_dir, moltype, Nprecursors, sim_description, O_Si_ratios,\
							porosities, Nconditions, Ntrials_per_cond, sim_type)
	
		# copy the hostfile to the simulation dir
		CopyFileToNewFolder('database/hostfile', sim_dir)

		# copy the potential to the simulation dir
		CopyFileToNewFolder('database/OCS.sw', sim_dir)

#-------------------------------------------------------------------------------------------------#
def CreateMixedSimulation(sim_dir, precursor_filename1, moltype1, precursor_filename2, moltype2,\
							fname_dihedrals, Nprecursors1, Nprecursors2, sim_type, sim_params):
	O_Si_ratios = getOSiRatioList(sim_type, sim_params)
	porosities = getPorosityList(sim_type, sim_params)

	O_Si_ratios, porosities = NormalizeParamsList(O_Si_ratios, porosities)

	# load the fist precursor
	data1, atom_metrics1 = LoadData(precursor_filename1)
	NSi_prec1, NO_prec1, NC_prec1 = getPrecursorAtoms(atom_metrics1, Nprecursors1)
	PrecursorStructure1, BondedPotentials = AnalyzePrecursor(moltype1, data1, fname_dihedrals)

	# load the second precursor
	data2, atom_metrics2 = LoadData(precursor_filename2)
	NSi_prec2, NO_prec2, NC_prec2 = getPrecursorAtoms(atom_metrics2, Nprecursors2)
	PrecursorStructure2, BondedPotentials = AnalyzePrecursor(moltype2, data2, fname_dihedrals,\
																BondedPotentials)

	# load the porogen template (Note: the bond coeffs may change)
	PorogenTemplates, BondedPotentials = LoadPorogen(BondedPotentials)
	
	# sort the bonded potentials 
	BondedPotentials = sortBondedPotentials(BondedPotentials)

	# get total number of simulations to run
	Nconditions = sim_params['Nconditions']
	Ntrials_per_cond = sim_params['Ntrials_per_cond']

	inputfile_params = []
	# genereate and write the LAMMPS data file for each simulation
	for i in range(Nconditions):
		for j in range(Ntrials_per_cond):
			trial = '%d.%d' %(i+1, j+1)
			
			NO_free = getNoxygen(O_Si_ratios[i], NSi_prec1+NSi_prec2)
			NO = NO_prec1 + NO_prec2 + NO_free
			Nporogen = getNporogen(NSi_prec1+NSi_prec2, NO, NC_prec1+NC_prec2, porosities[i])

			# replicate the basic structures (precursor, free O, porogen)
			atoms, bonds, angles,\
			dihedrals, cell_length = ReplicateMixedStructure(PrecursorStructure1,\
																PrecursorStructure2,\
																PorogenTemplates,\
																Nprecursors1,\
																Nprecursors2,\
																NO, Nporogen) 	

			# write out the LAMMPS data file
			datafile_base = WriteMixedDataFile(sim_dir, trial, BondedPotentials,\
												Nprecursors1, moltype1, Nprecursors2,\
												moltype2, atoms, bonds, angles, dihedrals,\
												cell_length)

			# parameters from the data file that will be used to generate the input file
			# data file order: precursors, free O and then porogen molecules 
			# NOTE: for simulations with 0% porosity, using a start molecleID of Nprec+1 and an 
			# end moleculeID of Nprec (since Nporogen=0), leads to 0 porogen molecules in LAMMPS,
			# so current logic works! 
			params = {'precursor1_start': 1,
						'precursor1_end': Nprecursors1,
						'precursor2_start': Nprecursors1+1,
						'precursor2_end': Nprecursors1+Nprecursors2,
						'freeO_start': Nprecursors1*PrecursorStructure1['metrics']['Natoms']+\
										Nprecursors2*PrecursorStructure2['metrics']['Natoms']+1,
						'freeO_end': Nprecursors1*PrecursorStructure1['metrics']['Natoms']+\
										Nprecursors2*PrecursorStructure2['metrics']['Natoms']+NO_free,
						'porogen_start': Nprecursors1+Nprecursors2+1,
						'porogen_end': Nprecursors1+Nprecursors2+Nporogen,
						'trial': trial,
						'datafile_base': datafile_base}

			inputfile_params.append(params)

	# write the LAMMPS input file
	if sim_type != 1: # dont write an input file when just generating a general data file 
		sim_description = GetSimulationDescriptions(sim_type)
		WriteMixedInputFile(sim_dir, moltype1, moltype2, sim_type, sim_description, inputfile_params)

		# write the simulation description file
		WriteMixedSimDescription(sim_dir, moltype1, Nprecursors1, moltype2, Nprecursors2,\
									sim_description, O_Si_ratios, porosities,\
									Nconditions, Ntrials_per_cond, sim_type)
	
		# copy the hostfile to the simulation dir
		CopyFileToNewFolder('database/hostfile', sim_dir)

		# copy the potential to the simulation dir
		CopyFileToNewFolder('database/OCS.sw', sim_dir)


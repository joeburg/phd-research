import math
import numpy
import operator
import sys
import time

from StructureAnalysisFuctions import ComputeAtoms, ComputeBonds, ComputeAngles,\
										ComputeDihedrals, LoadDihedrals

from utils import GetDihedralFile

#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
# utility functions

def Distance(dx,dy,dz):
	return (dx*dx + dy*dy + dz*dz)**0.5

def Angle(atom1,vertex,atom2):
	# use law of cosines
	D1V = Distance(atom1[1]-vertex[1],atom1[2]-vertex[2],atom1[3]-vertex[3])
	D2V = Distance(atom2[1]-vertex[1],atom2[2]-vertex[2],atom2[3]-vertex[3])
	D12 = Distance(atom1[1]-atom2[1],atom1[2]-atom2[2],atom1[3]-atom2[3])
	# return angle in degrees	
	return math.acos((D1V*D1V + D2V*D2V - D12*D12)/(2*D1V*D2V))*(180.0/math.pi)

def BondLength(atom1, atom2):
	# this function computes the bond length between atoms in lammps format
	dx = abs(atom1[4] - atom2[4])
	dy = abs(atom1[5] - atom2[5])
	dz = abs(atom1[6] - atom2[6])
	return Distance(dx,dy,dz)

def getBondType(bond_length, bond_coeffs):
	# only call this function if the bond length is know to be 
	# within 0.01 A of a given bond coeff
	for bond in bond_coeffs:
		# use same bond coeff if within 0.01 (Angstrom) tolerance
		if abs(bond_length - bond) <= 0.01:
			return bond_coeffs[bond]

def getAngleType(angle_deg, angle_coeffs):
	# only call this function if the angle is know to be 
	# within 0.1 degrees of a given angle coeff
	for angle in angle_coeffs:
		# use same bond coeff if within 0.01 (Angstrom) tolerance
		if abs(angle_deg - angle) <= 0.1:
			return angle_coeffs[angle]

def updateBondType(bond_length, bond_coeffs):
	for bond in bond_coeffs:
		# use same bond coeff if within 0.01 (Angstrom) tolerance
		if abs(bond_length - bond) <= 0.01:
			return bond_coeffs
	# if no match found, updated the bond_coeffs
	# the bond coeffs are from 1 to N
	bond_coeffs[bond_length] = len(bond_coeffs)+1
	return bond_coeffs


def updateAngleType(angle_deg, angle_coeffs):
	for angle in angle_coeffs:
		# use same bond coeff if within 0.1 (degree) tolerance
		if abs(angle_deg - angle) <= 0.1:
			return angle_coeffs
	# if no match found, return none
	# the angle coeffs are from 1 to N
	angle_coeffs[angle_deg] = len(angle_coeffs)+1
	return angle_coeffs

def updateDihedralType(dihedral, dihedral_coeffs):
	for fit_dihedral in dihedral_coeffs:
		# dihedrals are accessed by name e.g 'SiCCC'
		if dihedral == fit_dihedral:
			return dihedral_coeffs
	# if no match found, raise exception telling the user that the dihedral must be parameterized
	raise RuntimeError, "The '%s' dihedral was not found in the database. You must parameterize it!" %dihedral


def updatePotentialDispatcher(icoeffs, saved_coeffs, itype):
	# check and update the coeffs if necessary
	PotentailDispatcher = {'bonds' : updateBondType,
							'angles': updateAngleType,
							'dihedrals': updateDihedralType}
	for coeff in icoeffs:
		saved_coeffs = PotentailDispatcher[itype](coeff, saved_coeffs)
	return saved_coeffs 


def updateBondedPotentials(BondedPotentials, ibond_coeffs, iangle_coeffs, idihedral_coeffs):
	if ibond_coeffs:
		bond_coeffs = BondedPotentials['bond_coeffs']
		bond_coeffs = updatePotentialDispatcher(ibond_coeffs, bond_coeffs, 'bonds')
		BondedPotentials['bond_coeffs'] = bond_coeffs

	if iangle_coeffs: 
		angle_coeffs = BondedPotentials['angle_coeffs']
		angle_coeffs = updatePotentialDispatcher(iangle_coeffs, angle_coeffs, 'angles')
		BondedPotentials['angle_coeffs'] = angle_coeffs	

	if idihedral_coeffs:
		dihedral_coeffs = BondedPotentials['dihedral_coeffs']
		dihedral_coeffs = updatePotentialDispatcher(idihedral_coeffs, dihedral_coeffs, 'dihedrals')
		BondedPotentials['dihedral_coeffs'] = dihedral_coeffs

	return BondedPotentials


def sortPotential(coeffs, sortby):
	# sortby refers to sorting by the key (0) or value (1) assuming the 
	# coeffs are a dictionary
	# note: sorted() returns a list of tuples
	return sorted(coeffs.items(), key=operator.itemgetter(sortby))


def sortBondedPotentials(BondedPotentials):
	# sort bond, angle, dihedral coeffs by increasing ID number (returns list of tuples)
	BondedPotentials['bond_coeffs'] = sortPotential(BondedPotentials['bond_coeffs'], 1)
	BondedPotentials['angle_coeffs'] = sortPotential(BondedPotentials['angle_coeffs'], 1)
	BondedPotentials['dihedral_coeffs'] = sortPotential(BondedPotentials['dihedral_coeffs'], 0)
	return BondedPotentials	


def updateBonds(BondedPotentials, data, bonds):
	bond_coeffs = BondedPotentials['bond_coeffs']

	# bond structure is [ID, bondtype, atom1_idx, atom2_idx]
	# note: atom_idx = data_idx+1	
	for i in range(len(bonds)):
		atom1 = bonds[i][2]-1
		atom2 = bonds[i][3]-1
		# update the bond type
		dx = abs(data[atom1][1] - data[atom2][1])
		dy = abs(data[atom1][2] - data[atom2][2])
		dz = abs(data[atom1][3] - data[atom2][3])
		bonds[i][1] = getBondType(Distance(dx,dy,dz), bond_coeffs)
	return bonds


def updateAngles(BondedPotentials, data, angles):
	angle_coeffs = BondedPotentials['angle_coeffs']

	# angle structure is [ID, angletype, atom1_idx, atom2_idx (vertex), atom3_idx]
	# note: atom_idx = data_idx+1
	for i in range(len(angles)):
		# print angle
		# update the angle type
		angle = Angle(data[angles[i][2]-1], data[angles[i][3]-1], data[angles[i][4]-1])
		angles[i][1] = getAngleType(angle, angle_coeffs)
	return angles

#-------------------------------------------------------------------------------------------------#
# load the structure file

def ShiftToPositiveCoords(data):
	''' ensures all coordinate values are positive '''
	coord_shfit = abs(numpy.min(data))
	data[:,1:] = data[:,1:] + coord_shfit
	return data

def LoadData(inputfile):
	""" reads .xyz file """ 
	data = []
	f = open(inputfile)
	f.readline()

	NSi = 0
	NO = 0
	NC = 0

	while True:
		fields = f.readline().strip().split()
		if fields:
			atomtype = int(fields[0])
			xcoord = float(fields[1])
			ycoord = float(fields[2])
			zcoord = float(fields[3])

			# determine type of atom
			if atomtype == 1:
				NSi += 1
			elif atomtype == 2:
				NO += 1
			elif atomtype == 3:
				NC += 1

			# populate the data array
			data.append([atomtype,xcoord,ycoord,zcoord])
		else:
			break
	f.close()
	data = numpy.array(data)
	data = ShiftToPositiveCoords(data)

	atom_metrics = {'NSi_per_precursor': NSi,
					'NO_per_precursor': NO,
					'NC_per_precursor': NC}

	return data, atom_metrics

# load the porogen structure 
def LoadPorogen(BondedPotentials):
	# load the atom template 
	filename = 'database/porogen40_atoms.txt'
	atomtemplate = []
	f = open(filename)
	f.readline()

	while True:
		fields = f.readline().strip().split()
		if fields:
			ID = int(fields[0])
			molID = int(fields[1])
			atomtype = int(fields[2])
			charge = int(fields[3])
			xcoord = float(fields[4])
			ycoord = float(fields[5])
			zcoord = float(fields[6])

			atomtemplate.append([ID, molID, atomtype, charge, xcoord, ycoord, zcoord])
		else:
			break
	f.close()
	atomtemplate = numpy.array(atomtemplate)

	# compare the C-C bond length in the porogen and set the bond_coeff
	d = BondLength(atomtemplate[0], atomtemplate[1])
	bond_coeffs = BondedPotentials['bond_coeffs']
	bond_coeffs = updateBondType(d, bond_coeffs)
	BondedPotentials['bond_coeffs'] = bond_coeffs
	bond_type = getBondType(d, bond_coeffs)

	# load the bond template 
	filename = 'database/porogen40_bonds.txt'
	bondtemplate = []
	f = open(filename)
	f.readline()

	while True:
		fields = f.readline().strip().split()
		if fields:
			ID = int(fields[0])
			atom1 = int(fields[2])
			atom2 = int(fields[3])

			bondtemplate.append([ID, bond_type, atom1, atom2])
		else:
			break
	f.close()
	bondtemplate = numpy.array(bondtemplate)

	PorogenTemplates = {'atoms': atomtemplate,
						'bonds': bondtemplate}

	return (PorogenTemplates, BondedPotentials)

#-------------------------------------------------------------------------------------------------#
# analyze the precursor structure 

def getPrecursorMetrics(atoms, bonds, angles, dihedrals):
	precursor_metrics = {'Natoms': len(atoms),
							'Nbonds': len(bonds),
							'Nangles': len(angles),
							'Ndihedrals': len(dihedrals)}
	return precursor_metrics


def AnalyzePrecursor(moltype, data, fname_dihedrals, BondedPotentials=False):
	# get the atoms, bonds and angles
	atoms = ComputeAtoms(data)
	bonds, AM, bond_coeffs = ComputeBonds(data)
	angles, angle_coeffs = ComputeAngles(data, AM)

	# get the dihedrals
	if not fname_dihedrals:
		# compute and generate the dihedrals 
		dihedrals, dihedral_coeffs = ComputeDihedrals(data, AM)
	else:
		# if the database path is given, find the relevant precursor filename
		if fname_dihedrals == 'database/':
			fname_dihedrals = '%s%s' %(fname_dihedrals, GetDihedralFile(moltype))
		# otherwise, load the filename of the dihedrals provided by the user
		dihedrals, dihedral_coeffs = LoadDihedrals(fname_dihedrals)

	# update the bonded potentials
	if BondedPotentials:
		# update the existing bonded potentials
		BondedPotentials = updateBondedPotentials(BondedPotentials, bond_coeffs,\
													angle_coeffs, dihedral_coeffs)

		# update the bond-types and angle-types for the autogenerated bonds and angles
		bonds = updateBonds(BondedPotentials, data, bonds)
		angles = updateAngles(BondedPotentials, data, angles)

	else:
		# initialize and update the bonded potentials
		BondedPotentials = {}
		BondedPotentials['bond_coeffs'] = bond_coeffs
		BondedPotentials['angle_coeffs'] = angle_coeffs
		BondedPotentials['dihedral_coeffs'] = dihedral_coeffs

	# get the precursor metrics 
	precursor_metrics = getPrecursorMetrics(atoms, bonds, angles, dihedrals)
	PrecursorStructure = {'atoms': atoms,
							'bonds': bonds,
							'angles': angles,
							'dihedrals': dihedrals,
							'metrics': precursor_metrics}

	return (PrecursorStructure, BondedPotentials)

#-------------------------------------------------------------------------------------------------#
# functions to replicate the structure 

def ReplicatePrecursor(Nprecursors, PrecursorStructure, atoms, bonds, angles, dihedrals,\
						atom_idx, mol_idx, bond_idx, angle_idx, dihedral_idx, L):
	# get the precursor templates and metrics 
	Natoms_precursor = PrecursorStructure['metrics']['Natoms']
	Nbonds_precursor = PrecursorStructure['metrics']['Nbonds']
	Nangles_precursor = PrecursorStructure['metrics']['Nangles']
	Ndihrdrals_precursor = PrecursorStructure['metrics']['Ndihedrals']
	precursor_atom_template = PrecursorStructure['atoms']
	precursor_bond_template = PrecursorStructure['bonds']
	precursor_angle_template = PrecursorStructure['angles']
	precursor_dihedral_template = PrecursorStructure['dihedrals']

	# generate random positions [0, 1] for the atom positions 
	# and then scale by the simulation cell length
	precursor_rand = numpy.random.rand(Nprecursors, 3)*L

	# replicate the precursor
	for i in range(Nprecursors):
		# replicate the precursor atoms 
		# [atomID, moleculeID, atom_type, charge, x, y, z]
		if precursor_atom_template:
			atom_trans = numpy.array([atom_idx+Natoms_precursor*i,\
										mol_idx+i, 0, 0,\
										precursor_rand[i][0],\
										precursor_rand[i][1],\
										precursor_rand[i][2]])
			# constant transformation for the whole precursor
			addatoms = numpy.tile(atom_trans, (Natoms_precursor,1))

			atoms[atom_idx+Natoms_precursor*i:\
					atom_idx+Natoms_precursor*(i+1), :] = precursor_atom_template + addatoms

		# replicate the precursor bonds
		# [bondID, bond_type, atom1, atom2]
		if precursor_bond_template:
			bond_trans = numpy.array([bond_idx+Nbonds_precursor*i,\
										0,\
										atom_idx+Natoms_precursor*i,\
										atom_idx+Natoms_precursor*i])

			addbonds = numpy.tile(bond_trans, (Nbonds_precursor,1))

			bonds[bond_idx+Nbonds_precursor*i:\
					bond_idx+Nbonds_precursor*(i+1), :] = precursor_bond_template + addbonds

		# replicate the precursor angles
		# [angleID, angle_type, atom1, atom2, atom3] 
		if precursor_angle_template:
			angle_trans = numpy.array([angle_idx+Nangles_precursor*i,\
										0,\
										atom_idx+Natoms_precursor*i,\
										atom_idx+Natoms_precursor*i,\
										atom_idx+Natoms_precursor*i])

			addangles = numpy.tile(angle_trans, (Nangles_precursor,1))

			angles[angle_idx+Nangles_precursor*i:\
					angle_idx+Nangles_precursor*(i+1), :] = precursor_angle_template + addangles

		# replicate the precursor dihedrals 
		# [dihedralID, dihedral_type, atom1, atom2, atom3, atom4]
		if precursor_dihedral_template:
			dihedral_trans = numpy.array([dihedral_idx+Ndihrdrals_precursor*i,\
											0,\
											atom_idx+Natoms_precursor*i,\
											atom_idx+Natoms_precursor*i,\
											atom_idx+Natoms_precursor*i,\
											atom_idx+Natoms_precursor*i])

			adddihedrals = numpy.tile(dihedral_trans, (Ndihrdrals_precursor,1))

			dihedrals[dihedral_idx+Ndihrdrals_precursor*i:\
						dihedral_idx+Ndihrdrals_precursor*(i+1), :] = precursor_dihedral_template + adddihedrals

	return (atoms, bonds, angles, dihedrals)

def ReplicatePorogen(Nporogen, PorogenTemplates, atom_idx, mol_idx, bond_idx, atoms, bonds, L):
	# get the porogen template data
	porogen_atom_template = PorogenTemplates['atoms']
	porogen_bond_template = PorogenTemplates['bonds']
	Natoms_porogen = len(porogen_atom_template)
	Nbonds_porogen = len(porogen_bond_template)

	# generate random positions [0, 1] for the atom positions 
	# and then scale by the simulation cell length
	porogen_rand = numpy.random.rand(Nporogen, 3)*L

	# replicate the porogen molecules
	for i in range(Nporogen):
		# replicate the porogen atoms 
		atom_trans = numpy.array([atom_idx+Natoms_porogen*i,\
									mol_idx+i,\
									0, 0,\
									porogen_rand[i][0],\
									porogen_rand[i][1],\
									porogen_rand[i][2]])

		addatoms = numpy.tile(atom_trans, (Natoms_porogen,1))

		atoms[atom_idx+Natoms_porogen*i:\
				atom_idx+Natoms_porogen*(i+1), :] = porogen_atom_template + addatoms

		# repliate the porogen bonds
		bond_trans = numpy.array([bond_idx+Nbonds_porogen*i,\
									0,\
									atom_idx+Natoms_porogen*i,\
									atom_idx+Natoms_porogen*i])

		addbonds = numpy.tile(bond_trans, (Nbonds_porogen,1))

		bonds[bond_idx+Nbonds_porogen*i:\
				bond_idx+Nbonds_porogen*(i+1), :] = porogen_bond_template + addbonds

	return (atoms, bonds)

def ReplicateFreeO(NO, atom_idx, atoms, L):
	# generate random positions [0, 1] for the atom positions 
	# and then scale by the simulation cell length
	Oatoms_rand = numpy.random.rand(NO, 3)*L

	# replicate the free O atoms 
	for i in range(NO):
		atoms[atom_idx+i,:] = numpy.array([atom_idx+i+1,\
											0, 2, -2,\
											Oatoms_rand[i][0],\
											Oatoms_rand[i][1],\
											Oatoms_rand[i][2]])
	return atoms 

def CreateSimulationBox(atoms):
	# ensure that the final simulation cell contains all the atoms 
	x_min = min(atoms[:,4])
	x_max = max(atoms[:,4])
	y_min = min(atoms[:,5])
	y_max = max(atoms[:,5]) 
	z_min = min(atoms[:,6])
	z_max = max(atoms[:,6])

	cell_min = min(x_min, y_min, z_min) - 5
	cell_max = max(x_max, y_max, z_max) + 5
	return (cell_min, cell_max)

#-------------------------------------------------------------------------------------------------#
# replicate precursor structure and porogen molecules 

def ReplicateStructure(PrecursorStructure, PorogenTemplates, Nprecursors, NO, Nporogen):
	# get the precursor templates and metrics 
	Natoms_precursor = PrecursorStructure['metrics']['Natoms']
	Nbonds_precursor = PrecursorStructure['metrics']['Nbonds']
	Nangles_precursor = PrecursorStructure['metrics']['Nangles']
	Ndihrdrals_precursor = PrecursorStructure['metrics']['Ndihedrals']

	# get the number of porogen atoms and bonds
	Natoms_porogen = len(PorogenTemplates['atoms'])
	Nbonds_porogen = len(PorogenTemplates['bonds'])

	# get the total number of atoms, bonds, angles, dihedrals
	Natoms = Nprecursors*Natoms_precursor + Nporogen*Natoms_porogen + NO
	Nbonds = Nprecursors*Nbonds_precursor + Nporogen*Nbonds_porogen
	Nangles = Nprecursors*Nangles_precursor
	Ndihedrals = Nprecursors*Ndihrdrals_precursor	

	# generate the dimensions of the cubic simulation cell 
	# we will allocate 120 A^3 / atom which is ~10*V_Si_atom
	L = (Natoms*120.0)**(1.0/3)

	# create matricies to store the atom, bond, angle, dihedral data
	atoms = numpy.zeros((Natoms, 7))
	bonds = numpy.zeros((Nbonds, 4))
	angles = numpy.zeros((Nangles, 5))
	dihedrals = numpy.zeros((Ndihedrals, 6))

	# replicate the precursor
	atoms, bonds,\
	angles, dihedrals = ReplicatePrecursor(Nprecursors, PrecursorStructure,\
											atoms, bonds, angles, dihedrals,\
											0, 0, 0, 0, 0, L)

	# replicate the free O atoms 
	atom_idx = Nprecursors*Natoms_precursor
	atoms = ReplicateFreeO(NO, atom_idx, atoms, L)

	# replicate the porogen molecules
	atom_idx = Nprecursors*Natoms_precursor + NO
	mol_idx = Nprecursors
	bond_idx = Nprecursors*Nbonds_precursor
	atoms, bonds = ReplicatePorogen(Nporogen, PorogenTemplates, atom_idx,\
										mol_idx, bond_idx, atoms, bonds, L)

	# get the dimensions of the simulation cell
	cell_size = CreateSimulationBox(atoms)

	return atoms, bonds, angles, dihedrals, cell_size

#-------------------------------------------------------------------------------------------------#
# replicate mixed precursor structures and porogen molecules 
def ReplicateMixedStructure(PrecursorStructure1, PrecursorStructure2, PorogenTemplates,\
								Nprecursors1, Nprecursors2, NO, Nporogen):
	# get the frist precursor templates and metrics 
	Natoms_precursor1 = PrecursorStructure1['metrics']['Natoms']
	Nbonds_precursor1 = PrecursorStructure1['metrics']['Nbonds']
	Nangles_precursor1 = PrecursorStructure1['metrics']['Nangles']
	Ndihrdrals_precursor1 = PrecursorStructure1['metrics']['Ndihedrals']

	# get the second precursor templates and metrics 
	Natoms_precursor2 = PrecursorStructure2['metrics']['Natoms']
	Nbonds_precursor2 = PrecursorStructure2['metrics']['Nbonds']
	Nangles_precursor2 = PrecursorStructure2['metrics']['Nangles']
	Ndihrdrals_precursor2 = PrecursorStructure2['metrics']['Ndihedrals']

	# get the porogen template
	Natoms_porogen = len(PorogenTemplates['atoms'])
	Nbonds_porogen = len(PorogenTemplates['bonds'])

	# get the total number of atoms, bonds, angles, dihedrals
	Natoms = Nprecursors1*Natoms_precursor1 + Nprecursors2*Natoms_precursor2 +\
				Nporogen*Natoms_porogen + NO
	Nbonds = Nprecursors1*Nbonds_precursor1 + Nprecursors2*Nbonds_precursor2 +\
				Nporogen*Nbonds_porogen
	Nangles = Nprecursors1*Nangles_precursor1 + Nprecursors2*Nangles_precursor2
	Ndihedrals = Nprecursors1*Ndihrdrals_precursor1 + Nprecursors2*Ndihrdrals_precursor2

	# generate the dimensions of the cubic simulation cell 
	# we will allocate 120 A^3 / atom which is ~10*V_Si_atom
	L = (Natoms*120.0)**(1.0/3)

	# create matricies to store the atom, bond, angle, dihedral data
	atoms = numpy.zeros((Natoms, 7))
	bonds = numpy.zeros((Nbonds, 4))
	angles = numpy.zeros((Nangles, 5))
	dihedrals = numpy.zeros((Ndihedrals, 6))

	# replicate the first precursor
	atoms, bonds,\
	angles, dihedrals = ReplicatePrecursor(Nprecursors1, PrecursorStructure1,\
											atoms, bonds, angles, dihedrals,\
											0, 0, 0, 0, 0, L)

	# replicate the second precursor
	atom_idx = Natoms_precursor1*Nprecursors1
	mol_idx = Nprecursors1
	bond_idx = Nbonds_precursor1*Nprecursors1
	angle_idx = Nangles_precursor1*Nprecursors1
	dihedral_idx = Ndihrdrals_precursor1*Nprecursors1
	atoms, bonds,\
	angles, dihedrals = ReplicatePrecursor(Nprecursors2, PrecursorStructure2,\
											atoms, bonds, angles, dihedrals,\
											atom_idx, mol_idx, bond_idx, angle_idx,\
											dihedral_idx, L)

	# replicate the free O atoms 
	atom_idx = Nprecursors1*Natoms_precursor1 + Nprecursors2*Natoms_precursor2
	atoms = ReplicateFreeO(NO, atom_idx, atoms, L)

	# replicate the porogen molecules
	atom_idx = Nprecursors1*Natoms_precursor1 + Nprecursors2*Natoms_precursor2 + NO
	mol_idx = Nprecursors1 + Nprecursors2
	bond_idx = Nprecursors1*Nbonds_precursor1 + Nprecursors2*Nbonds_precursor2
	atoms, bonds = ReplicatePorogen(Nporogen, PorogenTemplates, atom_idx,\
										mol_idx, bond_idx, atoms, bonds, L)

	# get the dimensions of the simulation cell
	cell_size = CreateSimulationBox(atoms)

	return atoms, bonds, angles, dihedrals, cell_size


#-------------------------------------------------------------------------------------------------#
# functions to write the data files

def WriteSimulationParams(f, BondedPotentials, atoms, bonds, angles, dihedrals, cell_length):
	bond_coeffs = BondedPotentials['bond_coeffs']
	angle_coeffs = BondedPotentials['angle_coeffs']
	dihedral_coeffs = BondedPotentials['dihedral_coeffs']
	cell_min, cell_max = cell_length

	f.write('%d %d xlo xhi\n' %(cell_min, cell_max))
	f.write('%d %d ylo yhi\n' %(cell_min, cell_max))
	f.write('%d %d zlo zhi\n\n' %(cell_min, cell_max))

	f.write('%d atoms\n' %len(atoms))
	f.write('%d bonds\n' %len(bonds))
	f.write('%d angles\n' %len(angles))
	f.write('%d dihedrals\n\n' %len(dihedrals))

	f.write('3 atom types\n')
	f.write('%d bond types\n' %len(bond_coeffs))
	f.write('%d angle types\n' %len(angle_coeffs))
	f.write('%d dihedral types\n\n' %len(dihedral_coeffs))

	f.write('Masses\n\n')
	f.write('\t1\t28\n')
	f.write('\t2\t16\n')
	f.write('\t3\t12\n\n')


def WriteAngleCoeffs(f, BondedPotentials):
	angle_coeffs = BondedPotentials['angle_coeffs']
	f.write('Angle Coeffs\n\n')
	for angle_coeff in angle_coeffs:
		f.write('\t%d\t3.38\t%.4f\n' %(angle_coeff[1], angle_coeff[0]))


def WriteBondCoeffs(f, BondedPotentials):
	bond_coeffs = BondedPotentials['bond_coeffs']
	f.write('\nBond Coeffs\n\n')
	for bond_coeff in bond_coeffs:
		f.write('\t%d\t20\t%.4f\n' %(bond_coeff[1], bond_coeff[0]))


def WriteDihedralCoeff(f, BondedPotentials): 
	dihedral_coeffs = BondedPotentials['dihedral_coeffs']
	# we use a hybrid opls and compass potential to describe the dihedrals
	f.write('\nDihedral Coeffs\n\n')
	for dihedral_coeff in dihedral_coeffs:
		dihedralID = dihedral_coeff[0]
		coeffs = dihedral_coeff[1]
		if coeffs[0] == 'opls':
			f.write('\t%d\topls\t%.4f\t%.4f\t%.4f\t%.4f\n'\
					%(dihedralID, coeffs[1], coeffs[2],\
					 	coeffs[3], coeffs[4]))
		elif coeffs[0] == 'compass':
			f.write('\t%d\tclass2\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n'\
					%(dihedralID, coeffs[1], coeffs[2],\
						coeffs[3], coeffs[4],\
						coeffs[5], coeffs[6]))

	# for the compass potential, we only use the dihedral component, but 
	# we must still specify the mbt, ebt, at, aat, and bb13; we just set
	# the rest of the coeffs to 0 in these cases 
	f.write('\nMiddleBondTorsion Coeffs\n\n')
	for dihedral_coeff in dihedral_coeffs:
		dihedralID = dihedral_coeff[0]
		coeffs = dihedral_coeff[1]
		if coeffs[0] == 'opls':
			f.write('\t%d\topls\n' %dihedralID)
		elif coeffs[0] == 'compass':
			f.write('\t%d\tclass2\t0\t0\t0\t0\n' %dihedralID)

	f.write('\nEndBondTorsion Coeffs\n\n')
	for dihedral_coeff in dihedral_coeffs:
		dihedralID = dihedral_coeff[0]
		coeffs = dihedral_coeff[1]
		if coeffs[0] == 'opls':
			f.write('\t%d\topls\n' %dihedralID)
		elif coeffs[0] == 'compass':
			f.write('\t%d\tclass2\t0\t0\t0\t0\t0\t0\t0\t0\n' %dihedralID)

	f.write('\nAngleTorsion Coeffs\n\n')
	for dihedral_coeff in dihedral_coeffs:
		dihedralID = dihedral_coeff[0]
		coeffs = dihedral_coeff[1]
		if coeffs[0] == 'opls':
			f.write('\t%d\topls\n' %dihedralID)
		elif coeffs[0] == 'compass':
			f.write('\t%d\tclass2\t0\t0\t0\t0\t0\t0\t0\t0\n' %dihedralID)

	f.write('\nAngleAngleTorsion Coeffs\n\n')
	for dihedral_coeff in dihedral_coeffs:
		dihedralID = dihedral_coeff[0]
		coeffs = dihedral_coeff[1]
		if coeffs[0] == 'opls':
			f.write('\t%d\topls\n' %dihedralID)
		elif coeffs[0] == 'compass':
			f.write('\t%d\tclass2\t0\t0\t0\n' %dihedralID)

	f.write('\nBondBond13 Coeffs\n\n')
	for dihedral_coeff in dihedral_coeffs:
		dihedralID = dihedral_coeff[0]
		coeffs = dihedral_coeff[1]
		if coeffs[0] == 'opls':
			f.write('\t%d\topls\n' %dihedralID)
		elif coeffs[0] == 'compass':
			f.write('\t%d\tclass2\t0\t0\t0\n' %dihedralID)

def WriteAtoms(f, atoms):
	f.write('\nAtoms\n\n')
	charge = 0
	for i in range(len(atoms)):
		f.write('%d\t%d\t%d\t%d\t%.4f\t%.4f\t%.4f\n' %(atoms[i][0], atoms[i][1], 
														atoms[i][2], atoms[i][3], 
														atoms[i][4], atoms[i][5],
														atoms[i][6]))

def WriteBonds(f, bonds):
	f.write('\nBonds\n\n')
	for bond in bonds:
		f.write('%d\t%d\t%d\t%d\n' %(bond[0], bond[1], bond[2], bond[3]))

def WriteAngles(f, angles):
	f.write('\nAngles\n\n')
	for angle in angles:
		f.write('%d\t%d\t%d\t%d\t%d\n' %(angle[0], angle[1], angle[2], 
											angle[3], angle[4]))

def WriteDihedrals(f, dihedrals):
	f.write('\nDihedrals\n\n')
	for dihedral in dihedrals:
		f.write('%d\t%d\t%d\t%d\t%d\t%d\n' %(dihedral[0], dihedral[1], dihedral[2], 
											dihedral[3], dihedral[4], dihedral[5]))

#-------------------------------------------------------------------------------------------------#
# write out the data file

def WriteDataFile(sim_dir, Nprecursors, moltype, trial, BondedPotentials, atoms, bonds,\
					angles, dihedrals, cell_length):
	filename = '%sdata.%s_%d_%s' %(sim_dir, moltype, Nprecursors, trial)
	f = open(filename, 'w')

	f.write('# data file for %d %s precursors with porosity capabilities \n' %(Nprecursors, moltype))

	WriteSimulationParams(f, BondedPotentials, atoms, bonds, angles, dihedrals, cell_length)
	WriteAngleCoeffs(f, BondedPotentials)
	WriteBondCoeffs(f, BondedPotentials)
	WriteDihedralCoeff(f, BondedPotentials)
	WriteAtoms(f, atoms)
	WriteBonds(f, bonds)
	WriteAngles(f, angles)	
	WriteDihedrals(f, dihedrals)
	f.close()

	datafile_base = 'data.%s_%d_' %(moltype, Nprecursors)
	return datafile_base

#-------------------------------------------------------------------------------------------------#
# write out the mixed precursor data file
def WriteMixedDataFile(sim_dir, trial, BondedPotentials, Nprecursors1, moltype1, Nprecursors2,\
						moltype2, atoms, bonds, angles, dihedrals, cell_length):
	filename = '%sdata.%s_%d_%s_%d_%s' %(sim_dir, moltype1, Nprecursors1, moltype2, Nprecursors2, trial)
	f = open(filename, 'w')

	f.write('# data file for %d %s and %d %s precursors with porosity capabilities \n' %(Nprecursors1, moltype1,\
																						Nprecursors2, moltype2))

	WriteSimulationParams(f, BondedPotentials, atoms, bonds, angles, dihedrals, cell_length)
	WriteAngleCoeffs(f, BondedPotentials)
	WriteBondCoeffs(f, BondedPotentials)
	WriteDihedralCoeff(f, BondedPotentials)
	WriteAtoms(f, atoms)
	WriteBonds(f, bonds)
	WriteAngles(f, angles)	
	WriteDihedrals(f, dihedrals)
	f.close()

	datafile_base = 'data.%s_%d_%s_%d_' %(moltype1, Nprecursors1, moltype2, Nprecursors2)
	return datafile_base

#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#



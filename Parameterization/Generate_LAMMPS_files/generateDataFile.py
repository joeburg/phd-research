import numpy
import operator
import sys
import time

from StructureAnalysisFuctions import ComputeAtoms, ComputeBonds, ComputeAngles, ComputeDihedrals

#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
# utility functions

def Distance(dx,dy,dz):
	return (dx*dx + dy*dy + dz*dz)**0.5

def BondLength(atom1, atom2):
	dx = abs(atom1[4] - atom2[4])
	dy = abs(atom1[5] - atom2[5])
	dz = abs(atom1[6] - atom2[6])
	return Distance(dx,dy,dz)

def getBondType(bond_length, bond_types):
	for bond in bond_types:
		# use same bond coeff if within 0.01 (Angstrom) tolerance
		if abs(bond_length - bond) <= 0.01:
			return bond_types[bond]
	# if no match found, return none
	return None

def getBondCoeffs(bond_length, bond_coeffs):
	# bond_coeffs is a list of tuples [(1.31, 1), (1.45, 2), ...]
	# (bond_length, bond_type_idx)
	bond_type = None
	bond_type_found = False

	for item in bond_coeffs:
		bond, bond_type_idx = item
		# use same bond coeff if within 0.01 (Angstrom) tolerance
		if abs(bond_length - bond) <= 0.01:
			bond_type = bond_type_idx
			bond_type_found = True
			break

	if not bond_type_found:
		# bond types go from 1 to N 
		bond_type = len(bond_coeffs)+1
		bond_coeffs.append((bond_length, bond_type))

	return bond_type, bond_coeffs

#-------------------------------------------------------------------------------------------------#

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

def LoadSavedPrecursor(moltype):
	filename = moltype + '.xyz'
	data, atom_metrics = LoadData(filename)
	return data, atom_metrics

#-------------------------------------------------------------------------------------------------#
# analyze the precursor structure 

def getPrecursorMetrics(atoms, bonds, angles, dihedrals, \
						bond_coeffs, angle_coeffs, dihedral_coeffs):
	precursor_metrics = {'Natoms': len(atoms),
							'Nbonds': len(bonds),
							'Nangles': len(angles),
							'Ndihedrals': len(dihedrals),
							'Nbond_coeffs': len(bond_coeffs),
							'Nangle_coeffs': len(angle_coeffs),
							'Ndihedral_coeffs': len(dihedral_coeffs)}
	return precursor_metrics


def AnalyzePrecursor(data):
	atoms = ComputeAtoms(data)
	bonds, AM, bond_coeffs = ComputeBonds(data)
	angles, angle_coeffs = ComputeAngles(data, AM)
	dihedrals, dihedral_coeffs = ComputeDihedrals(data, AM)

	# sort bond, angle, dihedral coeffs by increasing ID number (returns list of tuples)
	bond_coeffs = sorted(bond_coeffs.items(), key=operator.itemgetter(1))
	angle_coeffs = sorted(angle_coeffs.items(), key=operator.itemgetter(1))
	dihedral_coeffs = sorted(dihedral_coeffs.items(), key=operator.itemgetter(0))

	# get the precursor metrics 
	precursor_metrics = getPrecursorMetrics(atoms, bonds, angles, dihedrals, \
											bond_coeffs, angle_coeffs, dihedral_coeffs)

	structure = {'atoms': atoms,
					'bonds': bonds,
					'angles': angles,
					'dihedrals': dihedrals,
					'bond_coeffs': bond_coeffs,
					'angle_coeffs': angle_coeffs,
					'dihedral_coeffs': dihedral_coeffs,
					'metrics': precursor_metrics}
	return structure

#-------------------------------------------------------------------------------------------------#
# load the porogen structure 

def LoadPorogen(bond_coeffs):
	# load the atom template 
	filename = 'porogen40_atoms.txt'
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
	bond_type, bond_coeffs = getBondCoeffs(d, bond_coeffs)

	# load the bond template 
	filename = 'porogen40_bonds.txt'
	bondtemplate = []
	f = open(filename)
	f.readline()

	while True:
		fields = f.readline().strip().split()
		if fields:
			ID = int(fields[0])
			# bondtype = int(fields[1])
			atom1 = int(fields[2])
			atom2 = int(fields[3])

			bondtemplate.append([ID, bond_type, atom1, atom2])
		else:
			break
	f.close()
	bondtemplate = numpy.array(bondtemplate)
	return (atomtemplate, bondtemplate), bond_coeffs


#-------------------------------------------------------------------------------------------------#
# replicate precursor structure and porogen molecules 

def ReplicateStructure(atom_metrics, precursor_structure, porogen_templates, Nprecursors):
	# get the precursor templates and metrics 
	Natoms_precursor = precursor_structure['metrics']['Natoms']
	Nbonds_precursor = precursor_structure['metrics']['Nbonds']
	Nangles_precursor = precursor_structure['metrics']['Nangles']
	Ndihrdrals_precursor = precursor_structure['metrics']['Ndihedrals']
	precursor_atom_template = precursor_structure['atoms']
	precursor_bond_template = precursor_structure['bonds']
	precursor_angle_template = precursor_structure['angles']
	precursor_dihedral_template = precursor_structure['dihedrals']

	# take the number of free oxygen to be 3 times the number of Si 
	NSi_per_precursor = atom_metrics['NSi_per_precursor']
	NO = 3*NSi_per_precursor*Nprecursors

	# generate 1 porogen molecule for every 2 precursors
	Nporogen = int(0.5*Nprecursors)
	porogen_atom_template = porogen_templates[0]
	porogen_bond_template = porogen_templates[1]
	Natoms_porogen = len(porogen_atom_template)
	Nbonds_porogen = len(porogen_bond_template)

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

	# generate random positions [0, 1] for the atom positions 
	# and then scale by the simulation cell length
	precursor_rand = numpy.random.rand(Nprecursors, 3)*L
	Oatoms_rand = numpy.random.rand(NO, 3)*L
	porogen_rand = numpy.random.rand(Nporogen, 3)*L

	# replicate the precursor
	for i in range(Nprecursors):
		# replicate the precursor atoms 
		# [atomID, moleculeID, atom_type, charge, x, y, z]
		atom_trans = numpy.array([Natoms_precursor*i,\
									i, 0, 0,\
									precursor_rand[i][0],\
									precursor_rand[i][1],\
									precursor_rand[i][2]])
		# constant transformation for the whole precursor
		addatoms = numpy.tile(atom_trans, (Natoms_precursor,1))

		atoms[Natoms_precursor*i:\
				Natoms_precursor*(i+1), :] = precursor_atom_template + addatoms

		# replicate the precursor bonds
		# [bondID, bond_type, atom1, atom2]
		bond_trans = numpy.array([Nbonds_precursor*i,\
									0,\
									Natoms_precursor*i,\
									Natoms_precursor*i])

		addbonds = numpy.tile(bond_trans, (Nbonds_precursor,1))

		bonds[Nbonds_precursor*i:\
				Nbonds_precursor*(i+1), :] = precursor_bond_template + addbonds

		# replicate the precursor angles
		# [angleID, angle_type, atom1, atom2, atom3] 
		angle_trans = numpy.array([Nangles_precursor*i,\
									0,\
									Natoms_precursor*i,\
									Natoms_precursor*i,\
									Natoms_precursor*i])

		addangles = numpy.tile(angle_trans, (Nangles_precursor,1))

		angles[Nangles_precursor*i:\
				Nangles_precursor*(i+1), :] = precursor_angle_template + addangles

		# replicate the precursor dihedrals 
		# [dihedralID, dihedral_type, atom1, atom2, atom3, atom4]
		dihedral_trans = numpy.array([Ndihrdrals_precursor*i,\
										0,\
										Natoms_precursor*i,\
										Natoms_precursor*i,\
										Natoms_precursor*i,\
										Natoms_precursor*i])

		adddihedrals = numpy.tile(dihedral_trans, (Ndihrdrals_precursor,1))

		dihedrals[Ndihrdrals_precursor*i:\
					Ndihrdrals_precursor*(i+1), :] = precursor_dihedral_template + adddihedrals


	# replicate the free O atoms 
	atom_idx = Nprecursors*Natoms_precursor
	for i in range(NO):
		atoms[atom_idx+i,:] = numpy.array([atom_idx+i,\
											0, 2, -2,\
											Oatoms_rand[i][0],\
											Oatoms_rand[i][1],\
											Oatoms_rand[i][2]])

	# replicate the porogen molecules
	atom_idx = Nprecursors*Natoms_precursor + NO
	bond_idx = Nprecursors*Nbonds_precursor
	for i in range(Nporogen):
		# replicate the porogen atoms 
		atom_trans = numpy.array([atom_idx+Natoms_porogen*i,\
									Nprecursors+i,\
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

	# ensure that the final simulation cell contains all the atoms 
	x_min = min(atoms[:,4])
	x_max = max(atoms[:,4])
	y_min = min(atoms[:,5])
	y_max = max(atoms[:,5]) 
	z_min = min(atoms[:,6])
	z_max = max(atoms[:,6])

	cell_min = min(x_min, y_min, z_min) - 5
	cell_max = max(x_max, y_max, z_max) + 5

	return atoms, bonds, angles, dihedrals, (cell_min, cell_max)

#-------------------------------------------------------------------------------------------------#
# write out the data file

def WriteDataFile(Nprecursors, moltype, precursor_structure, atoms, bonds, angles, dihedrals, cell_length):
	bond_coeffs = precursor_structure['bond_coeffs']
	angle_coeffs = precursor_structure['angle_coeffs']
	dihedral_coeffs = precursor_structure['dihedral_coeffs']

	filename = 'data.%s_%d' %(moltype, Nprecursors)
	f = open(filename, 'w')

	f.write('# data file for %d %s precursors with porosity capabilities \n' %(Nprecursors, moltype))
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

	f.write('Angle Coeffs\n\n')
	for angle_coeff in angle_coeffs:
		f.write('\t%d\t3.38\t%.4f\n' %(angle_coeff[1], angle_coeff[0]))

	f.write('\nBond Coeffs\n\n')
	for bond_coeff in bond_coeffs:
		f.write('\t%d\t20\t%.4f\n' %(bond_coeff[1], bond_coeff[0]))

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
			f.write('\t%d\tskip\n' %dihedralID)
		elif coeffs[0] == 'compass':
			f.write('\t%d\tclass2\t0\t0\t0\t0\n' %dihedralID)

	f.write('\nEndBondTorsion Coeffs\n\n')
	for dihedral_coeff in dihedral_coeffs:
		dihedralID = dihedral_coeff[0]
		coeffs = dihedral_coeff[1]
		if coeffs[0] == 'opls':
			f.write('\t%d\tskip\n' %dihedralID)
		elif coeffs[0] == 'compass':
			f.write('\t%d\tclass2\t0\t0\t0\t0\t0\t0\t0\t0\n' %dihedralID)

	f.write('\nAngleTorsion Coeffs\n\n')
	for dihedral_coeff in dihedral_coeffs:
		dihedralID = dihedral_coeff[0]
		coeffs = dihedral_coeff[1]
		if coeffs[0] == 'opls':
			f.write('\t%d\tskip\n' %dihedralID)
		elif coeffs[0] == 'compass':
			f.write('\t%d\tclass2\t0\t0\t0\t0\t0\t0\t0\t0\n' %dihedralID)

	f.write('\nAngleAngleTorsion Coeffs\n\n')
	for dihedral_coeff in dihedral_coeffs:
		dihedralID = dihedral_coeff[0]
		coeffs = dihedral_coeff[1]
		if coeffs[0] == 'opls':
			f.write('\t%d\tskip\n' %dihedralID)
		elif coeffs[0] == 'compass':
			f.write('\t%d\tclass2\t0\t0\t0\n' %dihedralID)

	f.write('\nBondBond13 Coeffs\n\n')
	for dihedral_coeff in dihedral_coeffs:
		dihedralID = dihedral_coeff[0]
		coeffs = dihedral_coeff[1]
		if coeffs[0] == 'opls':
			f.write('\t%d\tskip\n' %dihedralID)
		elif coeffs[0] == 'compass':
			f.write('\t%d\tclass2\t0\t0\t0\n' %dihedralID)

	f.write('\nAtoms\n\n')
	charge = 0
	for i in range(len(atoms)):
		f.write('%d\t%d\t%d\t%d\t%.4f\t%.4f\t%.4f\n' %(atoms[i][0], atoms[i][1], 
														atoms[i][2], atoms[i][3], 
														atoms[i][4], atoms[i][5],
														atoms[i][6]))

	f.write('\nBonds\n\n')
	for bond in bonds:
		f.write('%d\t%d\t%d\t%d\n' %(bond[0], bond[1], bond[2], bond[3]))

	f.write('\nAngles\n\n')
	for angle in angles:
		f.write('%d\t%d\t%d\t%d\t%d\n' %(angle[0], angle[1], angle[2], 
											angle[3], angle[4]))

	f.write('\nDihedrals\n\n')
	for dihedral in dihedrals:
		f.write('%d\t%d\t%d\t%d\t%d\t%d\n' %(dihedral[0], dihedral[1], dihedral[2], 
											dihedral[3], dihedral[4], dihedral[5]))
	f.close()


#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
if len(sys.argv) < 3:
	print 'Usage:'
	print '  python %s <N precursors> <new precursor [type "new"] or saved precursor [type "saved"]>' %sys.argv[0]
	exit()

Nprecursors = int(sys.argv[1])
precursor_type = sys.argv[2]

filename = ''
if precursor_type == 'new':
	filename = raw_input("Give the new precursor's filename: ")
	moltype = filename[:-4]
else:
	moltype = raw_input('Type the name of the saved precursor:'+\
						'\n\t\t\t"EtOCS"\n\t\t\t"EtOCSMethyl"\n\t\t\t'+\
						'"EtOCSVinyl"\n\t\t\t"EtOCSPhenyl"\n\t\t\t'+\
						'"135Benz"\n\t\t\t"SiO2"\n')
t0 = time.time()
# load the new file or saved files
if filename:
	data, atom_metrics = LoadData(filename)
else:
	data, atom_metrics = LoadSavedPrecursor(moltype)

# analyze the precursor structure
precursor_structure = AnalyzePrecursor(data)
# load the porogen template (Note: the bond coeffs may change)
porogen_templates, bond_coeffs = LoadPorogen(precursor_structure['bond_coeffs'])
precursor_structure['bond_coeffs'] = bond_coeffs

# replicate the precursor
atoms, bonds, angles, dihedrals, cell_length = ReplicateStructure(atom_metrics,\
														precursor_structure,\
														porogen_templates,\
														Nprecursors) 

# write out the LAMMPS data file
WriteDataFile(Nprecursors, moltype, precursor_structure, atoms, bonds, angles, dihedrals, cell_length)

print 'Analyzed the precursor and generated the data file in %.4f seconds.' %(time.time()-t0)


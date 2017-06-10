import math
import numpy
import operator

#-------------------------------------------------------------------------------------------------#
def createDataFile(SimulationData, Structure, Potentials):
	# generate the structure data
	StructureData = ReplicateStructures(Structure, SimulationData)

	# create the data file in lammps format 
	WriteDataFile(SimulationData, Potentials, StructureData)


#-------------------------------------------------------------------------------------------------#
# replicate precursor structure and porogen molecules 

def ReplicateStructures(Structure, SimulationData):
	# get the precursor templates and metrics 
	Natoms_precursor = len(Structure['atomstemplate'])
	Nbonds_precursor = len(Structure['bondstemplate'])
	Nangles_precursor = len(Structure['anglestemplate'])
	Ndihrdrals_precursor = len(Structure['dihedralstemplate'])
	Nprecursors = SimulationData['Nprecursors']
	NfreeAtoms = SimulationData['NfreeAtoms']

	# get the total number of atoms, bonds, angles, dihedrals
	Natoms = Nprecursors*Natoms_precursor + NfreeAtoms
	Nbonds = Nprecursors*Nbonds_precursor
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
	angles, dihedrals = ReplicatePrecursor(Nprecursors, Structure,\
											atoms, bonds, angles, dihedrals,\
											0, 0, 0, 0, 0, L)

	# replicate the free atoms 
	atom_idx = Nprecursors*Natoms_precursor
	free_atom_data = Structure['free_atom_data']
	atoms = ReplicateFreeAtoms(NfreeAtoms, free_atom_data, atom_idx, atoms, L)

	# get the dimensions of the simulation cell
	cell_size = CreateSimulationBox(atoms)

	StructureData = {
		'atoms' : atoms,
		'bonds' : bonds,
		'angles' : angles,
		'dihedrals' : dihedrals,
		'cell_size' : cell_size,
	}

	print '\nCreated simulation data:'
	print 'Number of %s precursors = %d' %(SimulationData['mol_name'], Nprecursors)
	print 'Number of free %s atoms = %d' %(SimulationData['atom_name'], NfreeAtoms)

	return StructureData


#-------------------------------------------------------------------------------------------------#

def ReplicatePrecursor(Nprecursors, Structure, atoms, bonds, angles, dihedrals,\
						atom_idx, mol_idx, bond_idx, angle_idx, dihedral_idx, L):
	# get the precursor templates and metrics 
	Natoms_precursor = len(Structure['atomstemplate'])
	Nbonds_precursor = len(Structure['bondstemplate'])
	Nangles_precursor = len(Structure['anglestemplate'])
	Ndihrdrals_precursor = len(Structure['dihedralstemplate'])
	precursor_atom_template = Structure['atomstemplate']
	precursor_bond_template = Structure['bondstemplate']
	precursor_angle_template = Structure['anglestemplate']
	precursor_dihedral_template = Structure['dihedralstemplate']

	# generate random positions [0, 1] for the atom positions 
	# and then scale by the simulation cell length
	precursor_rand = numpy.random.rand(Nprecursors, 3)*L

	# replicate the precursor
	for i in range(Nprecursors):
		# replicate the precursor atoms 
		# [atomID, moleculeID, atom_type, charge, x, y, z]
		if len(precursor_atom_template):
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
		if len(precursor_bond_template):
			bond_trans = numpy.array([bond_idx+Nbonds_precursor*i,\
										0,\
										atom_idx+Natoms_precursor*i,\
										atom_idx+Natoms_precursor*i])

			addbonds = numpy.tile(bond_trans, (Nbonds_precursor,1))

			bonds[bond_idx+Nbonds_precursor*i:\
					bond_idx+Nbonds_precursor*(i+1), :] = precursor_bond_template + addbonds

		# replicate the precursor angles
		# [angleID, angle_type, atom1, atom2, atom3] 
		if len(precursor_angle_template):
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
		if len(precursor_dihedral_template):
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


def ReplicateFreeAtoms(NfreeAtoms, free_atom_data, atom_idx, atoms, L):
	free_atom_ID, charge = free_atom_data
	# generate random positions [0, 1] for the atom positions 
	# and then scale by the simulation cell length
	atoms_rand = numpy.random.rand(NfreeAtoms, 3)*L

	# replicate the free atoms
	# [atomID, moleculeID, atom_type, charge, x, y, z] 
	for i in range(NfreeAtoms):
		atoms[atom_idx+i,:] = numpy.array([atom_idx+i+1, 0,\
											free_atom_ID,\
											charge,\
											atoms_rand[i][0],\
											atoms_rand[i][1],\
											atoms_rand[i][2]])
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
#-------------------------------------------------------------------------------------------------#
# functions to write the data files

def WriteSimulationParams(f, Potentials, atoms, bonds, angles, dihedrals, cell_size):
	Nmasses = len(Potentials['masses'])
	Nbond_coeffs = len(Potentials['bond_coeffs'])
	Nangle_coeffs = len(Potentials['angle_coeffs'])
	Ndihedral_coeffs = len(Potentials['dihedral_coeffs'])
	cell_min, cell_max = cell_size

	f.write('%d %d xlo xhi\n' %(cell_min, cell_max))
	f.write('%d %d ylo yhi\n' %(cell_min, cell_max))
	f.write('%d %d zlo zhi\n\n' %(cell_min, cell_max))

	f.write('%d atoms\n' %len(atoms))
	f.write('%d bonds\n' %len(bonds))
	f.write('%d angles\n' %len(angles))
	f.write('%d dihedrals\n\n' %len(dihedrals))

	f.write('%d atom types\n' %Nmasses)
	f.write('%d bond types\n' %Nbond_coeffs)
	f.write('%d angle types\n' %Nangle_coeffs)
	f.write('%d dihedral types\n\n' %Ndihedral_coeffs)


def WriteMasses(f, masses):
	f.write('Masses\n\n')
	for ID, mass in masses:
		f.write('\t%d\t%.6f\n' %(ID, mass))


def WriteAngleCoeffs(f, angle_coeffs):
	f.write('\nAngle Coeffs\n\n')
	for ID, K, angle in angle_coeffs:
		f.write('\t%d\t%.8f\t%.8f\n' %(ID, K, angle))


def WriteBondCoeffs(f, bond_coeffs):
	f.write('\nBond Coeffs\n\n')
	for ID, K, bond in bond_coeffs:
		f.write('\t%d\t%.8f\t%.8f\n' %(ID, K, bond))


def WriteDihedralCoeff(f, dihedral_coeffs): 
	f.write('\nDihedral Coeffs\n\n')
	for ID, K1, phi1, K2, phi2, K3, phi3 in dihedral_coeffs:
		f.write('\t%d\t%.8f\t%.2f\t%.8f\t%.2f\t%.8f\t%.2f\n'\
				%(ID, K1, phi1, K2, phi2, K3, phi3))

	# for the compass potential, we only use the dihedral component, but 
	# we must still specify the mbt, ebt, at, aat, and bb13; we just set
	# the rest of the coeffs to 0 in these cases 
	f.write('\nMiddleBondTorsion Coeffs\n\n')
	for dihedral in dihedral_coeffs:
		ID = dihedral[0]
		f.write('\t%d\t0\t0\t0\t0\n' %ID)


	f.write('\nEndBondTorsion Coeffs\n\n')
	for dihedral in dihedral_coeffs:
		ID = dihedral[0]
		f.write('\t%d\t0\t0\t0\t0\t0\t0\t0\t0\n' %ID)


	f.write('\nAngleTorsion Coeffs\n\n')
	for dihedral in dihedral_coeffs:
		ID = dihedral[0]
		f.write('\t%d\t0\t0\t0\t0\t0\t0\t0\t0\n' %ID)


	f.write('\nAngleAngleTorsion Coeffs\n\n')
	for dihedral in dihedral_coeffs:
		ID = dihedral[0]
		f.write('\t%d\t0\t0\t0\n' %ID)	


	f.write('\nBondBond13 Coeffs\n\n')
	for dihedral in dihedral_coeffs:
		ID = dihedral[0]
		f.write('\t%d\t0\t0\t0\n' %ID)


def WritePairPotential(f, pair_potential):
	f.write('\nPair Coeffs\n\n')
	for ID, epsilon, sigma, Rc_vwd, Rc_coul in pair_potential:
		f.write('\t%d\t%.8f\t%.8f\t%.2f\t%.2f\n'\
					%(ID, epsilon, sigma, Rc_vwd, Rc_coul))


def WriteAtoms(f, atoms):
	f.write('\nAtoms\n\n')
	for atomID, moleculeID, atom_type, charge, x, y, z in atoms:
		f.write('%d\t%d\t%d\t%.4f\t%.4f\t%.4f\t%.4f\n'\
				%(atomID, moleculeID, atom_type,\
					charge, x, y, z))


def WriteBonds(f, bonds):
	f.write('\nBonds\n\n')
	for bond in bonds:
		f.write('%d\t%d\t%d\t%d\n' %(bond[0], bond[1],\
										bond[2], bond[3]))

def WriteAngles(f, angles):
	f.write('\nAngles\n\n')
	for angle in angles:
		f.write('%d\t%d\t%d\t%d\t%d\n'\
					%(angle[0], angle[1], angle[2],\
						angle[3], angle[4]))

def WriteDihedrals(f, dihedrals):
	f.write('\nDihedrals\n\n')
	for dihedral in dihedrals:
		f.write('%d\t%d\t%d\t%d\t%d\t%d\n'\
					%(dihedral[0], dihedral[1], dihedral[2],\
						dihedral[3], dihedral[4], dihedral[5]))

#-------------------------------------------------------------------------------------------------#
# write out the data file

def WriteDataFile(SimulationData, Potentials, StructureData):
	# simulation data
	sim_dir = SimulationData['sim_dir']
	mol_name = SimulationData['mol_name']
	atom_name = SimulationData['atom_name']
	Nprecursors = SimulationData['Nprecursors']
	NfreeAtoms = SimulationData['NfreeAtoms']

	# structure data
	atoms = StructureData['atoms']
	bonds = StructureData['bonds']
	angles = StructureData['angles']
	dihedrals = StructureData['dihedrals']
	cell_size = StructureData['cell_size']

	# potential data 
	masses = Potentials['masses']
	angle_coeffs = Potentials['angle_coeffs']
	bond_coeffs = Potentials['bond_coeffs']
	dihedral_coeffs = Potentials['dihedral_coeffs']
	pair_potential = Potentials['pair_potential']

	# create the filename 
	filename = '%sdata.%s_%d_%s_%d' %(sim_dir, mol_name, Nprecursors, atom_name, NfreeAtoms)
	f = open(filename, 'w')

	f.write('# data file for %d %s precursors and %d %s atoms.\n'\
			%(Nprecursors, mol_name, NfreeAtoms, atom_name))

	WriteSimulationParams(f, Potentials, atoms, bonds, angles, dihedrals, cell_size)
	WriteMasses(f, masses)
	WriteAngleCoeffs(f, angle_coeffs)
	WriteBondCoeffs(f, bond_coeffs)
	WriteDihedralCoeff(f, dihedral_coeffs)
	# need to write pair potential in the input script 
	# WritePairPotential(f, pair_potential)
	WriteAtoms(f, atoms)
	WriteBonds(f, bonds)
	WriteAngles(f, angles)	
	WriteDihedrals(f, dihedrals)
	f.close()


#-------------------------------------------------------------------------------------------------#


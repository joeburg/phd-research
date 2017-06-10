import copy
import math
import numpy

#-------------------------------------------------------------------------------------------------#
# utility functions 

def ShiftToPositiveCoords(data):
	''' ensures all coordinate values are positive '''
	coord_shfit = abs(numpy.min(data))
	data[:,1:] = data[:,1:] + coord_shfit
	return data


def Distance(dx,dy,dz):
	return (dx*dx + dy*dy + dz*dz)**0.5


def Angle(atom1,vertex,atom2):
	# use law of cosines
	D1V = Distance(atom1[1]-vertex[1],atom1[2]-vertex[2],atom1[3]-vertex[3])
	D2V = Distance(atom2[1]-vertex[1],atom2[2]-vertex[2],atom2[3]-vertex[3])
	D12 = Distance(atom1[1]-atom2[1],atom1[2]-atom2[2],atom1[3]-atom2[3])

	# return angle in degrees	
	return math.acos((D1V*D1V + D2V*D2V - D12*D12)/(2*D1V*D2V))*(180.0/math.pi)

def AtomType(atomID):
	if atomID == 1:
		return 'Si'
	elif atomID == 2:
		return 'O'
	elif atomID == 3:
		return 'C'

#-------------------------------------------------------------------------------------------------#
# bond analysis functions

def getBondType(bond_length, bond_types):
	for bond in bond_types:
		# use same bond coeff if within 0.01 (Angstrom) tolerance
		if abs(bond_length - bond) <= 0.01:
			return bond_types[bond]
	# if no match found, return none
	return None


def ComputeBonds(data):
	
	# define the adjacency matrix to help compute angles and dihedrals
	AM = numpy.zeros((len(data),len(data)))

	# bond template will be in the lammps format [ID, type, atomID1, atomID2]
	bondstemplate = []

	# lammps uses index 1 to N
	ID = 1

	SiC_bond = 1 # ~1.86 A
	CC_single_bond = 2 # ~1.53 A
	CC_double_bond = 3 # ~1.31 A

	bond_types = {}
	Nbond_types = 0

	# precursor molecules are small so we will compute all distances
	for i in range(len(data)):
		atom1 = data[i][0]
		for j in range(i+1,len(data)):
			atom2 = data[j][0]			
	
			dx = abs(data[j][1] - data[i][1])
			dy = abs(data[j][2] - data[i][2])
			dz = abs(data[j][3] - data[i][3])

			d = Distance(dx,dy,dz)
			
			# for initial paramaterization, we set all bond stiffness coeffs 
			# to the same value; so just search for unique bond lengths
			if d < 2.0:
				bond_type = getBondType(d, bond_types)
				if bond_type:
					bondstemplate.append([ID, bond_type, i+1, j+1])		

				else:
					Nbond_types += 1
					bond_types[d] = Nbond_types
					bondstemplate.append([ID, bond_types[d], i+1, j+1])

				# update the adjacency matrix
				AM[i][j] = d
				AM[j][i] = d	
				ID += 1			

			# # Si-C bond (bond type 1)
			# if (atom1==1 and atom2==3) or (atom1==3 and atom2==1):
			# 	if d < 2.0:
			# 		# lammps uses index 1 to N
			# 		bondstemplate.append([ID, SiC_bond, i+1, j+1])
			# 		ID += 1

			# 		# update the adjacency matrix
			# 		AM[i][j] = d
			# 		AM[j][i] = d

			# # C-C bond (bond type 2 and 3)
			# if atom1==3 and atom2==3:
			# 	if d < 1.7:
			# 		if d < 1.4:
			# 			bondstemplate.append([ID, CC_double_bond, i+1, j+1])
			# 		else:
			# 			bondstemplate.append([ID, CC_single_bond, i+1, j+1])

			# 		ID += 1

			# 		# update the adjacency matrix
			# 		AM[i][j] = d
			# 		AM[j][i] = d

	return bondstemplate, AM, bond_types

#-------------------------------------------------------------------------------------------------#
# angle analysis functions

def getAngleType(angle_deg, angle_types):
	for angle in angle_types:
		# use same bond coeff if within 0.1 (degree) tolerance
		if abs(angle_deg - angle) <= 0.1:
			return angle_types[angle]
	# if no match found, return none
	return None


def ComputeAngles(data, AM):

	# angle template will be in the lammps format [ID, type, atomID1, atomID2 (vertex), atomID3]
	anglestemplate = []
	ID = 1

	angle_tet = 1 # tetradedral (109.5)
	angle_trig = 2 # trigonal (120.0)

	angle_types = {}
	Nangle_types = 0

	# use the AM to find which verticies are relevant
	for i in range(len(AM)):

		# indicies of atoms bonded to reference atom
		bonds = numpy.nonzero(AM[i])[0]
		
		# if there is more than 1 atom bonded to the reference atom, then it's a vertex
		if len(bonds) > 1:

			vertex = data[i]

			# compute the all angles for the given vertex
			# given the bonds to a vertex, iterate over pairs of other 2 atoms with no repeats
			for j in range(len(bonds)):
				atom1 = data[bonds[j]]
				for k in range(j+1,len(bonds)):
					atom2 = data[bonds[k]]

					angle = Angle(atom1,vertex,atom2)

					angle_type = getAngleType(angle, angle_types)
					if angle_type:
						anglestemplate.append([ID, angle_type, bonds[j]+1, i+1, bonds[k]+1])	

					else:
						Nangle_types += 1
						angle_types[angle] = Nangle_types
						anglestemplate.append([ID, angle_types[angle], bonds[j]+1, i+1, bonds[k]+1])

					ID += 1
					
					# if angle >= 106 and angle <= 112:
					# 	# lammps uses index 1 to N
					# 	anglestemplate.append([ID, angle_tet, bonds[j]+1, i+1, bonds[k]+1])
					# 	ID += 1

					# elif angle >= 113 and angle <=125:
					# 	anglestemplate.append([ID, angle_trig, bonds[j]+1, i+1, bonds[k]+1])
					# 	ID += 1

	return anglestemplate, angle_types

#-------------------------------------------------------------------------------------------------#
# dihedral analysis functions

def CheckPreviousStep(step_previous, next_steps):
	previous_idx = numpy.where(next_steps==step_previous)[0]

	if len(previous_idx):
		next_steps = numpy.delete(next_steps, previous_idx[0])
		
	return next_steps


def CleanPaths(paths):
	# remove any paths that are reverse paths
	paths_clean = copy.deepcopy(paths)

	for i in range(len(paths)):
		path1 = paths[i] 

		for j in range(i+1,len(paths)):
			path2 = paths[j]

			if 	path1 == list(reversed(path2)):
				paths_clean.remove(path2)

	return paths_clean

def ComputePaths(AM):
	# idea is to take 3 steps (joining 4 atoms) and record the path
	# must move to different atom with each step
	# previous steps are not allowed (e.g. if you move 1->2, then you cannot move 2->1)
	# finally check that paths are not reverse of other paths
	# this algo is horribly written; probably can make it much more elegant with recursion 
	paths = []

	# start with each atom and look for paths of length 3 following the rules above
	for i in range(len(AM)):
		step0 = i
		
		# start with each atom as a starting position (frist step)
		next_steps1 = numpy.nonzero(AM[i])[0]

		for j in range(len(next_steps1)):
			step1 = next_steps1[j]

			# second step
			next_steps2 = numpy.nonzero(AM[step1])[0]
			next_steps2 = CheckPreviousStep(step0, next_steps2)

			if len(next_steps2):
				for k in range(len(next_steps2)):
					step2 = next_steps2[k]

					# third step
					next_steps3 = numpy.nonzero(AM[step2])[0]
					next_steps3 = CheckPreviousStep(step1, next_steps3)

					if len(next_steps3):
						for l in range(len(next_steps3)):
							step3 = next_steps3[l]

							paths.append([step0, step1, step2, step3]) 

	# remove reverse paths
	paths = CleanPaths(paths)
	return paths

def DihedralAngle(p0,p1,p2,p3):
    """Praxeolitic formula: 1 sqrt, 1 cross product"""

    b0 = -1.0*(p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2

    # normalize b1 so that it does not influence magnitude of vector
    # rejections that come next
    b1 /= numpy.linalg.norm(b1)

    # vector rejections
    # v = projection of b0 onto plane perpendicular to b1
    #   = b0 minus component that aligns with b1
    # w = projection of b2 onto plane perpendicular to b1
    #   = b2 minus component that aligns with b1
    v = b0 - numpy.dot(b0, b1)*b1
    w = b2 - numpy.dot(b2, b1)*b1

    # angle between v and w in a plane is the torsion angle
    # v and w may not be normalized but that's fine since tan is y/x
    x = numpy.dot(v, w)
    y = numpy.dot(numpy.cross(b1, v), w)
    return numpy.degrees(numpy.arctan2(y, x))


def DihedralPotentialDict(parameterization):
	if parameterization:
		return {1: (0,0,0,0)}	

	else:
		dihedral_types = {	1 : ('opls', 0.375, -0.2, 0.625, -0.1),
							2 : ('opls', 0.01787, -0.007509, 0.06132, -0.002474),
							3 : ('opls', 0.01689, -0.003413, 0.05959, -0.002919),
							4 : ('opls', 0.006063, -0.006386, 0.06797, -0.003199),
							5 : ('opls', 5.1866, -0.7573, 0.1955, -0.04839),
							6 : ('compass', 0.003522, 180, 0.002722, 180, 0.02927, 0),
							7 : ('compass', 0, 0, 0.01072, 180, 0, 0),
							8 : ('compass', 5.5263, 0, 0.05301, -180, -0.4549, 0)}
		return dihedral_types

def DihedralTypesDict(parameterization):
	if parameterization:
		return 1

	else:
		dihedral_types = {	'SiCCSi': 		1,
							'CCSiCH3': 		2,
							'CCSiC2': 		3,
							'CCSiCar': 		4,
							'SiCarCarCar': 	5,
							'CSiC2C2': 		6,
							'CSiCarCar': 	7,
							'CarCarCarCar': 8}
		return dihedral_types

def GetDihedralType(path, data, AM):
	# this function will need to be extended if more dihedrals are added 
	# as of now it's matches given the assumption that only a finite set of dihedrals exist
	dihedral_type = ''
	for atom in path:
		dihedral_type += AtomType(data[atom][0])

	# distinguish between CCSiCH3, CCSiC2, CCSiCar
	if dihedral_type == 'CCSiC' or dihedral_type == 'CSiCC':
		if dihedral_type == 'CCSiC':
			atom = path[3]
		else:
			atom = path[0]

		Nbonds = numpy.nonzero(AM[atom])[0]

		# if the terminal C only has 1 bond, then it's a methyl group
		if len(Nbonds) == 1:
			dihedral_type = 'CCSiCH3'

		elif len(Nbonds) == 2:
			# find the C-C bond and check if it's a single or double bond
			for atom2 in Nbonds:
				atom_type2 = data[atom2][0]

				# only look for C atoms 
				if atom_type2 == 3:
					# if the C-C bond is < 1.4 then it's C=C 
					if AM[atom][atom2] < 1.4:
						dihedral_type = 'CCSiC2'
					else:
						dihedral_type = 'CCSiCar'

	elif dihedral_type == 'SiCCC' or dihedral_type == 'CCCSi':
		dihedral_type = 'SiCarCarCar'

	# distinguish between CSiC2C2 and CSiCarCar
	elif dihedral_type == 'CSiCC' or dihedral_type == 'CCSiC':
		if dihedral_type == 'CSiCC':
			atom1 = path[2]
			atom2 = path[3]
		else:
			atom1 = path[0]
			atom2 = path[1]

		# check to see if C-C bond is a double or single bond
		if AM[atom1][atom2] < 1.4:
			dihedral_type = 'CSiC2C2'
		else:
			dihedral_type = 'CSiCarCar'

	# only have CCCC in aromatic rings as of now
	elif dihedral_type == 'CCCC':
		dihedral_type = 'CarCarCarCar'

	return DihedralTypesDict(dihedral_type)



def ComputeDihedrals(data, AM, parameterization=False):
	# this program computes the dihedral template
	dihedralstemplate = []

	# get the atom index lists for the dihedrals
	paths = ComputePaths(AM)

	# comptue the dihedral angle
	ID = 0 
	for path in paths:
		ID += 1

		atom0 = numpy.array(data[path[0]][1:])
		atom1 = numpy.array(data[path[1]][1:])
		atom2 = numpy.array(data[path[2]][1:])
		atom3 = numpy.array(data[path[3]][1:])

		dihedral = DihedralAngle(atom0,atom1,atom2,atom3)

		if parameterization:
			# for parameterization, set all dihedrals to 0
			dihedraltype = 1

		else:
			dihedraltype = GetDihedralType(path, data, AM)

		dihedralstemplate.append([ID, dihedraltype, path[0]+1, 
									path[1]+1, path[2]+1, path[3]+1])

	dihedral_potentials = DihedralPotentialDict(parameterization)

	return dihedralstemplate, dihedral_potentials






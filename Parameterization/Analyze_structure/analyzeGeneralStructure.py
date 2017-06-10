import copy
import glob
import math
import numpy
import os
import sys
import time
import yaml

from masses import Masses

#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
''' 
This progam analyzes a precursor structure and creates templates that include 
the atom positions, bonds, angles, and dihedrals
'''


def LoadData(atomfile, filenames):
	""" reads .xyz file """ 
	atomTypesFile, bondTypesFile,\
		angleTypesFile, dihedralTypesFile = filenames

	# load the atom types 
	atomTypes = {}
	f = open(atomTypesFile)
	f.readline()
	while True:
		fields = f.readline().strip().split()
		if fields:
			atomID, atom_type, charge = fields
			atomID = int(atomID)
			charge = float(charge)
			atomTypes[atomID] = (atom_type, charge)
		else:
			break
	f.close()

	# load the bond types 
	bondTypes = {}
	f = open(bondTypesFile)
	f.readline()
	while True:
		fields = f.readline().strip().split()
		if fields:
			bond_type, bond_type_ID, bond_coeff = fields
			bondTypes[bond_type] = (bond_type_ID, bond_coeff)
		else:
			break
	f.close()

	# load the angle types
	angleTypes = {}
	f = open(angleTypesFile)
	f.readline()
	while True:
		fields = f.readline().strip().split()
		if fields:
			angle_type, angle_type_ID, angle_coeff = fields
			angleTypes[angle_type] = (angle_type_ID, angle_coeff)
		else:
			break
	f.close()

	# load the dihedral types
	dihedralTypes = {}
	f = open(dihedralTypesFile)
	f.readline()
	while True:
		fields = f.readline().strip().split()
		if fields:
			dihedral_type, dihedral_type_ID,\
					K1, phi1, K2, phi2, K3, phi3 = fields
			dihedralTypes[dihedral_type] = (dihedral_type_ID, K1,\
											phi1, K2, phi2, K3, phi3)
		else:
			break
	f.close()

	# data structures to store atom coords and lammps atom template 
	# atom template: [atomID, moleculeID, atom_type, charge, x, y, z]
	data = []
	atomstemplate = []

	f = open(atomfile)
	f.readline()

	while True:
		fields = f.readline().strip().split()
		if fields:
			atom_type = int(fields[0])
			xcoord = float(fields[1])
			ycoord = float(fields[2])
			zcoord = float(fields[3])

			# populate the data array
			data.append([atom_type, xcoord, ycoord, zcoord])
		else:
			break
	f.close()

	data = numpy.array(data)
	data = ShiftToPositiveCoords(data)

	# create atom template with shifted data 
	# lammps uses index 1 to N
	ID = 1
	molID = 1 
	for atom in data:
		atom_type, xcoord, ycoord, zcoord = atom
		charge = atomTypes[atom_type][1]

		# populate the atom template 
		atomstemplate.append([ID, molID, atom_type,\
								charge, xcoord, ycoord, zcoord])
		ID += 1

	return data, atomstemplate, atomTypes, bondTypes,\
			angleTypes, dihedralTypes

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

#-------------------------------------------------------------------------------------------------#
def getBondType(bond_length, bond_types):
	for bond in bond_types:
		# use same bond coeff if within 0.01 (Angstrom) tolerance
		if abs(bond_length - bond) <= 0.01:
			return bond_types[bond]
	# if no match found, return None
	return None

def ComputeBonds(data, atomTypes, bondTypes):
	
	# define the adjacency matrix to help compute angles and dihedrals
	AM = numpy.zeros((len(data),len(data)))

	# bond template will be in the lammps format [ID, type, atomID1, atomID2]
	bondstemplate = []

	# lammps uses index 1 to N
	ID = 1

	bond_types = {}
	Nbond_types = 0
	bonds = []

	# precursor molecules are small so we will compute all distances
	for i in range(len(data)):
		atom1 = data[i][0]
		for j in range(i+1,len(data)):
			atom2 = data[j][0]			
	
			dx = abs(data[j][1] - data[i][1])
			dy = abs(data[j][2] - data[i][2])
			dz = abs(data[j][3] - data[i][3])

			d = Distance(dx,dy,dz)
			
			# for initial paramaterization, get the atoms ID 
			# and bond lengths  
			if d < 2.0:
				# get the atom types and atom type ID 
				atom1_type = atomTypes[atom1][0]
				atom2_type = atomTypes[atom2][0]
				bond_type = '%s-%s' %(atom1_type, atom2_type)

				# ensure that the atoms can bond
				if bond_type in bondTypes:		
					# update the bond data 
					bonds.append((atom1_type, atom2_type, d))

					# update the adjacency matrix
					AM[i][j] = d
					AM[j][i] = d	
						
					# get the bond type ID 
					bond_type_ID = bondTypes[bond_type][0]
					bondstemplate.append([ID, bond_type_ID, i+1, j+1])
					ID += 1

	return bonds, AM, bondstemplate

#-------------------------------------------------------------------------------------------------#

def ComputeAngles(data, AM, atomTypes, angleTypes):

	# angle template will be in the lammps format [ID, type, atomID1, atomID2 (vertex), atomID3]
	anglestemplate = []
	angles = []
	ID = 1

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

					# get the angle atom types and angle type ID 
					atom1_type = atomTypes[atom1[0]][0]
					vertex_type = atomTypes[vertex[0]][0]
					atom2_type = atomTypes[atom2[0]][0]
					angle_type = '%s-%s-%s' %(atom1_type, vertex_type, atom2_type)
					angle_type_ID = angleTypes[angle_type][0]

					# for the initial parameterization, store the angle data 
					angles.append((atom1_type, vertex_type, atom2_type, angle))

					# update the angle template 
					# lammps uses index 1 to N
					anglestemplate.append([ID, angle_type_ID, bonds[j]+1, i+1, bonds[k]+1])
					ID += 1

	return angles, anglestemplate

#-------------------------------------------------------------------------------------------------#

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


def ComputeDihedrals(data, AM, atomTypes, dihedralTypes):
	# this program computes the dihedral template
	# lammps format: [ID, dihedral_type_ID, idx1, idx2, idx3, idx4]
	dihedralstemplate = []
	dihedrals = [] 

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
		
		# get the dihedral type 
		atom1_type = atomTypes[data[path[0]][0]][0]
		atom2_type = atomTypes[data[path[1]][0]][0]
		atom3_type = atomTypes[data[path[2]][0]][0]
		atom4_type = atomTypes[data[path[3]][0]][0]
		dihedral_type = '%s-%s-%s-%s' %(atom1_type, atom2_type, atom3_type, atom4_type)
		dihedral_type_ID = dihedralTypes[dihedral_type][0]

		# create the dihedral template 
		dihedralstemplate.append([ID, dihedral_type_ID, path[0]+1, 
									path[1]+1, path[2]+1, path[3]+1])
	
		# update the dihedral data
		dihedrals.append((atom1_type, atom2_type, atom3_type, atom4_type, dihedral))

	return dihedrals, dihedralstemplate

#-------------------------------------------------------------------------------------------------#

def getMasses(atomTypes):
	massData = []
	for atomID in atomTypes:
		atom_name = atomTypes[atomID][0]
		mass = Masses[atom_name]
		massData.append((atomID, mass))

	massData = sorted(massData, key=lambda x: x[0])
	return massData


def getBondCoeffs(bonds, bondTypes):
	# get the bond lengths 
	uniqueBonds = {}
	for bond in bonds:
		atom1, atom2, bond_length = bond
		bond_type, bond_coeff = bondTypes['%s-%s' %(atom1,atom2)]
		if bond_type in uniqueBonds:
			current_bond_length = uniqueBonds[bond_type][1]

			if abs(bond_length - current_bond_length) > 0.01:
				raise RuntimeError, "The %s-%s bonds have differences > 0.01 A. Check bond data!"\
											%(atom1, atom2)
		else:
			uniqueBonds[bond_type] = (bond_coeff, bond_length)

	# create the bond coeffs
	bondCoeffs = []
	for bond_type in uniqueBonds:
		bond_coeff, bond_length = uniqueBonds[bond_type]
		bondCoeffs.append((bond_type, bond_coeff, bond_length))

	bondCoeffs = sorted(bondCoeffs, key=lambda x: x[0])
	return bondCoeffs


def getAngleCoeffs(angles, angleTypes):
	# get the angles 
	uniqueAngles = {}
	for angle in angles:
		atom1, vertex, atom2, angle_size = angle
		angle_type, angle_coeff = angleTypes['%s-%s-%s' %(atom1, vertex, atom2)]
		if angle_type in uniqueAngles:
			current_angle_size = uniqueAngles[angle_type][1]

			if abs(current_angle_size - angle_size) > 1.0:
				raise RuntimeError, "The %s-%s-%s angles have differences > 1.0 deg. Check angle data!"\
											%(atom1, vertex, atom2)
		else:
			uniqueAngles[angle_type] = (angle_coeff, angle_size)

	# create the angle coeffs
	angleCoeffs = []
	for angle_type in uniqueAngles:
		angle_coeff, angle_size = uniqueAngles[angle_type]
		angleCoeffs.append((angle_type, angle_coeff, angle_size))

	angleCoeffs = sorted(angleCoeffs, key=lambda x: x[0])
	return angleCoeffs


def getDihedralCoeffs(dihedralTypes):
	uniqueDihedrals = {}
	for dihedral in dihedralTypes:
		dihedral_type, K1, phi1, K2, phi2, K3, phi3 = dihedralTypes[dihedral]
		if not dihedral_type in uniqueDihedrals:
			uniqueDihedrals[dihedral_type] = (K1, phi1, K2, phi2, K3, phi3)

	dihedralCoeffs = []
	for dihedral_type in uniqueDihedrals:
		K1, phi1, K2, phi2, K3, phi3 = uniqueDihedrals[dihedral_type]
		dihedralCoeffs.append((dihedral_type, K1, phi1, K2, phi2, K3, phi3))

	dihedralCoeffs = sorted(dihedralCoeffs, key=lambda x:x[0])
	return dihedralCoeffs


#-------------------------------------------------------------------------------------------------#
def WriteBonds(bonds, dir_name, filename):
	filename = '%s%s_bonds.txt' %(dir_name, filename[11:-4])
	f = open(filename, 'w')

	for bond in bonds:
		atom1, atom2, d = bond
		f.write('%s\t%s\t%.4f\n' %(atom1, atom2, d))
	f.close()


def WriteAngles(angles, dir_name, filename):
	filename = '%s%s_angles.txt' %(dir_name, filename[11:-4])
	f = open(filename, 'w')

	for angle in angles:
		atom1, vertex, atom2, theta = angle
		f.write('%s\t%s\t%s\t%.4f\n' %(atom1, vertex, atom2, theta))
	f.close()


def WriteDihdrals(dihedrals, dir_name, filename):
	filename = '%s%s_dihedrals.txt' %(dir_name, filename[11:-4])
	f = open(filename, 'w')

	for dihedral in dihedrals:
		atom1, atom2, atom3, atom4, phi = dihedral
		f.write('%s\t%s\t%s\t%s\t%.4f\n' %(atom1, atom2, atom3, atom4, phi))
	f.close()


def WriteAtomTemplate(atomstemplate, dir_name, filename):
	filename = '%s%s_atomTemplate.txt' %(dir_name, filename[11:-4])
	f = open(filename, 'w')

	for atom in atomstemplate:
		ID, molID, atom_type, charge, xcoord, ycoord, zcoord = atom
		f.write('%d\t%d\t%d\t%.4f\t\t%.6f\t\t%.6f\t\t%.6f\n'\
				%(ID, molID, atom_type, charge, xcoord, ycoord, zcoord))
	f.close()


def WriteBondTemplate(bondstemplate, dir_name, filename):
	filename = '%s%s_bondTemplate.txt' %(dir_name, filename[11:-4])
	f = open(filename, 'w')

	for bond in bondstemplate:
		ID, bond_type, idx1, idx2 = bond
		f.write('%d\t%s\t%d\t%d\n' %(ID, bond_type, idx1, idx2))
	f.close()


def WriteAngleTemplate(anglestemplate, dir_name, filename):
	filename = '%s%s_angleTemplate.txt' %(dir_name, filename[11:-4])
	f = open(filename, 'w')

	for angle in anglestemplate:
		ID, angle_type, idx1, idx2, idx3 = angle
		f.write('%d\t%s\t%d\t%d\t%d\n' %(ID, angle_type, idx1, idx2, idx3))
	f.close()


def WriteDihedralTemplate(dihedralstemplate, dir_name, filename):
	filename = '%s%s_dihedralTemplate.txt' %(dir_name, filename[11:-4])
	f = open(filename, 'w')

	for dihedral in dihedralstemplate:
		ID, dihedral_type, idx1, idx2, idx3, idx4 = dihedral
		f.write('%d\t%s\t%d\t%d\t%d\t%d\n' %(ID, dihedral_type, idx1, idx2, idx3, idx4))
	f.close()

#-------------------------------------------------------------------------------------------------#
def GetFileTypes(filename):
	atomTypesFile = '%s_atomTypes.txt' %filename[:-4]
	bondTypesFile = '%s_bondTypes.txt' %filename[:-4]
	angleTypesFile = '%s_angleTypes.txt' %filename[:-4]
	dihedralTypesFile = '%s_dihedralTypes.txt' %filename[:-4]

	return (atomTypesFile, bondTypesFile,\
			angleTypesFile, dihedralTypesFile)

def MakeResultsDirectory():
	# make nice path slug if spaces are given
	dir_name = 'structure_data/'
	# make the results directory if it dosent exit
	if not os.path.exists(dir_name):
		os.makedirs(dir_name)
	return dir_name

#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
# main progam

if len(sys.argv) < 2:
	print 'Usage:'
	print '  python %s <Usage Info: True/False>' %sys.argv[0]
	exit()

t0 = time.time()

usageInfo = sys.argv[1]

if usageInfo == 'True':
	print '''\nProvide the "atom file" with the xyz position of 
	each atom. The atom types, bond types, angle types, 
	and dihedral types must be defined in separate files.
	The naming scheme for the files are as follows:

	atom file: precursorName.xyz (Note: let precursorName = *)
	atom types: *_atomTypes.txt
	bond types: *_bondTypes.txt
	angle types: *_angleTypes.txt
	dihedral types: *_dihedralTypes.txt

	IMPORTANT: place all files in the 'inputfiles' directory.
	'''
	exit()

# get the files in the inputfiles directory 
atomfiles = glob.glob('inputfiles/*.xyz')

for atomfile in atomfiles:
	print '\nWorking with %s...' %atomfile
	filenames = GetFileTypes(atomfile)

	data, atomstemplate, atomTypes,\
		bondTypes, angleTypes, dihedralTypes = LoadData(atomfile, filenames)
	# print data
	# print atomTypes

	# compute bonds, angles, dihedrals
	bonds, AM, bondstemplate = ComputeBonds(data, atomTypes, bondTypes)
	angles, anglestemplate = ComputeAngles(data, AM, atomTypes, angleTypes)
	dihedrals, dihedralstemplate = ComputeDihedrals(data, AM, atomTypes, dihedralTypes)

	# get the bond, angle and dihedral coeffs
	masses = getMasses(atomTypes)
	bond_coeffs = getBondCoeffs(bonds, bondTypes)
	angle_coeffs = getAngleCoeffs(angles, angleTypes)
	dihedral_coeffs = getDihedralCoeffs(dihedralTypes)

	# write out data 
	dir_name = MakeResultsDirectory()
	WriteBonds(bonds, dir_name, atomfile)
	WriteAngles(angles, dir_name, atomfile)
	WriteDihdrals(dihedrals, dir_name, atomfile)
	WriteAtomTemplate(atomstemplate, dir_name, atomfile)
	WriteBondTemplate(bondstemplate, dir_name, atomfile)
	WriteAngleTemplate(anglestemplate, dir_name, atomfile)
	WriteDihedralTemplate(dihedralstemplate, dir_name, atomfile)

	print 'Number of atoms = %d' %len(data)
	print 'Number of bonds = %d' %len(bonds)
	print 'Number of angles = %d' %len(angles)
	print 'Number of dihedrals = %d' %len(dihedrals)

print '\nAnalyzed precursor structure in %.4f seconds.\n' %(time.time()-t0)


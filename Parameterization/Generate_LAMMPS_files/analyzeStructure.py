import copy
import math
import numpy
import sys
import time

#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
''' 
This progam analyzes a precursor structure and creates templates that include 
the atom positions, bonds, angles, and dihedrals
'''


def LoadData(inputfile):
	""" reads .xyz file """ 

	data = []
	NSi = 0
	NO = 0
	NC = 0

	f = open(inputfile)

	f.readline()

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
			else: 
				raise RuntimeError, "Incorrect atom type."

			# populate the data array
			data.append([atomtype,xcoord,ycoord,zcoord])

		else:
			break
	f.close()

	return numpy.array(data)

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

def ComputeBonds(data):
	
	# define the adjacency matrix to help compute angles and dihedrals
	AM = numpy.zeros((len(data),len(data)))

	# bond template will be in the lammps format [ID, type, atomID1, atomID2]
	bondstemplate = []
	ID = 1

	SiC_bond = 1 # ~1.86 A
	CC_single_bond = 2 # ~1.53 A
	CC_double_bond = 3 # ~1.31 A

	# precursor molecules are small so we will compute all distances
	for i in range(len(data)):
		atom1 = data[i][0]
		for j in range(i+1,len(data)):
			atom2 = data[j][0]			
	
			dx = abs(data[j][1] - data[i][1])
			dy = abs(data[j][2] - data[i][2])
			dz = abs(data[j][3] - data[i][3])

			d = Distance(dx,dy,dz)
			
			# Si-C bond (bond type 1)
			if (atom1==1 and atom2==3) or (atom1==3 and atom2==1):
				if d < 2.0:
					# lammps uses index 1 to N
					bondstemplate.append([ID, SiC_bond, i+1, j+1])
					ID += 1

					# update the adjacency matrix
					AM[i][j] = d
					AM[j][i] = d

			# C-C bond (bond type 2 and 3)
			if atom1==3 and atom2==3:
				if d < 1.7:
					if d < 1.4:
						bondstemplate.append([ID, CC_double_bond, i+1, j+1])
					else:
						bondstemplate.append([ID, CC_single_bond, i+1, j+1])

					ID += 1

					# update the adjacency matrix
					AM[i][j] = d
					AM[j][i] = d

	return bondstemplate, AM

#-------------------------------------------------------------------------------------------------#

def ComputeAngles(data, AM):

	# angle template will be in the lammps format [ID, type, atomID1, atomID2 (vertex), atomID3]
	anglestemplate = []
	ID = 1

	angle_tet = 1 # tetradedral (109.5)
	angle_trig = 2 # trigonal (120.0)

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
					print angle

					if angle >= 106 and angle <= 112:
						# lammps uses index 1 to N
						anglestemplate.append([ID, angle_tet, bonds[j]+1, i+1, bonds[k]+1])
						ID += 1

					elif angle >= 113 and angle <=125:
						anglestemplate.append([ID, angle_trig, bonds[j]+1, i+1, bonds[k]+1])
						ID += 1

	return anglestemplate

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


def GetDihedralType(path, dihedral):
	return 1

def ComputeDihedrals(data, AM):
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
		# print dihedral, path
		dihedraltype = GetDihedralType(path, dihedral)

		dihedralstemplate.append([ID, dihedraltype, path[0]+1, 
									path[1]+1, path[2]+1, path[3]+1])


	return dihedralstemplate

	# for i in range(len(AM)):
	# 	paths = []
	# 	while True:	
	# 		path = [i]
	# 		Nsteps = 0

	# 		while Nsteps < path_length:
	# 			# Nsteps0 = Nsteps

	# 			next_steps = numpy.nonzero(AM[i])[0]

	# 			if len(next_steps) > 0:
	# 				Nsteps += 1

	# 				for step in next_steps:

	# 					''' think of branching problem; need to add steps and whatever steps are generated next'''
	# 					path.append(step)
	# 			else:
	# 				break

	# 			# # if no steps could be make, stop looking
	# 			# if Nsteps0 == Nsteps:
	# 			# 	break 

#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#


if len(sys.argv) < 2:
	print 'Usage:'
	print '  python %s <inputfile>' %sys.argv[0]
	exit()

inputfile = sys.argv[1]

t0 = time.time()

data = LoadData(inputfile)
data = ShiftToPositiveCoords(data)
print data

bondstemplate, AM = ComputeBonds(data)
print bondstemplate
print AM

angletemplate = ComputeAngles(data, AM)
print angletemplate

dihedralstemplate = ComputeDihedrals(data, AM)
print dihedralstemplate

print '\n\n'
print 'Number of atoms = %d' %len(data)
print 'Number of bonds = %d' %len(bondstemplate)
print 'Number of angles = %d' %len(angletemplate)


print 'Analyzed network in %.4f seconds.' %(time.time()-t0)
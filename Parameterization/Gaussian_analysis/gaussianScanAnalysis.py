import copy
import math
import numpy
import operator
import sys
import time

from StructureAnalysisFuctions import ComputeBonds, ComputeAngles, ComputeDihedrals

#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#

def getCoordinates(line):
	atomtype = 'NotSupported'
	fields = line.split()
	atomic_number = fields[1]
	if atomic_number == '14':
		atomtype = '1'
	elif atomic_number == '8':
		atomtype = '2'
	elif atomic_number == '6':
		atomtype = '3'
	elif atomic_number == '1':
		atomtype = '0'

	xcoord = fields[3]
	ycoord = fields[4]
	zcoord = fields[5]
	return (atomtype, xcoord, ycoord, zcoord)


def LoadData(inputfile):
	""" reads Gaussian .log file """ 


	# end of document flag
	stop_flag = 'termination of Gaussian'

	# flags used to find the number of atoms
	atom_flag = 'Title Card Required'
	atom_active = False	

	# flags to get the atom coordinates 
	opt_flag = 'Optimized Parameters'
	coord_flag = 'Standard orientation'
	coord_active = False
	coord_read_active = False
	all_opt_coords = []
	atom_coords = []

	f = open(inputfile)

	Nline = 0
	Natoms = 0
	Nskip_lines = 0
	Nlines_read = 0
	Natoms_read = 0

	while True:
		line = f.readline().strip()
		if stop_flag in line:
			break

		Nline += 1

		# look for the title card to get the number of atoms
		if atom_flag in line:
			atom_active = True
			# ensure that atom flag is not used again
			atom_flag = 'Title Card Required Already Found!!!!!!'
			continue

		if atom_active:
			Nskip_lines += 1
			if line:
				# skip 4 lines to get to the atom data
				if Nskip_lines > 3:
					Natoms += 1
			else:
				atom_active = False
				continue


		# if the stop_flag has not been reached, search 
		# for the optimization flag and then the atom coords
		if opt_flag in line:
			coord_active = True
			continue

		if coord_active:
			if coord_flag in line:
				coord_active = False
				coord_read_active = True
				continue

		if coord_read_active:
			Nlines_read += 1
			if Nlines_read > 4:
				Natoms_read += 1
				if Natoms_read > Natoms:
					coord_read_active = False
					Nlines_read = 0
					Natoms_read = 0
					all_opt_coords.append(atom_coords)
					atom_coords = []
				else:
					coords = getCoordinates(line)
					atom_coords.append(coords)
	f.close()

	return (Nline, Natoms, all_opt_coords)


def ComputeBox(data):
	x_min = min(data[:,1]) - 5
	x_max = max(data[:,1]) + 5 

	y_min = min(data[:,2]) - 5 
	y_max = max(data[:,2]) + 5 

	z_min = min(data[:,3]) - 5
	z_max = max(data[:,3]) + 5 

	return (x_min, x_max, y_min, y_max, z_min, z_max)


def CleanAtoms(coords):
	''' Do not inlcude H atoms in LAMMPS files (modeled implicitly) '''
	atoms = []
	for atom in coords:
		if atom[0] != '0':
			atoms.append([int(atom[0]), float(atom[1]), 
						float(atom[2]), float(atom[3])])

	return numpy.array(atoms)


def GenerateLAMMPSDataFile(angle_start, angle_step, dihedral, precursor, all_opt_coords):
	''' this function prepares a LAMMPS data file '''

	dihedral_angles = []
	filename_base = 'data.%s_%s_{}_deg' %(precursor, dihedral)
	filename_base_lmp = 'data.%s_%s_$a_deg' %(precursor, dihedral)

	dihedral_angle = angle_start

	for coords in all_opt_coords:
		filename = filename_base.format(dihedral_angle)
		dihedral_angles.append(dihedral_angle)
		dihedral_angle += angle_step
		

		atoms = CleanAtoms(coords)
		bonds, AM, bond_types = ComputeBonds(atoms)
		angles, angle_types = ComputeAngles(atoms, AM)
		dihedrals, dihedral_potentials = ComputeDihedrals(atoms, AM, parameterization=True)

		# sort bond and angle types by increasing ID number (returns list of tuples)
		bond_types = sorted(bond_types.items(), key=operator.itemgetter(1))
		angle_types = sorted(angle_types.items(), key=operator.itemgetter(1))


		WriteLAMMPSDataFile(filename, atoms, bonds, bond_types, \
							angles, angle_types, dihedrals)

	return dihedral_angles, filename_base_lmp



def WriteCoords(angle_start, angle_step, dihedral, precursor, all_opt_coords):
	filename_base = '%s_%s_opt_structure_{}_deg.xyz' %(precursor, dihedral)
	angle = angle_start

	for coords in all_opt_coords:
		filename = filename_base.format(angle)
		f = open(filename, 'w')

		for atom in coords:
			f.write('%s\t%s\t%s\t%s\n' %(atom[0],atom[1],atom[2],atom[3]))
		f.close()
		angle += angle_step


def WriteLAMMPSDataFile(filename, atoms, bonds, bond_types, angles, angle_types, dihedrals):
	dimensions = ComputeBox(atoms)

	f = open(filename, 'w')
	f.write('# optimized configuration from Gaussian\n')
	f.write('%d %d xlo xhi\n' %(dimensions[0], dimensions[1]))
	f.write('%d %d ylo yhi\n' %(dimensions[2], dimensions[3]))
	f.write('%d %d zlo zhi\n\n' %(dimensions[4], dimensions[5]))

	f.write('%d atoms\n' %len(atoms))
	f.write('%d bonds\n' %len(bonds))
	f.write('%d angles\n' %len(angles))
	f.write('%d dihedrals\n\n' %len(dihedrals))

	f.write('3 atom types\n')
	f.write('%d bond types\n' %len(bond_types))
	f.write('%d angle types\n' %len(angle_types))
	f.write('1 dihedral types\n\n')

	f.write('Masses\n\n')
	f.write('\t1\t28\n')
	f.write('\t2\t16\n')
	f.write('\t3\t12\n\n')

	f.write('Angle Coeffs\n\n')
	for angle_type in angle_types:
		f.write('\t%d\t3.38\t%.4f\n' %(angle_type[1], angle_type[0]))
	# f.write('\t1\t3.38\t109.5\n')
	# f.write('\t2\t3.38\t120\n\n')

	f.write('\nBond Coeffs\n\n')
	for bond_type in bond_types:
		f.write('\t%d\t20\t%.4f\n' %(bond_type[1], bond_type[0]))
	# f.write('\t1\t20\t1.86\n')
	# f.write('\t2\t20\t1.53\n')
	# f.write('\t3\t20\t1.31\n\n')

	# dihedrals are set to 0 initially to calibrate the opls potential
	f.write('\nDihedral Coeffs\n\n')
	f.write('\t1\t0\t0\t0\t0\n\n\n')

	f.write('Atoms\n\n')
	charge = 0
	for i in range(len(atoms)):
		if atoms[i][0] == 1:
			charge = 4
		elif atoms[i][0] == 2:
			charge = -2
		elif atoms[i][0] == 3:
			charge = 0

		f.write('%d\t%d\t%d\t%d\t%.4f\t%.4f\t%.4f\n' %(i+1, 1, atoms[i][0], charge, 
													atoms[i][1], atoms[i][2], 
													atoms[i][3]))

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


def WriteLAMMPSInputFile(dihedral_angles, dihedral, precursor, filename_base_lmp):
	filename = 'in.%s_%s' %(precursor, dihedral)
	f = open(filename, 'w')

	f.write('# Computes point energy of torsion scans from Gaussian\n\n')
	f.write('variable a index ')
	for angle in dihedral_angles:
		f.write('%d ' %angle)

	# put all the point energies in the smae log file
	f.write('\nlog log_$a.lammps\n\n')

	f.write('dimension 3\n')
	f.write('atom_style full\n')

	# bonded potentials
	f.write('bond_style harmonic\n')
	f.write('angle_style harmonic\n')
	f.write('dihedral_style opls\n\n')

	f.write('boundary p p p\n')
	f.write('units metal\n')
	f.write('neighbor 2 bin\n')
	f.write('neigh_modify check yes\n\n')

	f.write('timestep 0.001\n\n')
	f.write('read_data %s\n\n' %filename_base_lmp)

	# non-bonded potential
	f.write('pair_style sw\n')
	f.write('pair_coeff * * OCS.sw Si O C\n\n')

	f.write('thermo_style custom step temp tave press pave lx ly lz vol epair ke pe etotal\n')
	f.write('thermo 1\n')
	f.write('minimize 1e-4 0 0\n\n')
	f.write('clear\n')
	f.write('next a\n')
	f.write('jump %s' %filename)
	f.close()



#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#

if len(sys.argv) < 6:
	print 'Usage:'
	print '''  python %s 
				<inputfile> 
				<precursor [e.g. EtOCS]> 
				<dihedral name [e.g. Si-C-C-Si]> 
				<angle start> 
				<angle step size>
			'''%sys.argv[0]
	exit()

inputfile = sys.argv[1]
precursor = sys.argv[2]
dihedral = sys.argv[3]
angle_start = int(sys.argv[4])
angle_step = int(sys.argv[5])


t0 = time.time()

Nlines, Natoms, all_opt_coords = LoadData(inputfile)
dihedral_angles, filename_base_lmp = GenerateLAMMPSDataFile(angle_start, angle_step, 
													dihedral, precursor, 
													all_opt_coords)
WriteLAMMPSInputFile(dihedral_angles, dihedral, precursor, filename_base_lmp)
# WriteCoords(angle_start, angle_step, dihedral, precursor, all_opt_coords)


print 'Number of atoms = %d' %Natoms
print "Number of Optimizations = %d" %len(all_opt_coords)
print 'Number of lines analyzed = %d' %Nlines
print 'Analyzed network in %.4f seconds.' %(time.time()-t0)
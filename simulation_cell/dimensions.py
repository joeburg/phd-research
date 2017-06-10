import numpy
import sys


#------------------------------------------------------------------------------#

def LoadData(inputfile):
	""" reads LAMMPS file, sorts data, and computes cell dimensions """ 
	f = open(inputfile)
	Natoms = int(f.readline())

	f.readline()

	data = numpy.zeros((Natoms,4))

	NSi = 0
	NO = 0
	NC = 0
	for i in range(Natoms):
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
			data[i] = [atomtype,xcoord,ycoord,zcoord]
	f.close()
	return data

def LoadDataMolFile(inputfile):
	f = open(inputfile)
	f.readline()
	data = []
	while True:
		fields = f.readline().strip().split()
		if fields:
			atomtype = int(fields[0])
			xcoord = float(fields[1])
			ycoord = float(fields[2])
			zcoord = float(fields[3])

			data.append([atomtype,xcoord,ycoord,zcoord])
		else:
			break
	f.close()
	return numpy.array(data)

def Distance(i, j, data):
	dx = abs(data[j][1] - data[i][1])
	dy = abs(data[j][2] - data[i][2])
	dz = abs(data[j][3] - data[i][3])
	return (dx*dx + dy*dy + dz*dz)**0.5


def CellDimensions(data):
	# find the dimensions of the simulation cell 
	x_min = min(data[:,1])
	x_max = max(data[:,1])
	Lx = x_max - x_min

	z_min = min(data[:,3])
	z_max = max(data[:,3])
	Lz = z_max - z_min

	y_min = min(data[:,2])
	y_max = max(data[:,2])
	Ly = y_max - y_min
	y_avg = 0.5*(y_max + y_min)

	return (Lx, Ly, Lz)

#------------------------------------------------------------------------------#
# main program 
if len(sys.argv) < 3:
	print 'Usage:'
	print '  python %s <inputfile> <file type: mol, lammps>' %sys.argv[0]
	exit()

inputfile = sys.argv[1]
filetype = sys.argv[2]

if filetype == 'mol':
	data = LoadDataMolFile(inputfile)
	distance = input('\nCompute distance between atoms?\n'+\
						'\t1) yes\n'+\
						'\t2) no\n'+\
						'\tchoice: ')
	if distance == 1:
		atom1 = input('Give the id of atom 1: ')
		atom2 = input('Give the id of atom 2: ')
		d = Distance(atom1-1, atom2-1, data)
		print data[atom1-1]
		print data[atom2-1]
		print 'Distance between %d and %d: %.4f' %(atom1, atom2, d)

elif filetype == 'lammps':
	data = LoadData(inputfile)


Lx, Ly, Lz = CellDimensions(data)
print 'Cell dimensions: %.2f, %.2f, %.2f' %(Lx,Ly,Lz)

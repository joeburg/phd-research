

def LoadData(inputfile):
	""" reads LAMMPS .xyz file, sorts data, and computes cell dimensions """ 
	f = open(inputfile)
	Natoms = int(f.readline())

	f.readline()

	data = set([])
	while True:
		fields = f.readline().strip().split()
		if fields:
			atomtype = fields[0]
			xcoord = float(fields[1])
			ycoord = float(fields[2])
			zcoord = float(fields[3])

			# populate the data array
			data.add((atomtype,xcoord,ycoord,zcoord))
		else:
			break
	f.close()
	return data

def GetUniqueData(data):
	clean_data = set([])
	for atom in data:
		if atom not in clean_data:
			clean_data.add(atom)
	return clean_data

def WriteAtoms(data, inputfile):
	Natoms = len(data)
	filename = '%s_clean.xyz' %inputfile[:-4]
	f = open(filename, 'w')
	f.write('%d\n' %Natoms)
	f.write('Terminal Groups atoms\n')	
	for atom in data:
		atomtype, x, y, z = atom
		f.write('%s  %.4f  %.4f  %.4f\n' %(atomtype, x, y, z))		
	f.close()	


#-----------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------#

import glob
import sys
import time


if len(sys.argv) < 2:
	print 'Usage:'
	print '  python %s <filname>' %sys.argv[0]
	exit()

filename = sys.argv[1]

t0 = time.time()

data = LoadData(filename)
clean_data = GetUniqueData(data)
WriteAtoms(data, filename)

print 'Analyzed network in %.4f seconds.' %(time.time()-t0)
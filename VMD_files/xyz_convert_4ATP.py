#The purpose of this program is to convert .xyz files from LAMMPS
# to a form in which VMD can read the .xyz file
#We must convert the species identification to its atomic number
#OCS model

from Tkinter import Tk
from tkFileDialog import askopenfilename
Tk().withdraw()

#main program

filename = askopenfilename()
print "Working with file:", filename

#----------------------------------------------------------------------------#
data = []
f = open(filename)
Natoms = f.readline().strip()
description = f.readline().strip()

while True:
	fields = f.readline().strip().split()
	if fields:
		atom_type, x, y, z = fields
		# if atom_type == '9' or atom_type == '10':
		# 	data.append([atom_type, x, y, z])
		data.append([atom_type, x, y, z])
	else:
		break
f.close()

#----------------------------------------------------------------------------#
atomTypes = {
	'1' : 'H',
	'2' : 'H',
	'3' : 'H',
	'4' : 'C',
	'5' : 'C',
	'6' : 'C',
	'7' : 'C',
	'8' : 'N',
	'9' : 'S',
	'10': 'Cu',
}

for i in range(len(data)):
	data[i][0] = atomTypes[data[i][0]]

#----------------------------------------------------------------------------#

#write to a text file
#the .join() method takes an array, i, and concantenates all the elements together
# with a space " " between each element.  Then a newline "\n" is added to make sure
# your output is broken up into separate lines

datafile = '%s_VMD.xyz' %filename[:-4]
# datafile = '%s_Cu_only_VMD.xyz' %filename[:-4]
f = open(datafile, 'w')
f.write('%s\n' %len(data))
f.write('%s\n' %description)
for atom in data:
	atom_type, x, y, z = atom
	f.write('%s\t%s\t%s\t%s\n' %(atom_type, x, y, z))
f.close()

print "All done!" 

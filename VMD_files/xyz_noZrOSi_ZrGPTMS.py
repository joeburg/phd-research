#remove Zr,O,Si from .xyz simulation files 
from Tkinter import Tk
from tkFileDialog import askopenfilename
Tk().withdraw()

import numpy as np

#main program

filename = askopenfilename()
print "Working with file:", filename

atom_data = []
with open(filename) as inputfile:
    for line in inputfile:
        atom_data.append(line.strip().split( ))

#access number of atoms from file and remove first 2 lines 
N_atoms = int(atom_data[0][0])
atom_data = atom_data[2:]

#transpose data into columns
atom_data = [list(x) for x in zip(*atom_data)]

atom_type = atom_data[0]
x_coord = atom_data[1]
y_coord = atom_data[2]
z_coord = atom_data[3]

no_ZrOSi_atom = []
for i in range(len(atom_type)):
    if atom_type[i]=='3' or atom_type[i]=='5':
        no_ZrOSi_atom.append([i,atom_type[i]])


N_atoms = int(len(no_ZrOSi_atom))
no_ZrOSi_atom = [list(x) for x in zip(*no_ZrOSi_atom)]

atom_type = no_ZrOSi_atom[1]
index = no_ZrOSi_atom[0]

for i in range(N_atoms):
    if atom_type[i]=='3':
        atom_type[i]='C'
    elif atom_type[i]=='5':
        atom_type[i]='O'

new_x_coord = []
for i in range(len(x_coord)):
    if i in index:
        new_x_coord.append(x_coord[i])

new_y_coord = []
for i in range(len(y_coord)):
    if i in index:
        new_y_coord.append(y_coord[i])

new_z_coord = []
for i in range(len(z_coord)):
    if i in index:
        new_z_coord.append(z_coord[i])

data = [atom_type] + [new_x_coord] + [new_y_coord] + [new_z_coord]
data = [list(x) for x in zip(*data)]

data = [[str(N_atoms)]]+[['Atoms']]+data

dataFile = open(filename[:-4]+"_noZrOSi_VMD.xyz", 'w')
for eachitem in data:
    dataFile.write("\t".join(eachitem)+'\n')

dataFile.close()


print 'All done!'

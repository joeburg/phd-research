#remove C from .xyz simulation files 
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

no_Si_C_atoms = []
for i in range(len(atom_type)):
    if atom_type[i]=='2' or atom_type[i]=='4':
        no_Si_C_atoms.append([i,atom_type[i]])


N_atoms = int(len(no_Si_C_atoms))
no_Si_C_atoms = [list(x) for x in zip(*no_Si_C_atoms)]

atom_type = no_Si_C_atoms[1]
index = no_Si_C_atoms[0]

for i in range(N_atoms):
    if atom_type[i]=='2':
        atom_type[i]='O'
    elif atom_type[i]=='4':
        atom_type[i]='Zr'

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

dataFile = open(filename[:-4]+"_noSiC_VMD.xyz", 'w')
for eachitem in data:
    dataFile.write("\t".join(eachitem)+'\n')

dataFile.close()


print 'All done!'

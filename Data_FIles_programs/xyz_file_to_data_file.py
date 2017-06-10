#Program to turn .xyz file into data file

from Tkinter import Tk
from tkFileDialog import askopenfilename
Tk().withdraw()

import numpy as np

#main program

filename = askopenfilename()
print "Working with file:", filename

N_silane = input('How many silane molecules are in the simulation? ')
N_Zr = input('How many Zr atoms are in the simulation? ')
N_O = input('How many O atoms are in the simulation? ')

datafile_name = raw_input('Provide output data file name: ')

#split data into array
data = []
with open(filename) as inputfile:
    for line in inputfile:
        data.append(line.strip().split( ))

#from .xyz file get number of atoms in simulation
N_atoms = int(data[0][0])

#remove the first two rows of .xyz file to obtain relevant data
relevant_xyz_data = []
for i in range(len(data)):
    if i > 1:
        relevant_xyz_data.append(data[i])

#sort the relevant data into columns
relevant_xyz_data = [list(x) for x in zip(*relevant_xyz_data)]

#Put data into correct LAMMPS input format
###############################################################
#Atom portion of data file
#column 1: index
atom_index = []
for i in range(N_atoms):
    atom_index = atom_index + [str(i+1)]
    
#column 2: molecule-tag
#assume 
molecule_tag = []
for i in range(N_silane):
    molecule_tag = molecule_tag + [str(i+1)]*8

for i in range(N_Zr):
    molecule_tag = molecule_tag + [str(0)]

for i in range(N_O):
    molecule_tag = molecule_tag + [str(0)]

#column 3: atom-type
atom_type = relevant_xyz_data[0]

#column 4: charge
charge = []
for i in range(len(atom_type)):
    if atom_type[i] == '1':
        charge = charge + [str(4)]
    elif atom_type[i] == '2':
        charge = charge + [str(-2)]
    elif atom_type[i] == '3':
        charge = charge + [str(0)]
    elif atom_type[i] == '4':
        charge = charge + [str(4)]
    elif atom_type[i] == '5':
        charge = charge + [str(0)]

#column 5: x-coord
x_coord = relevant_xyz_data[1]

#column 6: y-coord
y_coord = relevant_xyz_data[2]

#column 7: z-coord
z_coord = relevant_xyz_data[3]

atom_portion = [atom_index] + [molecule_tag] + [atom_type] + [charge] + [x_coord] + \
               [y_coord] + [z_coord]

atom_portion = [list(x) for x in zip(*atom_portion)]
   

###############################################################
#Bond portion of data file
#column 1: index
N_bonds = N_silane*7
N_atoms_w_bonds = N_silane*8

bond_index = []
for i in range(N_bonds):
    bond_index = bond_index + [str(i+1)]
    
#column 2: bond-type
bond_type = []
for i in range(N_bonds):
    if (i+1) % 7 == 1:
        bond_type = bond_type + [str(1)]
    else:
        bond_type = bond_type + [str(2)]

#column 3: atom-1
#column 4: atom-2
atom_1 = []
atom_2 = []
i = 1
while i <= N_atoms_w_bonds:
    for j in range(7):
        atom_1 = atom_1 + [str(i+j)]
        atom_2 = atom_2 + [str(i+j+1)]
    i = i + 8

bond_portion = [bond_index] + [bond_type] + [atom_1] + [atom_2]

bond_portion = [list(x) for x in zip(*bond_portion)]


###############################################################
#Angle portion of data file

#column 1: index
N_angles = N_silane*6

angle_index = []
for i in range(N_angles):
    angle_index = angle_index + [str(i+1)]

#column 2: angle-type
angle_type = []
for i in range(N_angles):
    angle_type = angle_type + [str(1)]

#column 3: atom-1
#column 4: atom-2
#column 5: atom-3
atom1 = []
atom2 = []
atom3 = []
i = 1
while i <= N_atoms_w_bonds:
    for j in range(6):
        atom1 = atom1 + [str(i+j)]
        atom2 = atom2 + [str(i+j+1)]
        atom3 = atom3 + [str(i+j+2)]
    i = i + 8

angle_portion = [angle_index] + [angle_type] + [atom1] + [atom2] + [atom3]

angle_portion = [list(x) for x in zip(*angle_portion)]


###############################################################
#put together atom, bond and angle portions 

data_array = [['Atoms']]+[['']]+atom_portion+[['']]+[['']]+[['Bonds']]\
             +[['']]+bond_portion+[['']]+[['']]+[['Angles']]+[['']]+angle_portion


###############################################################

dataFile = open(datafile_name, 'w')
for eachitem in data_array:
    dataFile.write(" ".join(eachitem)+'\n')

dataFile.close()


print "All done!" 

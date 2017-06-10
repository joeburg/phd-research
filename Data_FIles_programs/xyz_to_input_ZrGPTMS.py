#Program to turn .xyz file into data file for Zr/GPTMS model

from Tkinter import Tk
from tkFileDialog import askopenfilename
Tk().withdraw()

import numpy as np

###############################################################
#import atom data
atom_data = []
with open('ZrGPTMS.66000.xyz') as inputfile:
    for line in inputfile:
        atom_data.append(line.strip().split( ))

#import bond data
bond_data = []
with open('bonds.txt') as inputfile:
    for line in inputfile:
        bond_data.append(line.strip().split( ))


#access number of atoms from file and remove first 2 lines 
N_atoms = int(atom_data[0][0])
atom_data = atom_data[2:]

#access number of bonds from file and remove first 2 lines 
bond_data = bond_data[2:]
N_bonds = len(bond_data)

###############################################################################
#Find the dimensions of the cell and create cutoff in z-direction
#transpose data into columns
atom_data = [list(x) for x in zip(*atom_data)]

#access column to get specific data
atom_type = atom_data[0]
x_coord = atom_data[1]
y_coord = atom_data[2]
z_coord = atom_data[3]

x_coord = map(float, x_coord)
y_coord = map(float, y_coord)
z_coord = map(float, z_coord)

x_min = int(min(x_coord)-1)
x_max = int(max(x_coord)+1)

y_min = int(min(y_coord)-1)
y_max = int(max(y_coord)+1)

z_min = int(min(z_coord)-1)
z_max = int(max(z_coord)+1)

###############################################################
#Put data into correct LAMMPS input format
#Atom portion of data file
#column 1 - index
atom_index = []
for i in range(N_atoms):
    atom_index = atom_index + [str(i+1)]
    
#column 2 - molecule-tag
#assume 
molecule_tag = []
for i in range(N_atoms):
    molecule_tag.append(str(0))

#column 3 - atom-type

#column 4 - charge
charge = []
for i in range(len(atom_type)):
    if atom_type[i] == '1':
        charge.append(str(4))
    elif atom_type[i] == '2':
        charge.append(str(-2))
    elif atom_type[i] == '3':
        charge.append(str(0))
    elif atom_type[i] == '4':
        charge.append(str(4))
    elif atom_type[i] == '5':
        charge.append(str(0))

#column 5 - x-coord
x_coord = map(str, x_coord)

#column 6 - y-coord
y_coord = map(str, y_coord)

#column 7 - z-coord
z_coord = map(str, z_coord)

atom_portion = [atom_index]+[molecule_tag]+[atom_type]+[charge]+\
                      [x_coord]+[y_coord]+[z_coord]

atom_portion = [list(x) for x in zip(*atom_portion)]


###############################################################################
#Make LAMMPS data file without bonds 
                                        
N_bond_atoms_angles = [[str(N_atoms)+' atoms'],['0 bonds'],['0 angles'],\
                       ['0 dihedrals'],['0 impropers']]

atom_bond_angle_type = [['5 atom types'],['0 bond types'],['0 angle types'],\
                        ['0 dihedral types'],['0 improper types']]

cell_dim_x=[[str(x_min)+' '+str(x_max)+' xlo xhi']]
cell_dim_y=[[str(y_min)+' '+str(y_max)+' ylo yhi']]
cell_dim_z=[[str(z_min)+' '+str(z_max)+' zlo zhi']]

masses = [['Masses'],[''],['\t1 \t28.085'],['\t2 \t15.999'],['\t3 \t12.01'],\
          ['\t4 \t91.22'],['\t5 \t12.01']]

data_array = [['LAMMPS data file']]+[['']]+N_bond_atoms_angles+[['']]+\
             atom_bond_angle_type+[['']]+cell_dim_x+cell_dim_y+cell_dim_z+\
             [['']]+masses+[['']]+[['Atoms']]+[['']]+atom_portion

dataFile = open('data.ZrGPTMS_nobonds', 'w')
for eachitem in data_array:
    dataFile.write(" ".join(eachitem)+'\n')

dataFile.close()


print "All done!" 


#################################################################
###Bond portion of data file
###index, bond-type, atom-1, atom-2
##
##bond_data = [list(x) for x in zip(*bond_data)]
##
###increase indicies by 1 so index goes from 1 to N
##atom1 = bond_data[0]
##for i in range(len(atom1)):
##    atom1[i] = str(int(atom1[i])+1)
##    
##atom2 = bond_data[1]
##for i in range(len(atom2)):
##    atom2[i] = str(int(atom2[i])+1)
##
##bond_data = [list(x) for x in zip(*bond_data)]
##
##bond_portion = []
##for i in range(N_bonds):
##    #C-C bonds are bond-type 1
##    if bond_data[i][2]=='3.0' and bond_data[i][3]=='3.0':
##        bond_portion.append([str(i+1),str(1),atom1[i],atom2[i]])
##
##    #C-Q bonds are bond-type 1
##    elif bond_data[i][2]=='3.0' and bond_data[i][3]=='5.0':
##        bond_portion.append([str(i+1),str(1),atom1[i],atom2[i]])
##
##    #Q-Q bonds are bond-type 1
##    elif bond_data[i][2]=='5.0' and bond_data[i][3]=='5.0':
##        bond_portion.append([str(i+1),str(1),atom1[i],atom2[i]])
##
##    #Si-C bonds are bond-type 2
##    elif bond_data[i][2]=='1.0' and bond_data[i][3]=='3.0':
##        bond_portion.append([str(i+1),str(2),atom1[i],atom2[i]])
##
##    #Si-O bonds are bond-type 3
##    elif bond_data[i][2]=='1.0' and bond_data[i][3]=='2.0':
##        bond_portion.append([str(i+1),str(3),atom1[i],atom2[i]])
##
##    #O-Zr bonds are bond-type 4
##    elif bond_data[i][2]=='2.0' and bond_data[i][3]=='4.0':
##        bond_portion.append([str(i+1),str(4),atom1[i],atom2[i]])
##
#################################################################################
###Make LAMMPS data file with bonds 
##                                        
##N_bond_atoms_angles = [[str(N_atoms)+' atoms'],[str(N_bonds)+' bonds'],['0 angles'],\
##                       ['0 dihedrals'],['0 impropers']]
##
##atom_bond_angle_type = [['5 atom types'],['4 bond types'],['0 angle types'],\
##                        ['0 dihedral types'],['0 improper types']]
##
##cell_dim_x=[[str(x_min)+' '+str(x_max)+' xlo xhi']]
##cell_dim_y=[[str(y_min)+' '+str(y_max)+' ylo yhi']]
##cell_dim_z=[[str(z_min)+' '+str(z_max)+' zlo zhi']]
##
##masses = [['Masses'],[''],['\t1 \t28.085'],['\t2 \t15.999'],['\t3 \t12.01'],\
##          ['\t4 \t91.22'],['\t5 \t12.01']]
##
##bond_coeff = [['Bond Coeffs'],[''],['\t1 \t500 \t1.54'],['\t2 \t500 \t1.9'],\
##              ['\t3 \t250 \t1.63'],['\t4 \t250 \t2.1']]
##
##data_array = [['LAMMPS data file']]+[['']]+N_bond_atoms_angles+[['']]+\
##             atom_bond_angle_type+[['']]+cell_dim_x+cell_dim_y+cell_dim_z+\
##             [['']]+masses+[['']]+bond_coeff+[['']]+[['Atoms']]+[['']]+\
##             atom_portion+[['']]+[['Bonds']]+[['']]+bond_portion
##
#################################################################
##
##dataFile = open('data.ZrGPTMS', 'w')
##for eachitem in data_array:
##    dataFile.write(" ".join(eachitem)+'\n')
##
##dataFile.close()
##
##
##print "All done!" 
##

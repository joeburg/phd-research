#Purpose: remove nonbonded O from substrate simulation

from Tkinter import Tk
from tkFileDialog import askopenfilename
Tk().withdraw()

from pylab import *
from scipy import *
from numpy import *
import numpy as np
import math


################################################################################
################################################################################
def distance(x1,x2,y1,y2,z1,z2):
    d = ((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)**0.5
    return d

def dist(dx,dy,dz):
    d = (dx**2+dy**2+dz**2)**0.5
    return d

################################################################################
#################################################################################import data
filename = askopenfilename()
print "Working with file:", filename

atom_data = []
with open(filename) as inputfile:
    for line in inputfile:
        atom_data.append(line.strip().split( ))

#access number of atoms from file and remove first 2 lines 
N_atoms = int(atom_data[0][0])
atom_data = atom_data[2:]

#get atom type indicies
O_atom_index = []
Si_atom_index = []
for i in range(len(atom_data)):
    if atom_data[i][0]=='1':
        O_atom_index.append(int(i))
    elif atom_data[i][0]=='6':
        Si_atom_index.append(int(i))

N_O = len(O_atom_index)
N_Si = len(Si_atom_index)

print N_O
print N_Si

#transpose data into columns
atom_data = [list(x) for x in zip(*atom_data)]

for i in range(len(atom_data)):
    for j in range(N_atoms):
        atom_data[i][j] = float(atom_data[i][j])

#create first column to index the atoms (0 to N_atoms-1)
atom_index = []
for i in range(N_atoms):
    atom_index.append(i)

#access column to get specific data
atom_type = atom_data[0]
x_coord = atom_data[1]
y_coord = atom_data[2]
z_coord = atom_data[3]

x_min = min(x_coord)
x_max = max(x_coord)
Lx = x_max - x_min

y_min = min(y_coord)
y_max = max(y_coord)
Ly = y_max - y_min

z_min = min(z_coord)
z_max = max(z_coord)
Lz = z_max - z_min

atom_data = [list(x) for x in zip(*atom_data)]
print atom_data[:2]

print ''
print 'Data imported!'

################################################################################
#Find bonds in data; 1=O, 6=Si(S); rows in data will be [index1, index2, atom1, atom2]

print ''
print 'Finding bonds...'


SiO_bonds = []
SiO_bonds_across_BC = []
bond_index = []

for i in range(N_atoms):
    index1 = atom_index[i]
    atom1 = atom_type[i]
    x1 = x_coord[i]
    y1 = y_coord[i]
    z1 = z_coord[i]

    #to account for the PBC's, check if the distance between atoms
    # is greater than Lx/2. If so, translate the larger coord by Lx.
        
    for j in range(i+1,N_atoms):
        index2 = atom_index[j]
        atom2 = atom_type[j]
        x2 = x_coord[j]
        y2 = y_coord[j]
        z2 = z_coord[j]

        dx = abs(x1-x2)
        dy = abs(y1-y2)
        dz = abs(z1-z2)

        #determine Si-O bonds; cutoff at 2.0A; as a slight overestimate, we can say
        # that if the difference in the x-coords is <= 2.0 (we assume that y and z
        # coord are 0) then we will calculate the distance 
        if (atom1==1 and atom2==6) or (atom1==6 and atom2==1):
            d = dist(dx,dy,dz)
            if d<=2.0:
                atoms = sorted([atom1,atom2])
                indicies = sorted([index1,index2])
                index_atom = indicies + atoms + [d]
                if indicies not in bond_index:
                    bond_index.append(indicies)
                    SiO_bonds.append(index_atom)

            #account for PBCs
            if dx > Lx-2.0:
                d = dist(dx+Lx,dy,dz)
                if d <= 2.0:
                    atoms = sorted([atom1,atom2])
                    indicies = sorted([index1,index2])
                    index_atom = indicies + atoms + [d]
                    if indicies not in bond_index:
                        bond_index.append(indicies)
                        SiO_bonds_across_BC.append(index_atom)

                d = dist(dx-Lx,dy,dz)
                if d <= 2.0:                       
                    atoms = sorted([atom1,atom2])
                    indicies = sorted([index1,index2])
                    index_atom = indicies + atoms + [d]
                    if indicies not in bond_index:
                        bond_index.append(indicies)
                        SiO_bonds_across_BC.append(index_atom)

            if dy > Ly-2.0:
                d = dist(dx,dy+Ly,dz)
                if d <= 2.0:
                    atoms = sorted([atom1,atom2])
                    indicies = sorted([index1,index2])
                    index_atom = indicies + atoms + [d]
                    if indicies not in bond_index:
                        bond_index.append(indicies)
                        SiO_bonds_across_BC.append(index_atom)

                d = dist(dx,dy-Ly,dz)
                if d <= 2.0:
                    atoms = sorted([atom1,atom2])
                    indicies = sorted([index1,index2])
                    index_atom = indicies + atoms + [d]
                    if indicies not in bond_index:
                        bond_index.append(indicies)
                        SiO_bonds_across_BC.append(index_atom)

            if dz > Lz-2.0:
                d = dist(dx,dy,dz+Lz)
                if d <= 2.0:
                    atoms = sorted([atom1,atom2])
                    indicies = sorted([index1,index2])
                    index_atom = indicies + atoms + [d]
                    if indicies not in bond_index:
                        bond_index.append(indicies)
                        SiO_bonds_across_BC.append(index_atom)

                d = dist(dx,dy,dz-Lz)
                if d <= 2.0:
                    atoms = sorted([atom1,atom2])
                    indicies = sorted([index1,index2])
                    index_atom = indicies + atoms + [d]
                    if indicies not in bond_index:
                        bond_index.append(indicies)
                        SiO_bonds_across_BC.append(index_atom)


print 'bonds found!'

print ''
print 'Writing data...'

N_SiO_bonds = len(SiO_bonds)
N_SiO_bonds_across_BC = len(SiO_bonds_across_BC)
N_SiO_tot = N_SiO_bonds+N_SiO_bonds_across_BC
SiO_bonds_tot = SiO_bonds + SiO_bonds_across_BC
SiO_bonds_tot = [list(x) for x in zip(*SiO_bonds_tot)]

O_bond_index = SiO_bonds_tot[1]

#if an O atom is bonded, then keep it in file
atom_file = []
for i in range(N_atoms):
    if atom_data[i][0]==6:
        atom_file.append(atom_data[i])
    elif atom_data[i][0]==1:
        if i in O_bond_index:
            atom_file.append(atom_data[i])

#write out bond data
N_atoms = len(atom_file)

for i in range(N_atoms):
    for j in range(4):
        atom_file[i][j]=str(atom_file[i][j])

data = [[str(N_atoms)]]+[['Atoms']]+atom_file

dataFile = open("SiO2_substrate.xyz", 'w')
for eachitem in data:
    dataFile.write("\t".join(eachitem)+'\n')

dataFile.close()

print 'All done!'

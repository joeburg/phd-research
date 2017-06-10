#Purpose: find bonds in unit molecule with atoms in xyz format
# and then put in LAMMPS bond format; use file from replicate_molecule.py
# as the input

from Tkinter import Tk
from tkFileDialog import askopenfilename
Tk().withdraw()

import numpy as np
import math


def distance(x1,x2,y1,y2,z1,z2):
    d = ((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)**0.5
    return d

filename = askopenfilename()
print "Working with file:", filename

data = []
with open(filename) as inputfile:
    for line in inputfile:
        data.append(line.strip().split( ))

data = data[2:]

N_atoms = input('How many atoms are there in the base molecule? ')
bond_type = input('What is the bond type in the LAMMPS simulation? ')

L = len(data)
data = [list(x) for x in zip(*data)]

for i in range(len(data)):
    for j in range(L):
        if i != 0:
            data[i][j] = float(data[i][j])

x_coord = data[1]
y_coord = data[2]
z_coord = data[3]

print x_coord

#90 bonds in C60, 7 bonds from chain connecting backbone and C60,
#N_backbone-1 bonds in the backbone (+1 to link to new unit)
#Then multiply by number of replicated units 
#the bond length of C-chains is 1.54

bond_list = []
k=1
for i in range(N_atoms):
    x1 = x_coord[i]
    y1 = y_coord[i]
    z1 = z_coord[i]
    for j in range(i+1,N_atoms):
        x2 = x_coord[j]
        y2 = y_coord[j]
        z2 = z_coord[j]
        d = distance(x1,x2,y1,y2,z1,z2)
        if d <= 1.65:
            bond_list.append([str(k),str(bond_type),str(i+1),str(j+1)])
            k=k+1
            
dataFile = open(filename[:-4]+"_bonds.txt", 'w')
for eachitem in bond_list:
    dataFile.write("\t".join(eachitem)+'\n')

dataFile.close()


print "All done!" 
          



        

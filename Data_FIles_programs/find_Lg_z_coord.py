#Purpose: find the largest z-coord in the SiO2 substrate 

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

for i in range(N_atoms):
    for j in range(4):
        atom_data[i][j] = float(atom_data[i][j])

#access column to get specific data
atom_data = [list(x) for x in zip(*atom_data)]

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

print z_max

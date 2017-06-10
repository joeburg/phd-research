#Purpose: find the bonds, angles in Zr/GPTMS .xyz outpuf file from simulation 

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

def magnitude(x,y,z):
    mag = (x**2 + y**2 + z**2)**0.5
    return mag

def angle(x1,x2,y1,y2,z1,z2):
    mag1 = magnitude(x1,y1,z1)
    mag2 = magnitude(x2,y2,z2)
    theta = math.acos((x1*x2+y1*y2+z1*z2)/(mag1*mag2))*(180.0/math.pi)
    return theta

def vector(x1,x2,y1,y2,z1,z2):
    vx = x1-x2
    vy = y1-y2
    vz = z1-z2
    return [vx,vy,vz]

def unique_atoms(atoms_involved):
    atoms_unique = []
    for x in atoms_involved:
        if x not in atoms_unique:
            atoms_unique.append(x)
    return atoms_unique

def find_repeat(index1,index2,index3,index4,atom1,atom2,atom3,atom4):
    if index1==index2:
        return [atom1,atom3,atom4] #first index is center 
    elif index1==index3:
        return [atom1,atom2,atom4]
    elif index1==index4:
        return [atom1,atom2,atom3]
    elif index2==index3:
        return [atom2,atom1,atom4]
    elif index2==index4:
        return [atom2,atom1,atom3]
    elif index3==index4:
        return [atom3,atom1,atom2]

#bin the data using the histogram function from numpy

def get_histogram_data(data_list,N_bins):
    bins = linspace(min(data_list)-1,max(data_list)+1,N_bins+1)
    events, edges = histogram(data_list,bins)
    lower = resize(edges, len(edges)-1)
    tmid = lower + 0.5*diff(edges)
    most_prob_event = max(events)
    for i in range(len(events)):
        if most_prob_event == events[i]:
            index = i
    most_prob_tmid = tmid[index]
    return [tmid,events,most_prob_tmid]

def lists_overlap(a,b):
    overlap = []
    for i in a:
        if i in b:
            overlap.append(i)
    return overlap

def remove_overlap(intersection,a):
    difference = []
    for i in a:
        if i not in intersection:
            difference.append(i)
    return difference

################################################################################
################################################################################
#import data
filename = askopenfilename()
print "Working with file:", filename

atom_data = []
with open(filename) as inputfile:
    for line in inputfile:
        atom_data.append(line.strip().split( ))

N_bins_bonds = input('How many bins do you want for the bond histrogram? Note: on average,\
 10 bins will correspond to 0.1 Angstroms) ')

N_bins_angles = input('How many bins do you want for the angle histogram? Note: on average,\
 25 bins will correspond to 0.1 degrees) ')

####Lx = input('What is the final x-dimension of the simulation cell? ')
##
##N_O = input('How many oxygen atoms are in the simulation? ')
##N_Zr = input('How many Zr atoms are in the simulation? ')
##N_silane = input('How many silane molecules are in the simulation? ')
##
##N_C = N_silane*6
##N_Si = N_silane
##N_Q = N_silane

#access number of atoms from file and remove first 2 lines 
N_atoms = int(atom_data[0][0])
atom_data = atom_data[2:]

#get atom type indicies
O_atom_index = []
Zr_atom_index = []
Si_atom_index = []
C_atom_index = []
Q_atom_index = []
for i in range(len(atom_data)):
    if atom_data[i][0]=='2':
        O_atom_index.append(int(i))
    elif atom_data[i][0]=='1':
        Si_atom_index.append(int(i))
    elif atom_data[i][0]=='3':
        C_atom_index.append(int(i))
    elif atom_data[i][0]=='5':
        Q_atom_index.append(int(i))
    elif atom_data[i][0]=='4':
        Zr_atom_index.append(int(i))

N_O = len(O_atom_index)
N_Zr = len(Zr_atom_index)
N_silane = len(Si_atom_index)
N_Si = len(Si_atom_index)
N_Q = len(Q_atom_index)
N_C = len(C_atom_index)

print N_O
print N_Zr
print N_silane
print N_Si
print N_Q
print N_C


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

##print Lx,Ly,Lz

print ''
print 'Data imported!'

################################################################################
#Find bonds in data; 1=Si, 2=O, 3=C, 4=Zr, 5=Q; rows in data will be
# [index1, index2, atom1, atom2]

print ''
print 'Finding bonds...'

CC_bonds = []
CC_bonds_across_BC = []
CQ_bonds = []
CQ_bonds_across_BC = []
CSi_bonds = []
CSi_bonds_across_BC = []
SiO_bonds = []
SiO_bonds_across_BC = []
QQ_bonds = []
QQ_bonds_across_BC = []
ZrO_bonds = []
ZrO_bonds_across_BC = []
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

        #determine C-C bonds (only bonded on chains); cutoff 2.0A
        if atom1==3 and atom2==3:
##            if abs(index1-index2)<10: #only for non-parallel computation
            d = dist(dx,dy,dz)
            if d<=2.0:
                atoms = sorted([atom1,atom2])
                indicies = sorted([index1,index2])
                index_atom = indicies + atoms + [d]
                if indicies not in bond_index:
                    bond_index.append(indicies)
                    CC_bonds.append(index_atom)

            #account for PBCs
            if dx > Lx-2.0:
                d = dist(dx+Lx,dy,dz)
                if d <= 2.0:
                    atoms = sorted([atom1,atom2])
                    indicies = sorted([index1,index2])
                    index_atom = indicies + atoms + [d]
                    if indicies not in bond_index:
                        bond_index.append(indicies)
                        CC_bonds_across_BC.append(index_atom)

                d = dist(dx-Lx,dy,dz)
                if d <= 2.0:                       
                    atoms = sorted([atom1,atom2])
                    indicies = sorted([index1,index2])
                    index_atom = indicies + atoms + [d]
                    if indicies not in bond_index:
                        bond_index.append(indicies)
                        CC_bonds_across_BC.append(index_atom)

            if dy > Ly-2.0:
                d = dist(dx,dy+Ly,dz)
                if d <= 2.0:
                    atoms = sorted([atom1,atom2])
                    indicies = sorted([index1,index2])
                    index_atom = indicies + atoms + [d]
                    if indicies not in bond_index:
                        bond_index.append(indicies)
                        CC_bonds_across_BC.append(index_atom)

                d = dist(dx,dy-Ly,dz)
                if d <= 2.0:
                    atoms = sorted([atom1,atom2])
                    indicies = sorted([index1,index2])
                    index_atom = indicies + atoms + [d]
                    if indicies not in bond_index:
                        bond_index.append(indicies)
                        CC_bonds_across_BC.append(index_atom)

            if dz > Lz-2.0:
                d = dist(dx,dy,dz+Lz)
                if d <= 2.0:
                    atoms = sorted([atom1,atom2])
                    indicies = sorted([index1,index2])
                    index_atom = indicies + atoms + [d]
                    if indicies not in bond_index:
                        bond_index.append(indicies)
                        CC_bonds_across_BC.append(index_atom)

                d = dist(dx,dy,dz-Lz)
                if d <= 2.0:
                    atoms = sorted([atom1,atom2])
                    indicies = sorted([index1,index2])
                    index_atom = indicies + atoms + [d]
                    if indicies not in bond_index:
                        bond_index.append(indicies)
                        CC_bonds_across_BC.append(index_atom)

##                d = distance(x1,x2,y1,y2,z1,z2)
##                if d<=2.0:
##                    atoms = sorted([atom1,atom2])
##                    indicies = sorted([index1,index2])
##                    index_atom = indicies + atoms + [d]
##                    if index_atom not in CC_bonds:
##                        CC_bonds.append(index_atom)
##                        
##                elif d>Lx/2:
##                    if abs(x1-x2)>Lx/2:
##                        if x1>x2:
##                            x1=x1-Lx
##                        else:
##                            x2=x2-Lx
##                    if abs(y1-y2)>Lx/2:
##                        if y1>y2:
##                            y1=y1-Lx
##                        else:
##                            y2=y2-Lx
##                    if abs(z1-z2)>Lx/2:
##                        if z1>z2:
##                            z1=z1-Lx
##                        else:
##                            z2=z2-Lx
##                            
##                    d = distance(x1,x2,y1,y2,z1,z2)
##                    if d<=2.0:
##                        atoms = sorted([atom1,atom2])
##                        indicies = sorted([index1,index2])
##                        index_atom = indicies + atoms + [d]
##                        if index_atom not in CC_bonds_across_BC:
##                            CC_bonds_across_BC.append(index_atom)                    

        #determine C-Q bonds (next to each other in chain); cutoff 2.0A
        elif (atom1==3 and atom2==5) or (atom1==5 and atom2==3):
##            if abs(index1-index2)==1:
            d = dist(dx,dy,dz)
            if d<=2.0:
                atoms = sorted([atom1,atom2])
                indicies = sorted([index1,index2])
                index_atom = indicies + atoms + [d]
                if indicies not in bond_index:
                    bond_index.append(indicies)
                    CQ_bonds.append(index_atom)

            #account for PBCs
            if dx > Lx-2.0:
                d = dist(dx+Lx,dy,dz)
                if d <= 2.0:
                    atoms = sorted([atom1,atom2])
                    indicies = sorted([index1,index2])
                    index_atom = indicies + atoms + [d]
                    if indicies not in bond_index:
                        bond_index.append(indicies)
                        CQ_bonds_across_BC.append(index_atom)

                d = dist(dx-Lx,dy,dz)
                if d <= 2.0:                       
                    atoms = sorted([atom1,atom2])
                    indicies = sorted([index1,index2])
                    index_atom = indicies + atoms + [d]
                    if indicies not in bond_index:
                        bond_index.append(indicies)
                        CQ_bonds_across_BC.append(index_atom)

            if dy > Ly-2.0:
                d = dist(dx,dy+Ly,dz)
                if d <= 2.0:
                    atoms = sorted([atom1,atom2])
                    indicies = sorted([index1,index2])
                    index_atom = indicies + atoms + [d]
                    if indicies not in bond_index:
                        bond_index.append(indicies)
                        CQ_bonds_across_BC.append(index_atom)

                d = dist(dx,dy-Ly,dz)
                if d <= 2.0:
                    atoms = sorted([atom1,atom2])
                    indicies = sorted([index1,index2])
                    index_atom = indicies + atoms + [d]
                    if indicies not in bond_index:
                        bond_index.append(indicies)
                        CQ_bonds_across_BC.append(index_atom)

            if dz > Lz-2.0:
                d = dist(dx,dy,dz+Lz)
                if d <= 2.0:
                    atoms = sorted([atom1,atom2])
                    indicies = sorted([index1,index2])
                    index_atom = indicies + atoms + [d]
                    if indicies not in bond_index:
                        bond_index.append(indicies)
                        CQ_bonds_across_BC.append(index_atom)

                d = dist(dx,dy,dz-Lz)
                if d <= 2.0:
                    atoms = sorted([atom1,atom2])
                    indicies = sorted([index1,index2])
                    index_atom = indicies + atoms + [d]
                    if indicies not in bond_index:
                        bond_index.append(indicies)
                        CQ_bonds_across_BC.append(index_atom)

                
##                d = distance(x1,x2,y1,y2,z1,z2)
##                if d<=2.0:
##                    atoms = sorted([atom1,atom2])
##                    indicies = sorted([index1,index2])
##                    index_atom = indicies + atoms + [d]
##                    if index_atom not in CQ_bonds:
##                        CQ_bonds.append(index_atom)
##                        
##                elif d>Lx/2:
##                    if abs(x1-x2)>Lx/2:
##                        if x1>x2:
##                            x1=x1-Lx
##                        else:
##                            x2=x2-Lx
##                    if abs(y1-y2)>Lx/2:
##                        if y1>y2:
##                            y1=y1-Lx
##                        else:
##                            y2=y2-Lx
##                    if abs(z1-z2)>Lx/2:
##                        if z1>z2:
##                            z1=z1-Lx
##                        else:
##                            z2=z2-Lx
##                            
##                    d = distance(x1,x2,y1,y2,z1,z2)
##                    if d<=2.0:
##                        atoms = sorted([atom1,atom2])
##                        indicies = sorted([index1,index2])
##                        index_atom = indicies + atoms + [d]
##                        if index_atom not in CQ_bonds_across_BC:
##                            CQ_bonds_across_BC.append(index_atom)               

        #determine C-Si bonds (next to each other in chain); cutoff 2.3A
        elif (atom1==3 and atom2==1) or (atom1==1 and atom2==3):
##            if abs(index1-index2)==1:
            d = dist(dx,dy,dz)
            if d<=2.3:
                atoms = sorted([atom1,atom2])
                indicies = sorted([index1,index2])
                index_atom = indicies + atoms + [d]
                if indicies not in bond_index:
                    bond_index.append(indicies)
                    CSi_bonds.append(index_atom)

            #account for PBCs
            if dx > Lx-2.3:
                d = dist(dx+Lx,dy,dz)
                if d <= 2.3:
                    atoms = sorted([atom1,atom2])
                    indicies = sorted([index1,index2])
                    index_atom = indicies + atoms + [d]
                    if indicies not in bond_index:
                        bond_index.append(indicies)
                        CSi_bonds_across_BC.append(index_atom)

                d = dist(dx-Lx,dy,dz)
                if d <= 2.3:                       
                    atoms = sorted([atom1,atom2])
                    indicies = sorted([index1,index2])
                    index_atom = indicies + atoms + [d]
                    if indicies not in bond_index:
                        bond_index.append(indicies)
                        CSi_bonds_across_BC.append(index_atom)

            if dy > Ly-2.3:
                d = dist(dx,dy+Ly,dz)
                if d <= 2.3:
                    atoms = sorted([atom1,atom2])
                    indicies = sorted([index1,index2])
                    index_atom = indicies + atoms + [d]
                    if indicies not in bond_index:
                        bond_index.append(indicies)
                        CSi_bonds_across_BC.append(index_atom)

                d = dist(dx,dy-Ly,dz)
                if d <= 2.3:
                    atoms = sorted([atom1,atom2])
                    indicies = sorted([index1,index2])
                    index_atom = indicies + atoms + [d]
                    if indicies not in bond_index:
                        bond_index.append(indicies)
                        CSi_bonds_across_BC.append(index_atom)

            if dz > Lz-2.3:
                d = dist(dx,dy,dz+Lz)
                if d <= 2.3:
                    atoms = sorted([atom1,atom2])
                    indicies = sorted([index1,index2])
                    index_atom = indicies + atoms + [d]
                    if indicies not in bond_index:
                        bond_index.append(indicies)
                        CSi_bonds_across_BC.append(index_atom)

                d = dist(dx,dy,dz-Lz)
                if d <= 2.3:
                    atoms = sorted([atom1,atom2])
                    indicies = sorted([index1,index2])
                    index_atom = indicies + atoms + [d]
                    if indicies not in bond_index:
                        bond_index.append(indicies)
                        CSi_bonds_across_BC.append(index_atom)

               
##                d = distance(x1,x2,y1,y2,z1,z2)
##                if d<=2.3:
##                    atoms = sorted([atom1,atom2])
##                    indicies = sorted([index1,index2])
##                    index_atom = indicies + atoms + [d]
##                    if index_atom not in CSi_bonds:
##                        CSi_bonds.append(index_atom)
##                        
##                elif d>Lx/2:
##                    if abs(x1-x2)>Lx/2:
##                        if x1>x2:
##                            x1=x1-Lx
##                        else:
##                            x2=x2-Lx
##                    if abs(y1-y2)>Lx/2:
##                        if y1>y2:
##                            y1=y1-Lx
##                        else:
##                            y2=y2-Lx
##                    if abs(z1-z2)>Lx/2:
##                        if z1>z2:
##                            z1=z1-Lx
##                        else:
##                            z2=z2-Lx
##                            
##                    d = distance(x1,x2,y1,y2,z1,z2)
##                    if d<=2.3:
##                        atoms = sorted([atom1,atom2])
##                        indicies = sorted([index1,index2])
##                        index_atom = indicies + atoms + [d]
##                        if index_atom not in CSi_bonds_across_BC:
##                            CSi_bonds_across_BC.append(index_atom) 


        #determine Si-O bonds; cutoff at 2.0A; as a slight overestimate, we can say
        # that if the difference in the x-coords is <= 2.0 (we assume that y and z
        # coord are 0) then we will calculate the distance 
        elif (atom1==1 and atom2==2) or (atom1==2 and atom2==1):
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


        #determine Q-Q bonds; cuttoff at 2.0A; take difference in x-coord <=2.0
        elif (atom1==5 and atom2==5):
            d = dist(dx,dy,dz)
            if d<=2.0:
                atoms = sorted([atom1,atom2])
                indicies = sorted([index1,index2])
                index_atom = indicies + atoms + [d]
                if indicies not in bond_index:
                    bond_index.append(indicies)
                    QQ_bonds.append(index_atom)

            #account for PBCs
            if dx > Lx-2.0:
                d = dist(dx+Lx,dy,dz)
                if d <= 2.0:
                    atoms = sorted([atom1,atom2])
                    indicies = sorted([index1,index2])
                    index_atom = indicies + atoms + [d]
                    if indicies not in bond_index:
                        bond_index.append(indicies)
                        QQ_bonds_across_BC.append(index_atom)

                d = dist(dx-Lx,dy,dz)
                if d <= 2.0:                       
                    atoms = sorted([atom1,atom2])
                    indicies = sorted([index1,index2])
                    index_atom = indicies + atoms + [d]
                    if indicies not in bond_index:
                        bond_index.append(indicies)
                        QQ_bonds_across_BC.append(index_atom)

            if dy > Ly-2.0:
                d = dist(dx,dy+Ly,dz)
                if d <= 2.0:
                    atoms = sorted([atom1,atom2])
                    indicies = sorted([index1,index2])
                    index_atom = indicies + atoms + [d]
                    if indicies not in bond_index:
                        bond_index.append(indicies)
                        QQ_bonds_across_BC.append(index_atom)

                d = dist(dx,dy-Ly,dz)
                if d <= 2.0:
                    atoms = sorted([atom1,atom2])
                    indicies = sorted([index1,index2])
                    index_atom = indicies + atoms + [d]
                    if indicies not in bond_index:
                        bond_index.append(indicies)
                        QQ_bonds_across_BC.append(index_atom)

            if dz > Lz-2.0:
                d = dist(dx,dy,dz+Lz)
                if d <= 2.0:
                    atoms = sorted([atom1,atom2])
                    indicies = sorted([index1,index2])
                    index_atom = indicies + atoms + [d]
                    if indicies not in bond_index:
                        bond_index.append(indicies)
                        QQ_bonds_across_BC.append(index_atom)

                d = dist(dx,dy,dz-Lz)
                if d <= 2.0:
                    atoms = sorted([atom1,atom2])
                    indicies = sorted([index1,index2])
                    index_atom = indicies + atoms + [d]
                    if indicies not in bond_index:
                        bond_index.append(indicies)
                        QQ_bonds_across_BC.append(index_atom)            


        #determine Zr-O bonds; cutoff at 2.3A (Note: clustering affect will
        # change the amount of calculated bonds - some atoms within cutoff that
        # arnt actually bonded - source of error)
        elif (atom1==2 and atom2==4) or (atom1==4 and atom2==2):
            d = dist(dx,dy,dz)
            if d<=2.3:
                atoms = sorted([atom1,atom2])
                indicies = sorted([index1,index2])
                index_atom = indicies + atoms + [d]
                if indicies not in bond_index:
                    bond_index.append(indicies)
                    ZrO_bonds.append(index_atom)

            #account for PBCs
            if dx > Lx-2.3:
                d = dist(dx+Lx,dy,dz)
                if d<=2.3:
                    atoms = sorted([atom1,atom2])
                    indicies = sorted([index1,index2])
                    index_atom = indicies + atoms + [d]
                    if indicies not in bond_index:
                        bond_index.append(indicies)
                        ZrO_bonds_across_BC.append(index_atom)

                d = dist(dx-Lx,dy,dz)
                if d<=2.3:                       
                    atoms = sorted([atom1,atom2])
                    indicies = sorted([index1,index2])
                    index_atom = indicies + atoms + [d]
                    if indicies not in bond_index:
                        bond_index.append(indicies)
                        ZrO_bonds_across_BC.append(index_atom)

            if dy > Ly-2.3:
                d = dist(dx,dy+Ly,dz)
                if d<=2.3:
                    atoms = sorted([atom1,atom2])
                    indicies = sorted([index1,index2])
                    index_atom = indicies + atoms + [d]
                    if indicies not in bond_index:
                        bond_index.append(indicies)
                        ZrO_bonds_across_BC.append(index_atom)

                d = dist(dx,dy-Ly,dz)
                if d<=2.3:
                    atoms = sorted([atom1,atom2])
                    indicies = sorted([index1,index2])
                    index_atom = indicies + atoms + [d]
                    if indicies not in bond_index:
                        bond_index.append(indicies)
                        ZrO_bonds_across_BC.append(index_atom)

            if dz > Lz-2.3:
                d = dist(dx,dy,dz+Lz)
                if d<=2.3:
                    atoms = sorted([atom1,atom2])
                    indicies = sorted([index1,index2])
                    index_atom = indicies + atoms + [d]
                    if indicies not in bond_index:
                        bond_index.append(indicies)
                        ZrO_bonds_across_BC.append(index_atom)

                d = dist(dx,dy,dz-Lz)
                if d<=2.3:
                    atoms = sorted([atom1,atom2])
                    indicies = sorted([index1,index2])
                    index_atom = indicies + atoms + [d]
                    if indicies not in bond_index:
                        bond_index.append(indicies)
                        ZrO_bonds_across_BC.append(index_atom)
                        
                    
bonds_data = CC_bonds+CC_bonds_across_BC+CQ_bonds+CQ_bonds_across_BC+\
             CSi_bonds+CSi_bonds_across_BC+SiO_bonds+SiO_bonds_across_BC+\
             QQ_bonds+QQ_bonds_across_BC+ZrO_bonds+ZrO_bonds_across_BC

bonds_data_no_BC = CC_bonds+CQ_bonds+CSi_bonds+SiO_bonds+QQ_bonds+ZrO_bonds

N_CC_bonds = len(CC_bonds)
N_CC_bonds_across_BC = len(CC_bonds_across_BC)
N_CC_tot = N_CC_bonds+N_CC_bonds_across_BC

N_CQ_bonds = len(CQ_bonds)
N_CQ_bonds_across_BC = len(CQ_bonds_across_BC)
N_CQ_tot = N_CQ_bonds+N_CQ_bonds_across_BC

N_CSi_bonds = len(CSi_bonds)
N_CSi_bonds_across_BC = len(CSi_bonds_across_BC)
N_CSi_tot = N_CSi_bonds+N_CSi_bonds_across_BC

N_SiO_bonds = len(SiO_bonds)
N_SiO_bonds_across_BC = len(SiO_bonds_across_BC)
N_SiO_tot = N_SiO_bonds+N_SiO_bonds_across_BC

N_QQ_bonds = len(QQ_bonds)
N_QQ_bonds_across_BC = len(QQ_bonds_across_BC)
N_QQ_tot = N_QQ_bonds+N_QQ_bonds_across_BC

N_ZrO_bonds = len(ZrO_bonds)
N_ZrO_bonds_across_BC = len(ZrO_bonds_across_BC)
N_ZrO_tot = N_ZrO_bonds+N_ZrO_bonds_across_BC

#for clustering algorithm, the indicies of the Zr-O bonds and Si-O bonds are needed
N_CC_CQ_CSi = N_CC_tot + N_CQ_tot + N_CSi_tot
SiO_index = [str(N_CC_CQ_CSi), str(N_CC_CQ_CSi+N_SiO_tot)]

N_CC_CQ_CSi_SiO_QQ = N_CC_tot + N_CQ_tot + N_CSi_tot + N_SiO_tot + N_QQ_tot
ZrO_index = [str(N_CC_CQ_CSi_SiO_QQ), str(N_CC_CQ_CSi_SiO_QQ+N_ZrO_tot)]

#write out bond data
for bond in bonds_data:
    bond[0]=str(bond[0])
    bond[1]=str(bond[1])
    bond[2]=str(bond[2])
    bond[3]=str(bond[3])
    bond[4]=str(bond[4])

all_bond_data = [SiO_index]+[ZrO_index]+bonds_data

dataFile = open("bonds.txt", 'w')
for eachitem in all_bond_data:
    dataFile.write("\t".join(eachitem)+'\n')

dataFile.close()

print 'Bond data written to file!'

#change bond data back to int/float
for bond in bonds_data:
    bond[0]=int(bond[0])
    bond[1]=int(bond[1])
    bond[2]=int(float(bond[2]))
    bond[3]=int(float(bond[3]))
    bond[4]=float(bond[4])

                
#find number of types of bonds and transpose data to columns

CC_bonds_tot = CC_bonds+CC_bonds_across_BC
CC_bonds_tot = [list(x) for x in zip(*CC_bonds_tot)]
CC_bonds = [list(x) for x in zip(*CC_bonds)]

CQ_bonds_tot = CQ_bonds+CQ_bonds_across_BC
CQ_bonds_tot = [list(x) for x in zip(*CQ_bonds_tot)]
CQ_bonds = [list(x) for x in zip(*CQ_bonds)]

CSi_bonds_tot = CSi_bonds+CSi_bonds_across_BC
CSi_bonds_tot = [list(x) for x in zip(*CSi_bonds_tot)]
CSi_bonds = [list(x) for x in zip(*CSi_bonds)]

SiO_bonds_tot = SiO_bonds+SiO_bonds_across_BC
SiO_bonds_tot = [list(x) for x in zip(*SiO_bonds_tot)]
SiO_bonds = [list(x) for x in zip(*SiO_bonds)]

QQ_bonds_tot = QQ_bonds+QQ_bonds_across_BC
QQ_bonds_tot = [list(x) for x in zip(*QQ_bonds_tot)]
QQ_bonds = [list(x) for x in zip(*QQ_bonds)]

ZrO_bonds_tot = ZrO_bonds+ZrO_bonds_across_BC
ZrO_bonds_tot = [list(x) for x in zip(*ZrO_bonds_tot)]
ZrO_bonds = [list(x) for x in zip(*ZrO_bonds)]

print 'bonds found!'

#find the number of non-bonded O, non-bridging O (only connected to 1 ion), and connected O
print ''
print 'Finding number of non-bonded O, non-bridging O, and bridging O'

N_bonds = len(bonds_data)
N_bonds_no_BC = len(bonds_data_no_BC)
bonds_data = [list(x) for x in zip(*bonds_data)]
bonds_data_no_BC = [list(x) for x in zip(*bonds_data_no_BC)]

nonbonded_O = 0
for index in O_atom_index:
    if index not in bonds_data[0] and index not in bonds_data[1]:
        nonbonded_O=nonbonded_O+1

nonbridging_O=0
nonbridging_O_Si=0
nonbridging_O_Zr=0

bridging_O=0
bridging_O_Si=0
bridging_O_Zr=0

O_coord_3=0
O_coord_3_Si=0
O_coord_3_Zr=0

O_coord_4=0
O_coord_4_Si=0
O_coord_4_Zr=0

O_coord_5=0
O_coord_5_Si=0
O_coord_5_Zr=0

SiOZr_bridging=0
SiOZr_bridging_coord_3=0

index1=bonds_data[0]
index2=bonds_data[1]
atom1=bonds_data[2]
atom2=bonds_data[3]
for i in range(len(O_atom_index)):
    k=0
    l=0
    m=0
    for j in range(len(index1)):
        if O_atom_index[i]==index1[j] or O_atom_index[i]==index2[j]:
            m=m+1
            if atom1[j]==2: #O before Zr from sorted()
                k=k+1
            if atom2[j]==2: #Si before O from sorted()
                l=l+1

    if m==1:
        nonbridging_O=nonbridging_O+1
    elif m==2:
        bridging_O=bridging_O+1
    elif m==3:
        O_coord_3=O_coord_3+1
    elif m==4:
        O_coord_4=O_coord_4+1
    elif m==5:
        O_coord_5=O_coord_5+1
        
    if l==1 and k==0:
        nonbridging_O_Si=nonbridging_O_Si+1
    elif l==2 and k==0:
        bridging_O_Si=bridging_O_Si+1
    elif l==3 and k==0:
        O_coord_3_Si=O_coord_3_Si+1
    elif l==4 and k==0:
        O_coord_4_Si=O_coord_4_Si+1
    elif l==5 and k==0:
        O_coord_5_Si=O_coord_5_Si+1
        
    elif k==1 and l==0:
        nonbridging_O_Zr=nonbridging_O_Zr+1
    elif k==2 and l==0:
        bridging_O_Zr=bridging_O_Zr+1
    elif k==3 and l==0:
        O_coord_3_Zr=O_coord_3_Zr+1
    elif k==4 and l==0:
        O_coord_4_Zr=O_coord_4_Zr+1
    elif k==5 and l==0:
        O_coord_5_Zr=O_coord_5_Zr+1

    elif l==1 and k==1:
        SiOZr_bridging=SiOZr_bridging+1
    elif l==1 and k==2 or l==2 and k==1:
        SiOZr_bridging_coord_3=SiOZr_bridging_coord_3+1
        

percent_nonbonded_O = float(nonbonded_O)/float(N_O)
percent_nonbridging_O = float(nonbridging_O)/float(N_O)
percent_bridging_O = float(bridging_O)/float(N_O)

percent_interface = float(SiOZr_bridging+SiOZr_bridging_coord_3)/float(N_O)

percent_silane_to_ZrO = float(SiOZr_bridging+SiOZr_bridging_coord_3)/float(N_silane)


#Find percent of oxygen atoms that formed bonds
#find unique atoms in Si-O and Zr-O bond list; Note: for SiO index[0]=Si, index[1]=O
# and for ZrO, index[0]=O index[1]=Zr

O_atoms_in_SiO = SiO_bonds_tot[1]
O_atoms_in_SiO = unique_atoms(O_atoms_in_SiO)

##O_atoms_in_ZrO = ZrO_bonds_tot[0]
##O_atoms_in_ZrO = unique_atoms(O_atoms_in_ZrO)
##
###find intersection of SiO and ZrO lists    
##O_atoms_at_interface = lists_overlap(O_atoms_in_SiO,O_atoms_in_ZrO)
##
###remove intersection from SiO and ZrO lists
##O_atoms_in_SiO = remove_overlap(O_atoms_at_interface,O_atoms_in_SiO)
##O_atoms_in_ZrO = remove_overlap(O_atoms_at_interface,O_atoms_in_ZrO)

N_O_atoms_in_SiO = float(len(O_atoms_in_SiO))
##N_O_atoms_in_ZrO = float(len(O_atoms_in_ZrO))

percent_SiO = N_O_atoms_in_SiO/float(N_O)
##percent_ZrO = N_O_atoms_in_ZrO/float(N_O)

#find percent of silane molecules that bonded via Q-Q bond
Q_atoms_in_QQbonds = []
for i in range(N_QQ_tot):
    index1 = QQ_bonds_tot[0][i]
    index2 = QQ_bonds_tot[1][i]
    if index1 not in Q_atoms_in_QQbonds:
        Q_atoms_in_QQbonds.append(index1)
    if index2 not in Q_atoms_in_QQbonds:
        Q_atoms_in_QQbonds.append(index2)

N_Q_atoms_in_QQbonds = float(len(Q_atoms_in_QQbonds))
percent_silane_via_QQ = N_Q_atoms_in_QQbonds/float(N_silane)

#find average size of ZrO cluster




print 'Oxygen atoms statistics finished!'

################################################################################
#create bond bins and then calculate the most probable bonds

print ''
print 'Caclulating bond bins and most probable bond lengths...'

CC_hist_data = get_histogram_data(CC_bonds_tot[4],N_bins_bonds)
CQ_hist_data = get_histogram_data(CQ_bonds_tot[4],N_bins_bonds)
CSi_hist_data = get_histogram_data(CSi_bonds_tot[4],N_bins_bonds)
SiO_hist_data = get_histogram_data(SiO_bonds_tot[4],N_bins_bonds)
QQ_hist_data = get_histogram_data(QQ_bonds_tot[4],N_bins_bonds)
##ZrO_hist_data = get_histogram_data(ZrO_bonds_tot[4],N_bins_bonds)

high_prob_CC = CC_hist_data[2]
high_prob_CQ = CQ_hist_data[2]
high_prob_CSi = CSi_hist_data[2]
high_prob_SiO = SiO_hist_data[2]
high_prob_QQ = QQ_hist_data[2]
##high_prob_ZrO = ZrO_hist_data[2]


print 'Most probable bond lengths obtained!'

################################################################################
#plot bond RDFs

print ''
print 'Plotting bond distribution functions...'

plt.figure(8)
plt.plot(CC_hist_data[0],CC_hist_data[1],'k')
plt.ylabel('g(r)')
plt.xlabel('Distance, r (A)')
plt.title('C-C RDF')
plt.savefig('CC_bonds.png')

plt.figure(9)
plt.plot(CQ_hist_data[0],CQ_hist_data[1],'k')
plt.ylabel('g(r)')
plt.xlabel('Distance, r (A)')
plt.title('C-Q RDF')
plt.savefig('CQ_bonds.png')

plt.figure(10)
plt.plot(CSi_hist_data[0],CSi_hist_data[1],'g')
plt.ylabel('g(r)')
plt.xlabel('Distance, r (A)')
plt.title('C-Si RDF')
plt.savefig('CSi_bonds.png')

plt.figure(11)
plt.plot(SiO_hist_data[0],SiO_hist_data[1],'b')
plt.ylabel('g(r)')
plt.xlabel('Distance, r (A)')
plt.title('Si-O RDF')
plt.savefig('SiO_bonds.png')

plt.figure(12)
plt.plot(QQ_hist_data[0],QQ_hist_data[1],'m')
plt.ylabel('g(r)')
plt.xlabel('Distance, r (A)')
plt.title('Q-Q RDF')
plt.savefig('QQ_bonds.png')

##plt.figure(13)
##plt.plot(ZrO_hist_data[0],ZrO_hist_data[1],'r')
##plt.ylabel('g(r)')
##plt.xlabel('Distance, r (A)')
##plt.title('Zr-O RDF')
##plt.savefig('ZrO_bonds.png')


print 'RDFs obtained!'

################################################################################        
#from bond data, represent bonds as vectors; then for vectors that share a
# common atom, calculate the angle between them

print  ''
print 'Calculating bond vectors...'
        
vectors = []
for i in range(len(bonds_data_no_BC[0])):
    index1=bonds_data_no_BC[0][i]
    index2=bonds_data_no_BC[1][i]

    atom1=bonds_data_no_BC[2][i]
    atom2=bonds_data_no_BC[3][i]
    
    x1=x_coord[index1]
    y1=y_coord[index1]
    z1=z_coord[index1]

    x2=x_coord[index2]
    y2=y_coord[index2]
    z2=z_coord[index2]

    v=vector(x1,x2,y1,y2,z1,z2)

    v_index = v+[index1,index2,atom1,atom2]

    vectors.append(v_index)    
                
N_vectors = len(vectors)
N_v_index = len(zip(*vectors))

vectors = [list(x) for x in zip(*vectors)]

print 'Bond vectors obtained!'

print N_vectors
print N_bonds_no_BC
print N_bonds

################################################################################  
#calculate the angles; vectors must share a common atom to calculate angle

print ''
print 'Calculating the angles... '

CCC_angles = []
CQQ_angles = []
CCQ_angles = []
CCSi_angles = []
CSiO_angles = []
SiOZr_angles = []
OZrO_angles = []
SiOSi_angles = []
OSiO_angles = []
ZrOZr_angles = []
for i in range(N_vectors):
    x1 = vectors[0][i]
    y1 = vectors[1][i]
    z1 = vectors[2][i]
    index1=vectors[3][i]
    index2=vectors[4][i]
    atom1=vectors[5][i]
    atom2=vectors[6][i]
    for j in range(i+1,N_vectors):
        x2 = vectors[0][j]
        y2 = vectors[1][j]
        z2 = vectors[2][j]
        index3=vectors[3][j]
        index4=vectors[4][j]
        atom3=vectors[5][j]
        atom4=vectors[6][j]
        
        if index1==index3 or index1==index4 or index2==index3 or index2==index4:
            atoms_involved = unique_atoms([index1,index2,index3,index4])
            if len(atoms_involved)==3:
                theta=angle(x1,x2,y1,y2,z1,z2)
                if theta < 100:
                    theta = 180.0-theta
                
                #determine the center atom by finding the repeat atom
                nonrepeat_list = find_repeat(index1,index2,index3,index4,\
                                             atom1,atom2,atom3,atom4)

                center_atom = nonrepeat_list[0]
                Latom = nonrepeat_list[1]
                Ratom = nonrepeat_list[2]

                if center_atom == 3:
                    #C-C-C angles
                    if Latom==3 and Ratom==3:
                        CCC_angles.append(theta)

                    #C-C-Q angles
                    elif (Latom==3 or Latom==5) and (Ratom==5 or Ratom==3):
                        CCQ_angles.append(theta)
                        
                    #C-C-Si angles
                    elif (Latom==3 or Latom==1) and (Ratom==1 or Ratom==3):
                        CCSi_angles.append(theta)

                #C-Q-Q angles
                elif center_atom == 5:
                    if (Latom==3 or Latom==5) and (Ratom==5 or Ratom==3):
                        CQQ_angles.append(theta)

                elif center_atom == 1:
                    #C-Si-O angles
                    if (Latom==3 or Latom==2) and (Ratom==3 or Ratom==2):
                        CSiO_angles.append(theta)

                    #O-Si-O angles
                    if (Latom==2 and Ratom==2):
                        OSiO_angles.append(theta)

                elif center_atom == 2:
##                    #Si-O-Zr angles
##                    if (Latom==1 or Latom==4) and (Ratom==4 or Ratom==1):
##                        SiOZr_angles.append(theta)

                    #Si-O-Si angles
                    if (Latom==1 and Ratom==1):
                        SiOSi_angles.append(theta)

##                    #Zr-O-Zr angles
##                    if (Latom==4 and Ratom==4):
##                        ZrOZr_angles.append(theta)
##    
##                #O-Zr-O angles
##                elif center_atom == 4:
##                    if Latom==2 and Ratom==2:
##                        OZrO_angles.append(theta)


N_CCC_angles = len(CCC_angles)
N_CQQ_angles = len(CQQ_angles)
N_CCQ_angles = len(CCQ_angles)
N_CCSi_angles = len(CCSi_angles)
N_CSiO_angles = len(CSiO_angles)
##N_SiOZr_angles = len(SiOZr_angles)
##N_OZrO_angles = len(OZrO_angles)
N_SiOSi_angles = len(SiOSi_angles)
N_OSiO_angles = len(OSiO_angles)
##N_ZrOZr_angles = len(ZrOZr_angles)

print 'Angles found!'

################################################################################
#Create bond distribution functions and plot them

print ''
print 'Cacluating the angle distributions and most probable angles...'

CCC_hist_data = get_histogram_data(CCC_angles,N_bins_angles)
CQQ_hist_data = get_histogram_data(CQQ_angles,N_bins_angles)
CCQ_hist_data = get_histogram_data(CCQ_angles,N_bins_angles)
CCSi_hist_data = get_histogram_data(CCSi_angles,N_bins_angles)
CSiO_hist_data = get_histogram_data(CSiO_angles,N_bins_angles)
##SiOZr_hist_data = get_histogram_data(SiOZr_angles,N_bins_angles)
##OZrO_hist_data = get_histogram_data(OZrO_angles,N_bins_angles)
SiOSi_hist_data = get_histogram_data(SiOSi_angles,N_bins_angles)
OSiO_hist_data = get_histogram_data(OSiO_angles,N_bins_angles)
##ZrOZr_hist_data = get_histogram_data(ZrOZr_angles,N_bins_angles)

high_prob_CCC = CCC_hist_data[2]
high_prob_CQQ = CQQ_hist_data[2]
high_prob_CCQ = CCQ_hist_data[2]
high_prob_CCSi = CCSi_hist_data[2]
high_prob_CSiO = CSiO_hist_data[2]
##high_prob_SiOZr = SiOZr_hist_data[2]
##high_prob_OZrO = OZrO_hist_data[2]
high_prob_SiOSi = SiOSi_hist_data[2]
high_prob_OSiO = OSiO_hist_data[2]
##high_prob_ZrOZr = ZrOZr_hist_data[2]

print 'Angle distributions and most probable angles finished!'

print ''
print 'Plotting angle distributions...'

plt.figure(1)
plt.plot(CCC_hist_data[0],CCC_hist_data[1],'k')
plt.ylabel('g($\\theta$)')
plt.xlabel('Angle, $\\theta$ (degrees)')
plt.title('C-C-C Angle Distribution')
plt.savefig('CCC_angles.png')

plt.figure(2)
plt.plot(CQQ_hist_data[0],CQQ_hist_data[1],'k')
plt.ylabel('g($\\theta$)')
plt.xlabel('Angle, $\\theta$ (degrees)')
plt.title('C-Q-Q Angle Distribution')
plt.savefig('CQQ_angles.png')

plt.figure(3)
plt.plot(CCQ_hist_data[0],CCQ_hist_data[1],'k')
plt.ylabel('g($\\theta$)')
plt.xlabel('Angle, $\\theta$ (degrees)')
plt.title('C-C-Q Angle Distribution')
plt.savefig('CCQ_angles.png')

plt.figure(4)
plt.plot(CCSi_hist_data[0],CCSi_hist_data[1],'g')
plt.ylabel('g($\\theta$)')
plt.xlabel('Angle, $\\theta$ (degrees)')
plt.title('C-Si-O Angle Distribution')
plt.savefig('CCSi_angles.png')

plt.figure(5)
plt.plot(CSiO_hist_data[0],CSiO_hist_data[1],'m')
plt.ylabel('g($\\theta$)')
plt.xlabel('Angle, $\\theta$ (degrees)')
plt.title('C-Si-O Angle Distribution')
plt.savefig('CSiO_angles.png')

##plt.figure(6)
##plt.plot(SiOZr_hist_data[0],SiOZr_hist_data[1],'b')
##plt.ylabel('g($\\theta$)')
##plt.xlabel('Angle, $\\theta$ (degrees)')
##plt.title('Si-O-Zr Angle Distribution')
##plt.savefig('SiOZr_angles.png')
##
##plt.figure(7)
##plt.plot(OZrO_hist_data[0],OZrO_hist_data[1],'r')
##plt.ylabel('$\\theta$')
##plt.xlabel('Angle, $\\theta$ (degrees)')
##plt.title('O-Zr-O Angle Distribution')
##plt.savefig('OZrO_angles.png')

plt.figure(14)
plt.plot(SiOSi_hist_data[0],SiOSi_hist_data[1],'g')
plt.ylabel('g($\\theta$)')
plt.xlabel('Angle, $\\theta$ (degrees)')
plt.title('Si-O-Si Angle Distribution')
plt.savefig('SiOSi_angles.png')

plt.figure(15)
plt.plot(OSiO_hist_data[0],OSiO_hist_data[1],'g')
plt.ylabel('g($\\theta$)')
plt.xlabel('Angle, $\\theta$ (degrees)')
plt.title('O-Si-O Angle Distribution')
plt.savefig('OSiO_angles.png')

##plt.figure(16)
##plt.plot(ZrOZr_hist_data[0],ZrOZr_hist_data[1],'b')
##plt.ylabel('$\\theta$')
##plt.xlabel('Angle, $\\theta$ (degrees)')
##plt.title('Zr-O-Zr Angle Distribution')
##plt.savefig('ZrOZr_angles.png')

print 'Plots finished!'

###############################################################################
#Create summary of results and then write to text file

print ''
print 'Creating summary file...'

atoms_in_sim = [['Number of atoms in the simulation: '+str(N_atoms)]]

bond_results = [['Bond Results:'],\
                [''],\
                ['Total number of bonds: '+str(N_bonds)],\
                [''],\
                ['Number of C-C bonds: '+str(N_CC_tot)],\
                ['Most probable C-C bond length: '+str(high_prob_CC)],\
                ['Number of C-C bonds across BCs: '+str(N_CC_bonds_across_BC)],\
                [''],\
                ['Number of C-Q bonds: '+str(N_CQ_tot)],\
                ['Most probable C-Q bond length: '+str(high_prob_CQ)],\
                ['Number of C-Q bonds across BCs: '+str(N_CQ_bonds_across_BC)],\
                [''],\
                ['Number of C-Si bonds: '+str(N_CSi_tot)],\
                ['Most probable C-Si bond length: '+str(high_prob_CSi)],\
                ['Number of C-Si bonds across BCs: '+str(N_CSi_bonds_across_BC)],\
                [''],\
                ['Number of Si-O bonds: '+str(N_SiO_tot)],\
                ['Most probable Si-O bond length: '+str(high_prob_SiO)],\
                ['Number of Si-O bonds across BCs: '+str(N_SiO_bonds_across_BC)],\
                [''],\
                ['Number of Q-Q bonds: '+str(N_QQ_tot)],\
                ['Most probable Q-Q bond length: '+str(high_prob_QQ)],\
                ['Number of Q-Q bonds across BCs: '+str(N_QQ_bonds_across_BC)]]
##                [''],\
##                ['Number of Zr-O bonds: '+str(N_ZrO_tot)],\
##                ['Most probable Zr-O bond length: '+str(high_prob_ZrO)],\
##                ['Number of Zr-O bonds across BCs: '+str(N_ZrO_bonds_across_BC)]]


O_bonding_results = [['Oxygen atoms reacting statistics: '],\
                     ['Percent of Oxygen atoms that formed Si-O bonds: '+\
                      str(percent_SiO)],\
##                     ['Percent of Oxygen atoms that formed Zr-O bonds: '+\
##                      str(percent_ZrO)],\
##                     ['Percent of Oxygen atoms that formed Si-O-Zr (interface) bonds: '\
##                      +str(percent_interface)],\
                     [''],\
                     ['Percent of non-bonded O: '+str(percent_nonbonded_O)],\
                     ['Number of non-bonded O: '+str(nonbonded_O)],\
                     [''],\
                     ['Percent of non-bridging O: '+str(percent_nonbridging_O)],\
                     ['Number of non-bridging O: '+str(nonbridging_O)],\
                     ['Number of non-bridging O with Si: '+str(nonbridging_O_Si)],\
##                     ['Number of non-bridging O with Zr: '+str(nonbridging_O_Zr)],\
                     [''],\
                     ['Percent of bridging O: '+str(percent_bridging_O)],\
                     ['Number of bridging O: '+str(bridging_O)],\
                     ['Number of bridging O with Si: '+str(bridging_O_Si)],\
##                     ['Number of bridging O with Zr: '+str(bridging_O_Zr)],\
##                     ['Number of bridging O with Si and Zr: '+str(SiOZr_bridging)],\
##                     ['Number of bridging O with Si and Zr (coord 3): '+str(SiOZr_bridging_coord_3)],\
                     [''],\
                     ['Number of O with coordination 3: '+str(O_coord_3)],\
                     ['Number of O with coordination 3 with Si: '+str(O_coord_3_Si)],\
##                     ['Number of O with coordination 3 with Zr: '+str(O_coord_3_Zr)],\
                     [''],\
                     ['Number of O with coordination 4: '+str(O_coord_4)],\
                     ['Number of O with coordination 4 with Si: '+str(O_coord_4_Si)],\
##                     ['Number of O with coordination 4 with Zr: '+str(O_coord_4_Zr)],\
                     [''],\
                     ['Number of O with coordination 5: '+str(O_coord_5)],\
                     ['Number of O with coordination 5 with Si: '+str(O_coord_5_Si)]]
##                     ['Number of O with coordination 5 with Zr: '+str(O_coord_5_Zr)]]

                     
Silane_bonding_results = [['Silane molecule reacting statistics: '],\
                          ['Percent of silane molecules that bonded via Q-Q bonds: '\
                           +str(percent_silane_via_QQ)]]
##                          ['Percent of silane molecules that bonded to ZrO cluster '\
##                           +'via Si-O bonds: '+str(percent_silane_to_ZrO)]]



angle_results = [['Angle Results: '],\
                 [''],\
                 ['Number of C-C-C angles: '+str(N_CCC_angles)],\
                 ['Most probable C-C-C angle: '+str(high_prob_CCC)],\
                 [''],\
                 ['Number of C-Q-Q angles: '+str(N_CQQ_angles)],\
                 ['Most probable C-Q-Q angle: '+str(high_prob_CQQ)],\
                 [''],\
                 ['Number of C-C-Q angles: '+str(N_CCQ_angles)],\
                 ['Most probable C-C-Q angle: '+str(high_prob_CCQ)],\
                 [''],\
                 ['Number of C-C-Si angles: '+str(N_CCSi_angles)],\
                 ['Most probable C-C-Si angle: '+str(high_prob_CCSi)],\
                 [''],\
                 ['Number of C-Si-O angles: '+str(N_CSiO_angles)],\
                 ['Most probable C-Si-O angle: '+str(high_prob_CSiO)],\
##                 [''],\
##                 ['Number of Si-O-Zr angles: '+str(N_SiOZr_angles)],\
##                 ['Most probable Si-O-Zr angle: '+str(high_prob_SiOZr)],\
##                 [''],\
##                 ['Number of O-Zr-O angles: '+str(N_OZrO_angles)],\
##                 ['Most probable O-Zr-O angle: '+str(high_prob_OZrO)],\
                 [''],\
                 ['Number of Si-O-Si angles: '+str(N_SiOSi_angles)],\
                 ['Most probable Si-O-Si angle: '+str(high_prob_SiOSi)],\
                 [''],\
                 ['Number of O-Si-O angles: '+str(N_OSiO_angles)],\
                 ['Most probable O-Si-O angle: '+str(high_prob_OSiO)]]
##                 [''],\
##                 ['Number of Zr-O-Zr angles: '+str(N_ZrOZr_angles)],\
##                 ['Most probable Zr-O-Zr angles: '+str(high_prob_ZrOZr)]]

                 
results = atoms_in_sim+[['']]+[['']]+bond_results+[['']]+O_bonding_results+\
          [['']]+Silane_bonding_results+[['']]+[['']]+angle_results

################################################################################  
#write to a text file
#the .join() method takes an array, i, and concantenates all the elements together
# with a space " " between each element.  Then a newline "\n" is added to make sure
# your output is broken up into separate lines

dataFile = open(filename[:-4]+"_results.txt", 'w')
for eachitem in results:
    dataFile.write("\t".join(eachitem)+'\n')

dataFile.close()


print "All done!" 
          




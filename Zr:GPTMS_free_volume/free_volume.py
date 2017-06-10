#Purpose: find volume of Zr/O clusters and the free volume of organic/inorganic 

from pylab import *
import pylab
from scipy import *
from numpy import *
import numpy as np
import math
from scipy.spatial import ConvexHull
from scipy.spatial import Delaunay

################################################################################

def distance(dx,dy,dz):
    d = (dx**2+dy**2+dz**2)**0.5
    return d

def region(orig_cluster,ZrO_atoms_guage,cutoff,Lx,Ly,Lz):
    new_cluster = []
    for atom1 in orig_cluster:
        for atom2 in ZrO_atoms_guage:
            dx = abs(atom1[2]-atom2[2])
            dy = abs(atom1[3]-atom2[3])
            dz = abs(atom1[4]-atom2[4])
            d = distance(dx,dy,dz)
            
            if d <= cutoff:
                new_cluster.append(atom2)
                ZrO_atoms_guage.remove(atom2)
                
##            #account for PBCs
##            if dx > Lx-cutoff:
##                d = distance(dx+Lx,dy,dz)
##                if d <= cutoff:
##                    new_cluster.append(atom2)
##                    ZrO_atoms_guage.remove(atom2)
##
##                d = distance(dx-Lx,dy,dz)
##                if d <= cutoff:
##                    new_cluster.append(atom2)
##                    ZrO_atoms_guage.remove(atom2)
##
##            if dy > Ly-cutoff:
##                d = distance(dx,dy+Ly,dz)
##                if d <= cutoff:
##                    new_cluster.append(atom2)
##                    ZrO_atoms_guage.remove(atom2)
##
##                d = distance(dx,dy-Ly,dz)
##                if d <= cutoff:
##                    new_cluster.append(atom2)
##                    ZrO_atoms_guage.remove(atom2)
##
##            if dz > Lz-cutoff:
##                d = distance(dx,dy,dz+Lz)
##                if d <= cutoff:
##                    new_cluster.append(atom2)
##                    ZrO_atoms_guage.remove(atom2)
##
##                d = distance(dx,dy,dz-Lz)
##                if d <= cutoff:
##                    new_cluster.append(atom2)
##                    ZrO_atoms_guage.remove(atom2)

    new_cluster = new_cluster+orig_cluster
    return [new_cluster,ZrO_atoms_guage]

def average(cluster_sizes):
    avg = sum(cluster_sizes)/len(cluster_sizes)
    return avg
        

################################################################################
#import bond data [index1, index2, atom1, atom2, bond_length] and .xyz file to
# correlate indicies with positions

xyz_file = raw_input('What is the .xyz filename? ')
##bond_file = raw_input('What is the bond filename? ')
cutoff = input('What cutoff do you want for the region around each atom? ')

##atom_data = []
##with open('ZrGPTMS.66000.xyz') as inputfile:
##    for line in inputfile:
##        atom_data.append(line.strip().split())

atom_data = []
with open(xyz_file) as inputfile:
    for line in inputfile:
        atom_data.append(line.strip().split())

N_atoms = int(atom_data[0][0])
atom_data=atom_data[2:]

bonds_data = []
with open('bonds.txt') as inputfile:
    for line in inputfile:
        bonds_data.append(line.strip().split())

#find number of Si,O,Zr, and C atoms; input volume of each atom (treat each as sphere)
N_O=0
N_Si=0
N_Zr=0
N_C=0
for i in range(N_atoms):
    if atom_data[i][0]=='1':
        N_Si=N_Si+1
    elif atom_data[i][0]=='2':
        N_O=N_O+1
    elif atom_data[i][0]=='3':
        N_C=N_C+1
    elif atom_data[i][0]=='4':
        N_Zr=N_Zr+1
    elif atom_data[i][0]=='5':
        N_C=N_C+1

V_C=1.417542
V_Si=0.4946877
V_O=6.284367
V_Zr=1.9982289

#get the bond indicies for Si-O and Zr-O
SiO_index = bonds_data[0]
SiO_start = int(SiO_index[0])
SiO_end = int(SiO_index[1])

ZrO_index = bonds_data[1]
ZrO_start = int(ZrO_index[0])
ZrO_end = int(ZrO_index[1])

bonds_data = bonds_data[2:]

atom_index = []
for i in range(len(atom_data)):
    atom_index.append(i)

#transpose to columns
atom_data = [list(x) for x in zip(*atom_data)]
atom_data = [atom_index]+atom_data

for i in range(len(atom_data)):
    for j in range(N_atoms):
        atom_data[i][j] = float(atom_data[i][j])

#find dimensions of box
x_coord = atom_data[2]
y_coord = atom_data[3]
z_coord = atom_data[4]

x_min = min(x_coord)
x_max = max(x_coord)
Lx = x_max - x_min

y_min = min(y_coord)
y_max = max(y_coord)
Ly = y_max - y_min

z_min = min(z_coord)
z_max = max(z_coord)
Lz = z_max - z_min

volume = Lx*Ly*Lz

#transpose back to rows 
atom_data = [list(x) for x in zip(*atom_data)]


################################################################################
#begin with list of all Zr and O (in Zr-O bond) atoms in a list; for each atom in
# list, place it in a cluster (based off spatial location within 8A sphere); remove
# atom from original list; continue until no atoms left

#change str to ints/floats
for atom in atom_data:
    atom[0]=int(atom[0])
    atom[1]=int(atom[1])
    atom[2]=float(atom[2]) #x-coord 
    atom[3]=float(atom[3]) #y-coord 
    atom[4]=float(atom[4]) #z-coord

for bond in bonds_data:
    bond[0]=int(bond[0])
    bond[1]=int(bond[1])
    bond[2]=int(float(bond[2]))
    bond[3]=int(float(bond[3]))
    bond[4]=float(bond[4])

    
#only consider Zr-O and Si-O bonds and Zr/O/Si atoms
print ''
print 'Finding Zr, Si, and O atoms from the bond list...'

SiO_atoms = []
O_atoms_in_SiO = []
for i in range(SiO_start,SiO_end):
    if bonds_data[i][0] not in SiO_atoms: #dont overcount Si atoms
        SiO_atoms.append(bonds_data[i][0])
        O_atoms_in_SiO.append(bonds_data[i][0])

    if bonds_data[i][1] not in SiO_atoms: #dont overcount O atoms
        SiO_atoms.append(bonds_data[i][1])




print 'SiO atoms, ',len(SiO_atoms) 
print 'O in SiO, ',len(O_atoms_in_SiO)

ZrO_atoms = []
O_atoms_in_ZrO = []
for i in range(ZrO_start,ZrO_end):  
    if atom_data[bonds_data[i][0]] not in ZrO_atoms: #dont overcount O atoms
        ZrO_atoms.append(atom_data[bonds_data[i][0]])
        O_atoms_in_ZrO.append(atom_data[bonds_data[i][0]])        

    if atom_data[bonds_data[i][1]] not in ZrO_atoms: #dont overcount Zr atoms
        ZrO_atoms.append(atom_data[bonds_data[i][1]])
        
print 'ZrO atoms, ',len(ZrO_atoms)
print 'O in ZrO, ',len(O_atoms_in_ZrO)


ZrO_atoms_guage = ZrO_atoms
print 'Zr and O atoms obtained!'
print 'Number of Zr and O atoms: ',len(ZrO_atoms)

#compare each atom in ZrO_atoms to others and find ones close to it (8A)
# must consider periodic boundary conditions (in region function)
print ''
print 'Finding clusters...'

all_clusters = []
t=0
m=0
while len(ZrO_atoms_guage)>0:
    t=t+1
    print t
    
    l=True
    
    cluster = [ZrO_atoms_guage[0]]
    ZrO_atoms_guage.remove(ZrO_atoms_guage[0])
    
    while l==True:
        size_orig_cluster = len(cluster)
        get_region = region(cluster,ZrO_atoms_guage,cutoff,Lx,Ly,Lz)
        cluster = get_region[0]
        ZrO_atoms_guage = get_region[1]
        size_new_cluster = len(cluster)
        
        if size_new_cluster==size_orig_cluster:
            l=False
            all_clusters.append(cluster)

N_clusters = len(all_clusters)

print 'Found clusters!'
print 'Number of clusters: ',N_clusters

#################################################################################
#use the convex Hull algorithm to calculate the perimeter points of the clusters
# and then the Delaunay triangulation to find the volume of each cluster; turn off
# the PCBs for the clustering algorithm 
def convex_hull_volume(pts):
    ch = ConvexHull(pts)
    dt = Delaunay(pts[ch.vertices])
    tets = dt.points[dt.simplices]
    return np.sum(tetrahedron_volume(tets[:, 0], tets[:, 1],
                                     tets[:, 2], tets[:, 3]))


def tetrahedron_volume(a, b, c, d):
    return np.abs(np.einsum('ij,ij->i', a-d, np.cross(b-d, c-d))) / 6


print ''
print 'Finding the volume of each cluster...'

V_clusters=[]
cluster_atoms=[]
for i in range(N_clusters):
    print i
    NO=0
    NZr=0
    atom_positions=[]
    for j in range(len(all_clusters[i])):
        if all_clusters[i][j][1]==2:
            N0=NO+1
        elif all_clusters[i][j][1]==4:
            NZr=NZr+1
        atom_positions.append([all_clusters[i][j][2],all_clusters[i][j][3],
                               all_clusters[i][j][4]])
        
    cluster_atoms.append([NO,NZr])
     
    atom_positions=np.array(atom_positions)
    
    #to compute the volume, there must be 4 or more Hull points to perform
    # a Delaunay triangulation; for clusters with 3 or less atoms, merely append
    # the volume of the atoms

    if len(atom_positions) > 3:
        V_clusters.append(convex_hull_volume(atom_positions))
    else:
        V_O_atoms = NO*V_O
        V_Zr_atoms = NZr*V_Zr
        V_clusters.append(V_O_atoms+V_Zr_atoms)


#calculate free volume of organic and inorgamic components
V_inorganic = sum(V_clusters)

V_organic = volume - V_inorganic

#find free volume in inorganic clusters
V_free_inorganic=0
for i in range(N_clusters):
    V_free_inorganic = V_free_inorganic + (V_clusters[i]-cluster_atoms[i][0]*V_O-
                                           cluster_atoms[i][1]*V_Zr)


#find free volume in organic regions
V_free_organic = V_organic - N_C*V_C - N_Si*V_Si  #think about O in organic portion


###############################################################################
#Create summary of results and then write to text file
#find largest, smallest, and average cluster size
print ''
print 'Finding cluster statistics...'

cluster_sizes = []
for cluster in all_clusters:
    cluster_sizes.append(len(cluster))

cluster_sizes = sorted(cluster_sizes,reverse=True)

print 'Largest cluster: ',cluster_sizes[0]
print 'Second largest cluster: ',cluster_sizes[1]
print 'Average cluster: ',average(cluster_sizes)

print 'Free organic volume: ',V_free_organic
print 'Free inorganic volume: ',V_free_inorganic

cluster_sizes = map(str,cluster_sizes) 

num_clusters = [['Number of clusters: '+str(N_clusters)],\
                [''],\
                [''],\
                ['Free volume of organic network: '+str(V_free_organic)],\
                ['Free volume of inorganic network: '+str(V_free_inorganic)],\
                [''],\
                [''],\
                ['Cluster sizes: ']]

results = num_clusters+cluster_sizes


################################################################################  
#write to a text file
#the .join() method takes an array, i, and concantenates all the elements together
# with a space " " between each element.  Then a newline "\n" is added to make sure
# your output is broken up into separate lines

#cluster statistics 
dataFile = open("cluster_size_with_"+str(cutoff)+"A_cutoff.txt", 'w')
for eachitem in results:
    dataFile.write("".join(eachitem)+'\n')

dataFile.close()


print 'All done!'



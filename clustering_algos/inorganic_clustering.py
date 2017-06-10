#Purpose: find extent of Zr/O clusters

from pylab import *
from scipy import *
from numpy import *
import numpy as np
import math

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
    
##ZrO_bonds = []
##ZrO_atoms = []
##ZrO_index = []
##for bond in bonds_data:
##    if bond[3]==4:
##        ZrO_bonds.append(bond)
##        if bond[0] not in ZrO_index:
##            ZrO_index.append(bond[0])
##            l=True
##            for atom in atom_data:
##                if l==True:
##                    if atom[0]==bond[0]:
##                        l=False
##                        ZrO_atoms.append(atom)
##        if bond[1] not in ZrO_index:
##            ZrO_index.append(bond[1])
##            k=True
##            for atom in atom_data:
##                if k==True:
##                    if atom[0]==bond[1]:
##                        k==False
##                        ZrO_atoms.append(atom)
##
##print len(ZrO_atoms)

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
##print 'Thrid largest cluster: ',cluster_sizes[2]
##print 'Fourth largest cluster: ',cluster_sizes[3]
##print 'Fifth largest cluster: ',cluster_sizes[4]
##print 'Smallest cluster: ',min(cluster_sizes)
print 'Average cluster: ',average(cluster_sizes)

cluster_sizes = map(str,cluster_sizes) 

num_clusters = [['Number of clusters: '+str(N_clusters)],\
                [''],\
                [''],\
                ['Cluster sizes: ']]

results = num_clusters+cluster_sizes

###############################################################################
###To confirm algorithm, create VMD file to visualize clusters
##
##for cluster in all_clusters:
##    for atom in cluster:
##        atom[1]=str(atom[1])
##        atom[2]=str(atom[2])
##        atom[3]=str(atom[3])
##        atom[4]=str(atom[4])
##        atom.remove(atom[0])
##
##clusters_separate = []
##clusters_together = []
##for cluster in all_clusters:
##    clusters_separate = clusters_separate+[[str(len(cluster))]]+[['Atoms']]+cluster
##
##    clusters_together = clusters_together+cluster
##    
##    if len(cluster)==int(cluster_sizes[0]):
##        largest_cluster = cluster
##        
##    elif len(cluster)==int(cluster_sizes[1]):
##        second_largest_cluster = cluster
##
##    elif len(cluster)==int(cluster_sizes[2]):
##        third_largest_cluster = cluster       
##
##    elif len(cluster)==int(cluster_sizes[3]):
##        fourth_largest_cluster = cluster
##
##
##clusters_together = [[str(len(clusters_together))]]+[['Atoms']]+clusters_together
##
##clusters_xyz = clusters_separate+clusters_together
## 
##largest_cluster = [[str(len(largest_cluster))]]+[['Atoms']]+largest_cluster
##second_largest_cluster = [[str(len(second_largest_cluster))]]+[['Atoms']]+second_largest_cluster
##third_largest_cluster = [[str(len(third_largest_cluster))]]+[['Atoms']]+third_largest_cluster
##fourth_largest_cluster = [[str(len(fourth_largest_cluster))]]+[['Atoms']]+fourth_largest_cluster
##

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

###all clusters
##dataFile = open('clusters_'+str(cutoff)+'A_cutoff.xyz', 'w')
##for eachitem in clusters_xyz:
##    dataFile.write("\t".join(eachitem)+'\n')
##
##dataFile.close()
##
#largest cluster
##dataFile = open('largest_cluster_'+str(cutoff)+'A_cutoff.xyz', 'w')
##for eachitem in largest_cluster:
##    dataFile.write("\t".join(eachitem)+'\n')
##
##dataFile.close()
##
###second largest cluster
##dataFile = open('2nd_largest_cluster_'+str(cutoff)+'A_cutoff.xyz', 'w')
##for eachitem in second_largest_cluster:
##    dataFile.write("\t".join(eachitem)+'\n')
##
##dataFile.close()

###third largest cluster 
##dataFile = open('3rd_largest_cluster_'+str(cutoff)+'A_cutoff.xyz', 'w')
##for eachitem in third_largest_cluster:
##    dataFile.write("\t".join(eachitem)+'\n')
##
##dataFile.close()

###fourth largest cluster
##dataFile = open('4th_largest_cluster_'+str(cutoff)+'A_cutoff.xyz', 'w')
##for eachitem in fourth_largest_cluster:
##    dataFile.write("\t".join(eachitem)+'\n')
##
##dataFile.close()


print 'All done!'



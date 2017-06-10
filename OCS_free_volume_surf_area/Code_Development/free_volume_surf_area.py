#Purpose: find volume of progen clusters and the free volume of OCS 

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

def load_data(filename):
    data = []
    f = open(filename)
    #first 2 lines are not relvant 
    f.readline()
    f.readline()
    for line in f:
        line = line.strip().split()
        data.append([int(line[0]),float(line[1]),float(line[2]),float(line[3])])
    f.close()
    return data

def region(orig_cluster,porogen_guage,cutoff):
    new_cluster = []
    for atom1 in orig_cluster:
        for atom2 in porogen_guage:
            dx = abs(atom1[1]-atom2[1])
            dy = abs(atom1[2]-atom2[2])
            dz = abs(atom1[3]-atom2[3])
            d = distance(dx,dy,dz)
            
            if d <= cutoff:
                new_cluster.append(atom2)
                porogen_guage.remove(atom2)

            #account for PBCs
            if dx > Lx-cutoff:
                d = distance(dx+Lx,dy,dz)
                if d <= cutoff:
                    new_cluster.append(atom2)
                    porogen_guage.remove(atom2)

                d = distance(dx-Lx,dy,dz)
                if d <= cutoff:
                    new_cluster.append(atom2)
                    porogen_guage.remove(atom2)

            if dy > Ly-cutoff:
                d = distance(dx,dy+Ly,dz)
                if d <= cutoff:
                    new_cluster.append(atom2)
                    porogen_guage.remove(atom2)

                d = distance(dx,dy-Ly,dz)
                if d <= cutoff:
                    new_cluster.append(atom2)
                    porogen_guage.remove(atom2)

            if dz > Lz-cutoff:
                d = distance(dx,dy,dz+Lz)
                if d <= cutoff:
                    new_cluster.append(atom2)
                    porogen_guage.remove(atom2)

                d = distance(dx,dy,dz-Lz)
                if d <= cutoff:
                    new_cluster.append(atom2)
                    porogen_guage.remove(atom2)                
                
    new_cluster = new_cluster+orig_cluster
    return new_cluster,porogen_guage

def average(cluster_sizes):
    avg = sum(cluster_sizes)/len(cluster_sizes)
    return avg
        

################################################################################
#import porogen .xyz file to correlate indicies with positions

porogen_file = raw_input('What is the porogen.xyz filename? ')
atom_file = raw_input('What is the .xyz file name of the total simulation? ')
cutoff = input('What cutoff do you want for the region around each atom? ')
##volume = input('What is the V0 you want to normalize to? ')

porogen_data = load_data(porogen_file)
N_P = len(porogen_data)

atom_data = load_data(atom_file)
N_atoms = len(atom_data)

#find number of Si,O, and C atoms; input volume of each atom (treat each as sphere)
N_O=0
N_Si=0
N_C=0
for i in range(N_atoms):
    if atom_data[i][0]==1:
        N_Si += 1
    elif atom_data[i][0]==2:
        N_O += 1
    elif atom_data[i][0]==3:
        N_C += 1

V_C=1.417542
V_Si=0.4946877
V_O=6.284367
V_P=1.417542

N_C = N_C - N_P

#find dimensions of box; transpose array
atom_data = [list(x) for x in zip(*atom_data)]

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

volume = Lx*Ly*Lz
##print volume

##volume = 232638.632649

atom_data = [list(x) for x in zip(*atom_data)]

################################################################################
#begin with list of all porogen atoms; for each atom in list, place it in a
# cluster (based off spatial location within cutoff sphere); remove
# atom from original list; continue until no atoms left

porogen_guage = porogen_data

print '\nFinding clusters...'

all_clusters = []
t=0
while len(porogen_guage)>0:
    t += 1
    print t
    
    cluster = [porogen_guage[0]]
    porogen_guage.remove(porogen_guage[0])
    
    while True:
        size_orig_cluster = len(cluster)
        get_region = region(cluster,porogen_guage,cutoff)
        cluster, porogen_guage = get_region
        size_new_cluster = len(cluster)
        
        if size_new_cluster==size_orig_cluster:
            all_clusters.append(cluster)
            break

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

def chull_V_Asurf(pts):
    ch = ConvexHull(pts)
    dt = Delaunay(pts[ch.verticies])
    tets = dt.points[dt.simplices]
    volume = np.sum(tetrahedron_volume(tets[:, 0], tets[:, 1],
                                     tets[:, 2], tets[:, 3]))
    #find surface area
    


print ''
print 'Finding the volume of each cluster...'

V_clusters=[]
cluster_atoms=[]
for i in range(N_clusters):
    print i
    atom_positions=[]
    N = len(all_clusters[i])
    for j in range(N):
        atom_positions.append([all_clusters[i][j][1],all_clusters[i][j][2],
                               all_clusters[i][j][3]])
        
    cluster_atoms.append(N)
    atom_positions=np.array(atom_positions)
    
    #to compute the volume, there must be 4 or more Hull points to perform
    # a Delaunay triangulation; for clusters with 3 or less atoms, merely append
    # the volume of the atoms

    if N > 3:
        V_clusters.append(convex_hull_volume(atom_positions))
    else:
        V_clusters.append(V_P*N)


#calculate free volume of the porogen (pores)
V_pores = sum(V_clusters)
V_OCS = volume - V_pores

#find free volume in OCS 
V_free_OCS = V_OCS - N_C*V_C - N_Si*V_Si - N_O*V_O

#############################################################################
#Create summary of results and then write to text file
#find largest, smallest, and average cluster size
print ''
print 'Finding cluster statistics...'

cluster_sizes = []
for cluster in all_clusters:
    cluster_sizes.append(len(cluster))

cluster_sizes = sorted(cluster_sizes,reverse=True)

print 'Largest cluster: ',cluster_sizes[0]
##print 'Second largest cluster: ',cluster_sizes[1]
##print 'Average cluster: ',average(cluster_sizes)

print 'Volume of pores per unit V0: ',V_pores/volume
print 'Free volume in OCS per unit V0: ',V_free_OCS/volume

cluster_sizes = map(str,cluster_sizes) 

num_clusters = [['Number of clusters: '+str(N_clusters)],\
                [''],\
                [''],\
                ['Volume of pores per unit V0: '+str(V_pores/volume)],\
                ['Free volume in OCS per unit V0: '+str(V_free_OCS/volume)],\
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
dataFile = open("cluster_size_with_"+str(cutoff)+"A_cutoff"+\
                porogen_file[-5]+".txt", 'w')
for eachitem in results:
    dataFile.write("".join(eachitem)+'\n')

dataFile.close()


print 'All done!'



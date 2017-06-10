"""
Purpose: find the free volume of OCS with no porogen; this algorithm 1) triangulates the
space using the delaunay triangulation algorithm; 2) finds the center of the circumscribed
sphere around the delaunay tetrahedrons which corresponds to a voronoi vertex (mathematical
duality to voronoi tesselation); 3) each voronoi vertex is the center of an interstitial sphere
of radius Di = D0 - d, where D0 is the diameter of the circumscribed delanay sphere and d is
the diameter of the atom; 4) overlapping voronoi vertices are calculated via a = Ri + Rj - Rij
where Ri and Rj are radii of interstitial spheres and Rij is the distance between them 
"""
from circumcenter import circumcenter
import copy
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy
import numpy.linalg
import pylab
import scipy
from scipy.spatial import ConvexHull
from scipy.spatial import Delaunay
from scipy.spatial import Voronoi
import time

################################################################################
def LoadData(filename):
    data = []
    f = open(filename)
    #first 2 lines are not relvant 
    N_atoms = int(f.readline())
    f.readline()
    for line in f:
        line = line.strip().split()
        data.append([int(line[0]),float(line[1]),float(line[2]),float(line[3])])
    f.close()
    data = numpy.array(data)

    # make arrays of atom types and coordinates
    atom_type = numpy.array(data[:,0],dtype=numpy.int32)
    atom_pos = numpy.delete(data,0,1)

    # create dictionary with {atom_index : atom_type}
    # 1 = Si, 2 = O, 3 = C
    AtomType = {}
    for i in range(N_atoms):
        AtomType[i] = atom_type[i]
        
    return N_atoms,atom_pos,atom_type,AtomType

def CellVolume(atom_type,atom_pos):
    x_coord = atom_pos[:,0]
    y_coord = atom_pos[:,1]
    z_coord = atom_pos[:,2]
    xmin = min(x_coord)
    xmax = max(x_coord)
    ymin = min(y_coord)
    ymax = max(y_coord)
    zmin = min(z_coord)
    zmax = max(z_coord)
    Lx = xmax - xmin
    Ly = ymax - ymin
    Lz = zmax - zmin
    volume = Lx*Ly*Lz

    print Lx,Ly,Lz

    N_O=0
    N_Si=0
    N_C=0
    for i in range(N_atoms):
        if atom_type[i]==1:
            N_Si += 1
        elif atom_type[i]==2:
            N_O += 1
        elif atom_type[i]==3:
            N_C += 1
            
    R_Si = 2.1 #van der waals radius in Angstroms 
    R_C = 1.7
    R_O = 1.52
    atom_volume = SphereVol(R_Si)*N_Si + SphereVol(R_O)*N_O +\
                  SphereVol(R_C)*N_C
    return volume, atom_volume, (xmin,xmax,ymin,ymax,zmin,zmax)

def SphereVol(R):
    return 4/3*math.pi*R**3

def SphereOverlap(R1,R2,d):
    if d == 0:
        return 0
    else:
        return (numpy.pi/12/d)*((R1 + R2 - d)**2)*\
           (d**2 + 2*d*(R1 + R2) - 3*((R1 - R2)**2))

def Distance(dx,dy,dz):
    d = (dx**2+dy**2+dz**2)**0.5
    return d

def VoronoiNetwork(tets_coords,boundary):
    N_tets = len(tets_coords)
    vor_vertices = []
    vor_R = []
    for i in range(N_tets):
        # cirumcenter() computes the center and radius of a circumscribing
        # sphere around a given simplex
        vertex,R = circumcenter(tets_coords[i])
        vor_vertices.append(vertex)
        vor_R.append(R)

    # remove atoms that go beyond the boundaries
    vor_vertices_BCs = []
    vor_R_BCs = []
    for i in range(len(vor_vertices)):
        if (boundary[0] < vor_vertices[i][0] < boundary[1]) and \
           (boundary[2] < vor_vertices[i][1] < boundary[3]) and \
           (boundary[4] < vor_vertices[i][2] < boundary[5]):
            vor_vertices_BCs.append(vor_vertices[i])
            vor_R_BCs.append(vor_R[i])
            
    return numpy.array(vor_vertices_BCs),numpy.array(vor_R_BCs)

def VoronoiRadii(tets_coords):
    N_tets = len(tets_coords)
    vor_R = []
    for i in range(N_tets):
        vor_R.append(circumcenter(tets_coords[i])[1])
    return numpy.array(vor_R)

def InterstitialRadii(tets_idx,vor_R,AtomType):
    R_Si = 2.1 #van der waals radius in Angstroms 
    R_C = 1.7
    R_O = 1.52
    Ri = []
    for i in range(len(vor_R)):
        #consider radius of largest atom
        idx_types = [AtomType[tets_idx[i][0]],AtomType[tets_idx[i][1]],\
                     AtomType[tets_idx[i][2]],AtomType[tets_idx[i][3]]]
        if 1 in idx_types:
            Ri.append(vor_R[i] - R_Si)
        elif 3 in idx_types:
            Ri.append(vor_R[i] - R_C)
        else:
            Ri.append(vor_R[i] - R_O)
    return numpy.array(Ri)

##def Sphereicity(vor_R,Ri):
##    sphereicity = numpy.zeros(len(Ri))
##    for i in range(len(Ri)):
##        sphereicity[i] = (Ri[i]**3)/(vor_R[i]**3)
##    return sphereicity

def region(orig_cluster,vertex_guage,vertex_idx,Ri):
##    new_cluster = set([])
##    vertex_guage2 = copy.deepcopy(vertex_guage)
    new_cluster = []
    for vertex1 in orig_cluster:
##        R1 = Ri[vertex_idx[(vertex1[0],vertex1[1],vertex1[2])]]
        R1 = Ri[vertex1]
        for vertex2 in vertex_guage:
##            R2 = Ri[vertex_idx[(vertex2[0],vertex2[1],vertex2[2])]]
            R2 = Ri[vertex2]
            
            dx = abs(vertex_idx[vertex1][0]-vertex_idx[vertex2][0])
            if dx > R1 + R2: continue
            dy = abs(vertex_idx[vertex1][1]-vertex_idx[vertex2][1])
            if dy > R1 + R2: continue
            dz = abs(vertex_idx[vertex1][2]-vertex_idx[vertex2][2])
            if dz > R1 + R2: continue
            d = Distance(dx,dy,dz)

            Dij = R1 + R2 - d
            if Dij >= 0:
                new_cluster.append(vertex2)
                vertex_guage.remove(vertex2)               
##                if vertex2 in vertex_guage:
##                    vertex_guage.remove(vertex2)
                
    new_cluster = new_cluster + orig_cluster
##    new_cluster = new_cluster | orig_cluster
    return new_cluster,vertex_guage
        
def ComputeHoles(Ri,vor_vertices):
##    vertex_idx = {}
##    idx = 0
##    for vertex in vor_vertices:
##        vertex_idx[(vertex[0],vertex[1],vertex[2])] = idx
##        idx += 1

    idx_vertex = {}
    idx = 0
    for vertex in vor_vertices:
        idx_vertex[idx] = vertex
        idx += 1
    
##    vertex_guage = set(range(len(vor_vertices)))
    vertex_guage = range(len(vor_vertices))
    vor_vertices = vor_vertices.tolist()
    all_clusters = []
    t=0
    while len(vertex_guage):
        t += 1
        print t

##        vertex0 = next(iter(vertex_guage))
##        cluster = set([vertex0])
##        vertex_guage.remove(vertex0)
        cluster = [vertex_guage[0]]
        vertex_guage.remove(vertex_guage[0])

        t1 = 0
        while True:
            t1 += 1
            print 'inner loop iteration %d' %t1
            print 'current cluster size %d' %len(cluster)
            size_orig_cluster = len(cluster)
            get_region = region(cluster,vertex_guage,idx_vertex,Ri)
            cluster, vertex_guage = get_region
            size_new_cluster = len(cluster)
            
            if size_new_cluster == size_orig_cluster:
                all_clusters.append(cluster)
                break
    return all_clusters        


def ComputeDist(Ri,vor_vertices):
    N_tets = 0
    N_octs = 0
    for i in range(len(Ri)):
        if 0.85 < Ri[i] < 0.87:
            N_octs += 1
        elif 0.45 < Ri[i] < 0.48:
            N_tets += 1

    vor_dist = numpy.zeros((len(vor_vertices),len(vor_vertices)))
    vor_overlap = numpy.zeros((len(vor_vertices),len(vor_vertices)))
    n = 0
    t = 0
    Voids_row = []
    Voids = []
    void0 = set([0])
    free_volume = 0
    VO = 0
    for i in range(len(vor_vertices)):
        void = [i]
##        free_volume += SphereVol(Ri[i])
        for j in range(i+1,len(vor_vertices)):
            dx = vor_vertices[i][0] - vor_vertices[j][0]
            dy = vor_vertices[i][1] - vor_vertices[j][1]
            dz = vor_vertices[i][2] - vor_vertices[j][2]
            d = Distance(dx,dy,dz)
            
            vor_dist[i][j] = d

            R1 = Ri[i]
            R2 = Ri[j]
            if R1 + R2 >= d:
                vor_overlap[i][j] = 1
                n += 1
                void.append(j)
##                free_volume -= SphereOverlap(R1,R2,d)
                free_volume += SphereVol(R1) + SphereVol(R2) - SphereOverlap(R1,R2,d)
                VO += SphereOverlap(R1,R2,d)

            t += 1               
            if t % 1000000 == 0:
                print "Iteration number: %d" %t

        void = set(void)
        Voids_row.append(void)

        if void & void0:
            void0 = void0 | void

    print 'Volume of interstitial overlap: %.4f' %VO           

##    Voids = []
##    for void1 in Voids_row:
##        for void2 in Voids_row:
##            if void1 & void2:
##                void1 = void1 | void2
            
        
##    Voids = []
##    for row in vor_overlap:
##        void = [row]
##        for idx in row:
##            if idx:
##                void.append(idx)
##        Voids.append(void)
    
    return vor_dist,vor_overlap,n,Voids,void0,free_volume,N_octs,N_tets
            
def ComputeVoids(Ri,vor_vertices):
    clusters = []
    cluster = [Ri[0],vor_vertices[0]]
    vor_vertices.remove(vor_vertices[0])
    clusters.append(cluster)
    for i in range(len(vor_vertices)):
        for j in range(len(clusters)):
            print i
    


####    overlap = numpy.zeros((len(Ri),len(Ri)))
##    # sparse matrix in COO format 
##    col_idx = []
##    row_idx = []
####    val_d = []
##    val_V = []
##    holes = []
##    for i in range(len(Ri)-1):
##        R1 = Ri[i]
##        hole = set([i])
##        holes.append(hole)
####        # find which hole interstitial sphere is in
####        notinhole = 1
####        for hole in holes:
####            if i in hole:
####                current_hole = hole
####                notinhole = 0
####                break
####        if notinhole:
####            # first element of list is hole volume
####            current_hole = [4/3*numpy.pi*R1**3,i]
####            holes.append(current_hole)
##
##        # compute overlap between all other interstitials
##        for j in range(i+1,len(Ri)):
##            R2 = Ri[j]
##            dx = abs(vor_vertices[i][0] - vor_vertices[j][0])
##            dy = abs(vor_vertices[i][1] - vor_vertices[j][1])
##            dz = abs(vor_vertices[i][2] - vor_vertices[j][2])
##            d = Distance(dx,dy,dz)
##            
##            Dij = R1 + R2 - d
##            if Dij >= 0:
####                overlap[i][j] = Dij
##                row_idx.append(i)
##                col_idx.append(j)
####                val_d.append(d)
##                hole.add(j)
##                
##                V_overlap = (numpy.pi/12/d)*((R1 + R2 - d)**2)*\
##                            (d**2 + 2*d*(R1 + R2) - 3*((R1 - R2)**2))
##                val_V.append(V_overlap)
####                overlap[i][j] = V_overlap
####
####                # adjust hole volume
####                current_hole[0] += (4/3*numpy.pi*R2**3) - V_overlap
####                current_hole.append(j)
####
####    free_volume = 0
####    for i in range(len(holes)):
####        free_volume += holes[i][0]
##
##    # find the holes and free volume
##    # add all verticies that overlap with 0; check rows of all those that overlap and
##    # continue to check until all that overlap are in same hole; continue to add
##    # volumes and subtract overlaps
##
####    idx = 0
####    holes = []
####    hole = set([0])
####    holes.append(hole)
####    for i in range(len(row_idx)-1):
####        if row_idx[i] == row_idx[i+1]:
####            hole.add(col_idx[i])
####        else:
####            hole = set([row_idx[i+1]])
####            holes.append(hole)
##
####    holes = set([])
####    for i in range(len(row_idx)):
####        # search to see if row index is already in hole
####        for hole in holes:
####            if i in hole:
####                current_hole = hole
####            
####        for j in range(len(col_idx)):
####            overlap[i][j] 
####        
##
####    holes_guage = range(len(Ri))
####    Holes = []
####    idx = 0
####    while len(holes_guage):
####        hole1 = holes[idx]
####        idx += 1
####        for hole2 in Holes:
####            if hole1 & hole2:
####                hole1 = hole1 | hole2
####
####        Holes.append(hole1)
####        for i in hole1:
####            if i in holes_guage:
####                print i
####                holes_guage.remove(i)
####
####    for i in range(len(Ri)):
####        bag = []
####        for hole in holes:
####            if i in hole:
####                bag.append(hole)
####
####    for hole in holes:
####        while bag & hole:
####            bag = bag | hole
##
####    print holes_guage
####    Holes = []
####    for i in range(len(holes)-1):
####        hole = holes[i]
####        for hole1 in holes:
####            if hole & hole1:
####                hole = hole | hole1
####        
####        for j in range(i+1,len(holes)-1):
####            if hole & holes[j]:
####                hole = hole | holes[j]
####        Holes.append(hole)
##                
####    for hole1 in holes:
####        for hole2 in holes:
####            # if there is an intersection, take the union
####            if hole1 & hole2:
####                hole1 = hole1 | hole2
####                print hole1                
            

################################################################################
#import .xyz file to correlate indicies with positions
filename = raw_input('What is the .xyz file name of the total simulation? ')

t0 = time.time()

N_atoms, atom_pos, atom_type, AtomType = LoadData(filename)
volume, atom_volume, boundary = CellVolume(atom_type,atom_pos)

print 'Packing Fraction: %.4f' %(atom_volume/volume)

# Perfrom delaunay triangulation; Delanay() creates a triangulation object where
# the x,y,z points can be accessed by .point and the inidices can be accessed
# by .simplices
dt = Delaunay(atom_pos)
tets_coords = dt.points[dt.simplices]
tets_idx = dt.simplices

print 'N_atoms = %d' %N_atoms
print 'Number Delaunay tetrahedrons: %d' %len(tets_coords)
##print tets_coords
##print tets_idx

##fig = plt.figure()
##ax = fig.add_subplot(111, projection='3d')
####ax.scatter(tets_coords[:,0],tets_coords[:,1],tets_coords[:,2],c='b',marker='o')
##ax.scatter(atom_pos[:,0],atom_pos[:,1],atom_pos[:,2],c='r',marker='o')
####ax.triplot(tets[:,0],tets[:,1],tets[:,2],tets1,'go-')
##plt.show()

# Perform Voronoi tessellation; x,y,z vertices are accessed by .verticies;
# each vertex corresponds to the circumcenter of the Delaunay tetrahedrons
# in the same order (can confirm with VeronoiVerts function)
##vor = Voronoi(atom_pos)
##vor_vertices = vor.vertices
##vor_R = VoronoiRadii(tets_coords)
vor_vertices, vor_R = VoronoiNetwork(tets_coords,boundary)
print 'Number of Voronoi vertices: %d' %len(vor_vertices)

# find radii of interstitial spheres centered at each voronoi vertex (Di = D0-d)
Ri = InterstitialRadii(tets_idx,vor_R,AtomType)

V = 0
for R in Ri:
  V += SphereVol(R)

print 'Volume of box: %.4f' %volume
print 'Volume of atoms in box: %.4f' %atom_volume
print 'Free volume of interstitials w/ no overlap: %.4f' %V

# find degree of overlap between interstitial spheres and compute free volume
vor_dist,vor_overlap,n,Voids,void0,free_volume,N_octs,N_tets = ComputeDist(Ri,vor_vertices)
####print vor_dist,vor_overlap
print 'Number of overlapping Voronoi vertices: %d' %n
print 'Number of octahedral sites: %d' %N_octs
print 'Number of tetrahedral sites: %d' %N_tets
##print Voids
##print void0
##print len(void0)


##clusters = ComputeHoles(Ri,vor_vertices)
##print row_idx
##print col_idx
##print val_V
##print len(clusters)

##print "number of holes = %d " %len(holes)
##print free_volume
##print atom_volume
##print volume
print "percent free volume = %.5f" %(free_volume/atom_volume)



### find sphereicity (Vvoid / Vsphere)
##sphereicity = Sphereicity(vor_R,Ri)

t1 = time.time()
print 'Elapsed time: %.3f seconds' % (t1-t0)


##################################################################################  
### write Voronoi vertices to file for VMD
##vor_vertices_BCs = numpy.array(vor_vertices_BCs)    
##
##vor_vertices_BCs = numpy.append(numpy.ones([len(vor_vertices_BCs),1]),vor_vertices_BCs,1)
##numpy.savetxt("voidsVMD.xyz",vor_vertices_BCs,fmt='%.4f',delimiter='\t',newline='\n')


###write to a text file
##results = "\nFree vol at Timestep %s = %.5f" %(filename[-12:-4],perc_free_vol)
##
##dataFile = open("free_volume_"+filename[-5]+".txt", 'a')
##dataFile.write(results)
##dataFile.close()
##
##
##print 'All done!'


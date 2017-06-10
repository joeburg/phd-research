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

##    # compute distances between all atoms
##    AtomDist = numpy.zeros((len(atom_pos),len(atom_pos)))
##    t = 0
##    for i in range(len(atom_pos)):
##        for j in range(i+1,len(atom_pos)):
##            t += 1
##            if t % 1000000 == 0:
##                print "iteration: %d" %t
##            d = Distance(atom_pos[i][0]-atom_pos[j][0],
##                         atom_pos[i][1]-atom_pos[j][1],
##                         atom_pos[i][2]-atom_pos[j][2])
##            AtomDist[i][j]=d
    
    return N_atoms,atom_pos,AtomType

def AtomVolume(AtomType,atom_pos):
    N_O=0
    N_Si=0
    N_C=0
    for i in range(N_atoms):
        if AtomType[i]==1:
            N_Si += 1
        elif AtomType[i]==2:
            N_O += 1
        elif AtomType[i]==3:
            N_C += 1
    
    R_Si = 2.1 #van der waals radius in Angstroms 
    R_C = 1.7
    R_O = 1.52
    atom_volume = SphereVol(R_Si)*N_Si + SphereVol(R_O)*N_O +\
                  SphereVol(R_C)*N_C
    return atom_volume

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

def TetVolume(a, b, c, d):
    return numpy.abs(numpy.dot(a-d,numpy.cross(b-d, c-d))) / 6

def Angle(vec1,vec2):
    ab = numpy.dot(vec1,vec2)
    a = numpy.linalg.norm(vec1)
    b = numpy.linalg.norm(vec2)

##    print ab,a,b,ab/a/b
    return numpy.arccos(ab/a/b)

def PointH(vec1,vec2,vertexA,vertexB):
    """
    This function finds the H point, which is the height from vertex
    A(x1,y1,z1), so AH vector is normal to BCD plane, and H lies in the plane BCD.
    H = (x1+ka,y1+kb,z1+kc), where ax+by+cz-d=0 is the equation of the plane BCD.
    To ensure the correct plane orientation, choose point (B) in BCD plane and find
    vector BA and check if each spatial component has the same sign as the normal vector.
    """
    normal = numpy.cross(vec1,vec2)
    a = normal[0]
    b = normal[1]
    c = normal[2]
    
##    vecBA = vertexA - vertexB 
##    if numpy.array_equal(numpy.sign(normal),numpy.sign(vecBA)): 
##        a = normal[0]
##        b = normal[1]
##        c = normal[2]
##    else:
##        normal = -normal
##        a = normal[0]
##        b = normal[1]
##        c = normal[2]
        
    x1 = vertexA[0]
    y1 = vertexA[1]
    z1 = vertexA[2]

    d = -(vertexB[0]*normal[0]+vertexB[1]*normal[1]+vertexB[2]*normal[2])

    k = -(d + a*x1 + b*y1 + c*z1)/(a**2 + b**2 + c**2)

    H = numpy.array([x1+k*a,y1+k*b,z1+k*c])
##    print 'H = ', H
    vecAH = vertexA - H
##    vecHA = H - vertexA
    if numpy.array_equal(numpy.sign(normal),numpy.sign(vecAH)):
##        print 'height of tet: %.4f' %Distance(H[0]-vertexA[0],H[1]-vertexA[1],H[2]-vertexA[2])
##        print 'Same direction'
        if numpy.array_equal(H,vertexA):
            # if the atoms are all coplanar, take V~0
            return 0,H
        else:
            return 1,H
    else:
##        print '\nNot same direction'
##        normal = -normal
##        a = normal[0]
##        b = normal[1]
##        c = normal[2]
##
##        x1 = vertexA[0]
##        y1 = vertexA[1]
##        z1 = vertexA[2]
##
##        d = -(vertexB[0]*normal[0]+vertexB[1]*normal[1]+vertexB[2]*normal[2])
##        k = -(d + a*x1 + b*y1 + c*z1)/(a**2 + b**2 + c**2)
##        H = numpy.array([x1+k*a,y1+k*b,z1+k*c])
##
##        vecAH = vertexA - H
##        if numpy.array_equal(numpy.sign(normal),numpy.sign(vecAH)):
##            print 'Fixed direction!'
##        else:
##            print H,normal,vertexA,vertexB,d,k

        # if they are not in the same direction, switch orientation of plane
        T,H = PointH(vec2,vec1,vertexA,vertexB)
        return T, H
##        if numpy.array_equal(H,vertexA):
##            return 0
##        else:
##            return H

##    print 'height of tet: %.4f' %Distance(H[0]-vertexA[0],H[1]-vertexA[1],H[2]-vertexA[2])
##    return H

def FracAtomVol(R,vertex1,vertex2,vertex3,vertex4):
    """
    vertex1, vertex2, and vertex3 are coplanar; vertex4 sits above the
    plane
    """
    # using angle1 as origin, find angle between vectors in base plane 
    vec1 = vertex2 - vertex1
    vec2 = vertex3 - vertex1

##    print vec1,vec2
    alpha = Angle(vec1,vec2)

    # find angle in from base plane towards 'z' direcion 
    # must find vector in base plane 
    vec3 = vertex4 - vertex1

    T,pointH = PointH(vec1,vec2,vertex4,vertex1)

##    print T,pointH
    if T:
        vec4 = pointH - vertex1
        
        beta = Angle(vec3,vec4)

    ##    print alpha*180/numpy.pi,beta*180/numpy.pi,SphereVol(R)
    ##    print alpha,beta,SphereVol(R)

        # compute fraction of volume between sphere and tetrahedron
        V_frac = (0.5)*SphereVol(R)*(beta/numpy.pi)*(alpha/2/numpy.pi)

    ##    print 'Fraction of atom: %.10f\n' %(V_frac/SphereVol(R))
        return V_frac
    else:
        return 0

def ComputeOverlap(atom1,atom2,tet_idx,idx1,idx2,tets_idx,AtomType):
    R_Si = 2.1 #van der waals radius in Angstroms 
    R_C = 1.7
    R_O = 1.52
    
    V_overlap = 0
    # get coords and size of vdW radius for atom1
    x1 = atom1[0]
    y1 = atom1[1]
    z1 = atom1[2]
    if AtomType[tets_idx[tet_idx][idx1]] == 1:               
        R1 = R_Si
    elif AtomType[tets_idx[tet_idx][idx1]] == 2:               
        R1 = R_O
    elif AtomType[tets_idx[tet_idx][idx1]] == 3:               
        R1 = R_C

    # get coords and size of vdW radius for atom2
    x2 = atom2[0]
    y2 = atom2[1]
    z2 = atom2[2]
    if AtomType[tets_idx[tet_idx][idx2]] == 1:               
        R2 = R_Si
    elif AtomType[tets_idx[tet_idx][idx2]] == 2:               
        R2 = R_O
    elif AtomType[tets_idx[tet_idx][idx2]] == 3:               
        R2 = R_C

    d = Distance(x1-x2,y1-y2,z1-z2)
##    print 'Distance: %.4f' %d

    # determine overlap
    if R1 > R2:
        RLg = R1
        Rsm = R2
    else:
        RLg = R2
        Rsm = R1
        
    if 0 <= d <= RLg - Rsm:
        V_overlap += SphereVol(Rsm)
    elif RLg - Rsm < d <= RLg + Rsm:
        V_overlap += SphereOverlap(RLg,Rsm,d)
    else:
        V_overlap = 0

    return V_overlap

def FreeVolume(tets_coords,tets_idx,AtomType):
    R_Si = 2.1 #van der waals radius in Angstroms 
    R_C = 1.7
    R_O = 1.52

    V_free = 0
    V_Tets = 0
    N_Sm_holes = 0
    for i in range(len(tets_coords)):
        if i % 10000 == 0:
            print 'Iteration: %d' %i
        tet = tets_coords[i]
        V_tet = TetVolume(tet[0],tet[1],tet[2],tet[3])
        V_Tets += V_tet

        # compute overlap of tet vertices (atoms w/ vdW radii)
        V_overlap = 0
        N_overlap = 0
        for j in range(len(tet)):
            for k in range(j+1,len(tet)):
                V_overlap += ComputeOverlap(tet[j],tet[k],i,j,k,tets_idx,AtomType)
                if V_overlap:
                    N_overlap += 1

        # if 3 or more overlapping atoms then space is
        # too small for diffusing particles 
        if N_overlap >= 3:
            N_Sm_holes += 1
                    
        # compute volume of atoms that overlap with tetrahedron
        V_atoms = 0
        R = []
        for j in range(len(tet)):
            if AtomType[tets_idx[i][j]] == 1:
                R = R_Si
            elif AtomType[tets_idx[i][j]] == 2:
                R = R_O
            elif AtomType[tets_idx[i][j]] == 3:
                R = R_C

            V_atoms += FracAtomVol(R,tet[j%4],tet[(j+1)%4],tet[(j+2)%4],tet[(j+3)%4])
            
##        V_atoms += FracAtomVol(R[1],tet[3],tet[0],tet[1],tet[2])
##        V_atoms += FracAtomVol(R[2],tet[2],tet[3],tet[0],tet[1])
##        V_atoms += FracAtomVol(R[3],tet[1],tet[2],tet[3],tet[0])
    
        # add to free volume
##        print "V_tet, V_atom, V_overlap = %.4f  %.4f  %.4f" %(V_tet,V_atoms,V_overlap)

        V_free += V_tet - V_atoms + V_overlap
        
##        # only consider + volumes (- volumes means that the tetrahedron is
##        # too small for any free volume
##        N_bad_tets = 0
##        if V_tet - V_atoms + V_overlap > 0:
##            V_free += V_tet - V_atoms + V_overlap
##        else:
##            N_bad_tets += 1
            

    return V_free,V_Tets,N_Sm_holes

        

################################################################################
#import .xyz file to correlate indicies with positions
filename = raw_input('What is the .xyz file name of the total simulation? ')

t0 = time.time()

N_atoms,atom_pos,AtomType = LoadData(filename)
atom_volume = AtomVolume(AtomType,atom_pos)
dimension = 3

print '\nAtom volume: %.4f' %atom_volume
print 'N_atoms = %d' %N_atoms






# Perfrom delaunay triangulation; Delanay() creates a triangulation object where
# the x,y,z points can be accessed by .point and the inidices can be accessed
# by .simplices
dt = Delaunay(atom_pos,incremental=False,qhull_options="QJ")
tets_coords = dt.points[dt.simplices]
tets_idx = dt.simplices

print dt.coplanar
print 'Number of atoms not used: %d' %(len(dt.coplanar))

print 'Number Delaunay tetrahedrons: %d\n' %len(tets_coords)
##print tets_coords[:15]

##print type(AtomType)
##print tets_idx[:15]


##idx_used = set([])
##for tet in tets_idx:
##    for idx in tet:
##        if idx not in idx_used:
##            idx_used.add(idx)
##
##
##print len(idx_used)


##V_free,V_Tets,N_Sm_holes = FreeVolume(tets_coords,tets_idx,AtomType)
####print V_free
####print N_sm_holes
####print N_bad_tets
##
##print 'Tetrahedron volume: %.4f' %V_Tets
##print 'Vol_atom + Vol_free = %.4f' %(atom_volume+V_free)
##print 'Packing Fraction: %.4f' %(atom_volume/(atom_volume+V_free))
##print 'Free volume: %.4f' %(V_free/(atom_volume+V_free))

##tets_coords = tets_coords[:200]
##fig = plt.figure()
##ax = fig.add_subplot(111, projection='3d')
####ax.scatter(tets_coords[:,0],tets_coords[:,1],tets_coords[:,2],c='b',marker='o')
##ax.scatter(atom_pos[:,0],atom_pos[:,1],atom_pos[:,2],c='r',marker='o')
####ax.triplot(tets[:,0],tets[:,1],tets[:,2],tets1,'go-')
##plt.show()

##print free_volume
##print atom_volume
##print "percent free volume = %.5f" %(free_volume/atom_volume)

t1 = time.time()
print 'Elapsed time: %.3f seconds' % (t1-t0)


##################################################################################  
### write Voronoi vertices to file for VMD
##vor_vertices_BCs = numpy.array(vor_vertices_BCs)    
##
##vor_vertices_BCs = numpy.append(numpy.ones([len(vor_vertices_BCs),1]),vor_vertices_BCs,1)
##numpy.savetxt("voidsVMD.xyz",vor_vertices_BCs,fmt='%.4f',delimiter='\t',newline='\n')

### write tetrahedrons to VMD file
##tetFile = open("tets_VMD.xyz",'w')
##tetFile.write('%d\n'%len(tets_coords))
##tetFile.write('Tetrahedrals\n')
##for i in range(len(tets_coords)):
####    tetFile.write('4\n')
####    tetFile.write('Tetrahedral_'+str(i)+'\n')
##    for j in range(len(tets_coords[i])):        
##        if AtomType[tets_idx[i][j]] == 1:
##            tp = 'Si'
##        elif AtomType[tets_idx[i][j]] == 2:
##            tp = 'O'
##        elif AtomType[tets_idx[i][j]] == 3:
##            tp = 'C'
##
##        x = tets_coords[i][j][0]
##        y = tets_coords[i][j][1]
##        z = tets_coords[i][j][2]
##        
##        line = '%s  %10f %10f %10f\n' %(tp,x,y,z)
##
##        tetFile.write(line)
##tetFile.close()

###write to a text file
##results = "\nFree vol at Timestep %s = %.5f" %(filename[-12:-4],perc_free_vol)
##
##dataFile = open("free_volume_"+filename[-5]+".txt", 'a')
##dataFile.write(results)
##dataFile.close()
##
##
##print 'All done!'


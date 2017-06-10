from circumcenter import circumcenter
import copy
import math
import matplotlib.pyplot as plt
import numpy
import numpy.linalg
import pylab
import scipy
import time

#################################################################################################
""" Mathematical computations"""

def Distance(dx,dy,dz):
    d = (dx**2+dy**2+dz**2)**0.5
    return d

def trunc(f, n):
    '''Truncates/pads a float f to n decimal places without rounding'''
    slen = len('%.*f' % (n, f))
    return str(f)[:slen]

def AtomVolume(AtomType,N_atoms):
    N_O=0
    N_Si=0
    N_C=0
    for i in range(N_atoms):
        if AtomType[i]==2.1:
            N_Si += 1
        elif AtomType[i]==1.52:
            N_O += 1
        elif AtomType[i]==1.7:
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

def TetVolume(a, b, c, d):
    return numpy.abs(numpy.dot(a-d,numpy.cross(b-d, c-d))) / 6

def Angle(vec1,vec2):
    ab = numpy.dot(vec1,vec2)
    a = numpy.linalg.norm(vec1)
    b = numpy.linalg.norm(vec2)
    return numpy.arccos(ab/a/b)

#################################################################################################
""" Writing Data"""

def GetType(R):
    if 2.0 < R < 2.2:
        return 'Si'
    elif 1.4 < R < 1.55:
        return 'O'
    elif 1.6 < R < 1.8:
        return 'C'

def WriteLine(typeID,atom):
    f.write('%s  %.4f  %.4f  %.4f\n' %(typeID,atom[0],atom[1],atom[2])) 

#################################################################################################
""" Loading Data """

def LoadDataUnweighted(filenameNodes,filenameTets,filenameLAMMPS):
    # import the atom types from LAMMPS file
    atom_R = []
    f = open(filenameLAMMPS)
    f.readline()
    f.readline()
    for line in f:
        line = line.strip().split()
        if line[0] == '1':
            atom_R.append(2.1)
        elif line[0] == '2':
            atom_R.append(1.52)
        elif line[0] == '3':
            atom_R.append(1.7)
    f.close()
    
    # import the atom coordinates
    atom_ID = []
    atom_pos = []
    f = open(filenameNodes)
    #first 2 lines are not relvant 
    f.readline()
    for line in f:
        line = line.strip().split()
        if not line[0][0] == '#':
            atom_ID.append(int(line[0]))
            atom_pos.append([float(line[1]),float(line[2]),float(line[3])])
    f.close()
    
    N_atoms = len(atom_ID)

    # make arrays of atom types and coordinates
    atom_R = numpy.array(atom_R)
    atom_ID = numpy.array(atom_ID,dtype=numpy.int32)
    atom_pos = numpy.array(atom_pos)

    # create dictionary with {atom_index : atom_type} and {atom_idex : [x,y,z]}
    # 2.1 = Si, 1.52 = O, 1.7 = C
    AtomType = {}
    AtomPos = {}
    for i in range(N_atoms):
        AtomType[i] = atom_R[i]
        AtomPos[i] = atom_pos[i]

    # import the Delaunay triangulations
    Tets = {} # dictionary {tetID : [atom1,atom2,atom3,atom4]}
    f = open(filenameTets)
    f.readline()
    for line in f:
        line = line.strip().split()
        if not line[0][0] == '#':
            Tets[int(line[0])] = [int(line[1]),int(line[2]),int(line[3]),int(line[4])]
    f.close()
    
    return N_atoms,AtomPos,AtomType,Tets


def LoadDataWeighted(filenameNodes,filenameTets):
    # import the atom coordinates
    atom_R = []
    atom_ID = []
    atom_pos = []
    f = open(filenameNodes)
    #first 2 lines are not relvant 
    f.readline()
    for line in f:
        line = line.strip().split()
        if not line[0][0] == '#':
            atom_ID.append(int(line[0]))
            atom_pos.append([float(line[1]),float(line[2]),float(line[3])])
            atom_R.append(float(trunc(float(line[4]),2)))
    f.close()
    
    N_atoms = len(atom_ID)

    # make arrays of atom types and coordinates
    atom_R = numpy.array(atom_R)
    atom_ID = numpy.array(atom_ID,dtype=numpy.int32)
    atom_pos = numpy.array(atom_pos)

    # create dictionary with {atom_index : atom_type} and {atom_idex : [x,y,z]}
    # 2.1 = Si, 1.52 = O, 1.7 = C
    AtomType = {}
    AtomPos = {}
    for i in range(N_atoms):
        AtomType[i] = atom_R[i]
        AtomPos[i] = atom_pos[i]

    # import the Delaunay triangulations
    Tets = {} # dictionary {tetID : [atom1,atom2,atom3,atom4]}
    f = open(filenameTets)
    f.readline()
    for line in f:
        line = line.strip().split()
        if not line[0][0] == '#':
            Tets[int(line[0])] = [int(line[1]),int(line[2]),int(line[3]),int(line[4])]
    f.close()
    
    return N_atoms,AtomPos,AtomType,Tets

#################################################################################################
"""Voronoi Network and void computation (overlapping Voronoi vertices)"""

def VoronoiVertices(Tets,Tets_free_idx,AtomPos,AtomType):
    """
    Computes the circumsphere around each tet with free volume
    to obtain the voronoi vertex and radius accounting for the
    surrounding atoms size. NOTE: the voronoi index will be the same
    as the tets index for future calculations. 
    """
    t0v = time.time()
    VorVerts = {}
    for tetID in Tets_free_idx:
        Ratom = max([AtomType[Tets[tetID][0]],AtomType[Tets[tetID][1]],AtomType[Tets[tetID][2]]])
        
        x = AtomPos[Tets[tetID][0]]
        y = AtomPos[Tets[tetID][1]]
        z = AtomPos[Tets[tetID][2]]
    
        vertex,R0 = circumcenter(numpy.array([x,y,z]))
        Rvert = R0 - Ratom
        
        VorVerts[tetID] = (vertex,Rvert)
    print '\nComputed Voronoi vertices in %.4f seconds' %(time.time()-t0v)
    return VorVerts

def ComputeVoids(VorVerts,Tets_free_idx,min_deg_overlap):
    """
    This function finds all the clusters of vertices 
    """
    # compute distances; store tetID's as tuple key values in dictionary
    t0d = time.time()
    Dist = {} # {(tet1,tet2) : d}
    for i in range(len(VorVerts)):
        tetID1 = Tets_free_idx[i]
        x1 = VorVerts[tetID1][0][0]
        y1 = VorVerts[tetID1][0][1]
        z1 = VorVerts[tetID1][0][2]
        for j in range(i+1,len(VorVerts)):
            tetID2 = Tets_free_idx[j]
            x2 = VorVerts[tetID2][0][0]
            y2 = VorVerts[tetID2][0][1]
            z2 = VorVerts[tetID2][0][2]
            
            Dist[tuple(sorted([tetID1,tetID2]))] = Distance(x1-x2,y1-y2,z1-z2)
    print 'Computed distances between Voronoi vertices in %.4f seconds.' %(time.time()-t0d)
    
    # compute void clusters
    t0c = time.time()
    tets_guage = set([])
    TetsFreeID = set(Tets_free_idx)
    Voids = []
    while True:
        # get arbitrary tetID from TetsFreeID after taking
        # away the elements that are in the tets_guage 
        randID = next(iter(TetsFreeID.difference(tets_guage)))
        void = [randID]
        tets_guage.add(randID)
        
        while True:
            TetsLeft = TetsFreeID.difference(tets_guage)
            size_orig_void = len(void)
            void, tets_guage = region(void,tets_guage,TetsLeft,VorVerts,Dist,min_deg_overlap)
            size_new_void = len(void)
            
            if size_new_void == size_orig_void:
                Voids.append(void)
                break
            
        # if all indicies are used, break out of loop
        if not len(TetsFreeID.difference(tets_guage)): 
            break
    print 'Computed clusters in %.4f seconds.\n' %(time.time()-t0c)
    return Voids    
    
def region(orig_cluster,tets_guage,TetsLeft,VorVerts,Dist,min_deg_overlap):
    """
    This function searches through verticies outside of the
    current void and adds them to the void if they are overlapping
    with any vertex in the original void
    """
    new_cluster = []
    for tetID1 in orig_cluster:
        R1 = VorVerts[tetID1][1]
        for tetID2 in TetsLeft:
            R2 = VorVerts[tetID2][1]
            d = Dist[tuple(sorted([tetID1,tetID2]))]
            
            Dij = R1 + R2 - d
            if Dij >= min_deg_overlap:
                new_cluster.append(tetID2)
                tets_guage.add(tetID2)               
                
    new_cluster += orig_cluster
    
    return new_cluster, tets_guage

def VoidVolumes(Voids,V_free_ID):
    VolVoids = []
    for void in Voids:
        Vvoid = 0
        for idx in void:
            Vvoid += V_free_ID[idx]
        VolVoids.append(Vvoid)
    return VolVoids
            
            
#################################################################################################
"""Free volume computation within tetrahedrons"""

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
        
    x1 = vertexA[0]
    y1 = vertexA[1]
    z1 = vertexA[2]

    d = -(vertexB[0]*normal[0]+vertexB[1]*normal[1]+vertexB[2]*normal[2])

    k = -(d + a*x1 + b*y1 + c*z1)/(a**2 + b**2 + c**2)

    H = numpy.array([x1+k*a,y1+k*b,z1+k*c])
    vecAH = vertexA - H

    # check if the normal of the plane is in the same direction as the height vector
    if numpy.array_equal(numpy.sign(normal),numpy.sign(vecAH)):
        if numpy.array_equal(H,vertexA):
            # if the atoms are all coplanar, take V~0
            return 0,H
        else:
            return 1,H
    else:
        # if they are not in the same direction, switch orientation of plane
        T,H = PointH(vec2,vec1,vertexA,vertexB)
        return T, H

def FracAtomVol(R,vertex1,vertex2,vertex3,vertex4):
    """
    vertex1, vertex2, and vertex3 are coplanar; vertex4 sits above the
    plane
    """
    # using angle1 as origin, find angle between vectors in base plane 
    vec1 = vertex2 - vertex1
    vec2 = vertex3 - vertex1
    alpha = Angle(vec1,vec2)

    # find angle in from base plane towards 'z' direcion 
    # must find vector in base plane 
    vec3 = vertex4 - vertex1
    T,pointH = PointH(vec1,vec2,vertex4,vertex1)

    if T:
        vec4 = pointH - vertex1
        beta = Angle(vec3,vec4)

        # compute fraction of volume between sphere and tetrahedron
        V_frac = (0.5)*SphereVol(R)*(beta/numpy.pi)*(alpha/2/numpy.pi)
        return V_frac
    else:
        return 0

def ComputeOverlap(atom1,atom2,R1,R2):
    """Computes the volume overlap between two spheres """
    d = Distance(atom1[0]-atom2[0],atom1[1]-atom2[1],atom1[2]-atom2[2])
    
    V_overlap = 0
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

    return V_overlap

def FreeVolume(Tets,AtomPos,AtomType,min_free_vol):
    """Computes the free volume within a tetrahedron formed by 4 atoms"""
    R_Si = 2.1 #van der waals radius in Angstroms 
    R_C = 1.7
    R_O = 1.52

    t0 = time.time()
    V_free = [] # free volume in tets that meet criteria 
    V_Tets = [] # volume of all tets 
    Tets_free_idx = [] # index of tets with free volume
    V_free_ID = {} # {tetID: freeVol}
    N_Sm_holes = 0
    NSmTets = 0
    NLgOverlap = 0
    for i in range(len(Tets)):
        tet = Tets[i]
        
        # compute overlap of tet vertices (atoms w/ vdW radii)
        V_overlap = 0
        N_overlap = 0
        for j in range(len(tet)):
            for k in range(j+1,len(tet)):
                V_overlap += ComputeOverlap(AtomPos[tet[j]],AtomPos[tet[k]],
                                            AtomType[tet[j]],AtomType[tet[k]])
                if V_overlap:
                    N_overlap += 1

        # compute volume of tetrahedrons; if volume is < 3.5 A^3
        # then its too small and neglect it; also if overlap is
        # close to tet volume then neglect 
        V_tet = TetVolume(AtomPos[tet[0]],AtomPos[tet[1]],
                          AtomPos[tet[2]],AtomPos[tet[3]])

        V_Tets.append(V_tet)

        if V_tet < min_free_vol:
            NSmTets += 1
            continue
        elif V_tet - V_overlap < 3.5:
            NLgOverlap += 1

        # if 3 or more overlapping atoms then space is
        # too small for diffusing particles 
        if N_overlap >= 3:
            N_Sm_holes += 1
                    
        # compute volume of atoms that overlap with tetrahedron
        V_atoms = 0
        R = []
        for j in range(len(tet)):
            V_atoms += FracAtomVol(AtomType[tet[j]],AtomPos[tet[j%4]],AtomPos[tet[(j+1)%4]],
                                   AtomPos[tet[(j+2)%4]],AtomPos[tet[(j+3)%4]])
    
        # add to free volume; free volume must be large enough to occupy
        # a hydrogen atom and cannot exceed the volume of the tet
        Vfree = V_tet - V_atoms + V_overlap

        if min_free_vol < Vfree < V_tet:
            V_free.append(Vfree)
            Tets_free_idx.append(i)
            V_free_ID[i] = Vfree

    V_free = numpy.array(V_free)
    V_Tets = numpy.array(V_Tets)
    Tets_free_idx = numpy.array(Tets_free_idx)
    print 'Computed tetrahedron volumes in %.4f seconds.' %(time.time()-t0)

    return V_free,V_Tets,Tets_free_idx,V_free_ID,(N_Sm_holes,NSmTets,NLgOverlap)   
    
    
        
        

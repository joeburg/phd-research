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
import numpy
import numpy.linalg
import pylab
import scipy
import time

##from matplotlib import rc
##rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
##rc('text', usetex=True)


##font = {'family' : 'normal','weight' : 'bold','size' : 22}
##matplotlib.rc('font', **font)
################################################################################
def trunc(f, n):
    '''Truncates/pads a float f to n decimal places without rounding'''
    slen = len('%.*f' % (n, f))
    return str(f)[:slen]

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

def AtomVolume(AtomType):
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

def Distance(dx,dy,dz):
    d = (dx**2+dy**2+dz**2)**0.5
    return d      

def TetVolume(a, b, c, d):
    return numpy.abs(numpy.dot(a-d,numpy.cross(b-d, c-d))) / 6

def Angle(vec1,vec2):
    ab = numpy.dot(vec1,vec2)
    a = numpy.linalg.norm(vec1)
    b = numpy.linalg.norm(vec2)
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

def FreeVolume(Tets,AtomPos,AtomType):
    R_Si = 2.1 #van der waals radius in Angstroms 
    R_C = 1.7
    R_O = 1.52

    V_free = [] # free volume in tets that meet criteria 
    V_Tets = [] # volume of all tets 
    Tets_free_idx = [] # index of tets with free volume 
    N_Sm_holes = 0
    NSmTets = 0
    NLgOverlap = 0
    for i in range(len(Tets)):
##        if i % 10000 == 0:
##            print 'Iteration: %d' %i

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

        if V_tet < 5.4:
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

        if 5.4 < Vfree < V_tet:
            V_free.append(Vfree)
            Tets_free_idx.append(i)

    V_free = numpy.array(V_free)
    V_Tets = numpy.array(V_Tets)
    Tets_free_idx = numpy.array(Tets_free_idx)

    return V_free,V_Tets,Tets_free_idx,(N_Sm_holes,NSmTets,NLgOverlap)      

##def Plot(ax,xlabel,ylabel,fontsize):
##    ax.plot([1,2])
##    ax.locator_params(nbins=3)
##    ax.set_xlabel(xlabel,fontsize=fontsize)
##    ax.set_ylabel(ylabel,fontsize=fontsize)
##    ax.set_xlim((0,xmax))
##    ax.set_ylim((0,ymax))

################################################################################
#import .xyz file to correlate indicies with positions
##filenameNodes = raw_input('What is the .node filename? ')
##filenameTets = raw_input('What is the .ele filename? ')

filenameNodes = 'OCSEt_140000.1.node'
filenameTets = 'OCSEt_140000.1.ele'
filenameLAMMPS = 'OCSEt_140000.xyz'

t0 = time.time()

N_atoms,AtomPos,AtomType,Tets = LoadDataUnweighted(filenameNodes,filenameTets,filenameLAMMPS)
atom_volume = AtomVolume(AtomType)

print '\nAtom volume: %.4f' %atom_volume
print 'N_atoms = %d' %N_atoms

print 'Number Delaunay tetrahedrons: %d\n' %len(Tets)

V_free,V_Tets,Tets_free_idx,TetStats = FreeVolume(Tets,AtomPos,AtomType)
##print V_free
print 'Number of small holes: %d' %TetStats[0]
print 'Number of small tets: %d' %TetStats[1]
print 'Number of large overlaps: %d' %TetStats[2]
print 'Percent of tets with acceptable free volume: %.4f' %(float(len(V_free))/len(Tets))

##print V_free[:5]
##print V_Tets[:5]


pylab.figure(1)
n1,bins1,patches1 = pylab.hist(V_free,bins=200,normed=1,histtype='stepfilled')
pylab.setp(patches1,'facecolor','b','alpha',0.75)
pylab.xlabel('Free volume per tetrahedron (A^3)')
pylab.ylabel('Frequency')
pylab.xlim((0,25))
pylab.ylim((0,0.6))
pylab.savefig('FreeVolume.jpg')
pylab.close(1)

pylab.figure(2)
n2,bins2,patches2 = pylab.hist(V_Tets,bins=200,normed=1,histtype='stepfilled')
pylab.setp(patches2,'facecolor','g','alpha',0.75)
pylab.xlabel('Tetrahedron Volume (A^3)')
pylab.ylabel('Frequency')
pylab.xlim((0,25))
pylab.ylim((0,0.4))
pylab.savefig('TetsVolume.jpg')
pylab.close(2)

Vtets = numpy.sum(V_Tets)
VFree = numpy.sum(V_free)
Vtot = atom_volume + VFree

print 'Tetrahedron volume: %.4f' %Vtets
print 'Vol_atom + Vol_free = %.4f' %(atom_volume+VFree)
print '(Vatoms + Vfree)/Vtot = %.4f' %(Vtot/Vtets)
print 'Packing Fraction: %.4f' %(atom_volume/Vtot)
print 'Free volume: %.4f' %(VFree/Vtot)

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


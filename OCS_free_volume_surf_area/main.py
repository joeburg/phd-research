"""This main program invokes functions in FreeVolFunctions.py and
circumsphere.py to compute the free volume in an OCS structure"""

import FreeVolFunctions as FV
import matplotlib.pyplot as plt
import pylab
import time
import numpy
################################################################################
t0 = time.time()

filenameNodes = 'OCSEt_140000.1.node'
filenameTets = 'OCSEt_140000.1.ele'
filenameLAMMPS = 'OCSEt_140000.xyz'
N_atoms,AtomPos,AtomType,Tets = FV.LoadDataUnweighted(filenameNodes,filenameTets,filenameLAMMPS)

##filenameNodes = 'OCSEt_140000_w.1.node'
##filenameTets = 'OCSEt_140000_w.1.ele'
##N_atoms,AtomPos,AtomType,Tets = FV.LoadDataWeighted(filenameNodes,filenameTets)

min_free_vol = 5.0;
min_deg_overlap = 0;

atom_volume = FV.AtomVolume(AtomType,N_atoms)

print '\nAtom volume: %.4f' %atom_volume
print 'N_atoms = %d' %N_atoms
print 'Number Delaunay tetrahedrons: %d\n' %len(Tets)

V_free,V_Tets,Tets_free_idx,V_free_ID,TetStats = FV.FreeVolume(Tets,AtomPos,AtomType,min_free_vol)

print 'Number of small holes: %d' %TetStats[0]
print 'Number of small tets: %d' %TetStats[1]
print 'Number of large overlaps: %d' %TetStats[2]
print 'Percent of tets with acceptable free volume: %.4f' %(float(len(V_free))/len(Tets))

##def GetType(R):
##    if 2.0 < R < 2.2:
##        return 'Si'
##    elif 1.4 < R < 1.55:
##        return 'O'
##    elif 1.6 < R < 1.8:
##        return 'C'
##
##def WriteLine(typeID,atom):
##    f.write('%s  %.4f  %.4f  %.4f\n' %(typeID,atom[0],atom[1],atom[2])) 
##
##
##maxID = 0
##for idx in Tets:
##    for i in range(len(Tets[idx])):
##        if Tets[idx][i] > maxID:
##            maxID = Tets[idx][i]
##
##print maxID
##
##IDs0 = set([])
##for idx in AtomPos:
##    if idx not in IDs0:
##        IDs0.add(idx)
##
##IDs = set([])
##f = open('Atoms_in_delaunay.xyz','w')
####f.write('%d\n' %(len(Tets)*4))
##f.write('Atoms in the Delaunay triangulation\n')
##for idx in Tets:
##    for i in range(len(Tets[idx])):
##        if Tets[idx][i] not in IDs:
##            IDs.add(Tets[idx][i])
##            WriteLine(GetType(AtomType[Tets[idx][i]]),AtomPos[Tets[idx][i]])
##f.close()
##
##print IDs0.difference(IDs)

##f = open("Free_Volume_tets_VMD.xyz",'w')
##f.write('%d\n' %len(Tets_free_idx*4))
##f.write('Tets with relevant free volume\n')
##for idx in Tets_free_idx:
##    for i in range(len(Tets[idx])):
##        WriteLine(GetType(AtomType[i]),AtomPos[Tets[idx][i]])
##f.close()   

##Vtets = numpy.sum(V_Tets)
##VFree = numpy.sum(V_free)
##Vtot = atom_volume + VFree
##
##print 'Tetrahedron volume: %.4f' %Vtets
##print 'Vol_atom + Vol_free = %.4f' %(atom_volume+VFree)
##print '(Vatoms + Vfree)/Vtot = %.4f' %(Vtot/Vtets)
##print 'Packing Fraction: %.4f' %(atom_volume/Vtot)
##print 'Free volume: %.4f' %(VFree/Vtot)
##
### Compute Voids
##VorVerts = FV.VoronoiVertices(Tets,Tets_free_idx,AtomPos,AtomType)
##Voids = FV.ComputeVoids(VorVerts,Tets_free_idx,min_deg_overlap)
##VolVoids = FV.VoidVolumes(Voids,V_free_ID)
##
##print 'Number of voids: %d' %len(Voids)
##
##pylab.figure(1)
##weights1 = numpy.ones_like(V_free)/len(V_free)
##n1,bins1,patches1 = pylab.hist(V_free,bins=200,weights=weights1,histtype='stepfilled')
##pylab.setp(patches1,'facecolor','b','alpha',0.75)
##pylab.xlabel('Free volume per tetrahedron (A^3)')
##pylab.ylabel('Frequency')
##pylab.xlim((0,25))
##pylab.ylim((0,0.6))
##pylab.savefig('FreeVolume.jpg')
##pylab.close(1)
##
##pylab.figure(2)
##weights2 = numpy.ones_like(V_Tets)/len(V_Tets)
##n2,bins2,patches2 = pylab.hist(V_Tets,bins=200,weights=weights2,histtype='stepfilled')
##pylab.setp(patches2,'facecolor','g','alpha',0.75)
##pylab.xlabel('Tetrahedron Volume (A^3)')
##pylab.ylabel('Frequency')
##pylab.xlim((0,25))
##pylab.ylim((0,0.4))
##pylab.savefig('TetsVolume.jpg')
##pylab.close(2)
##
##pylab.figure(3)
##weights3 = numpy.ones_like(VolVoids)/len(VolVoids)
##n2,bins2,patches2 = pylab.hist(VolVoids,bins=200,weights=weights3,histtype='stepfilled')
##pylab.setp(patches2,'facecolor','b','alpha',0.75)
##pylab.xlabel('Void Volume (A^3)')
##pylab.ylabel('Frequency')
####pylab.xlim((0,25))
####pylab.ylim((0,0.4))
##pylab.savefig('VoidVolume.jpg')
##pylab.close(3)

t1 = time.time()
print 'Elapsed time: %.3f seconds' % (t1-t0)

""" This program analyzes the Voronoi cell network """

import math
import matplotlib.pyplot as plt
import pylab
import time
import numpy

################################################################################
def SphereVol(R):
    return 4/3*math.pi*R**3

def LoadDataLAMMPS(filename):
    atomType = {}
    f = open(filename)
    Natoms = int(f.readline())
    f.readline()
    for i in range(Natoms):
        fields = f.readline().strip().split()
        if fields[0] == '1':
            R = 2.1
        elif fields[0] == '2':
            R = 1.52
        elif fields[0] == '3':
            R = 1.7
        atomType[i] = R
    return atomType
        
def LoadDataUnWeighted(filename):
    data = {} # {atomID : [x,y,z,v]}
    f = open(filename)
    for line in f:
        fields = line.strip().split()
        atomID = int(fields[0])-1
        xcoord = float(fields[1])
        ycoord = float(fields[2])
        zcoord = float(fields[3])
        Vcell = float(fields[4])

        data[atomID] = [xcoord,ycoord,zcoord,Vcell]
    f.close()
    return data


def CellFreeVol(data,atomType):
    Vtot = 0
    Vfreetot = 0
    VorVol = []
    VCellsFree = []
    CVol = []
    CVolFree = []
    SiVol = []
    SiVolFree = []
    OVol = []
    OVolFree = []
    for idx in data:
        Vcell = data[idx][3]
        
        Vtot += Vcell
        
        Vfree = Vcell - SphereVol(atomType[idx])

##        if Vfree > 0:
        Vfreetot += Vfree
        VCellsFree.append(Vfree)

        VorVol.append(Vcell)

        # distinguish between atom types
        if atomType[idx] == 2.1:
            SiVol.append(Vcell)
            SiVolFree.append(Vfree)
        elif atomType[idx] == 1.52:
            OVol.append(Vcell)
            OVolFree.append(Vfree)
        elif atomType[idx] == 1.7:
            CVol.append(Vcell)
            CVolFree.append(Vfree)
        
    VCellsFree = numpy.array(VCellsFree)
    VorVol = numpy.array(VorVol)
    SiVol = numpy.array(SiVol)
    SiVolFree = numpy.array(SiVolFree)
    OVol = numpy.array(OVol)
    OVolFree = numpy.array(OVolFree)
    CVol = numpy.array(CVol)
    CVolFree = numpy.array(CVolFree)
    
    return VCellsFree,Vtot,Vfreetot,VorVol,SiVol,SiVolFree,OVol,OVolFree,CVol,CVolFree
        
    
################################################################################
t0 = time.time()

fileLAMMPS1 = "OCSEt_155000.xyz"
fileLAMMPS2 = "OCSEt_255000.xyz"
fileLAMMPS3 = "OCSEt_355000.xyz"
filename1 = "OCSEt_155000_voro.vol"
filename2 = "OCSEt_255000_voro.vol"
filename3 = "OCSEt_355000_voro.vol"

##atomType = LoadDataLAMMPS(fileLAMMPS)
##data = LoadDataUnWeighted(filename)
##VCellsFree,Vtot,Vfreetot,VorVol,SiVol,SiVolFree,OVol,OVolFree,CVol,CVolFree = CellFreeVol(data,atomType)

atomType1 = LoadDataLAMMPS(fileLAMMPS1)
data1 = LoadDataUnWeighted(filename1)
VCellsFree1,Vtot1,Vfreetot1,VorVol1,SiVol1,SiVolFree1,OVol1,OVolFree1,CVol1,CVolFree1 = CellFreeVol(data1,atomType1)

atomType2 = LoadDataLAMMPS(fileLAMMPS2)
data2 = LoadDataUnWeighted(filename2)
VCellsFree2,Vtot2,Vfreetot2,VorVol2,SiVol2,SiVolFree2,OVol2,OVolFree2,CVol2,CVolFree2 = CellFreeVol(data2,atomType1)

atomType3 = LoadDataLAMMPS(fileLAMMPS3)
data3 = LoadDataUnWeighted(filename3)
VCellsFree3,Vtot3,Vfreetot3,VorVol3,SiVol3,SiVolFree3,OVol3,OVolFree3,CVol3,CVolFree3 = CellFreeVol(data3,atomType1)

##print "\nWorking with %s\n" %filename
##print "Vtot = %.4f" %Vtot
##print "Vfreetot = %.4f" %Vfreetot
##print "Free volume fraction = %.4f" %(Vfreetot/Vtot)

print "\nFree volume fraction = %.4f" %(Vfreetot1/Vtot1)
print "Free volume fraction = %.4f" %(Vfreetot2/Vtot2)
print "Free volume fraction = %.4f" %(Vfreetot3/Vtot3)

##pylab.figure(1)
##pylab.hist((VCellsFree,VorVol),bins=300,normed=True,histtype='step',lw=1,color=('b','k'))
##pylab.xlabel('free volume per atom (A^3)')
##pylab.ylabel('frequency')
####pylab.xlim((0,60))
####pylab.ylim((0,0.20))
##pylab.savefig('FreeVol_'+filename[6:12]+'.jpg')
##pylab.close(1)

##pylab.figure(2)
##pylab.hist(SiVol,bins=300,normed=True,histtype='step',lw=1,color='k')
##pylab.xlabel('free volume per Si atom (A^3)')
##pylab.ylabel('frequency')
####pylab.xlim((0,60))
####pylab.ylim((0,0.20))
##pylab.savefig('SiFreeVol_'+filename[6:12]+'.jpg')
##pylab.close(2)
##
##pylab.figure(3)
##pylab.hist(OVol,bins=300,normed=True,histtype='step',lw=1,color='k')
##pylab.xlabel('free volume per O atom (A^3)')
##pylab.ylabel('frequency')
####pylab.xlim((0,60))
####pylab.ylim((0,0.20))
##pylab.savefig('OFreeVol_'+filename[6:12]+'.jpg')
##pylab.close(3)
##
##pylab.figure(4)
##pylab.hist(CVol,bins=300,normed=True,histtype='step',lw=1,color='k')
##pylab.xlabel('free volume per C atom (A^3)')
##pylab.ylabel('frequency')
####pylab.xlim((0,60))
####pylab.ylim((0,0.20))
##pylab.savefig('CFreeVol_'+filename[6:12]+'.jpg')
##pylab.close(4)

pylab.figure(5)
pylab.hist((VorVol1,VorVol2,VorVol3),bins=200,normed=True,histtype='step',lw=1,color=('b','g','r'))
pylab.xlabel('free volume per atom (A^3)')
pylab.ylabel('frequency')
##pylab.xlim((-10,60))
##pylab.ylim((0,0.20))
pylab.savefig('All_Cells_Peak_Shift.jpg')
pylab.close(5)

pylab.figure(6)
pylab.hist((SiVol1,SiVol2,SiVol3),bins=200,normed=True,histtype='step',lw=1,color=('b','g','r'))
pylab.xlabel('free volume per Si atom (A^3)')
pylab.ylabel('frequency')
##pylab.xlim((0,60))
##pylab.ylim((0,0.20))
pylab.savefig('Si_atom_peak_shift.jpg')
pylab.close(6)

pylab.figure(7)
pylab.hist((CVol1,CVol2,CVol3),bins=200,normed=True,histtype='step',lw=1,color=('b','g','r'))
pylab.xlabel('free volume per C atom (A^3)')
pylab.ylabel('frequency')
##pylab.xlim((0,60))
##pylab.ylim((0,0.20))
pylab.savefig('C_atom_peak_shift.jpg')
pylab.close(7)

pylab.figure(8)
pylab.hist((OVol1,OVol2,OVol3),bins=300,normed=True,histtype='step',lw=1,color=('b','g','r'))
pylab.xlabel('free volume per O atom (A^3)')
pylab.ylabel('frequency')
##pylab.xlim((15,40))
##pylab.ylim((0,0.04))
pylab.savefig('O_atom_peak_shift.jpg')
pylab.close(8)

pylab.figure(5)
pylab.hist((VCellsFree1,VCellsFree2,VCellsFree3),bins=200,normed=True,histtype='step',lw=1,color=('b','g','r'))
pylab.xlabel('free volume per atom (A^3)')
pylab.ylabel('frequency')
##pylab.xlim((-10,60))
##pylab.ylim((0,0.20))
pylab.savefig('All_Cells_Peak_Shift_free.jpg')
pylab.close(5)

t1 = time.time()
print 'Elapsed time: %.3f seconds' % (t1-t0)

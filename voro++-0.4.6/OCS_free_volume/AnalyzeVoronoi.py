""" This program analyzes the Voronoi cell network """

import math
import matplotlib.pyplot as plt
import pylab
import time
import numpy

################################################################################
def SphereVol(R):
    return 4/3*math.pi*R**3

def LoadDataWeighted(filename):
    data = {} # {atomID : [x,y,z,v,r]}
    f = open(filename)
    for line in f:
        fields = line.strip().split()
        atomID = int(fields[0])-1
        xcoord = float(fields[1])
        ycoord = float(fields[2])
        zcoord = float(fields[3])
        Vcell = float(fields[4])
        R = float(fields[5])
        Nverts = int(fields[6])
        Nedges = int(fields[7])
        Nfaces = int(fields[8])

        data[atomID] = [xcoord,ycoord,zcoord,Vcell,R,Nverts,Nedges,Nfaces]
    f.close()
    return data


def CellFreeVol(data):
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
        
        Vfree = Vcell - SphereVol(data[idx][4])

        if Vfree > 0:
            Vfreetot += Vfree
            VCellsFree.append(Vfree)

        VorVol.append(Vcell)

        # distinguish between atom types
        if data[idx][4] == 2.1:
            SiVol.append(Vcell)
            SiVolFree.append(Vfree)
        elif data[idx][4] == 1.52:
            OVol.append(Vcell)
            OVolFree.append(Vfree)
        elif data[idx][4] == 1.7:
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

def VoronoiStats(data):
    Verts = []
    O_verts = []
    C_verts = []
    Si_verts = []
    Edges = []
    O_edges = []
    C_edges = []
    Si_edges = []
    Faces = []
    O_faces = []
    C_faces = []
    Si_faces = []
    for idx in data:
        vert = data[idx][5]
        edge = data[idx][6]
        face = data[idx][7]

        Verts.append(vert)
        Edges.append(edge)
        Faces.append(face)

        if data[idx][4] == 2.1:
            Si_verts.append(vert)
            Si_edges.append(edge)
            Si_faces.append(face)
        elif data[idx][4] == 1.52:
            O_verts.append(vert)
            O_edges.append(edge)
            O_faces.append(face)
        elif data[idx][4] == 1.7:
            C_verts.append(vert)
            C_edges.append(edge)
            C_faces.append(face)

    return numpy.array(Verts),numpy.array(Edges),numpy.array(Faces),\
           numpy.array(Si_verts),numpy.array(Si_edges),numpy.array(Si_faces),\
           numpy.array(O_verts),numpy.array(O_edges),numpy.array(O_faces),\
           numpy.array(C_verts),numpy.array(C_edges),numpy.array(C_faces)

def Plot1(N,data,Nbins,filename,xlbl):
    pylab.figure(N)
    pylab.hist(data,bins=Nbins,normed=True,histtype='step',lw=1,color='k')
    pylab.xlabel(xlbl)
    pylab.ylabel('frequency')
##    pylab.xlim((0,60))
##    pylab.ylim((0,0.20))
    pylab.savefig(filename)
    pylab.close(N)

def Plot2(N,data1,data2,Nbins,filename,xlbl):
    pylab.figure(N)
    pylab.hist((data1,data2),bins=Nbins,normed=True,histtype='step',lw=1,color=('b','k'))
    pylab.xlabel(xlbl)
    pylab.ylabel('frequency')
##    pylab.xlim((0,60))
##    pylab.ylim((0,0.20))
    pylab.savefig(filename)
    pylab.close(N)

def Plot3(N,data1,data2,data3,Nbins,filename,xlbl):
    pylab.figure(N)
    pylab.hist((data1,data2,data3),bins=Nbins,normed=True,histtype='step',lw=1,color=('b','g','r'))
    pylab.xlabel(xlbl)
    pylab.ylabel('frequency')
##    pylab.xlim((0,60))
##    pylab.ylim((0,0.20))
    pylab.savefig(filename)
    pylab.close(N)

################################################################################
t0 = time.time()

filename1 = "OCSEt_155000_voro.vol"
filename2 = "OCsEt_255000_voro.vol"
filename3 = "OCSEt_355000_voro.vol"

data1 = LoadDataWeighted(filename1)
VCellsFree1,Vtot1,Vfreetot1,VorVol1,SiVol1,SiVolFree1,OVol1,OVolFree1,CVol1,CVolFree1 = CellFreeVol(data1)
Verts1,Edges1,Faces1,Si_verts1,Si_edges1,Si_faces1,O_verts1,O_edges1,O_faces1,C_verts1,C_edges1,C_faces1 = VoronoiStats(data1)

data2 = LoadDataWeighted(filename2)
VCellsFree2,Vtot2,Vfreetot2,VorVol2,SiVol2,SiVolFree2,OVol2,OVolFree2,CVol2,CVolFree2 = CellFreeVol(data2)
Verts2,Edges2,Faces2,Si_verts2,Si_edges2,Si_faces2,O_verts2,O_edges2,O_faces2,C_verts2,C_edges2,C_faces2 = VoronoiStats(data2)

data3 = LoadDataWeighted(filename3)
VCellsFree3,Vtot3,Vfreetot3,VorVol3,SiVol3,SiVolFree3,OVol3,OVolFree3,CVol3,CVolFree3 = CellFreeVol(data3)
Verts3,Edges3,Faces3,Si_verts3,Si_edges3,Si_faces3,O_verts3,O_edges3,O_faces3,C_verts3,C_edges3,C_faces3 = VoronoiStats(data3)


##print "\nWorking with %s\n" %filename
##
##print "Vtot = %.4f" %Vtot
##print "Vfreetot = %.4f" %Vfreetot
##print "Free volume fraction = %.4f" %(Vfreetot/Vtot)

##Plot2(1,VCellsFree1,VorVol1,300,'FreeVol_'+filename1[6:12]+'.jpg','free volume per atom (A^3)')
##Plot1(2,SiVol1,300,'Si_FreeVol_'+filename1[6:12]+'.jpg','free volume per Si atom (A^3)')
##Plot1(3,OVol1,300,'O_FreeVol_'+filename1[6:12]+'.jpg','free volume per O atom (A^3)')
##Plot1(4,CVol1,300,'C_FreeVol_'+filename1[6:12]+'.jpg','free volume per C atom (A^3)')
##Plot3(5,VorVol1,VorVol2,VorVol3,200,'All_Cells_Peak_Shift.jpg','free volume per atom (A^3)')
##Plot3(6,SiVol1,SiVol2,SiVol3,200,'Si_atom_peak_shift.jpg','free volume per Si atom (A^3)')
##Plot3(7,SiVol1,SiVol2,SiVol3,200,'C_atom_peak_shift.jpg','free volume per C atom (A^3)')
##Plot3(8,OVol1,OVol2,OVol3,200,'O_atom_peak_shift.jpg','free volume per O atom (A^3)')
##Plot3(9,VCellsFree1,VCellsFree2,VCellsFree3,200,'All_Cells_Peak_Shift_free.jpg','free volume per atom (A^3)')

##Plot3(10,Verts1,Verts2,Verts3,600,'All_voronoi_verts.jpg','Number of Vonoroi Vertices Per Cell')
##Plot3(11,Edges1,Edges2,Edges3,600,'All_voronoi_edges.jpg','Number of Vonoroi Edges Per Cell')
##Plot3(12,Faces1,Faces2,Faces3,600,'All_voronoi_faces.jpg','Number of Vonoroi Faces Per Cell')
##Plot3(13,Si_verts1,Si_verts2,Si_verts3,200,'Si_voronoi_verts.jpg','Number of Vonoroi Vertices Per Si Cell')
##Plot3(14,Si_edges1,Si_edges2,Si_edges3,200,'Si_voronoi_edges.jpg','Number of Vonoroi Edges Per Si Cell')
##Plot3(15,Si_faces1,Si_faces2,Si_faces3,200,'Si_voronoi_faces.jpg','Number of Vonoroi Faces Per Si Cell')
##Plot3(16,O_verts1,O_verts2,O_verts3,200,'O_voronoi_verts.jpg','Number of Vonoroi Vertices Per O Cell')
##Plot3(17,O_edges1,O_edges2,O_edges3,200,'O_voronoi_edges.jpg','Number of Vonoroi Edges Per O Cell')
##Plot3(18,O_faces1,O_faces2,O_faces3,200,'O_voronoi_faces.jpg','Number of Vonoroi Faces Per O Cell')
##Plot3(19,C_verts1,C_verts2,C_verts3,200,'C_voronoi_verts.jpg','Number of Vonoroi Vertices Per C Cell')
##Plot3(20,C_edges1,C_edges2,C_edges3,200,'C_voronoi_edges.jpg','Number of Vonoroi Edges Per C Cell')
##Plot3(21,C_faces1,C_faces2,C_faces3,200,'C_voronoi_faces.jpg','Number of Vonoroi Faces Per C Cell')


t1 = time.time()
print 'Elapsed time: %.3f seconds' % (t1-t0)

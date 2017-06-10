""" This program analyzes the Voronoi cell network """

import glob
import math
import matplotlib.pyplot as plt
import numpy
import pylab
import scipy.interpolate
import sys
import time

from scipy.stats import gaussian_kde


################################################################################
def SphereVol(R):
    return 4/3*math.pi*R**3

def Distance(dx,dy,dz):
    return (dx*dx + dy*dy + dz*dz)**0.5

def LoadDataWeighted(filename):
    NSi = 0
    NO = 0

    data = {} # {atomID : [x,y,z,v,r]}
    f = open(filename)
    for line in f:
        fields = line.strip().split()
        atomID = int(fields[0])-1 # convert to 0 based indexing
        xcoord = float(fields[1])
        ycoord = float(fields[2])
        zcoord = float(fields[3])
        Vcell = float(fields[4])
        R = float(fields[5])

        if R == 2.1:
            NSi += 1
        elif R == 1.52:
            NO += 1

        Nverts = int(fields[6])
        Nedges = int(fields[7])
        Nfaces = int(fields[8])

        data[atomID] = [atomID,xcoord,ycoord,zcoord,Vcell,R,Nverts,Nedges,Nfaces]
        # data[atomID] = [atomID,xcoord,ycoord,zcoord,Vcell,R]
    f.close()

    # get dimensions of simulation cell and create coord array
    data_xyz = numpy.zeros((NSi+NO,4)) # numpy array for computing terminal O
    idx_Si = 0
    idx_O = NSi 
    x_min = 20
    x_max = 20
    y_min = 20
    y_max = 20
    z_min = 20
    z_max = 20
    for atomID in data:
        x = data[atomID][1]
        y = data[atomID][2]
        z = data[atomID][3]

        if x < x_min:
            x_min = x
        if x > x_max:
            x_max = x

        if y < y_min:
            y_min = y
        if y > y_max:
            y_max = y

        if z < z_min:
            z_min = z
        if z > z_max:
            z_max = z

        if data[atomID][5] == 2.1:
            data_xyz[idx_Si][0] = atomID
            data_xyz[idx_Si][1] = x
            data_xyz[idx_Si][2] = y
            data_xyz[idx_Si][3] = z
            idx_Si += 1
        elif data[atomID][5] == 1.52:
            data_xyz[idx_O][0] = atomID
            data_xyz[idx_O][1] = x
            data_xyz[idx_O][2] = y
            data_xyz[idx_O][3] = z
            idx_O += 1

    Lx = x_max - x_min
    Ly = y_max - y_min
    Lz = z_max - z_min

    return data, data_xyz, NSi, NO, Lx, Ly, Lz


def CellVolStats(data):
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
        Vcell = data[idx][4]
        
        Vtot += Vcell
        
        Vfree = Vcell - SphereVol(data[idx][5])

        if Vfree > 0:
            Vfreetot += Vfree
            VCellsFree.append(Vfree)

        VorVol.append(Vcell)

        # distinguish between atom types
        if data[idx][5] == 2.1:
            SiVol.append(Vcell)
            SiVolFree.append(Vfree)
        elif data[idx][5] == 1.52:
            OVol.append(Vcell)
            OVolFree.append(Vfree)
        elif data[idx][5] == 1.7:
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
        vert = data[idx][6]
        edge = data[idx][7]
        face = data[idx][8]

        Verts.append(vert)
        Edges.append(edge)
        Faces.append(face)

        if data[idx][5] == 2.1:
            Si_verts.append(vert)
            Si_edges.append(edge)
            Si_faces.append(face)
        elif data[idx][5] == 1.52:
            O_verts.append(vert)
            O_edges.append(edge)
            O_faces.append(face)
        elif data[idx][5] == 1.7:
            C_verts.append(vert)
            C_edges.append(edge)
            C_faces.append(face)

    return [numpy.array(Verts),numpy.array(Edges),numpy.array(Faces),
            numpy.array(Si_verts),numpy.array(Si_edges),numpy.array(Si_faces),
            numpy.array(O_verts),numpy.array(O_edges),numpy.array(O_faces),
            numpy.array(C_verts),numpy.array(C_edges),numpy.array(C_faces)]

def GetTerminalO(data_all,data,NSi,NO,Lx,Ly,Lz):
    print '\nComputing terminal O ...'
    CM = numpy.zeros((NSi,NO))
    
    for i in range(NSi):
        for j in range(NSi,len(data)):

            dx = abs(data[j][1] - data[i][1])
            dy = abs(data[j][2] - data[i][2])
            dz = abs(data[j][3] - data[i][3])

            if Distance(dx,dy,dz) < 2.3:
                CM[i][j-NSi] = 1
            
            # consider PBCs
            if dx > Lx - 2.3:
                if Distance(dx + Lx, dy, dz) < 2.3:
                    CM[i][j-NSi] = 1

                if Distance(dx - Lx, dy, dz) < 2.3:
                    CM[i][j-NSi] = 1

            if dy > Ly - 2.3:
                if Distance(dx ,dy + Ly, dz) < 2.3:
                    CM[i][j-NSi] = 1

                if Distance(dx, dy - Ly, dz) < 2.3:
                    CM[i][j-NSi] = 1

            if dz > Lz - 2.3:
                if Distance(dx, dy, dz + Lz) < 2.3:
                    CM[i][j-NSi] = 1

                if Distance(dx, dy, dz - Lz) < 2.3:
                    CM[i][j-NSi] = 1

    terminalO = []
    bridgingO = []
    for i in range(NO):
        atomID = data[i+NSi][0]
        V = data_all[atomID][4]

        nO = numpy.sum(CM[:,i])

        if  nO == 0 or nO == 1:
            # free O and terminal O
            terminalO.append(V)
        elif nO >= 2:
            bridgingO.append(V)


    return numpy.array(terminalO), numpy.array(bridgingO)


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
##    pylab.xlim((0,100))
##    pylab.ylim((0,0.20))
    pylab.savefig(filename)
    pylab.close(N)

def GetDensity(data):
    density = gaussian_kde(data)
    density.covariance_factor = lambda : 0.075
    density._compute_covariance()
    return density     

def PlotDensity(data,color):
    density = GetDensity(data)
    xs = numpy.linspace(0,60,200)
    plt.plot(xs,density(xs),color)


def PlotDensity3(N,data1,data2,data3,filename,xlbl):
    density1 = GetDensity(data1)
    density2 = GetDensity(data2)
    density3 = GetDensity(data3)
    xs = numpy.linspace(0,60,200)
    
    plt.figure(N)
    plt.plot(xs,density1(xs),'b',xs,density2(xs),'g',xs,density3(xs),'r')
    plt.xlabel(xlbl)
    plt.ylabel('frequency')
##    plt.ylim((0,0.09))
##    plt.set_aspect(1)
    plt.savefig(filename)
    plt.close(N)

def PlotHistAndDensity(N,data,Nbins,filename):
    plt.figure(N)
    density = GetDensity(data)
    xs = numpy.linspace(0,60,200)
    
    plt.hist(data,bins=Nbins,normed=True,histtype='step',lw=1)
    plt.plot(xs,density(xs),'r')
    plt.xlim((0,50))
    plt.ylim((0,0.08))
    plt.savefig(filename)
    plt.close(N)
    

def WriteNormDensityStats(data,Nsimplicies,outfile):
    data = numpy.array(data)

    f = open(outfile,'w')
    bins = numpy.linspace(0,60,200)
    histdata = numpy.histogram(data,bins=bins,normed=False)
    counts = histdata[0]
    density = counts/float(Nsimplicies)

    for i in range(len(density)):
        f.write('%.8f,%.10f\n' % (bins[i],density[i]))
    f.close()


def WriteDensityStats(data,outfile):
    f = open(outfile,'w')
    bins = numpy.linspace(0,100,101)
    histdata = numpy.histogram(data,bins=bins,normed=False)
    density = histdata[0]

    for i in range(len(density)):
        if density[i] > 0.0000000:
            f.write('%.8f,%.10f\n' % (bins[i],density[i]))
    f.close()


def WriteStats(data,VtotEq,inputfile,outfile):
    VCellsFree,Vtot,Vfreetot,VorVol,SiVol,SiVolFree,OVol,OVolFree,CVol,CVolFree = CellVolStats(data)

    VolSi = sum(SiVol)
    VolO = sum(OVol)
    VolC = sum(CVol)

    f = open(outfile,'a')
    f.write('\n %s \n' %inputfile)
    f.write('Si Vol / V0 = %.4f\n' %(VolSi/VtotEq))
    f.write('O Vol / V0 = %.4f\n' %(VolO/VtotEq))
    f.write('C Vol / V0 = %.4f\n' %(VolC/VtotEq))
    f.write('Total Vol / V0 = %.4f\n' %(Vtot/VtotEq))
    f.close()



################################################################################
# analyze the command line arguments and setup corresponding parameters
if len(sys.argv) < 4:
    print 'Usage:'
    print '  python %s <input dir> <plot file name> <atom type (all, Si, O, C)>' %sys.argv[0]
    exit()

t0 = time.time()

inputdir = sys.argv[1]
filename = sys.argv[2]
atomtype = sys.argv[3]

filenames = glob.glob('%sOCSEt_*_voro.vol'%inputdir)
# filenames = glob.glob('%s135Benz_*_voro.vol'%inputdir)
# filenames = glob.glob('%sSiO2_*_voro.vol'%inputdir)

print filenames

if atomtype in ['total','Si','O','C']:
    colors = ['b','g','r']
    # colors = ['g','g','b','b','r','r']
elif atomtype == 'all':
    colors = ['b','g','r']
    # colors = ['g','g','g','b','b','b','r','r','r']
elif atomtype == 'ind+all':
    colors = []
else:
    raise RuntimeError, "Invalid atom type. Must be 'all', 'Si', 'O', or 'C'"

# colors = ['b','g','r']
# colors = ['g','b','b','r','r']
# colors = ['g','g','b','r']
# colors = ['g','g','r','b']
# colors = ['g','g','b','b','r','r']
# colors = ['g','g','b','b','r','r']
# colors = ['g','g','g','b','b','b','r','r','r']
# temps = numpy.array([995.65928,75.138305,225.05056,375.83718,499.07168,648.06385,845.4484])
# t = 0

#------------------------------------------------------------------------------#
# # set font type and axes thickness
# arialfont = {'fontname':'Arial'}
# plt.rcParams['axes.linewidth'] = 2
# plt.rcParams.update({'fontname':'Arial','font.size': 16})

# # analyze results directory
# f = open('results/results.txt','w')
# plt.figure(1)
# for fname in filenames:
#     terminalO = []
#     data, data_xyz, NSi, NO, Lx, Ly, Lz = LoadDataWeighted(fname)
#     VCellsFree,Vtot,Vfreetot,VorVol,SiVol,SiVolFree,OVol,OVolFree,CVol,CVolFree = CellVolStats(data)

#     terminalO, bridgingO = GetTerminalO(data,data_xyz,NSi,NO,Lx,Ly,Lz)

#     print "\nWorking with %s" %fname
#     print "VorVol = %.4f" %numpy.sum(VorVol)
#     print "SiVol = %.4f" %numpy.sum(SiVol)
#     print "OVol = %.4f" %numpy.sum(OVol)
#     print "CVol = %.4f" %numpy.sum(CVol)

#     f.write("Filename = %s\n"%fname)
#     f.write("VorVol = %.4f\n" %numpy.sum(VorVol))
#     f.write("SiVol = %.4f\n" %numpy.sum(SiVol))
#     f.write("OVol = %.4f\n" %numpy.sum(OVol))
#     f.write("CVol = %.4f\n\n" %numpy.sum(CVol))

#     if len(terminalO):
#         print 'Terminal O = %.4f' %numpy.sum(terminalO)
#         print 'Bridging O = %.4f' %numpy.sum(bridgingO)
#         f.write('Terminal O = %.4f\n' %numpy.sum(terminalO))
#         f.write('Bridging O = %.4f\n\n' %numpy.sum(bridgingO))



    # WriteNormDensityStats(OVol,numpy.sum(VorVol),'%s_%s.csv' % (fname[:-4],'O'))

    # names = ['all','Si','C','termO','bridgeO']
    # types = [VorVol,SiVol,CVol,terminalO,bridgingO]
    # Vnorm = numpy.sum(VorVol)

    # for i in range(len(names)):
    #     outfile = '%s_%s.csv' % (fname[:-4],names[i])

    #     WriteNormDensityStats(types[i],Vnorm,outfile)



#     if atomtype in ['total','Si','O','C']:
#         if atomtype == 'total':
#             data = VorVol
#         elif atomtype == 'Si':
#             data = SiVol
#         elif atomtype == 'O':
#             data = OVol
#         elif atomtype == 'C':
#             data = CVol

#         PlotDensity(data,colors[t])
#         t+=1

#     elif atomtype == 'all':
#         if len(filenames) > 3:
#             if len(SiVol):
#                 PlotDensity(SiVol,colors[t])
#             if len(OVol):
#                 PlotDensity(OVol,colors[t+1])
#             if len(CVol):
#                 PlotDensity(CVol,colors[t+2])
#             t+=3    
#         else:
#             if len(SiVol):
#                 PlotDensity(SiVol,colors[t])
#             if len(OVol):
#                 PlotDensity(OVol,colors[t])
#             if len(CVol):
#                 PlotDensity(CVol,colors[t])
#             t+=1

#     elif atomtype == 'ind+all':
#         PlotDensity(OVol,'r')
#         PlotDensity(CVol,'b')
#         PlotDensity(SiVol,'k')
#         PlotDensity(VorVol,'g')

# plt.xlabel('Volume per Atom, $\Omega$ ($\AA^{3}$)',size=20,**arialfont)
# plt.ylabel('Frequency',size=20,**arialfont)
# plt.xlim((0,60))
# plt.axes().set_aspect(1./plt.axes().get_data_ratio())
# plt.savefig('results/%s' %filename,bbox_inches='tight')

# plt.close(1)
# f.close()

#------------------------------------------------------------------------------#
for fname in filenames:
    data, data_xyz, NSi, NO, Lx, Ly, Lz = LoadDataWeighted(fname)
    attributes = VoronoiStats(data)

    names = ['vertices','edges','faces','Si_vertices','Si_edges','Si_faces',
            'O_vertices','O_edges','O_faces','C_vertices','C_edges','C_faces']

    for i in range(len(names)):
        outfile = '%s_%s.csv' % (fname[:-4],names[i])

        WriteDensityStats(attributes[i],outfile)




#------------------------------------------------------------------------------#
# VorVol_data = numpy.zeros(len(filenames))
# SiVol_data = numpy.zeros(len(filenames))
# OVol_data = numpy.zeros(len(filenames))
# CVol_data = numpy.zeros(len(filenames))

# n = 0
# for file in filenames:
#     data = LoadDataWeighted(file)
#     VCellsFree,Vtot,Vfreetot,VorVol,SiVol,SiVolFree,OVol,OVolFree,CVol,CVolFree = CellVolStats(data)

#     if n==0:
#         V0 = numpy.sum(VorVol)

#     print '\nWorking with file: %s' %file
#     VorVol_data[n] = numpy.sum(VorVol)/V0
#     SiVol_data[n] = numpy.sum(SiVol)/V0
#     OVol_data[n] = numpy.sum(OVol)/V0
#     CVol_data[n] = numpy.sum(CVol)/V0

#     n+=1

# plt.figure(2)
# plt.scatter(temps,VorVol_data,c='r')
# plt.scatter(temps,SiVol_data,c='g')
# plt.scatter(temps,OVol_data,c='c')
# plt.scatter(temps,CVol_data,c='b')
# plt.xlabel('Temperature, T (K)')
# plt.ylabel('Volume / V0')
# plt.ylim((0.22,0.3))
# plt.savefig('volume_change.pdf')
# plt.close(2)


#------------------------------------------------------------------------------#
# fdir = 'voro_data_out/'
# filename1 = "OCSEt_155000_voro.vol"
# filename2 = "OCSEt_255000_voro.vol"
# filename3 = "OCSEt_355000_voro.vol"
# ##fileEq = "OCSEt_255000_voro.vol"

# ##filename1 = "135Benz_185000_voro.vol"
# ##filename2 = "135Benz_285000_voro.vol"
# ##filename3 = "135Benz_385000_voro.vol"

# f1 = fdir + filename1
# f2 = fdir + filename2
# f3 = fdir + filename3
# ##feq = fdir + fileEq

# data1 = LoadDataWeighted(f1)
# VCellsFree1,Vtot1,Vfreetot1,VorVol1,SiVol1,SiVolFree1,OVol1,OVolFree1,CVol1,CVolFree1 = CellVolStats(data1)
# ##Verts1,Edges1,Faces1,Si_verts1,Si_edges1,Si_faces1,O_verts1,O_edges1,O_faces1,C_verts1,C_edges1,C_faces1 = VoronoiStats(data1)

# data2 = LoadDataWeighted(f2)
# VCellsFree2,Vtot2,Vfreetot2,VorVol2,SiVol2,SiVolFree2,OVol2,OVolFree2,CVol2,CVolFree2 = CellVolStats(data2)
# ##Verts2,Edges2,Faces2,Si_verts2,Si_edges2,Si_faces2,O_verts2,O_edges2,O_faces2,C_verts2,C_edges2,C_faces2 = VoronoiStats(data2)

# data3 = LoadDataWeighted(f3)
# VCellsFree3,Vtot3,Vfreetot3,VorVol3,SiVol3,SiVolFree3,OVol3,OVolFree3,CVol3,CVolFree3 = CellVolStats(data3)
# ##Verts3,Edges3,Faces3,Si_verts3,Si_edges3,Si_faces3,O_verts3,O_edges3,O_faces3,C_verts3,C_edges3,C_faces3 = VoronoiStats(data3)

# ##dataEq = LoadDataWeighted(feq)
# ##VCellsFreeEq,VtotEq,VfreetotEq,VorVolEq,SiVolEq,SiVolFreeEq,OVolEq,OVolFreeEq,CVolEq,CVolFreeEq = CellVolStats(dataEq)

# data1 = LoadDataWeighted(f1)
# WriteStats(data1,Vtot2,f1,'OCSEt_VorCellStats.txt')

# data2 = LoadDataWeighted(f2)
# WriteStats(data2,Vtot2,f2,'OCSEt_VorCellStats.txt')

# data3 = LoadDataWeighted(f3)
# WriteStats(data3,Vtot2,f3,'OCSEt_VorCellStats.txt')

# ##print "\nWorking with %s" %f1
# ##print "Vtot = %.4f" %Vtot1
# ##print "Vfreetot = %.4f" %Vfreetot1
# ##print "SiVol = %.4f" %SiVol1
# ##print "OVol = %.4f" %OVol1
# ##print "CVol1= %.4f" %CVol1
# ##print "Free volume fraction = %.4f" %(Vfreetot1/Vtot1)
# ##
# ##print "\nWorking with %s" %f2
# ##print "Vtot = %.4f" %Vtot2
# ##print "Vfreetot = %.4f" %Vfreetot2
# ##print "SiVol = %.4f" %SiVol2
# ##print "OVol = %.4f" %OVol2
# ##print "CVol1= %.4f" %CVol2
# ##print "Free volume fraction = %.4f" %(Vfreetot2/Vtot2)
# ##
# ##print "\nWorking with %s" %f3
# ##print "Vtot = %.4f" %Vtot3
# ##print "Vfreetot = %.4f" %Vfreetot3
# ##print "SiVol = %.4f" %SiVol3
# ##print "OVol = %.4f" %OVol3
# ##print "CVol1= %.4f" %CVol3
# ##print "Free volume fraction = %.4f" %(Vfreetot3/Vtot3)

# ##Plot2(1,VCellsFree1,VorVol1,300,'FreeVol_'+filename1[6:12]+'.jpg','free volume per atom (A^3)')
# ##Plot1(2,SiVol1,300,'Si_FreeVol_'+filename1[6:12]+'.jpg','free volume per Si atom (A^3)')
# ##Plot1(3,OVol1,300,'O_FreeVol_'+filename1[6:12]+'.jpg','free volume per O atom (A^3)')
# ##Plot1(4,CVol1,300,'C_FreeVol_'+filename1[6:12]+'.jpg','free volume per C atom (A^3)')
# ##Plot3(9,VCellsFree1,VCellsFree2,VCellsFree3,200,'All_Cells_Peak_Shift_free.jpg','free volume per atom (A^3)')

# ##Plot3(5,VorVol1,VorVol2,VorVol3,200,'All_Cells_Peak_Shift.jpg','volume per atom (A^3)')
# ##Plot3(6,SiVol1,SiVol2,SiVol3,200,'Si_atom_peak_shift.jpg','volume per Si atom (A^3)')
# ##Plot3(7,CVol1,CVol2,CVol3,200,'C_atom_peak_shift.jpg','volume per C atom (A^3)')
# ##Plot3(8,OVol1,OVol2,OVol3,200,'O_atom_peak_shift.jpg','volume per O atom (A^3)')
# ##
# ##Plot3(10,Verts1,Verts2,Verts3,600,'All_voronoi_verts.jpg','Number of Vonoroi Vertices Per Cell')
# ##Plot3(11,Edges1,Edges2,Edges3,600,'All_voronoi_edges.jpg','Number of Vonoroi Edges Per Cell')
# ##Plot3(12,Faces1,Faces2,Faces3,600,'All_voronoi_faces.jpg','Number of Vonoroi Faces Per Cell')

# ##Plot3(13,Si_verts1,Si_verts2,Si_verts3,200,'Si_voronoi_verts.jpg','Number of Vonoroi Vertices Per Si Cell')
# ##Plot3(14,Si_edges1,Si_edges2,Si_edges3,200,'Si_voronoi_edges.jpg','Number of Vonoroi Edges Per Si Cell')
# ##Plot3(15,Si_faces1,Si_faces2,Si_faces3,200,'Si_voronoi_faces.jpg','Number of Vonoroi Faces Per Si Cell')
# ##Plot3(16,O_verts1,O_verts2,O_verts3,200,'O_voronoi_verts.jpg','Number of Vonoroi Vertices Per O Cell')
# ##Plot3(17,O_edges1,O_edges2,O_edges3,200,'O_voronoi_edges.jpg','Number of Vonoroi Edges Per O Cell')
# ##Plot3(18,O_faces1,O_faces2,O_faces3,200,'O_voronoi_faces.jpg','Number of Vonoroi Faces Per O Cell')
# ##Plot3(19,C_verts1,C_verts2,C_verts3,200,'C_voronoi_verts.jpg','Number of Vonoroi Vertices Per C Cell')
# ##Plot3(20,C_edges1,C_edges2,C_edges3,200,'C_voronoi_edges.jpg','Number of Vonoroi Edges Per C Cell')
# ##Plot3(21,C_faces1,C_faces2,C_faces3,200,'C_voronoi_faces.jpg','Number of Vonoroi Faces Per C Cell')

# PlotDensity3(22,VorVol1,VorVol2,VorVol3,'All_Cells_Peak_Shift.jpg','volume per atom (A^3)')
# PlotDensity3(23,OVol1,OVol2,OVol3,'O_atom_peak_shift.jpg','volume per O atom (A^3)')
# PlotDensity3(24,CVol1,CVol2,CVol3,'C_atom_peak_shift.jpg','volume per C atom (A^3)')
# PlotDensity3(25,SiVol1,SiVol2,SiVol3,'Si_atom_peak_shift.jpg','volume per Si atom (A^3)')

# ##PlotHistAndDensity(24,VorVol1,200,'density_hist_test.jpg')
#------------------------------------------------------------------------------#
t1 = time.time()
print '\nElapsed time: %.3f seconds' % (t1-t0)


""" Generates a plane given 2 points """

import math
import numpy
import time

################################################################################
def LoadData(filename):
    f = open(filename)
    bounds = f.readline().strip().split()
    f.close()

    xmin = float(bounds[0])
    xmax = float(bounds[1])
    ymin = float(bounds[2])
    ymax = float(bounds[3])
    zmin = float(bounds[4])
    zmax = float(bounds[5])

    pointA = (xmin,ymin,zmin)
    pointB = (xmax,ymin,zmin)
    pointC = (xmin,ymax,zmin)
    pointD = (xmax,ymin,zmax)
    pointE = (xmax,ymax,zmax)
    pointF = (xmin,ymax,zmax)

    return pointA,pointB,pointC,pointD,pointE,pointF

def Plane(pointA,pointB,pointC):
    """ This funciton generates a plane given 3 points """
    vec1 = numpy.array([pointB[0] - pointA[0], pointB[1] - pointA[1], pointB[2] - pointA[2]])
    vec2 = numpy.array([pointC[0] - pointA[0], pointC[1] - pointA[1], pointC[2] - pointA[2]])
    
    normal = numpy.cross(vec1,vec2)
    a = normal[0]
    b = normal[1]
    c = normal[2]
        
    x = pointA[0]
    y = pointA[1]
    z = pointA[2]

    d = -(a*x + b*y + c*z)

    return a, b, c, d

def WriteData(pointA,pointB,pointC,pointD,pointE,pointF,outputfile):
    a1,b1,c1,d1 = Plane(pointA,pointB,pointC)
    a2,b2,c2,d2 = Plane(pointA,pointB,pointD)
    a3,b3,c3,d3 = Plane(pointC,pointE,pointF)
    a4,b4,c4,d4 = Plane(pointD,pointE,pointF)
    a5,b5,c5,d5 = Plane(pointF,pointC,pointA)
    a6,b6,c6,d6 = Plane(pointB,pointE,pointD)

    f = open(outputfile, 'w')
    f.write('%.4f  %.4f  %.4f  %.4f\n' %(a1,b1,c1,d1))
    f.write('%.4f  %.4f  %.4f  %.4f\n' %(a2,b2,c2,d2))
    f.write('%.4f  %.4f  %.4f  %.4f\n' %(a3,b3,c3,d3))
    f.write('%.4f  %.4f  %.4f  %.4f\n' %(a4,b4,c4,d4))
    f.write('%.4f  %.4f  %.4f  %.4f\n' %(a5,b5,c5,d5))
    f.write('%.4f  %.4f  %.4f  %.4f\n' %(a6,b6,c6,d6))
    f.close()
    

################################################################################
t0 = time.time()

BCdir = 'Boundary_conditions_data/'
BCfile = 'BCs_155000_voro.txt'
filename = BCdir + BCfile

outdir = 'Walls_data/'
ofile = 'Walls_'+BCfile[4:10]+'_voro.txt'
outputfile = outdir + ofile

pointA,pointB,pointC,pointD,pointE,pointF = LoadData(filename)

WriteData(pointA,pointB,pointC,pointD,pointE,pointF,outputfile)

print 'Elapsed time: %.4f seconds' %(time.time()-t0)

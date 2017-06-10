""" Convert .xyz files from LAMMPS to both weighted and unweighted
    formats that will be used in voro++
    weighted format: <ID: 1,...,N> <x> <y> <z> <radius>
    unweighed format: <ID: 1,...,N> <x> <y> <z>
"""
import sys
import numpy

# Analyze the command line arguments and setup the corresponding parameters
if len(sys.argv) < 7:
  print 'Usage:'
  print '  python %s <LAMMPS dir> <LAMMPS file> <in voro dir> <in voro file> <BC dir> <BC file> [weighted (default = True)]' % sys.argv[0]
  exit()

idir = sys.argv[1]
ifile = sys.argv[2]
inputfile = idir + ifile

odir = sys.argv[3]
ofile = sys.argv[4]
outputfile = odir + ofile

BCdir = sys.argv[5]
BCfile = sys.argv[6]
BCfile = BCdir + BCfile

if len(sys.argv) == 8:
  weighted = bool(sys.argv[7])
else:
  weighted = True

##datadir = 'LAMMPS_data/'
##ifile = 'OCSEt_335000.xyz'
##inputfile = datadir + ifile
##
##outdir = 'radical_plane_voronoi/voro_data_in/'
##ofile = ifile[:-4]+'_voro'
##outputfile = outdir + ofile
##
##BCdir = 'Boundary_conditions_data/'
##bcfile = 'BCs_'+ifile[6:12]+'_voro.txt'
##BCfile = BCdir + bcfile
##
##weighted = True

# load data
fin = open(inputfile)
Natoms = int(fin.readline())
atomR = numpy.zeros(Natoms)

data = []
atomID = 0
fin.readline()
while True:
  fields = fin.readline().strip().split()
  if fields:
    # get radius for each atom type
    if fields[0] == '1':
        R = 2.1
    elif fields[0] == '2':
        R = 1.52
    elif fields[0] == '3':
        R = 1.7
    atomR[atomID] = R
    atomID += 1
    data.append([int(atomID),float(fields[1]),float(fields[2]),float(fields[3])])
  else:
    break
fin.close()
data = numpy.array(data)

# convert to voro file (first line is xmin, xmax, ymin...)
f = open(BCfile, 'w')
f.write('%.4f  %.4f  %.4f  %.4f  %.4f  %.4f\n' %
           (min(data[:,1]),max(data[:,1]),min(data[:,2]),
           max(data[:,2]),min(data[:,3]),max(data[:,3])))
f.close()

fout = open(outputfile, 'w')
if weighted:
  for i in range(len(data)):
    fout.write('%d  %.4f  %.4f  %.4f  %.4f\n' %
               (data[i][0],data[i][1],data[i][2],data[i][3],atomR[i]))
else:
  for i in range(len(data)):
    fout.write('%d  %.4f  %.4f  %.4f\n' %
             (data[i][0],data[i][1],data[i][2],data[i][3]))
fout.close()

print "Successfully converted the file!"

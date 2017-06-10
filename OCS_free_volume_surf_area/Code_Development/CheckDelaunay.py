#The purpose of this program is to check if each of the vertices are
#used in the Delaunay triangulation
#OCS model

import numpy

#main program

filename = raw_input('What is the .ele filename: ')
print "Working with file: ", filename

data = []
with open(filename) as inputfile:
    for line in inputfile:
        data.append(line.strip().split( ))

data = data[2:] #remove first 2 lines

data = []
f = open(filename)
f.readline()
for line in f:
    line = line.strip().split()
    if not line[0][0] == '#':
        data.append([int(line[0]),int(line[1]),int(line[2]),int(line[3]),int(line[4])])
f.close()
data = numpy.array(data)

atom_idx = numpy.delete(data,0,1)

idx_used = set([])
for i in range(len(atom_idx)):
    for j in range(len(atom_idx[i])):
        if atom_idx[i][j] not in idx_used:
            idx_used.add(atom_idx[i][j])

print 'Number of vertices used: %d' %len(idx_used)

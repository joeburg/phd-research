#The purpose of this program is to convert .xyz files from LAMMPS
# to a weighted file format that will be used in TetGen
#We must convert the species identification to its atomic number
#OCS model

from Tkinter import Tk
from tkFileDialog import askopenfilename
Tk().withdraw()

import numpy as np

#main program

filename = askopenfilename()
print "Working with file:", filename

data = []
with open(filename) as inputfile:
    for line in inputfile:
        data.append(line.strip().split( ))

data = data[2:] #remove first 2 lines

#use the number of step in .xyz file to give number of rows, int(data[0][0])
#convert 1 to (Si)
#convert 2 to (O)
#convert 3 to (C)


# weighted Delaunay triangulation file    
##for i in range(len(data)):
##    if data[i][0] == '1':
##        R = '2.1'
##    elif data[i][0] == '2':
##        R = '1.52'
##    elif data[i][0] == '3':
##        R = '1.7'
##
##    data[i][0] = str(i)  #first column is index
##    data[i].append(R)    #last column is radius 

# regular file
for i in range(len(data)):
    data[i][0] = str(i)  #first column is index


#write to a text file
#the .join() method takes an array, i, and concantenates all the elements together
# with a space " " between each element.  Then a newline "\n" is added to make sure
# your output is broken up into separate lines

dataFile = open(filename[:-4]+".node", 'w')

### weighted
##dataFile.write('# Node count, 3 dim, 1 attribute, no boundary marker\n')
##dataFile.write('%d  %d  %d  %d\n' %(len(data),3,1,0))

# regular
dataFile.write('# Node count, 3 dim, 0 attribute, no boundary marker\n')
dataFile.write('%d  %d  %d  %d\n' %(len(data),3,0,0))

dataFile.write('# Node index, node coordinates, vdW radius\n')

for eachitem in data:
    dataFile.write(" ".join(eachitem)+'\n')

dataFile.close()


print "All done!" 

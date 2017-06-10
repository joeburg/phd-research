#The purpose of this program is to convert .xyz files from LAMMPS
# to a form in which VMD can read the .xyz file
#We must convert the species identification to its atomic number
#ZrGPTMS model


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

#use the number of step in .xyz file to give number of rows, int(data[0][0])
#convert 1 to 14 (Si)
#convert 2 to 8 (O)
#convert 3 to 6 (C)
#convert 4 to 40 (Zr)
#convert 5 to 8 (C)
    
for i in range(int(data[0][0]) + 2):
    for j in range(1):
        if i > 1:
            data[i][j]  = float(data[i][j])
            if data[i][j] == 1.0:
                data[i][j] = 'Si'
            elif data[i][j] == 2.0:
                data[i][j] = 'O'
            elif data[i][j] == 3.0:
                data[i][j] = 'C'
            elif data[i][j] == 4.0:
                data[i][j] = 'Zr'
            elif data[i][j] == 5.0:
                data[i][j] = 'C'
            data[i][j] = str(data[i][j])


#write to a text file
#the .join() method takes an array, i, and concantenates all the elements together
# with a space " " between each element.  Then a newline "\n" is added to make sure
# your output is broken up into separate lines

dataFile = open(filename[:-4]+"_VMD"+".xyz", 'w')
for eachitem in data:
    dataFile.write(" ".join(eachitem)+'\n')

dataFile.close()


print "All done!" 

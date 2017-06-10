#Purpose: Adjust the positions of the initial input files
#This program will ask for a number to scale the initial positions by
#don't forget to manually adjust the x,y,z range of the simulation box

from Tkinter import Tk
from tkFileDialog import askopenfilename
Tk().withdraw()

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
rc('mathtext', default='regular')


#function to change the sring vlaues to float numbers

def str_to_float(n_atoms, column):
    for i in range(n_atoms):
        column[i] = float(column[i])
    return column 


#function to scale the atoms positions
def scale(n_scale, column):
    for i in range(len(column)):
        column[i] = column[i]/n_scale
    return column    


#main program: 
filename = askopenfilename()
print "Working with file:", filename
n_atoms = input("Provide the number of atoms: ")
n_scale = input("What do you want to scale the initial positions by? \
NOTE: The atom positions will be multiplied by this number ")


    
#partition the data to access the atom intial position data 
data = []
with open(filename) as inputfile:
    for line in inputfile:
        data.append(line.strip(' \n\t\r').split())


initialization_data = []
relevant_data = []
bond_angle_data = []
for i in range(len(data)):
    if i < 39:
        initialization_data.append(data[i])
    elif i < 39+n_atoms:
        relevant_data.append(data[i])
    else:
        bond_angle_data.append(data[i])

#access each column in relevant_data and adjust by the input parameter
relevant_data = zip(*relevant_data)
for i in range(len(relevant_data)):
    relevant_data[i] = str_to_float(n_atoms,list(relevant_data[i]))
    for j in range(n_atoms):
        if i > 3:
            relevant_data[i][j] = relevant_data[i][j]*n_scale

relevant_data = [list(x) for x in zip(*relevant_data)]

for i in range(len(relevant_data)):
    for j in range(7):
        relevant_data[i][j] = str(relevant_data[i][j])

        
#reconstruct the final data array
data_final = initialization_data + relevant_data + bond_angle_data


#write to a text file
#the .join() method takes an array, i, and concantenates all the elements together
# with a space " " between each element.  Then a newline "\n" is added to make sure
# your output is broken up into separate lines

dataFile = open(filename[:-4]+"scale_by_"+str(n_scale)+".txt", 'w')
for eachitem in data_final:
    dataFile.write(" ".join(eachitem)+'\n')

dataFile.close()



        
print "All done!"
print "Don't forget to manually adjust the x,y,z range of the simulation box."

#Purpose: Analyze the RDF data
#This program will plot G(r) as well as the coordination number
#with respect to the interatomic distance
#This is a general program for any structure 

from Tkinter import Tk
from tkFileDialog import askopenfilename
Tk().withdraw()

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
rc('mathtext', default='regular')


#function to change the sring vlaues to float numbers

def str_to_float(n_t, column):
    for i in range(n_t):
        column[i] = float(column[i])
    return column 

#main program: 
filename = askopenfilename()
print "Working with file:", filename
n_t = input("Provide the number of timestep: ")


#separate data into separate blocks  
data = []
with open(filename) as inputfile:
    for line in inputfile:
        data.append(line.strip().split(' '))

data_clean = []
for i in data:
    if len(i) > 2:
        data_clean.append(i)

data_separate = zip(*[data_clean[i::n_t] for i in range(n_t)])

#plot the rdf and cn on separate plots:
for i in range(len(data_separate)):
    column = zip(*data_separate[i])
    num = i
    num1 = i+100
    for j in range(4):
        column[j] = str_to_float(n_t,list(column[j]))
        r = column[1]
        g_r = column[2]
        cn_r = column[3]

        plt.figure(num)
        plt.plot(r,g_r,'b')
        plt.ylabel("g(r)")
        plt.xlabel("Distance, r (A)")
        plt.savefig("RDF"+str(i)+".png")

        plt.figure(num1)
        plt.plot(r,cn_r,'r')
        plt.ylabel('g(r)')
        plt.xlabel('Distance, r (A)')
        plt.savefig("CN"+str(i)+".png")
        

print "All done! Don't forget to put the plots in a common folder."


        

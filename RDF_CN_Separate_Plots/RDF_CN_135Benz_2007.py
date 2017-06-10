#Purpose: Analyze the RDF data
#lammps-2007 version
#This program will plot G(r) as well as the coordination number
#with respect to the interatomic distance, r, on different plots
#the interactions considered (Zrgpmts model) are:
# Si-Si, Si-O, Si-Zr, O-Si, O-O, O-Zr, Zr-Si, Zr-O, Zr-Zr, Q-Q (methyl)
#If you don't want all plots, simple comment out undesired plot sections

from Tkinter import Tk
from tkFileDialog import askopenfilename
Tk().withdraw()

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
rc('mathtext', default='regular')
from matplotlib import cm


#function to change the sring vlaues to float numbers

def str_to_float(n_t, column):
    for i in range(n_t):
        column[i] = float(column[i])
    return column 

#main program: 
filename = askopenfilename()
print "Working with file:", filename
n_t = input("Provide the number of Nbins for the RDF compute: ")


#separate data into separate blocks  
data = []
with open(filename) as inputfile:
    for line in inputfile:
        data.append(line.strip().split(' '))

data_clean = []
for i in data:
    if len(i) > 2 and i[0] != 'r,':
        data_clean.append(i)
    

data_separate = zip(*[data_clean[i::n_t] for i in range(n_t)])
#plot all the RDF's on one plot with normalized values 
for i in range(4):
    column = zip(*data_separate[i])
    num = i
    for j in range(9):
        column[j] = str_to_float(n_t,list(column[j]))
        r = column[0] 
        g_SiSi = column[1]
        g_SiO = column[3]
        g_OO = column[7]
        
        plt.figure(num)
        plt.plot(r,g_SiSi, 'r', lw=0.05, label="Si-Si" if j==0 else "")
        plt.plot(r,g_SiO, 'g', lw=0.05, label="Si-O" if j==0 else "")
        plt.plot(r,g_OO, 'b', lw=0.05, label="O-O" if j==0 else "")
        leg = plt.legend(loc='best')
        for l in leg.legendHandles:
            l.set_linewidth(1)
        plt.title("RDF for all interactions")
        plt.ylabel("g(r)")
        plt.xlabel("Distance, r (A)")
        plt.axes().set_aspect(1./plt.axes().get_data_ratio())
        plt.savefig("RDF%d.png"%i,bbox_inches='tight')

#plot all the coordination numbers on one plot:
for i in range(4):
    column = zip(*data_separate[i])
    num1 = i+100
    for j in range(9):
        column[j] = str_to_float(n_t,list(column[j]))
        r = column[0]
        cn_SiSi = column[2]
        cn_SiO = column[4]
        cn_OO = column[8]
        
        plt.figure(num1)
        plt.plot(r,cn_SiSi, 'r', lw=0.8, label="Si-Si" if j==0 else "")
        plt.plot(r,cn_SiO, 'g', lw=0.8, label="Si-O" if j==0 else "")
        plt.plot(r,cn_OO, 'b', lw=0.8, label="O-O" if j==0 else "")     
        plt.legend(loc='best')
        plt.title("Coordination numbers for all interactions")
        plt.ylabel("CN")
        plt.xlabel("Distance, r (A)")
        plt.xlim((0,4))
        plt.ylim((0,8))
        plt.axes().set_aspect(1./plt.axes().get_data_ratio())
        plt.savefig("CN%d.png"%i,bbox_inches='tight')

        
print "All done! Don't forget to put the plots in a common folder."


        

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

###plot the RDFs on separate plots
##for i in range(len(data_separate)):
##    column = zip(*data_separate[i])
##    num = i+200
##    for j in range(21):
##        column[j] = str_to_float(n_t,list(column[j]))
##        r = column[0]
##        g_SiSi = column[1]
##        plt.figure(num)
##        plt.plot(r,g_SiSi, 'r',label="Si-Si" if j==0 else "")
##        plt.legend(loc='best')
##        plt.title("Si-Si Bonds")
##        plt.ylabel("g(r)")
##        plt.xlabel("Distance, r (A)")
##        plt.savefig("RDF_SiSi"+str(i)+".png")
##        
##for i in range(len(data_separate)):
##    column = zip(*data_separate[i])
##    num = i+300
##    for j in range(21):
##        column[j] = str_to_float(n_t,list(column[j]))
##        r = column[0]
##        g_SiO = column[3]
##        g_OSi = column[7]
##        plt.figure(num)
##        plt.plot(r,g_SiO, 'r',label="Si-O" if j==0 else "")
##        plt.plot(r,g_OSi, color=(0.0,0.7,0.0), label="O-Si" if j==0 else "")
##        plt.legend(loc='best')
##        plt.title("Si-O Bonds")
##        plt.ylabel("g(r)")
##        plt.xlabel("Distance, r (A)")
##        plt.savefig("RDF_SiO"+str(i)+".png")
##
##for i in range(len(data_separate)):
##    column = zip(*data_separate[i])
##    num = i+400
##    for j in range(21):
##        column[j] = str_to_float(n_t,list(column[j]))
##        r = column[0]
##        g_SiZr = column[5]
##        g_ZrSi = column[13]
##        plt.figure(num)
##        plt.plot(r,g_SiZr, 'b', label="Si-Zr" if j==0 else "")
##        plt.plot(r,g_ZrSi, color=(0,0,.8), label="Zr-Si" if j==0 else "")
##        plt.legend(loc='best')
##        plt.title("Si-Zr Bonds")
##        plt.ylabel("g(r)")
##        plt.xlabel("Distance, r (A)")
##        plt.savefig("RDF_SiZr"+str(i)+".png")
##
##for i in range(len(data_separate)):
##    column = zip(*data_separate[i])
##    num = i+500
##    for j in range(21):
##        column[j] = str_to_float(n_t,list(column[j]))
##        r = column[0]
##        g_OO = column[9]
##        plt.figure(num)
##        plt.plot(r,g_OO, 'c', label="O-O" if j==0 else "")
##        plt.legend(loc='best')
##        plt.title("O-O Bonds")
##        plt.ylabel("g(r)")
##        plt.xlabel("Distance, r (A)")
##        plt.savefig("RDF_OO"+str(i)+".png")
##
##for i in range(len(data_separate)):
##    column = zip(*data_separate[i])
##    num = i+600
##    for j in range(21):
##        column[j] = str_to_float(n_t,list(column[j]))
##        r = column[0]
##        g_OZr = column[11]
##        g_ZrO = column[15]
##        plt.figure(num)
##        plt.plot(r,g_OZr, 'm', label="O-Zr" if j==0 else "")
##        plt.plot(r,g_ZrO, color=(0.5,0,0.5), label="Zr-O" if j==0 else "")
##        plt.legend(loc='best')
##        plt.title("Zr-O Bonds")
##        plt.ylabel("g(r)")
##        plt.xlabel("Distance, r (A)")
##        plt.savefig("RDF_ZrO"+str(i)+".png")
##
##for i in range(len(data_separate)):
##    column = zip(*data_separate[i])
##    num = i+700
##    for j in range(21):
##        column[j] = str_to_float(n_t,list(column[j]))
##        r = column[0]
##        g_ZrZr = column[17]
##        plt.figure(num)
##        plt.plot(r,g_ZrZr, 'y', label="Zr-Zr" if j==0 else "")
##        plt.legend(loc='best')
##        plt.title("Zr-Zr Bonds")
##        plt.ylabel("g(r)")
##        plt.xlabel("Distance, r (A)")
##        plt.savefig("RDF_ZrZr"+str(i)+".png")
##
##for i in range(len(data_separate)):
##    column = zip(*data_separate[i])
##    num = i+700
##    for j in range(21):
##        column[j] = str_to_float(n_t,list(column[j]))
##        r = column[0]
##        g_QQ = column[19]
##        plt.figure(num)
##        plt.plot(r,g_QQ, 'k', label="Porogen-Porogen" if j==0 else "")
##        plt.legend(loc='best')
##        plt.title("Q-Q Bonds")
##        plt.ylabel("g(r)")
##        plt.xlabel("Distance, r (A)")
##        plt.savefig("RDF_PP"+str(i)+".png")
##
##        
#plot all the RDF's on one plot with normalized values 
for i in range(len(data_separate)):
    column = zip(*data_separate[i])
    num = i
    for j in range(21):
        column[j] = str_to_float(n_t,list(column[j]))
        r = column[0] 
        g_SiSi = column[1]
        g_SiO = column[3]
        g_SiZr = column[5]
        g_OSi = column[7]
        g_OO = column[9]
        g_OZr = column[11]
        g_ZrSi = column[13]
        g_ZrO = column[15]
        g_ZrZr = column[17]
        g_QQ = column[19]
        
        
        plt.figure(num)
        plt.plot(r,g_SiSi, 'r', lw=0.05, label="Si-Si" if j==0 else "")
        plt.plot(r,g_SiO, 'g', lw=0.05, label="Si-O" if j==0 else "")
        plt.plot(r,g_SiZr, 'b', lw=0.05, label="Si-Zr" if j==0 else "")
#        plt.plot(r,g_OSi, color=(0.0,0.7,0.0), lw=0.05, label="O-Si" if j==0 else "")
        plt.plot(r,g_OO, 'c', lw=0.05, label="O-O" if j==0 else "")
#        plt.plot(r,g_OZr, 'm', lw=0.05, label="O-Zr" if j==0 else "")
#        plt.plot(r,g_ZrSi, color=(0,0,.8), lw=0.05, label="Zr-Si" if j==0 else "")
        plt.plot(r,g_ZrO, color=(0.5,0,0.5), lw=0.05, label="Zr-O" if j==0 else "")
        plt.plot(r,g_ZrZr, 'y', lw=0.05, label="Zr-Zr" if j==0 else "")
        plt.plot(r,g_QQ, 'k', lw=0.05, label="Q-Q" if j==0 else "")
        leg = plt.legend(loc='best')
        for l in leg.legendHandles:
            l.set_linewidth(1)
        plt.title("RDF for all interactions")
        plt.ylabel("g(r)")
        plt.xlabel("Distance, r (A)")
        plt.savefig("RDF"+str(i)+".png")

#plot all the coordination numbers on one plot:
for i in range(len(data_separate)):
    column = zip(*data_separate[i])
    num1 = i+100
    for j in range(21):
        column[j] = str_to_float(n_t,list(column[j]))
        r = column[0]
        cn_SiSi = column[2]
        cn_SiO = column[4]
        cn_SiZr = column[6]
        cn_OSi = column[8]
        cn_OO = column[10]
        cn_OZr = column[12]
        cn_ZrSi = column[14]
        cn_ZrO = column[16]
        cn_ZrZr = column[18]
        cn_QQ = column[20]
        
        plt.figure(num1)
        plt.plot(r,cn_SiSi, 'r', lw=0.8, label="Si-Si" if j==0 else "")
        plt.plot(r,cn_SiO, 'g', lw=0.8, label="Si-O" if j==0 else "")
        plt.plot(r,cn_SiZr, 'b', lw=0.8, label="Si-Zr" if j==0 else "")
#        plt.plot(r,cn_OSi, color=(0.0,0.7,0.0), lw=0.8, label="O_Si" if j==0 else "")
        plt.plot(r,cn_OO, 'c', lw=0.8, label="O-O" if j==0 else "")
#        plt.plot(r,cn_OZr, 'm', lw=0.8, label="O-Zr" if j==0 else "")
#        plt.plot(r,cn_ZrSi, color=(0,0,.8), lw=0.8, label="Zr-Si" if j==0 else "")
        plt.plot(r,cn_ZrO, color=(0.5,0,0.5), lw=0.8, label="Zr-O" if j==0 else "")
        plt.plot(r,cn_ZrZr, 'y', lw=0.8, label="Zr-Zr" if j==0 else "")
        plt.plot(r,cn_QQ, 'k', lw=0.8, label="Q-Q" if j==0 else "")      
        plt.legend(loc='best')
        plt.title("Coordination numbers for all interactions")
        plt.ylabel("CN")
        plt.xlabel("Distance, r (A)")
        plt.xlim((0,4))
        plt.ylim((0,8))
        plt.savefig("CN"+str(i)+".png")

        
print "All done! Don't forget to put the plots in a common folder."


        

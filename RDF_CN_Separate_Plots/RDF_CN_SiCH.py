#Purpose: Analyze the RDF data
#This program will plot G(r) as well as the coordination number
#with respect to the interatomic distance, r, on different plots
#the interactions considered (SiCH model) are:
# Si-Si, Si-C, Si-H, C-C, C-H, H-H
#If you don't want all plots, simple comment out undesired plot sections

from Tkinter import Tk
from tkFileDialog import askopenfilename
Tk().withdraw()

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
rc('mathtext', default='regular')


#function to normalize the RDF
def normalize(g_r):
    for i in range(len(g_r)):
        g_r[i] = g_r[i]/g_sum
    return g_r

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
    if len(i) > 2:
        data_clean.append(i)

data_separate = zip(*[data_clean[i::n_t] for i in range(n_t)])

#add up all the probabilities so the RDFs can be normalized:
for i in range(len(data_separate)):
    column = zip(*data_separate[i])
    for j in range(14):
        column[j] = str_to_float(n_t, list(column[j]))
        g_SiSi_sum = sum(float(i) for i in column[2])
        g_SiC_sum = sum(float(i) for i in column[4])
        g_SiH_sum = sum(float(i) for i in column[6])
        g_CC_sum = sum(float(i) for i in column[8])
        g_CH_sum = sum(float(i) for i in column[10])
        g_HH_sum = sum(float(i) for i in column[12])
        #think about taking out of loop
        g_sum = g_SiSi_sum + g_SiC_sum + g_SiH_sum + g_CC_sum + g_CH_sum + g_HH_sum


###convert the tuples to lists and normalize the g(r) data; plot each interaction
##for i in range(len(data_separate)):
##    column = zip(*data_separate[i])
##    num1 = i
##    num2 = i+100
##    num3 = i+200
##    num4 = i+300
##    num5 = i+400
##    num6 = i+500
##    num7 = i+600
##    num8 = i+700
##    num9 = i+800
##    num10 = i+900
##    num11 = i+1000
##    num12 = i+1100
##    num13 = i+1200
##    for j in range(14):
##        column[j] = str_to_float(n_t,list(column[j]))
##                                    #LAMMPS claims that data is already normalized 
##        r = column[1]            #assign data to each interaction and plot each interaction
##
##        if j == 2:
##            g_SiSi = column[2]
###            g_SiSi = normalize(g_SiSi)
##            plt.figure(num1)
##            plt.title("RDF Si-Si") 
##            plt.ylabel("g(r)")
##            plt.xlabel("Distance, r (A)")
##            plt.plot(r,g_SiSi)
##            plt.savefig("Si_Si_RDF"+str(i)+".png")
##
##        if j == 3:
##            cn_SiSi = column[3]
##            plt.figure(num2)
##            plt.title("Si-Si Coordination") 
##            plt.ylabel("CN")
##            plt.xlabel("Distance, r (A)")
##            plt.plot(r,cn_SiSi)
##            plt.savefig("Si_Si_CN"+str(i)+".png")
##
##        if j == 4:
##            g_SiC = column[4]
###            g_SiC = normalize(g_SiC)
##            plt.figure(num3)
##            plt.title("RDF Si-C") 
##            plt.ylabel("g(r)")
##            plt.xlabel("Distance, r (A)")
##            plt.plot(r,g_SiC)
##            plt.savefig("Si_C_RDF"+str(i)+".png")
##
##        if j == 5:
##            cn_SiC = column[5]
##            plt.figure(num4)
##            plt.title("Si-C Coordination") 
##            plt.ylabel("CN")
##            plt.xlabel("Distance, r (A)")
##            plt.plot(r,cn_SiC)
##            plt.savefig("Si_C_CN"+str(i)+".png")
##
##        if j == 6:
##            g_SiH = column[6]
###            g_SiH = normalize(g_SiH)
##            plt.figure(num5)
##            plt.title("RDF Si-H") 
##            plt.ylabel("g(r)")
##            plt.xlabel("Distance, r (A)")
##            plt.plot(r,g_SiH)
##            plt.savefig("Si_H_RDF"+str(i)+".png")
##
##        if j == 7:
##            cn_SiH = column[7]
##            plt.figure(num6)
##            plt.title("Si-H Coordination") 
##            plt.ylabel("CN")
##            plt.xlabel("Distance, r (A)")
##            plt.plot(r,cn_SiH)
##            plt.savefig("Si_H_CN"+str(i)+".png")
##
##        if j == 8:
##            g_CC = column[8]
###            g_CC = normalize(g_CC)
##            plt.figure(num7)
##            plt.title("RDF C-C") 
##            plt.ylabel("g(r)")
##            plt.xlabel("Distance, r (A)")
##            plt.plot(r,g_CC)
##            plt.savefig("C_C_RDF"+str(i)+".png")
##
##        if j == 9:
##            cn_CC = column[9]
##            plt.figure(num8)
##            plt.title("C-C Coordination") 
##            plt.ylabel("CN")
##            plt.xlabel("Distance, r (A)")
##            plt.plot(r,cn_CC)
##            plt.savefig("C_C_CN"+str(i)+".png")
##
##        if j == 10:
##            g_CH = column[10]
###            g_CH = normalize(g_CH)
##            plt.figure(num9)
##            plt.title("RDF C-H") 
##            plt.ylabel("g(r)")
##            plt.xlabel("Distance, r (A)")
##            plt.plot(r,g_CH)
##            plt.savefig("C_H_RDF"+str(i)+".png")
##
##        if j == 11:
##            cn_CH = column[11]
##            plt.figure(num10)
##            plt.title("C-H Coordination") 
##            plt.ylabel("CN")
##            plt.xlabel("Distance, r (A)")
##            plt.plot(r,cn_CH)
##            plt.savefig("C_H_CN"+str(i)+".png")
##
##        if j == 12:
##            g_HH = column[12]
###            g_HH = normalize(g_HH)
##            plt.figure(num11)
##            plt.title("RDF H-H") 
##            plt.ylabel("g(r)")
##            plt.xlabel("Distance, r (A)")
##            plt.plot(r,g_HH)
##            plt.savefig("H_H_RDF"+str(i)+".png")
##
##        if j == 13:
##            cn_HH = column[13]
##            plt.figure(num12)
##            plt.title("H-H Coordination") 
##            plt.ylabel("CN")
##            plt.xlabel("Distance, r (A)")
##            plt.plot(r,cn_HH)
##            plt.savefig("H_H_CN"+str(i)+".png")

#plot all the RDF's on one plot with normalized values 
for i in range(len(data_separate)):
    column = zip(*data_separate[i])
    num13 = i+1200
    for j in range(14):
        column[j] = str_to_float(n_t,list(column[j]))
        r = column[1] 
        g_SiSi = column[2]
#        g_SiSi = normalize(g_SiSi)
        g_SiC = column[4]
#        g_SiC = normalize(g_SiC)
        g_SiH = column[6]
#        g_SiH = normalize(g_SiH)
        g_CC = column[8]
#        g_CC = normalize(g_CC)
        g_CH = column[10]
#        g_CH = normalize(g_CH)
        g_HH = column[12]
#        g_HH = normalize(g_HH)
        plt.figure(num13)
        plt.plot(r,g_SiSi, 'r', label="Si-Si" if j==0 else "")
        plt.plot(r,g_SiC, 'g', label="Si-C" if j==0 else "")
        plt.plot(r,g_SiH, 'b', label="Si-H" if j==0 else "")
        plt.plot(r,g_CC, 'y', label="C-C" if j==0 else "")
        plt.plot(r,g_CH, 'c', label="C-H" if j==0 else "")
        plt.plot(r,g_HH, 'm', label="H-H" if j==0 else "")
        plt.legend(loc='upper left')
        plt.title("RDF for all interactions")
        plt.ylabel("g(r)")
        plt.xlabel("Distance, r (A)")
        plt.savefig("RDF"+str(i)+".png")

#plot all the coordination numbers on one plot:
for i in range(len(data_separate)):
    column = zip(*data_separate[i])
    num14 = i+1300
    for j in range(14):
        column[j] = str_to_float(n_t,list(column[j]))
        r = column[1]
        cn_SiSi = column[3]
        cn_SiC = column[5]
        cn_SiH = column[7]
        cn_CC = column[9]
        cn_CH = column[11]
        cn_HH = column[13]
        plt.figure(num14)
        plt.plot(r,cn_SiSi, 'r', label="Si-Si" if j==0 else "")
        plt.plot(r,cn_SiC, 'g', label="Si-C" if j==0 else "")
        plt.plot(r,cn_SiH, 'b', label="Si-H" if j==0 else "")
        plt.plot(r,cn_CC, 'y', label="C-C" if j==0 else "")
        plt.plot(r,cn_CH, 'c', label="C-H" if j==0 else "")
        plt.plot(r,cn_HH, 'm', label="H-H" if j==0 else "")
        plt.legend(loc='upper left')
        plt.title("Coordination numbers for all interactions")
        plt.ylabel("CN")
        plt.xlabel("Distance, r (A)")
        plt.savefig("CN"+str(i)+".png")
        
print "All done! Don't forget to put the plots in a common folder."


        

#Purpose: Analyze the RDF data
#This program will plot G(r) as well as the coordination number
#with respect to the interatomic distance, r, on different plots
#the interactions considered (Zrgpmts model) are:
# Si-Si, Si-O, Si-Zr, O-Si, O-O, O-Zr, Zr-Si, Zr-O, Zr-Zr, Porogen-Porogen
#If you don't want all plots, simple comment out undesired plot sections

from Tkinter import Tk
from tkFileDialog import askopenfilename
Tk().withdraw()

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
rc('mathtext', default='regular')
from matplotlib import cm


#function to obtain the index with the maximum value of a list:
def maxindex(column):
    max_index = []
    for i in range(len(column)):
        if column[i] == max(column):
            max_index.append(i)
    max_index = max(max_index)
    return max_index

#function to choose the largest index in a list (could just use max(list))
def max_index(max_list):
    for i in range(len(max_list)):
        if len(max_list) == 1:
            return max_list
        elif max_list[i] < max(max_list):
            max_list.remove(max_list[i])

#function to set the cutoff for a peak
def cutoff(max_index,column):
    i = 0
    leading_zeros = []
    while column[i] == 0:
        leading_zeros.append(column[i])
        i = i+1

    half_peak_len = max_index - len(leading_zeros)
    cutoff_index = max_index + half_peak_len
    return cutoff_index

#function to set values of g(r) to zero after the cutoff
def clean_g_r(cutoff_index, column):
    for i in range(cutoff_index,len(column)):
        column[i] = 0
    return column
        
    

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

###plot the RDFs on separate plots
##for i in range(len(data_separate)):
##    column = zip(*data_separate[i])
##    num = i+200
##    for j in range(22):
##        column[j] = str_to_float(n_t,list(column[j]))
##        r = column[1]
##        g_SiSi = column[2]
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
##    for j in range(22):
##        column[j] = str_to_float(n_t,list(column[j]))
##        r = column[1]
##        g_SiO = column[4]
##        g_OSi = column[8]
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
##    for j in range(22):
##        column[j] = str_to_float(n_t,list(column[j]))
##        r = column[1]
##        g_SiZr = column[6]
##        g_ZrSi = column[14]
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
##    for j in range(22):
##        column[j] = str_to_float(n_t,list(column[j]))
##        r = column[1]
##        g_OO = column[10]
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
##    for j in range(22):
##        column[j] = str_to_float(n_t,list(column[j]))
##        r = column[1]
##        g_OZr = column[12]
##        g_ZrO = column[16]
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
##    for j in range(22):
##        column[j] = str_to_float(n_t,list(column[j]))
##        r = column[1]
##        g_ZrZr = column[18]
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
##    for j in range(22):
##        column[j] = str_to_float(n_t,list(column[j]))
##        r = column[1]
##        g_PP = column[20]
##        plt.figure(num)
##        plt.plot(r,g_PP, 'k', label="Porogen-Porogen" if j==0 else "")
##        plt.legend(loc='best')
##        plt.title("Porogen-Porogen Bonds")
##        plt.ylabel("g(r)")
##        plt.xlabel("Distance, r (A)")
##        plt.savefig("RDF_PP"+str(i)+".png")
##
##        
#plot all the RDF's on one plot with normalized values 
for i in range(len(data_separate)):
    column = zip(*data_separate[i])
    num = i
    for j in range(22):
        column[j] = str_to_float(n_t,list(column[j]))
#        column[j] = clean_g_r(cutoff(maxindex(column[j]),column[j]),column[j])
        r = column[1] 
        g_SiSi = column[2]
        g_SiO = column[4]
        g_SiZr = column[6]
        g_OSi = column[8]
        g_OO = column[10]
        g_OZr = column[12]
        g_ZrSi = column[14]
        g_ZrO = column[16]
        g_ZrZr = column[18]
        g_PP = column[20]
        
        
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
#        plt.plot(r,g_PP, 'k', lw=0.05, label="Porogen-Porogen" if j==0 else "")
        leg = plt.legend(loc='best')
        for l in leg.legendHandles:
            l.set_linewidth(1)
        plt.title("RDF for all interactions")
        plt.ylabel("g(r)")
        plt.xlabel("Distance, r (A)")
        plt.savefig("RDF_no_PP_ZrZr"+str(i)+".png")

#plot all the coordination numbers on one plot:
for i in range(len(data_separate)):
    column = zip(*data_separate[i])
    num1 = i+100
    for j in range(22):
        column[j] = str_to_float(n_t,list(column[j]))
        r = column[1]
        cn_SiSi = column[3]
        cn_SiO = column[5]
        cn_SiZr = column[7]
        cn_OSi = column[9]
        cn_OO = column[11]
        cn_OZr = column[13]
        cn_ZrSi = column[15]
        cn_ZrO = column[17]
        cn_ZrZr = column[19]
        cn_PP = column[21]
        
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
#        plt.plot(r,cn_PP, 'k', lw=0.8, label="P-P" if j==0 else "")      
        plt.legend(loc='best')
        plt.title("Coordination numbers for all interactions")
        plt.ylabel("CN")
        plt.xlabel("Distance, r (A)")
        plt.savefig("CN_no_PP_ZrZr"+str(i)+".png")

        
print "All done! Don't forget to put the plots in a common folder."


        

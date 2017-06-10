#Purpose: Analyze the RDF data
#This program will plot G(r) as well as the coordination number
#with respect to the interatomic distance, r, on the same plot
#the interactions considered are:
# Si-Si, Si-C, Si-H, C-C, C-H, H-H


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
        m = float(max(g_r))
        g_r[i] = g_r[i]/m
    return g_r

#function to change the sring vlaues to float numbers

def str_to_float(column):
    for i in range(n_t):
        column[i] = float(column[i])
    return column 
    
#functions for plotting the various interactions
def plot_SiSi(column):
    r = column[1]
    g_SiSi = column[2]
    g_SiSi = normalize(g_SiSi)
    cn_SiSi = column[3]

    plt.figure(1)
    ax = plt.gca()
    ax2 = ax.twinx()
    ax.plot(r,g_SiSi)
    ax2.plot(r,cn_SiSi)
    plt.savefig("Si_Si.png")

    return

def plot_SiC(column):
    r = column[1]
    g_SiC = column[4]
    g_SiC = normalize(g_SiC)
    cn_SiC = column[5]

    plt.figure(2)
    ax = plt.gca()
    ax2 = ax.twinx()
    ax.plot(r,g_SiC)
    ax2.plot(r,cn_SiC)
    plt.savefig("Si_C.png")


    return

def plot_SiH(column):
    r = column[1]
    g_SiH = column[6]
    g_SiH = normalize(g_SiH)
    cn_SiH = column[7]

    plt.figure(3)
    ax = plt.gca()
    ax2 = ax.twinx()
    ax.plot(r,g_SiH)
    ax2.plot(r,cn_SiH)
    plt.savefig("Si_H.png")

    return

def plot_CC(column):
    r = column[1]
    g_CC = column[8]
    g_CC = normalize(g_CC)
    cn_CC = column[9]

    plt.figure(4)
    ax = plt.gca()
    ax2 = ax.twinx()
    ax.plot(r,g_CC)
    ax2.plot(r,cn_CC)
    plt.savefig("C_C.png")

    return

def plot_CH(column):
    r = column[1]
    g_CH = column[10]
    g_CH = normalize(g_CH)
    cn_CH = column[11]

    plt.figure(5)
    ax = plt.gca()
    ax2 = ax.twinx()
    ax.plot(r,g_CH)
    ax2.plot(r,cn_CH)
    plt.savefig("C_H.png")
    return

def plot_HH(column):
    r = column[1]
    g_HH = column[12]
    g_HH = normalize(g_HH)
    cn_HH = column[13]

    plt.figure(6)
    ax = plt.gca()
    ax2 = ax.twinx()
    ax.plot(r,g_HH)
    ax2.plot(r,cn_HH)
    plt.savefig("H_H.png")    
    return


def plot(D):
    r = D[1]
    g_r = D[2]
    g_r = normalize(g_r)
    cn = D[3]
    
    plt.figure(1)
    ax = plt.gca()
    ax2 = ax.twinx()
    ax.plot(r,g_r)
    ax2.plot(r,cn)

    plt.savefig("plot.png")
    return 


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


#convert the tuples to lists and normalize the g(r) data 
for i in range(len(data_separate)):
    column = zip(*data_separate[i])
    for i in range(14):
        column[i] = str_to_float(list(column[i]))
        if i == 2 or 4 or 6 or 8 or 10 or 12:
            column[i] = normalize(column[i])

        r = column[1]        #assign data to each interaction

        g_SiSi = column[2]
        g_SiSi = normalize(g_SiSi)
        cn_SiSi = column[3]

        g_SiC = column[4]
        g_SiC = normalize(g_SiC)
        cn_SiC = column[5]

        g_SiH = column[6]
        g_SiH = normalize(g_SiH)
        cn_SiH = column[7]

        g_CC = column[8]
        g_CC = normalize(g_CC)
        cn_CC = column[9]

        g_CH = column[10]
        g_CH = normalize(g_CH)
        cn_CH = column[11]

        g_HH = column[12]
        g_HH = normalize(g_HH)
        cn_HH = column[13]
        
        plt.figure(1)            #plot each interaction
        ax = plt.gca()
        ax2 = ax.twinx()
        ax.plot(r,g_SiSi)
        ax2.plot(r,cn_SiSi)
        plt.savefig("Si_Si"+str(i)+".png")

        plt.figure(2)
        ax = plt.gca()
        ax2 = ax.twinx()
        ax.plot(r,g_SiC)
        ax2.plot(r,cn_SiC)
        plt.savefig("Si_C"+str(i)+".png")

        plt.figure(3)
        ax = plt.gca()
        ax2 = ax.twinx()
        ax.plot(r,g_SiH)
        ax2.plot(r,cn_SiH)
        plt.savefig("Si_H"+str(i)+".png")

        plt.figure(4)
        ax = plt.gca()
        ax2 = ax.twinx()
        ax.plot(r,g_CC)
        ax2.plot(r,cn_CC)
        plt.savefig("C_C"+str(i)+".png")

        plt.figure(5)
        ax = plt.gca()
        ax2 = ax.twinx()
        ax.plot(r,g_CH)
        ax2.plot(r,cn_CH)
        plt.savefig("C_H"+str(i)+".png")

        plt.figure(6)
        ax = plt.gca()
        ax2 = ax.twinx()
        ax.plot(r,g_HH)
        ax2.plot(r,cn_HH)
        plt.savefig("H_H"+str(i)+".png")
    

print "All done! Don't forget to put the plots in a common folder."


        

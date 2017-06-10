#Purpose: this program plots the pairwise energy versus temperature,
#pressure versus volume, temperature versus timesteps, and pressure versus timesteps
#and the pressure versus box length
#possibly pressure versus timestep
#all during the quench step
#Zrgptms 

#for the following thermo
#step temp press Lx Ly Lz E_pair volume 

from Tkinter import Tk
from tkFileDialog import askopenfilename
Tk().withdraw()

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
rc('mathtext', default='regular')




#main program: 
filename = askopenfilename()
print "Working with file:", filename

anneal_start = input("At what line in the log does the anneal step begin? ")
anneal_end = input("At what line does the anneal step end? ")

N_quench_steps = input("How many quench steps are there? ")

plot_type = raw_input("What type of plot do you want? 'points' or 'lines'? ")

for i in range(N_quench_steps):
    globals()['quench_start%s' % i] = input("At what line in the log does quench step "+str(i)+" begin? ")
    globals()['quench_end%s' % i] = input("At what line in the log does quench step "+str(i)+" end? ")

    
#partition the data to access the relevant data 
data = []
with open(filename) as inputfile:
    for line in inputfile:
        data.append(line.strip(' \n\t\r').split())

quench_data = []
anneal_schedule = []
for i in range(len(data)):
    if i > (anneal_start - 1) and i < (anneal_end - 1):
        anneal_schedule.append(data[i])
    for j in range(N_quench_steps):
        if j == 0:
            if i > (globals()['quench_start%s' % j] - 1) and i < (globals()['quench_end%s' % j] - 1):
                quench_data.append(data[i])
                if i > globals()['quench_start%s' % j]:
                    anneal_schedule.append(data[i])
        elif j > 1:
            if i > globals()['quench_start%s' % j] and i < (globals()['quench_end%s' % j] - 1):
                quench_data.append(data[i])
                anneal_schedule.append(data[i])
                


#assign representations to each column in quench_data
#plot E_pair vs T and P vs Volume

L = len(quench_data)        
quench_data = [list(x) for x in zip(*quench_data)]

for i in range(len(quench_data)):
    for j in range(len(quench_data[0])):
        quench_data[i][j] = float(quench_data[i][j])

step = quench_data[0]
temp = quench_data[1]
press = quench_data[2]
Lx = quench_data[3]
Ly = quench_data[4]
Lz = quench_data[5]
E_pair = quench_data[6]
E_mol = quench_data[7]
    
volume = []
for i in range(len(Lx)):
    x = [Lx[i]*Ly[i]*Lz[i]]
    volume = volume + x


#Plots with connnected lines 
if plot_type == 'lines':
    plt.figure(1)
    plt.plot(temp, E_pair,'r')
    plt.title("Quench Step")
    plt.xlabel("Temperature (K)")
    plt.ylabel("Pair Energy (eV)")
    plt.savefig("T_vs_Epair.png")


    plt.figure(2)
    plt.plot(volume,press,'b')
    plt.title("Quench Step")
    plt.ylabel("Pressure (bar)")
    plt.xlabel("Volume ($A^{3}$)")
    plt.savefig("P_vs_vol.png")

#Plots with data points only
if plot_type == 'points':
    plt.figure(1)
    plt.plot(temp, E_pair,'r', marker='o',linestyle='')
    plt.title("Quench Step")
    plt.xlabel("Temperature (K)")
    plt.ylabel("Pair Energy (eV)")
    plt.savefig("T_vs_Epair.png")


    plt.figure(2)
    plt.plot(volume,press,'b',marker='o',linestyle='')
    plt.title("Quench Step")
    plt.ylabel("Pressure (bar)")
    plt.xlabel("Volume ($A^{3}$)")
    plt.savefig("P_vs_vol.png")

#assign representations to each column in anneal_schedule
#plot T vs step and P vs step

anneal_schedule = [list(x) for x in zip(*anneal_schedule)]

step1 = anneal_schedule[0]
temp1 = anneal_schedule[1]
press1 = anneal_schedule[2]

       
for i in range(len(anneal_schedule)):
    for j in range(len(anneal_schedule[0])):
        anneal_schedule[i][j] = float(anneal_schedule[i][j])

Lx1 = anneal_schedule[3]
Ly1 = anneal_schedule[4]
Lz1 = anneal_schedule[5]

    
volume1 = []
for i in range(len(Lx1)):
    x = [Lx1[i]*Ly1[i]*Lz1[i]]
    volume1 = volume1+ x

#Plots with lines
if plot_type == 'lines':
    plt.figure(3)
    plt.plot(step1, temp1,'g')
    plt.title("Annealing Schedule")
    plt.ylabel('Temperature (K)')
    plt.xlabel('Timesteps (ps)')
    plt.savefig("T_vs_step.png")


    plt.figure(4)
    plt.plot(step1, press1,'m')
    plt.title('Annealing Schedule')
    plt.ylabel('Pressure (bar)')
    plt.xlabel('Timesteps (ps)')
    plt.savefig('P_vs_step.png')

#Plots with data points
if plot_type == 'points':
    plt.figure(3)
    plt.plot(step1, temp1,'g',marker='o',linestyle='')
    plt.title("Annealing Schedule")
    plt.ylabel('Temperature (K)')
    plt.xlabel('Timesteps (ps)')
    plt.savefig("T_vs_step.png")


    plt.figure(4)
    plt.plot(step1, press1,'m',marker='o',linestyle='')
    plt.title('Annealing Schedule')
    plt.ylabel('Pressure (bar)')
    plt.xlabel('Timesteps (ps)')
    plt.savefig('P_vs_step.png')


#Volume versus timestep for entire anneal scheule
if plot_type == 'lines':
    plt.figure(5)
    plt.plot(step1,volume1,'c')
    plt.title('Annealing Schedule')
    plt.ylabel('Volume ($A^{3}$)')
    plt.xlabel('Timesteps (ps)')
    plt.savefig('Vol_vs_step.png')

if plot_type == 'points':
    plt.figure(5)
    plt.plot(step1,volume1,'c',marker='o',linestyle='')
    plt.title('Annealing Schedule')
    plt.ylabel('Volume ($A^{3}$)')
    plt.xlabel('Timesteps (ps)')
    plt.savefig('Vol_vs_step.png')   



print "All done! Don't forget to add the files to the correct folder."


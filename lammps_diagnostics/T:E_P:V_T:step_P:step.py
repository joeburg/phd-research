#Purpose: this program plots the pairwise energy versus temperature,
#pressure versus volume, temperature versus timesteps, and pressure versus timesteps
#and the pressure versus box length
#possibly pressure versus timestep
#all during the quench step
#SiCH model - 7 quench steps 


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

quench_start1 = input("At what line in the log does the first quench step being? ")
quench_end1 = input("At what line does the first quench step end? ")

quench_start2 = input("At what line in the log does the second quench step being? ")
quench_end2 = input("At what line does the second quench step end? ")

quench_start3 = input("At what line in the log does the third quench step being? ")
quench_end3 = input("At what line does the third quench step end? ")

quench_start4 = input("At what line in the log does the fourth quench step being? ")
quench_end4 = input("At what line does the fourth quench step end? ")

quench_start5 = input("At what line in the log does the fifth quench step being? ")
quench_end5 = input("At what line does the fifth quench step end? ")

quench_start6 = input("At what line in the log does the sixth quench step being? ")
quench_end6 = input("At what line does the sixth quench step end? ")

quench_start7 = input("At what line in the log does the seventh quench step being? ")
quench_end7 = input("At what line does the fist seventh step end? ")

    
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
    elif i > (quench_start1 - 1) and i < (quench_end1 - 1):
        quench_data.append(data[i])
        if i > quench_start1:
            anneal_schedule.append(data[i])
    elif i > quench_start2 and i < (quench_end2 - 1):
        quench_data.append(data[i])
        anneal_schedule.append(data[i])
    elif i > quench_start3 and i < (quench_end3 - 1):
        quench_data.append(data[i])
        anneal_schedule.append(data[i])
    elif i > quench_start4 and i < (quench_end4 - 1):
        quench_data.append(data[i])
        anneal_schedule.append(data[i])
    elif i > quench_start5 and i < (quench_start5 - 1):
        quench_data.append(data[i])
        anneal_schedule.append(data[i])
    elif i > quench_start6 and i < (quench_start6 - 1):
        quench_data.append(data[i])
        anneal_schedule.append(data[i])
    elif i > quench_start7 and i < (quench_start7 - 1):
        quench_data.append(data[i])
        anneal_schedule.append(data[i])


##for i in range(len(data)):
##    if i > 271 and i < 324:
##        anneal_schedule.append(data[i])
##    elif i > 356 and i < 559:
##        quench_data.append(data[i])
##        if i > 357:
##            anneal_schedule.append(data[i])
##    elif i > 592 and i < 599:
##        quench_data.append(data[i])
##        anneal_schedule.append(data[i])
##    elif i > 631 and i < 638:
##        quench_data.append(data[i])
##        anneal_schedule.append(data[i])
##    elif i > 670 and i < 677:
##        quench_data.append(data[i])
##        anneal_schedule.append(data[i])
##    elif i > 709 and i < 716:
##        quench_data.append(data[i])
##        anneal_schedule.append(data[i])
##    elif i > 748 and i < 755:
##        quench_data.append(data[i])
##        anneal_schedule.append(data[i])
##    elif i > 787 and i < 794:
##        quench_data.append(data[i])
##        anneal_schedule.append(data[i])
##

#assign representations to each column in quench_data
#plot E_pair vs T and P vs Volume

L = len(quench_data)        
quench_data = [list(x) for x in zip(*quench_data)]

for i in range(len(quench_data)):
    for j in range(L):
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


#assign representations to each column in anneal_schedule
#plot T vs step and P vs step

anneal_schedule = [list(x) for x in zip(*anneal_schedule)]

step1 = anneal_schedule[0]
temp1 = anneal_schedule[1]
press1 = anneal_schedule[2]


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

print "All done! Don't forget to add the files to the correct folder."


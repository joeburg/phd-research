#Purpose: manipulate data for plotting of Pressure vs Dilitation

from Tkinter import Tk
from tkFileDialog import askopenfilename
Tk().withdraw()

###############################################################################

filename = askopenfilename()
print "Working with file:", filename

scale = input('What modulo do you want to manipulate the data? ')

data = []
with open(filename) as inputfile:
    for line in inputfile:
        data.append(line.strip().split())
        
###############################################################################
#Access data from columns 5 (avg pressure) and 9 (volume)

press_vol_data = []
for i in range(len(data)):
    press_vol_data.append([data[i][4],data[i][8]])

print press_vol_data
#Apply the scale to clean data

press_vol_data_scale = []
for i in range(len(press_vol_data)):
    if i%scale == 0:
        press_vol_data_scale.append(press_vol_data[i])
        
print len(press_vol_data_scale)
print press_vol_data_scale
#################################################################################

dataFile = open("data_scale"+str(scale)+".txt", 'w')
for i in range(len(press_vol_data_scale)):
    dataFile.write("\t".join(press_vol_data_scale[i])+'\n')

print "All done!"




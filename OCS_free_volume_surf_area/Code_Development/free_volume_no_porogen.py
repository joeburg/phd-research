"""
Purpose: find the free volume of OCS with no porogen; computes simulation cell volume
and substracts volume of atoms
"""
 
from Tkinter import Tk
from tkFileDialog import askopenfilename
Tk().withdraw()

################################################################################
def load_data(filename):
    data = []
    f = open(filename)
    #first 2 lines are not relvant 
    f.readline()
    f.readline()
    for line in f:
        line = line.strip().split()
        data.append([int(line[0]),float(line[1]),float(line[2]),float(line[3])])
    f.close()
    return data

################################################################################
#import porogen .xyz file to correlate indicies with positions
filename = askopenfilename()
print "Working with file:", filename

##atom_file = raw_input('What is the .xyz file name of the total simulation? ')
##cutoff = input('What cutoff do you want for the region around each atom? ')

atom_data = load_data(filename)
N_atoms = len(atom_data)

#find number of Si,O and C atoms; input volume of each atom (treat each as sphere)
N_O=0
N_Si=0
N_C=0
for i in range(N_atoms):
    if atom_data[i][0]==1:
        N_Si += 1
    elif atom_data[i][0]==2:
        N_O += 1
    elif atom_data[i][0]==3:
        N_C += 1

V_C=1.417542
V_Si=0.4946877
V_O=6.284367


#find dimensions of box; transpose array
atom_data = [list(x) for x in zip(*atom_data)]

x_coord = atom_data[1]
y_coord = atom_data[2]
z_coord = atom_data[3]

x_min = min(x_coord)
x_max = max(x_coord)
Lx = x_max - x_min

y_min = min(y_coord)
y_max = max(y_coord)
Ly = y_max - y_min

z_min = min(z_coord)
z_max = max(z_coord)
Lz = z_max - z_min

volume = Lx*Ly*Lz
print volume

atom_data = [list(x) for x in zip(*atom_data)]

#compute free volume by subtracting volume of each atom
free_volume = volume - N_C*V_C - N_O*V_O - N_Si*V_Si

perc_free_vol = free_volume/volume*100

print "percent free volume = %.5f" %perc_free_vol

################################################################################

################################################################################  
#write to a text file
results = "\nFree vol at Timestep %s = %.5f" %(filename[-12:-4],perc_free_vol)

dataFile = open("free_volume_"+filename[-5]+".txt", 'a')
dataFile.write(results)
dataFile.close()


print 'All done!'

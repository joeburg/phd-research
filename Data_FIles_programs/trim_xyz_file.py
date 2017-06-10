#Purpose: for visualization purposes (VMD) remove atoms beyond given cutoffs 

from Tkinter import Tk
from tkFileDialog import askopenfilename
Tk().withdraw()

################################################################################
#input atom data
filename = askopenfilename()

x_min = input('Give x-min: ')
x_max = input('Give x-max: ')
y_min = input('Give y-min: ')
y_max = input('Give y-max: ')
z_min = input('Give z-min: ')
z_max = input('Give z-max: ')

atom_data = []
with open(filename) as inputfile:
    for line in inputfile:
        atom_data.append(line.strip().split())

N_atoms = int(atom_data[0][0])
atom_data=atom_data[2:]

for i in range(N_atoms):
    for j in range(1,4):
        atom_data[i][j] = float(atom_data[i][j])

trim_atom_data = []
for i in range(N_atoms):
    if x_min < atom_data[i][1] < x_max:
        if y_min < atom_data[i][2] < y_max:
            if z_min < atom_data[i][3] < z_max:
                trim_atom_data.append(atom_data[i])

N_atoms = len(trim_atom_data)

###############################################################
#make output file
for i in range(N_atoms):
    for j in range(4):
        trim_atom_data[i][j]=str(trim_atom_data[i][j])

trim_atom_data = [[str(N_atoms)]]+[['Atoms']]+trim_atom_data
    
dataFile = open('ZrGPTMS_SiO2_substrate_VMD_trim.xyz', 'w')
for eachitem in trim_atom_data:
    dataFile.write(" ".join(eachitem)+'\n')

dataFile.close()


print "All done!" 
    


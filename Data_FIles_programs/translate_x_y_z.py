#Purpose: translate structure by given amount in x-dir, y-dir, z-dir (after
# structure has already been translated to origin)

from Tkinter import Tk
from tkFileDialog import askopenfilename
Tk().withdraw()

################################################################################
#input atom data
filename = askopenfilename()

translate_x = input('How much do you want to translate in the + x-direction? ')
translate_y = input('How much do you want to translate in the + y-direction? ')
translate_z = input('How much do you want to translate in the + z-direction? ')

atom_data = []
with open(filename) as inputfile:
    for line in inputfile:
        atom_data.append(line.strip().split())

N_atoms = int(atom_data[0][0])
atom_data=atom_data[2:]

for i in range(N_atoms):
    for j in range(4):
        atom_data[i][j] = float(atom_data[i][j])

for i in range(N_atoms):
    atom_data[i][1]=atom_data[i][1]+translate_x
    atom_data[i][2]=atom_data[i][2]+translate_y
    atom_data[i][3]=atom_data[i][3]+translate_z

###############################################################
#make output file
for i in range(N_atoms):
    for j in range(4):
        atom_data[i][j]=str(atom_data[i][j])

atom_data = [[str(N_atoms)]]+[['Atoms']]+atom_data
    
dataFile = open('ZrGPTMS_translate_xyz.xyz', 'w')
for eachitem in atom_data:
    dataFile.write(" ".join(eachitem)+'\n')

dataFile.close()


print "All done!" 
    


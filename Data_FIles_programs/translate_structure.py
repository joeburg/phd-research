#Purpose: translate structure to origin then translate to desired location

from Tkinter import Tk
from tkFileDialog import askopenfilename
Tk().withdraw()

################################################################################
#input atom data
filename = askopenfilename()

translate = input('Give the value of the z-coord for the bottom surface: ')

atom_data = []
with open(filename) as inputfile:
    for line in inputfile:
        atom_data.append(line.strip().split())

N_atoms = int(atom_data[0][0])
atom_data=atom_data[2:]

#transpose to columns
atom_data = [list(x) for x in zip(*atom_data)]

for i in range(len(atom_data)):
    for j in range(N_atoms):
        atom_data[i][j] = float(atom_data[i][j])

#find dimensions of box
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

#transpose back to rows 
atom_data = [list(x) for x in zip(*atom_data)]

################################################################################
#find the origin atom and translate to (0,0,0)
x_cut = x_min+5
y_cut = y_min+5
z_cut = z_min+5

min_points = []
for i in range(N_atoms):
    if atom_data[i][1] < x_cut:
        if atom_data[i][2] < y_cut:
            if atom_data[i][3] < z_cut:
                min_points.append(i)

##use min sum of x,y,z coords as origin
##sum_min_points = []
##for i in range(len(min_points)):
##    x=atom_data[min_points[i]][1]
##    y=atom_data[min_points[i]][2]
##    z=atom_data[min_points[i]][3]
##
##    sum_min_points.append(x+y+z)
##
##min_min_points = min(sum_min_points)
##
##for i in range(len(sum_min_points)):
##    if sum_min_points[i]==min_min_points:
##        origin = i
##
##origin = min_points[origin]
##
##print origin 

#use smallest z-coord as origin
min_atoms = []
for i in range(len(min_points)):
    min_atoms.append(atom_data[min_points[i]])

print min_atoms  
min_atoms = [list(x) for x in zip(*min_atoms)]
min_zcoord = min(min_atoms[3])

for i in range(len(min_points)):
    if atom_data[min_points[i]][3] == min_zcoord:
        origin = min_points[i]

#translate the structure so the origin is at (0,0,0)
x_translate = atom_data[origin][1]
y_translate = atom_data[origin][2]
z_translate = atom_data[origin][3]

for i in range(N_atoms):
    atom_data[i][1]=atom_data[i][1]-x_translate
    atom_data[i][2]=atom_data[i][2]-y_translate
    atom_data[i][3]=atom_data[i][3]-z_translate


#translate z-coord of structure as given in input

for i in range(N_atoms):
    atom_data[i][3]=atom_data[i][3]+translate

###############################################################
#make output file
for i in range(N_atoms):
    for j in range(4):
        atom_data[i][j]=str(atom_data[i][j])

atom_data = [[str(N_atoms)]]+[['Atoms']]+atom_data
    
dataFile = open('ZrGPTMS_translate.xyz', 'w')
for eachitem in atom_data:
    dataFile.write(" ".join(eachitem)+'\n')

dataFile.close()


print "All done!" 
    


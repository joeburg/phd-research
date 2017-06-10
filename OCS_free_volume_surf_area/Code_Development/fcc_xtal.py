#Purpose: This program will make a computational cell consisting of LxMxN periodic copies
# of a fcc crystal with the corresponding basis atoms 
# RSi = 2.1 A
# a = 2*R*sqrt(2)
# Z = 12

import math

##########################################################################################
L = input("Give the L dimension of a LxMxN simulation cell: ")
M = input("Give the M dimension of a LxMxN simulation cell: ")
N = input("Give the N dimension of a LxMxN simulation cell: ")

##########################################################################################
#for the LAMMPS simulation, Si=1
atom1 = [1,0.0,0.0,0.0]
atom2 = [1,0.5,0.5,0.0]
atom3 = [1,0.5,0.0,0.5]
atom4 = [1,0.0,0.5,0.5]

##basis_atoms = [atom1,atom2,atom3,atom4]

atom5 = [1,1,0,0]
atom6 = [1,0,1,0]
atom7 = [1,0,0,1]
atom8 = [1,1,1,0]
atom9 = [1,1,0,1]
atom10 = [1,0,1,1]
atom11 = [1,1,1,1]
atom12 = [1,0.5,0.5,1]
atom13 = [1,1,0.5,0.5]
atom14 = [1,0.5,1,0.5]

basis_atoms = [atom1,atom2,atom3,atom4,atom5,atom6,atom7,\
               atom8,atom9,atom10,atom11,atom12,atom13,atom14]

##basis_atoms = [atom1,atom2,atom3,atom4]

#lattice parameter
a = 2*2.1*math.sqrt(2)

#make copies of the unit cell with LxMxN dimension; displace each atom by some vector
# [i j k] in each loop; keep track of the number of atoms in the computational cell

N_atoms = 0
simulation_atoms = []
index = 0
for i in range(L):
    for j in range(M):
        for k in range(N):
            for l in range(len(basis_atoms)):
                #VMD
                #[atom-type, x, y, z]
                simulation_atoms.append([str(basis_atoms[l][0]),\
                                         str((basis_atoms[l][1]+i)*a),\
                                         str((basis_atoms[l][2]+j)*a),\
                                         str((basis_atoms[l][3]+k)*a)])
            N_atoms = N_atoms + len(basis_atoms)
            
print N_atoms

##########################################################################################
#VMD file
atom_data_xyz = [[str(N_atoms)],['fcc xtal']]+simulation_atoms

dataFile = open("fcc"+str(N_atoms)+".xyz", 'w')
for eachitem in atom_data_xyz:
    dataFile.write("\t".join(eachitem)+'\n')

dataFile.close()

print "All done!"

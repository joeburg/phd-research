#purpose: create initial data file
#SiCH model

atoms = input("How many atoms do you want to simulate? ")

#ask for the number of desired molecuels to be simulated

N_C = input("How many carbon atoms? ")
N_Me = input("How many methyl groups? ")
N_CH2 = input("How many CH2 molecules? ")
N_CH = input("How many CH molecules? ")
N_C_C = input("How many C-C molecules? ")
N_C_Me = input("How many C-Me molecules? ")
N_C_CH2 = input("How many C-CH2 molecules?")
N_C_CH = input("How many C-CH molecules? ")

N_Si = input("How many Si atoms? ")
N_SiH3 = input("How many SiH3 molecules? ")
N_SiH2 = input("How many SiH2 molecules? ")
N_SiH = input("How many SiH molecules? ")
N_Si_Si = input("How many Si-Si molecules? ")
N_Si_SiH3 = input("How many Si-SiH3 molecules? ")
N_Si_SiH2 = input("How many Si-SiH2 molecule? ")
N_Si_SiH = input("How may Si-SiH molecules? ")
N_SiH_SiH = input("How many SiH-SiH molecules? ")

#calculate total number of Si, C and H atoms 

N_Si_tot = N_Si+N_SiH3+N_SiH2+N_SiH+ \
           2*(N_Si_Si+N_Si_SiH3+N_Si_SiH2+N_Si_SiH+N_SiH_SiH)

N_C_tot = N_C+N_Me+N_CH2+N_CH+2*(N_C_C+N_C_Me+N_C_CH2+N_C_CH)

N_H_tot = 3*N_Me+2*N_CH2+N_CH+3*N_C_Me+2*N_C_CH2+N_C_CH+ \
          3*N_SiH3+2*N_SiH2+N_SiH+3*N_Si_SiH3+2*N_Si_SiH2+N_Si_SiH+2*N_SiH_SiH


#total number of atoms, bonds, angles 

N_atoms = N_Si_tot + N_C_tot + N_H_tot

N_bonds = 3*N_Me+2*N_CH2+N_CH+N_C_C+4*N_C_Me+3*N_C_CH2+2*N_C_CH+ \
          3*N_SiH3+2*N_SiH2+N_SiH+N_Si_Si+4*N_Si_SiH3+3*N_Si_SiH2+2*N_Si_SiH+3*N_SiH_SiH

N_angles = 3*N_Me+N_CH2+6*N_C_Me+3*N_C_CH2+N_C_CH \
           3*N_SiH3+N_SiH2+6*N_Si_SiH3+3*N_Si_SiH2+N_Si_SiH+2*N_SiH_SiH

#Simulation cell length
L = input("Provide the length of the simultion cell in terms of lattice constants: ")

#atoms per methyl %C, 3H%
atomsMe = 4

#atoms per CH2 
atomsCH2 = 3

#atoms per CH
atomsCH = 2

#atoms per C-Me
atomsC_Me = 5 

#atoms per C-CH2 
atomsC_CH2 = 4

#atoms per C-C
atomsC_C = 2

#atoms per C-CH
atomsC_CH = 3

#atoms per SiH 
atomsSiH = 2

#atoms per SiH2
atomsSiH2 = 3

#atoms per SiH3
atomsSiH3 = 4

#bonds per methyl 3*C-H
bondsMe = 3

#bonds per CH2 2*C-H
bondsCH2 = 2

#bond per CH
bondsCH = 1

#bonds per SiH 
bondsSiH = 1

#bonds per CH2 2*C-H
bondsSiH2 = 2

#bonds per SiH3
bondsSiH3 = 3

#angles per methyl
anglesMe = 3

#angle per CH2
anglesCH2 = 1

#angle per SiH
anglesSiH = 1

#angles per SiH2
anglesSiH2 = 1

#angles per SiH3
anglesSiH3 = 3





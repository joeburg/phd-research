# def getModelData(Nmols, moltype):

# 	if moltype == 'EtOCS':
# 		NSi_per_mol = 2
# 		Natoms_per_mol = 4
# 		Nbonds_per_mol = 3
# 		Nangles_per_mol = 2
# 		Ndihedrals_per_mol = 1

# 		# types of atoms, bonds, angles, dihedrals 
# 		# note: the porogen molecules have the same bond types as the C-C in the ethane bridge 
# 		Natomtypes = 3 # 1 = Si, 2 = O, 3 = C
# 		Nbondtypes = 2 # Si-C, C-C 
# 		Nangletypes = 1 # Si-C-C (C-C-Si)
# 		Ndihedraltypes = 1 # Si-C-C-Si

# 	elif moltype == 'EtOCSMe':
# 		NSi_per_mol = 2
# 		Natoms_per_mol = 5
# 		Nbonds_per_mol = 4
# 		Nangles_per_mol = 3
# 		Ndihedrals_per_mol = 2

# 		# types of atoms, bonds, angles, dihedrals 
# 		Natomtypes = 3 # 1 = Si, 2 = O, 3 = C
# 		Nbondtypes = 2 # Si-C, C-C 
# 		Nangletypes = 2 # Si-C-C (C-C-Si), C-Si-C
# 		Ndihedraltypes = 2 	# Si-C-C-Si, C-C-Si-C

# 	elif moltype == 'EtOCSVinyl':
# 		NSi_per_mol = 2
# 		Natoms_per_mol = 6
# 		Nbonds_per_mol = 5
# 		Nangles_per_mol = 4
# 		Ndihedrals_per_mol = 3

# 		# types of atoms, bonds, angles, dihedrals 
# 		Natomtypes = 3 # 1 = Si, 2 = O, 3 = C
# 		Nbondtypes = 3 # Si-C, C-C, C=C
# 		Nangletypes = 3 # Si-C-C (C-C-Si), C-Si-C, Si-C=C
# 		Ndihedraltypes = 3 # Si-C-C-Si, C-C-Si-C, C-Si-C=C


# 	elif moltype == 'EtOCSPhenyl':
# 		NSi_per_mol = 2
# 		Natoms_per_mol = 10
# 		Nbonds_per_mol = 10
# 		Nangles_per_mol = 11
# 		Ndihedrals_per_mol = 4 # or 10 if we add C-C-C-C

# 		# types of atoms, bonds, angles, dihedrals 
# 		Natomtypes = 3 # 1 = Si, 2 = O, 3 = C
# 		Nbondtypes = 2 # Si-C, C-C
# 		Nangletypes = 3 # Si-C-C (C-C-Si), C-Si-C, C-C-C
# 		Ndihedraltypes = 5 # Si-C-C-Si, C-C-Si-C, C-Si-C-C, Si-C-C-C, C-C-C-C 



# 	# get the molecule, bond, angle, dihedral templates based on precursor
# 	moltemplate = getMoleculeTemplate(moltype)
# 	bondtemplate = getBondTemplate(moltype)
# 	angletemplate = getAngleTemplate(moltype)
# 	dihedraltemplate = getDihedralTemplate(moltype)

# 	# take the number of free oxygen to be 3 times the number of Si 
# 	NO = 3*NSi_per_mol*Nmols

# 	# generate 1 porogen molecule for every 2 precursors
# 	Nporogen = int(0.5*Nmols)
# 	Natoms_per_porogen = 40
# 	Nbonds_per_porogen = 39
	
# 	porogentemplate = getMoleculeTemplate('porogen')
# 	porogenbondtemplate = getBondTemplate('porogen')


# 	# get total number of atoms, bonds, 
# 	Natoms = Nmols*Natoms_per_mol + Nporogen*Natoms_per_porogen + NO
# 	Nbonds = Nmols*Nbonds_per_mol + Nporogen*Nbonds_per_porogen
# 	Nangles = Nmols*Nangles_per_mol
# 	Ndihedrals = Nmols*Ndihedrals_per_mol

# 	modeldata = (Natoms,Natomtypes,Nbonds,Nbondtypes,
# 					Nangle,Nangletypes,Ndihedrals,Ndihedraltypes)

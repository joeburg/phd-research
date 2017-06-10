#Combine .xyz files into one (substrate + ZrGPTMS)


#input atom data

substrate = raw_input('What is the substrate filename? ')
ZrGPTMS = raw_input('What is the ZrGPTMS filename? ')

substrate_data = []
with open(substrate) as inputfile:
    for line in inputfile:
        substrate_data.append(line.strip().split())

ZrGPTMS_data = []
with open(ZrGPTMS) as inputfile:
    for line in inputfile:
        ZrGPTMS_data.append(line.strip().split())

N_substrate_atoms = int(substrate_data[0][0])
substrate_data = substrate_data[2:]

N_ZrGPTMS_atoms = int(ZrGPTMS_data[0][0])
ZrGPTMS_data = ZrGPTMS_data[2:]

data = []
for i in range(N_substrate_atoms):
    data.append(substrate_data[i])

for i in range(N_ZrGPTMS_atoms):
    data.append(ZrGPTMS_data[i])

N_atoms = str(N_substrate_atoms+N_ZrGPTMS_atoms)

data_file = [[N_atoms]]+[['Atoms']]+data

dataFile = open('ZrGPTMS_SiO2_substrate.xyz', 'w')
for eachitem in data_file:
    dataFile.write(" ".join(eachitem)+'\n')

dataFile.close()


print "All done!" 

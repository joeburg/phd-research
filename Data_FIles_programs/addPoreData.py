import sys
import numpy

#------------------------------------------------------------------------------#
def LoadAtomData(inputfile):
	""" reads raw atom data file and translates IDs from 0 to N-1 """ 
	data = []

	f = open(inputfile)
	while True:
		fields = f.readline().strip().split()
		if fields:
			atomID = int(fields[0])
			molID = int(fields[1])
			atomtype = int(fields[2])
			charge = int(fields[3])
			xcoord = float(fields[4])
			ycoord = float(fields[5])
			zcoord = float(fields[6])

			# populate the data array
			data.append([atomID,molID,atomtype,charge,xcoord,ycoord,zcoord])

		else:
			break
	f.close()

	return data

def TranslateAtomData0(data):
	# translate IDs (atom & molecule) from 0 to N-1
	Natoms = len(data)
	startAtomidx = data[0][0]
	startMolidx = data[0][1]

	for i in range(Natoms):
		data[i][0] -= startAtomidx
		data[i][1] -= startMolidx

	return data

def TranslateAtomDataN(poredata,LAMMPSdata):
	# translate IDs (atom & molecule) from some start idx to idx+N-1
	Natoms = len(poredata)
	startAtomidx = LAMMPSdata[-1][0] + 1

	idx = LAMMPSdata[-1][1]
	for i in range(len(LAMMPSbondData)):
		if LAMMPSdata[i][1] > idx:
			startMolidx = LAMMPSdata[i][1]
	startMolidx += 1

	for i in range(Natoms):
		poredata[i][0] += startAtomidx
		poredata[i][1] += startMolidx

	return poredata

def LoadBondData(inputfile):
	''' reads raw bond data file and translates IDs from 0 to N-1 '''
	data = []

	f = open(inputfile)
	while True:
		fields = f.readline().strip().split()
		if fields:
			bondID = int(fields[0])
			bondtype = int(fields[1])
			atomID1 = int(fields[2])
			atomID2 = int(fields[3])

			# populate the data array
			data.append([bondID,bondtype,atomID1,atomID2])

		else:
			break
	f.close()

	return data

def TranslateBondData0(bondData,atomData):
	# translate bond ID from 0 to N-1; atom IDs to IDs in translate atoms
	Nbonds = len(bondData)

	startBondidx = bondData[0][0]
	shift = abs(atomData[-1][0] - bondData[0][2]) - 1

	for i in range(Nbonds):
		bondData[i][0] -= startBondidx
		bondData[i][2] -= shift
		bondData[i][3] -= shift

	return bondData

def TranslateBondDataN(poredata,LAMMPSdata):
	# translate the bond ID from start idx to idx+N-1
	Nbonds = len(poredata)
	startBondidx = LAMMPSdata[-1][0] + 1

	for i in range(Nbonds):
		poredata[i][0] += startBondidx

	return poredata

def WriteAtomFile(poredata,LAMMPSdata):
	f = open('newAtomData.txt','w')
	for i in range(len(LAMMPSdata)):
		f.write('%d   %d   %d   %d   %.4f   %.4f   %.4f\n' 
			%(LAMMPSdata[i][0],LAMMPSdata[i][1],LAMMPSdata[i][2],
				LAMMPSdata[i][3],LAMMPSdata[i][4],LAMMPSdata[i][5],
				LAMMPSdata[i][6]))
	
	for i in range(len(poredata)):
		f.write('%d   %d   %d   %d   %.4f   %.4f   %.4f\n' 
			%(poredata[i][0],poredata[i][1],poredata[i][2],poredata[i][3],
				poredata[i][4],poredata[i][5],poredata[i][6]))
	
	print 'Natoms = %d' %(len(LAMMPSdata)+len(poredata))
	f.close()


def WriteBondFile(poredata,LAMMPSdata):
	f = open('newBondData.txt', 'w')
	for i in range(len(LAMMPSdata)):
		f.write('%d   %d   %d   %d\n' 
			%(LAMMPSdata[i][0],LAMMPSdata[i][1],LAMMPSdata[i][2],
				LAMMPSdata[i][3]))

	for i in range(len(poredata)):
		f.write('%d   %d   %d   %d\n' 
			%(poredata[i][0],poredata[i][1],poredata[i][2],poredata[i][3]))

	print 'Nbonds = %d' %(len(LAMMPSdata)+len(poredata))
	f.close()


#------------------------------------------------------------------------------#
''' main program '''

if len(sys.argv) < 5:
	print 'Usage:'
	print '  python %s <pore atom file> <pore bond file> <atom file> <bond file>' %sys.argv[0]
	exit()

poreAtomfile = sys.argv[1]
poreBondfile = sys.argv[2]
atomfile = sys.argv[3]
bondfile = sys.argv[4]

# load data files
poreAtomdata = LoadAtomData(poreAtomfile)
LAMMPSatomData = LoadAtomData(atomfile)

poreBonddata = LoadBondData(poreBondfile)
LAMMPSbondData = LoadBondData(bondfile)

# translate atom data
poreAtomdata = TranslateAtomData0(poreAtomdata)
poreAtomdata = TranslateAtomDataN(poreAtomdata,LAMMPSatomData)

# translate bond data
poreBonddata = TranslateBondData0(poreBonddata,LAMMPSatomData)
poreBonddata = TranslateBondDataN(poreBonddata,LAMMPSbondData)

# write out new data files
WriteAtomFile(poreAtomdata,LAMMPSatomData)
WriteBondFile(poreBonddata,LAMMPSbondData)

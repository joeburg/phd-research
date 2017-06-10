''' This program gives terminal O a new ID in the LAMMPS data file '''

import sys
import time

#--------------------------------------------------------------------------------#

def LoadLAMMPSfile(inputfile):
	NSi = 0
	NO = 0
	NC = 0
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

			# determine type of atom
			if atomtype == 1:
				NSi += 1
			elif atomtype == 2:
				NO += 1
			elif atomtype == 3:
				NC += 1
			else: 
				raise RuntimeError, "Incorrect atom type."

			# populate the data array
			data.append([atomID,molID,atomtype,charge,xcoord,ycoord,zcoord])

		else:
			break
	f.close()

	return NSi, NO, NC, data

def GetTerminalO(inputfile,NO):
	''' reads in the terminal O atomIDs and sets the bridging O IDs'''
	# just O atomIDs
	O = set(range(9001,9001 + NO))
	terminalO = set([])

	f = open(inputfile)
	while True:
		fields = f.readline().strip().split()
		if fields:
			terminalOidx = int(fields[0])
			terminalO.add(terminalOidx)
		else:
			break
	f.close()

	bridgingO = O - terminalO
	return terminalO, bridgingO  

def ReplaceTerminalO(data, terminalO, bridgingO, NSi, NC, NO, Nsilane):
	''' writes out the atom data in LAMMPS format with the terminal O as the first group of O atoms and 
	then the bridging O'''

	outfile = 'newAtomData.txt'
	f = open(outfile, 'w')

	# write out silane 
	for i in range(Nsilane*4):
		f.write('%d   %d   %d   %d   %.4f   %.4f   %.4f\n' %
			(data[i][0], data[i][1], data[i][2], data[i][3], data[i][4], data[i][5], data[i][6]))

	# write out terminal O first so they can be grouped together 
	atomID = 9001
	for i in range(Nsilane*4, Nsilane*4 + len(terminalO)+len(bridgingO)):
		if data[i][0] in terminalO:
			f.write('%d   %d   %d   %d   %.4f   %.4f   %.4f\n' % 
				(atomID, data[i][1], data[i][2], data[i][3], data[i][4], data[i][5], data[i][6]))

			atomID += 1 

	# write out bridging O 
	for i in range(Nsilane*4, Nsilane*4 + len(terminalO)+len(bridgingO)):
		if data[i][0] in bridgingO:
			f.write('%d   %d   %d   %d   %.4f   %.4f   %.4f\n' % 
				(atomID, data[i][1], data[i][2], data[i][3], data[i][4], data[i][5], data[i][6]))

			atomID += 1 

	# write out remaining O and porogen molecules
	for i in range(Nsilane*4 + len(terminalO)+len(bridgingO), len(data)):
		f.write('%d   %d   %d   %d   %.4f   %.4f   %.4f\n' %
			(data[i][0], data[i][1], data[i][2], data[i][3], data[i][4], data[i][5], data[i][6]))

	f.close()

#--------------------------------------------------------------------------------#

# main program
if len(sys.argv) < 4:
	print 'Usage:'
	print '  python %s <LAMMPS data file> <terminal O file> <number of silane>' %sys.argv[0]
	exit()

LAMMPSfile = sys.argv[1]
terminalOfile = sys.argv[2]
Nsilane = int(sys.argv[3])

t0 = time.time()

NSi, NO, NC, data = LoadLAMMPSfile(LAMMPSfile)
terminalO, bridgingO = GetTerminalO(terminalOfile, NO)
ReplaceTerminalO(data, terminalO, bridgingO, NSi, NC, NO, Nsilane)


print '\nNumber of Si: %d' % NSi
print 'Number of C: %d' % NC
print 'Number of O: %d' % NO
print 'Number of bridging O: %d' % len(bridgingO)
print 'Number of terminal O: %d\n' % len(terminalO)
print 'Generated new data file in %.4f seconds.' %(time.time()-t0)







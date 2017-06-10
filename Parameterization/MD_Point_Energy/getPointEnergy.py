import copy
import glob
import sys
import time


##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##

def ReadData(filename):
	# dihedral angle is in the filename 
	angle = int(filename[4:-7])

	# end of document flag
	stop_flag = 'jump'

	# flag to read energy 
	energy_flag = 'minimize'
	energy_read = False
	Nlines_read = 0

	f = open(filename)
	energy = 0
	while True:
		line = f.readline().strip()

		if stop_flag in line:
			break

		if energy_flag in line:
			energy_read = True 
			continue 

		if energy_read:
			Nlines_read += 1
			if Nlines_read == 3:
				fields = line.split()
				energy = fields[12]

	return angle, energy


def AnalyzeFiles(filenames):
	data = []
	for filename in filenames:
		angle, energy = ReadData(filename)
		data.append((angle,energy))

	# sort files in order of increasing dihedral angle
	data = sorted(data, key=lambda x: x[0])

	# use symmery of molecule to get 360 degree scan data
	# for 180 degree scan 
	data_sym = copy.deepcopy(data[:-1])
	data_sym = list(reversed(data_sym))

	angle = data[-1][0]
	step = data[1][0] - data[0][0]

	for i in range(len(data[:-1])):
		angle += step
		energy = data_sym[i][1]
		data.append((angle, energy)) 

	return data


def WriteData(data):
	filename = 'point_energies.txt'
	f = open(filename, 'w')

	for point in data:
		f.write('%s\t%s\n' %(point[0], point[1]))

	f.close()


##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##

if len(sys.argv) < 1:
	print 'Usage:'
	# print '  python %s <step size> <scan symmetry (C2, C3, C4, etc.)>' %sys.argv[0]
	print '  python %s' %sys.argv[0]
	exit()


t0 = time.time()
filenames = glob.glob('log_*.lammps')

data = AnalyzeFiles(filenames)
WriteData(data)

print 'Found all point energies in %.4f seconds.' %(time.time()-t0)

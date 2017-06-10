import copy
import glob
import lmfit
import math
import numpy
import sys
import time


##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##

def LoadDFTEnergy(filename):
	DFT_scan = []

	f = open(filename)

	# skip the first 4 lines
	for i in range(4):
		f.readline()

	while True:
		fields = f.readline().strip().split()
		if fields:
			dihedral_angle = int(float(fields[0]))
			# convert hartree to eV 
			energy = float(fields[1])*27.2144 
			DFT_scan.append((dihedral_angle, energy))
		else:
			break
	f.close()

	# use symmery of molecule to get 360 degree scan data
	# for 180 degree scan 
	data_sym = copy.deepcopy(DFT_scan[:-1])
	data_sym = list(reversed(data_sym))

	dihedral_angle = DFT_scan[-1][0]
	step = DFT_scan[1][0] - DFT_scan[0][0]

	for i in range(len(DFT_scan[:-1])):
		dihedral_angle += step
		energy = data_sym[i][1]
		DFT_scan.append((dihedral_angle, energy))

	# find the minimum energy and convert energy to change in energy
	# also convert hartree to eV 
	E_min = min(DFT_scan, key = lambda x: x[1])[1]
	DFT_scan = map(lambda x: (x[0], x[1]-E_min) , DFT_scan)

	WriteDFTScan(DFT_scan)

	return DFT_scan


def LoadMDEnergy(filename):
	MD_scan = []

	f = open(filename)

	while True:
		fields = f.readline().strip().split()
		if fields:
			dihedral_angle = int(fields[0])
			energy = float(fields[1])
			MD_scan.append((dihedral_angle, energy)) 
		else:
			break
	f.close()

	return MD_scan


def SubtractMDfromDFT(MD_scan, DFT_scan):
	return map(lambda x, y: (x[0], y[1]-x[1]), MD_scan, DFT_scan)


def opls(phi, K1, K2, K3, K4):
	# phi is given in degrees to convert to radians
	return (1.0/2)*(K1*(1 + numpy.cos(phi*numpy.pi/180.0)) + \
					K2*(1 - numpy.cos(2*phi*numpy.pi/180.0)) + \
					K3*(1 + numpy.cos(3*phi*numpy.pi/180.0)) + \
					K4*(1 - numpy.cos(4*phi*numpy.pi/180.0)))

def FitOPLS(data):
	data = numpy.array(data)
	phi = data[:,0] # dihedral angle
	energy = data[:,1] # energy in eV 

	opls_model = lmfit.Model(opls)

	# give initial guess for the fit 
	result = opls_model.fit(energy, phi=phi, K1=0.375, K2=-0.2, K3=0.625, K4=-0.1)

	return result.fit_report()


def WriteFitReport(fit_report):
	filename = 'opls_fit.txt'
	f = open(filename, 'w')
	f.write('%s' %fit_report)
	f.close()


def WriteDFTScan(DFT_scan):
	filename = 'DFT_scan_energies.txt'
	f = open(filename, 'w')
	for scan in DFT_scan:
		f.write('%d\t%.10f\n' %(scan[0], scan[1]))
	f.close()


##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##

if len(sys.argv) < 3:
	print 'Usage:'
	# print '  python %s <step size> <scan symmetry (C2, C3, C4, etc.)>' %sys.argv[0]
	print '  python %s <MD file> <DFT file>' %sys.argv[0]
	exit()


t0 = time.time()
# filenames = glob.glob('log_*.lammps')
filename_MD = sys.argv[1]
filename_DFT = sys.argv[2]

MD_scan = LoadMDEnergy(filename_MD)
DFT_scan = LoadDFTEnergy(filename_DFT)
DFTminusMD = SubtractMDfromDFT(MD_scan, DFT_scan)

fit_report = FitOPLS(DFTminusMD)
WriteFitReport(fit_report)


print 'Fit dihedral potential in %.4f seconds.' %(time.time()-t0)
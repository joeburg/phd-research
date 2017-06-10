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

def LoadData(filename):
	data = []

	f = open(filename)

	while True:
		fields = f.readline().strip().split()
		if fields:
			dihedral_angle = float(fields[0])
			energy = float(fields[1])
			data.append((dihedral_angle, energy)) 
		else:
			break
	f.close()
	return data

#---------------------------------------------------------------------------#
def opls1_1(phi, K1):
	# phi is given in degrees to convert to radians
	return (1.0/2)*(K1*(1 + numpy.cos(phi*numpy.pi/180.0)))

def opls1_2(phi, K2):
	# phi is given in degrees to convert to radians
	return (1.0/2)*(K2*(1 - numpy.cos((2*phi)*numpy.pi/180.0)))

def opls1_3(phi, K3):
	# phi is given in degrees to convert to radians
	return (1.0/2)*(K3*(1 + numpy.cos(3*phi*numpy.pi/180.0)))

def opls1_4(phi, K4):
	# phi is given in degrees to convert to radians
	return (1.0/2)*(K4*(1 - numpy.cos(4*phi*numpy.pi/180.0)))

def opls2(phi, K1, K2):
	# phi is given in degrees to convert to radians
	return (1.0/2)*(K1*(1 + numpy.cos(phi*numpy.pi/180.0)) + \
					K2*(1 - numpy.cos(2*phi*numpy.pi/180.0)))

def opls2_4(phi, K2, K4):
	# phi is given in degrees to convert to radians
	return (1.0/2)*(K2*(1 - numpy.cos(2*phi*numpy.pi/180.0)) + \
					K4*(1 - numpy.cos(4*phi*numpy.pi/180.0)))

def opls3(phi, K1, K2, K3):
	# phi is given in degrees to convert to radians
	return (1.0/2)*(K1*(1 + numpy.cos(phi*numpy.pi/180.0)) + \
					K2*(1 - numpy.cos(2*phi*numpy.pi/180.0)) + \
					K3*(1 + numpy.cos(3*phi*numpy.pi/180.0)))

def opls4(phi, K1, K2, K3, K4):
	# phi is given in degrees to convert to radians
	return (1.0/2)*(K1*(1 + numpy.cos(phi*numpy.pi/180.0)) + \
					K2*(1 - numpy.cos(2*phi*numpy.pi/180.0)) + \
					K3*(1 + numpy.cos(3*phi*numpy.pi/180.0)) + \
					K4*(1 - numpy.cos(4*phi*numpy.pi/180.0)))

def FitOPLS(data, Ncoeffs):
	data = numpy.array(data)
	phi = data[:,0] # dihedral angle
	energy = data[:,1] # energy in eV 

	if Ncoeffs == 1:
		opls_model = lmfit.Model(opls1_1)
		# give initial guess for the fit 
		result1 = opls_model.fit(energy, phi=phi, K1=0.375)

		opls_model = lmfit.Model(opls1_2)
		result2 = opls_model.fit(energy, phi=phi, K2=-0.2)

		opls_model = lmfit.Model(opls1_3)
		result3 = opls_model.fit(energy, phi=phi, K3=0.625)

		opls_model = lmfit.Model(opls1_4)
		result4 = opls_model.fit(energy, phi=phi, K4=-0.1)

		return result1.fit_report() + result2.fit_report() \
				+ result3.fit_report() + result4.fit_report()

	elif Ncoeffs == 2:
		opls_model = lmfit.Model(opls2)
		result = opls_model.fit(energy, phi=phi, K1=0.375, K2=-0.2)		

	elif Ncoeffs == 24:
		opls_model = lmfit.Model(opls2_4)
		result = opls_model.fit(energy, phi=phi, K2=-0.2,  K4=-0.1)		

	elif Ncoeffs == 3:
		opls_model = lmfit.Model(opls3)
		result = opls_model.fit(energy, phi=phi, K1=0.375, K2=-0.2, K3=0.625)

	elif Ncoeffs == 4:
		opls_model = lmfit.Model(opls4)
		result = opls_model.fit(energy, phi=phi, K1=0.375, K2=-0.2, K3=0.625, K4=-0.1)	

	return result.fit_report()
#---------------------------------------------------------------------------#

def charmm(phi, n, d, K):
	return K*(1 + numpy.cos((n*phi - d)*numpy.pi/180.0))

def FitCHARMM(data):
	data = numpy.array(data)
	phi = data[:,0] # dihedral angle
	energy = data[:,1] # energy in eV 

	charmm_model = lmfit.Model(charmm)
	# result = charmm_model.fit(energy, phi=phi, n=2, d=-3, K=0.0097)
	result = charmm_model.fit(energy, phi=phi, n=1, d=180, K=6)

	return result.fit_report()	

#---------------------------------------------------------------------------#
def compass(phi, K1, phi1, K2, phi2, K3, phi3):
	return K1*(1 - numpy.cos((phi - phi1)*numpy.pi/180.0)) +\
			K2*(1 - numpy.cos((2*phi - phi2)*numpy.pi/180.0)) +\
			K3*(1 - numpy.cos((3*phi - phi3)*numpy.pi/180.0))

def FitCOMPASS(data):
	data = numpy.array(data)
	phi = data[:,0] # dihedral angle
	energy = data[:,1] # energy in eV 

	compass_model = lmfit.Model(compass)
	# result = compass_model.fit(energy, phi=phi, K1=0, phi1=180, \
	# 									K2=3, phi2=0, K3=0, phi3=0)
	result = compass_model.fit(energy, phi=phi, K1=0, phi1=0, \
										K2=0.007, phi2=180, K3=0.005, phi3=-90)	
	return result.fit_report()	
#---------------------------------------------------------------------------#
def harmonic(r, r0, K):
	return K*(r-r0)**2

def FitHarmonic(data):
	data = numpy.array(data)
	r = data[:,0] # bond length
	energy = data[:,1] # energy in eV 

	harmonic_model = lmfit.Model(harmonic)
	# result = harmonic_model.fit(energy, r=r, r0=1.9, K=2)	
	result = harmonic_model.fit(energy, r=r, r0=120, K=2)

	return result.fit_report()	

#---------------------------------------------------------------------------#

def WriteFitReport(fit_report, filename):
	f = open(filename, 'w')
	f.write('%s' %fit_report)
	f.close()



##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##
##---------------------------------------------------------------------------##

if len(sys.argv) < 3:
	print 'Usage:'
	print '  python %s <analyzed data file> <opls, charmm, compass, harmonic>' %sys.argv[0]
	exit()


t0 = time.time()
filename = sys.argv[1]
fit_type = sys.argv[2]

data = LoadData(filename)

if fit_type == 'opls':
	Ncoeffs = input('<coeffs to fit (1, 2, 3, 4, 24)>? ')
	fit_report = FitOPLS(data, Ncoeffs)
	filename = 'opls_fit.txt'

elif fit_type == 'charmm':
	fit_report = FitCHARMM(data)
	filename = 'charmm_fit.txt'

elif fit_type == 'compass':
	fit_report = FitCOMPASS(data)
	filename = 'compass_fit.txt'

elif fit_type == 'harmonic':
	fit_report = FitHarmonic(data)
	filename = 'harmonic_fit.txt'

WriteFitReport(fit_report, filename)


print 'Fit dihedral potential in %.4f seconds.' %(time.time()-t0)




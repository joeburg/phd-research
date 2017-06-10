''' Estimate the Volume %% Porosity based on the number of Porogen molecules'''
import math
import numpy
import sys

#------------------------------------------------------------------------------#

def VolSphere(r):
	return 4.0/3*math.pi*r**3

def ComputePorosity(N_OCSEt, N_OCSEtMe, N_135Benz, NO, N_PP40, N_PP20):
	VSi = VolSphere(2.1)
	VC = VolSphere(1.7)
	VO = VolSphere(1.52)

	NSi = 2.0*N_OCSEt + 2*N_OCSEtMe + 3*N_135Benz
	NC = 2.0*N_OCSEt + 3*N_OCSEtMe + 6*N_135Benz + 40*N_PP40 + 20*N_PP20

	vol = VSi*NSi + VC*NC + VO*NO
	volPP = VC*(20.0*N_PP20 + 40*N_PP40)

	print 'Volume %% porosity = %.1f%%' %(100.0*volPP/vol)

#------------------------------------------------------------------------------#
if len(sys.argv) < 7:
	print 'Usage:'
	print '  python %s <N OCSEt> <N OCSEtMe> <N 135Benz> <N oxygen> <N porogen40> <N porogen20>' %sys.argv[0]
	exit()

N_OCSEt = int(sys.argv[1])
N_OCSEtMe = int(sys.argv[2])
N_135Benz = int(sys.argv[3])
NO = int(sys.argv[4])
N_PP40 = int(sys.argv[5])
N_PP20 = int(sys.argv[6])

ComputePorosity(N_OCSEt, N_OCSEtMe, N_135Benz, NO, N_PP40, N_PP20)
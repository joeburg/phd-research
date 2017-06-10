""" Computes the CTE of Et-OCS given a LAMMPS log file """

import numpy
import pylab
import sys
import time

#------------------------------------------------------------------------------------#
def LoadData(logfile,beginLine,Nsteps,Nramps):
	f = open(logfile)
	for i in range(beginLine-1):
		f.readline()

	# data structures: dataT = { T : (lx, ly, lz) } 
	dataT = {}
	T = numpy.zeros(Nsteps)
	Lx = numpy.zeros(Nsteps)
	Ly = numpy.zeros(Nsteps)
	Lz = numpy.zeros(Nsteps)

	#beginLine-1,beginLine+Nsteps-1
	for i in range(Nsteps):
		fields = f.readline().strip().split()
		if fields:
			Tavg = float(fields[2])
			lx = float(fields[5])
			ly = float(fields[6])
			lz = float(fields[7])

			dataT[Tavg] = (lx,ly,lz)

			T[i] = Tavg
			Lx[i] = lx
			Ly[i] = ly
			Lz[i] = lz

	if Nramps > 1:
		Ngarbage = 28
		for i in range(Ngarbage):
			f.readline()

		dataT2 = {}
		T2 = numpy.zeros(Nsteps)
		Lx2 = numpy.zeros(Nsteps)
		Ly2 = numpy.zeros(Nsteps)
		Lz2 = numpy.zeros(Nsteps)		
		for i in range(Nsteps):
			fields = f.readline().strip().split()
			if fields:
				Tavg = float(fields[2])
				lx = float(fields[5])
				ly = float(fields[6])
				lz = float(fields[7])

				dataT2[Tavg] = (lx,ly,lz)

				T2[i] = Tavg
				Lx2[i] = lx
				Ly2[i] = ly
				Lz2[i] = lz
		
		return dataT, T, Lx, Ly, Lz, dataT2, T2, Lx2, Ly2, Lz2
	else:
		return dataT, T, Lx, Ly, Lz


def Plot(Li,T,axis,Nfigure,clr,filename):
	if axis == 'x':
		ylbl = 'simulation cell length, $L_{x}$ ($\AA$)'
	elif axis == 'y':
		ylbl = 'simulation cell length, $L_{y}$ ($\AA$)'
	elif axis == 'z':
		ylbl = 'simulation cell length, $L_{z}$ ($\AA$)'

	pylab.figure(Nfigure)
	pylab.scatter(T,Li,marker='.',color=clr,edgecolors='none')
	pylab.xlabel('temperature, T (K)', size=20)
	pylab.ylabel(ylbl, size=20)
	pylab.savefig(filename)
	pylab.close(Nfigure)

def Plot2(Li1,T1,Li2,T2,axis,Nfigure,filename):
	if axis == 'x':
		ylbl = 'simulation cell length, $L_{x}$ ($\AA$)'
	elif axis == 'y':
		ylbl = 'simulation cell length, $L_{y}$ ($\AA$)'
	elif axis == 'z':
		ylbl = 'simulation cell length, $L_{z}$ ($\AA$)'

	pylab.figure(Nfigure)
	pylab.scatter(T1,Li1,marker='.',color='b',edgecolors='none')
	pylab.scatter(T2,Li2,marker='.',color='r',edgecolors='none')
	pylab.xlabel('temperature, T (K)', size=20)
	pylab.ylabel(ylbl, size=20)
	pylab.savefig(filename)
	pylab.close(Nfigure)


def SplitData(T,Lx,Ly,Lz,cutoff):
	tol = 1
	for i in range(len(T)):
		if abs(T[i] - cutoff) < tol:
			tempidx = i
			break

	# separate the two slopes
	T1 = T[:tempidx]
	T2 = T[tempidx:]

	Lx1 = Lx[:tempidx]
	Lx2 = Lx[tempidx:]

	Ly1 = Ly[:tempidx]
	Ly2 = Ly[tempidx:]

	Lz1 = Lz[:tempidx]
	Lz2 = Lz[tempidx:]

	return T1, T2, Lx1, Lx2, Ly1, Ly2, Lz1, Lz2

def PolyInterpolation(T,Li):
	degree = 1
	while True:
		interp = numpy.polyfit(T,Li,degree,full=True)

		Li_T = interp[0]

		if len(interp[1]) > 0:
			residual = interp[1][0]
		else: 
			raise RuntimeError, 'Connot compute polynomial interpolation.'

		if residual < 1e1:
			degree += 1 
		else:
			return Li_T

def PolyDerivative(Li_T):
	Li = numpy.poly1d(Li_T)
	Li_prime = numpy.polyder(Li)
	return Li_prime

def ComputeCTE(Lx_prime,Ly_prime,Lz_prime,dataT,temp):
	for t in dataT:
		if abs(temp - t) <= 1:
			temp = t
			break

	lx, ly, lz = dataT[temp]
	lx_prime_T = Lx_prime(temp)
	ly_prime_T = Ly_prime(temp)
	lz_prime_T = Lz_prime(temp)

	alpha = (1/3.0)*( (1/lx)*lx_prime_T + (1/ly)*ly_prime_T + (1/lz)*lz_prime_T )
	return alpha

def WritetoFile2(alpha1,alpha2,outfile):
	f = open(outfile,'w')
	f.write('CTE (low temp) = %.10f  K^-1 \n' %alpha1)
	f.write('CTE (high temp) = %.10f  K^-1' %alpha2)
	f.close()

def WritetoFile1(alpha,outfile):
	f = open(outfile,'w')
	f.write('CTE = %.10f  K^-1 \n' %alpha)
	f.close()

#------------------------------------------------------------------------------------#
# analyze the command line arguments and setup corresponding parameters
if len(sys.argv) < 6:
	print 'Usage:'
	print '  python %s <log file> <begin line> <N steps> <temp cutoff> <heating / cooling> [N ramps (default=1)]'%sys.argv[0]
	exit()

t0 = time.time()

logfile = sys.argv[1]
beginLine = int(sys.argv[2])
Nsteps = int(sys.argv[3])
cutoff = float(sys.argv[4])
direction = sys.argv[5]

if direction == 'heating':
	color = 'r'
elif direction == 'cooling':
	color = 'b'

if len(sys.argv) == 7:
	Nramps = int(sys.argv[6])
else:
	Nramps = 1

if Nramps > 1:
	dataT, T, Lx, Ly, Lz, dataT2, T_2, Lx_2, Ly_2, Lz_2 = LoadData(logfile,beginLine,Nsteps,Nramps)
else:
	dataT, T, Lx, Ly, Lz = LoadData(logfile,beginLine,Nsteps,Nramps)

# use for one CTE
Plot(Lx,T,'x',1,color,'T_Lx_%s.pdf' %direction)
# Plot2(Lx,T,Lx_2,T_2,'x',1,'T_Lx_cool_heat.pdf')

Lx_T = PolyInterpolation(T,Lx)
Ly_T = PolyInterpolation(T,Ly)
Lz_T = PolyInterpolation(T,Lz)
print Lx_T

Lx_prime = PolyDerivative(Lx_T)
Ly_prime = PolyDerivative(Ly_T)
Lz_prime = PolyDerivative(Lz_T)

alpha = ComputeCTE(Lx_prime,Ly_prime,Lz_prime,dataT,250)
WritetoFile1(alpha,'CTE_%s.txt' %direction)
print 'CTE = %.10f' %alpha

# Lx_T_2 = PolyInterpolation(T_2,Lx_2)
# Ly_T_2 = PolyInterpolation(T_2,Ly_2)
# Lz_T_2 = PolyInterpolation(T_2,Lz_2)
# print Lx_T_2

# Lx_prime_2 = PolyDerivative(Lx_T_2)
# Ly_prime_2 = PolyDerivative(Ly_T_2)
# Lz_prime_2 = PolyDerivative(Lz_T_2)

# alpha_2 = ComputeCTE(Lx_prime_2,Ly_prime_2,Lz_prime_2,dataT,250)
# print 'CTE = %.10f' %alpha_2

# WritetoFile2(alpha,alpha_2,'CTE_results.txt')




# for data with different CTE use the following code
# T1, T2, Lx1, Lx2, Ly1, Ly2, Lz1, Lz2 = SplitData(T,Lx,Ly,Lz,cutoff)

# # determine the CTE for low temp slope at 300 K
# Lx_T1 = PolyInterpolation(T1,Lx1)
# Ly_T1 = PolyInterpolation(T1,Ly1)
# Lz_T1 = PolyInterpolation(T1,Lz1)

# print Lx_T1

# Lx_prime1 = PolyDerivative(Lx_T1)
# Ly_prime1 = PolyDerivative(Ly_T1)
# Lz_prime1 = PolyDerivative(Lz_T1)

# alpha1 = ComputeCTE(Lx_prime1,Ly_prime1,Lz_prime1,dataT,300)
# print 'CTE = %.10f' %alpha1

# # determine the CTE for high temp slope at 600 K
# Lx_T2 = PolyInterpolation(T2,Lx2)
# Ly_T2 = PolyInterpolation(T2,Ly2)
# Lz_T2 = PolyInterpolation(T2,Lz2)

# print Lx_T2

# Lx_prime2 = PolyDerivative(Lx_T2)
# Ly_prime2 = PolyDerivative(Ly_T2)
# Lz_prime2 = PolyDerivative(Lz_T2)

# alpha2 = ComputeCTE(Lx_prime2,Ly_prime2,Lz_prime2,dataT,600)
# print 'CTE = %.10f' %alpha2

# WritetoFile2(alpha1,alpha2,'CTE_results.txt')



print 'Time elapsed %.4f seconds' %(time.time() - t0)


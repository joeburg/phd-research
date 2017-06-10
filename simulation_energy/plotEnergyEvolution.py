""" Computes the CTE of Et-OCS given a LAMMPS log file """

import numpy
import pylab
import sys
import time

#------------------------------------------------------------------------------------#
def LoadData(logfile,beginLine,Nsteps):
	timestep = 0.001

	f = open(logfile)
	for i in range(beginLine-1):
		f.readline()

	T = numpy.zeros(Nsteps)
	Epair = numpy.zeros(Nsteps)
	Etot = numpy.zeros(Nsteps)
	Time = numpy.zeros(Nsteps)

	#beginLine-1,beginLine+Nsteps-1
	for i in range(Nsteps):
		fields = f.readline().strip().split()
		if fields:
			step = float(fields[0])
			Tavg = float(fields[2])
			lx = float(fields[5])
			ly = float(fields[6])
			lz = float(fields[7])
			V = float(fields[8])
			Ep = float(fields[9])
			Et = float(fields[10])

			Time[i] = timestep*step
			T[i] = Tavg
			Epair[i] = Ep
			Etot[i] = Et

	return Time, T, Epair, Etot


def Plot(xdata,ydata,axis,Nfigure,clr,filename):
	ylbl = 'energy, E (eV)'

	if axis == 'time':
		xlbl = 'time, t (ps)'
	elif axis == 'temp':
		xlbl = 'temperature, T (K)'

	if clr == 'r':
		filename = filename+'_heating.pdf'
	elif clr == 'b':
		filename = filename+'_cooling.pdf'

	pylab.figure(Nfigure)
	pylab.scatter(xdata,ydata,marker='.',color=clr,edgecolors='none')
	pylab.xlabel(xlbl, size=20)
	pylab.ylabel(ylbl, size=20)
	pylab.savefig(filename,bbox_inches='tight')
	pylab.close(Nfigure)

def Plot2(xdata,ydata1,ydata2,axis,Nfigure,clr,filename):
	ylbl = 'energy, E (eV)'

	if axis == 'time':
		xlbl = 'time, t (ps)'
	elif axis == 'temp':
		xlbl = 'temperature, T (K)'

	if clr == 'r':
		filename = filename+'_heating.pdf'
	elif clr == 'b':
		filename = filename+'_cooling.pdf'

	pylab.figure(Nfigure)
	pylab.scatter(xdata,ydata1,marker='o',color=clr,edgecolors='none')
	pylab.scatter(xdata,ydata2,marker='s',color=clr,edgecolors='none')
	pylab.xlabel(xlbl, size=20)
	pylab.ylabel(ylbl, size=20)
	pylab.savefig(filename,bbox_inches='tight')
	pylab.close(Nfigure)

def Plot4(xdata1,xdata2,ydata1,ydata2,ydata3,ydata4,axis,Nfigure,clr,filename):
	ylbl = 'energy, E (eV)'

	if axis == 'time':
		xlbl = 'time, t (ps)'
	elif axis == 'temp':
		xlbl = 'temperature, T (K)'

	if clr == 'r':
		filename = filename+'_heating.pdf'
	elif clr == 'b':
		filename = filename+'_cooling.pdf'

	pylab.figure(Nfigure)
	pylab.scatter(xdata1,ydata1,marker='o',color='b',edgecolors='none')
	pylab.scatter(xdata1,ydata2,marker='s',color='b',edgecolors='none')
	pylab.scatter(xdata2,ydata3,marker='o',color='r',edgecolors='none')
	pylab.scatter(xdata2,ydata4,marker='s',color='r',edgecolors='none')
	pylab.xlabel(xlbl, size=20)
	pylab.ylabel(ylbl, size=20)
	pylab.savefig(filename,bbox_inches='tight')
	pylab.close(Nfigure)


#------------------------------------------------------------------------------------#
# analyze the command line arguments and setup corresponding parameters
if len(sys.argv) < 5:
	print 'Usage:'
	print '  python %s <log file> <begin line> <N steps> <heating / cooling> [log file 2]'%sys.argv[0]
	exit()

t0 = time.time()

logfile = sys.argv[1]
beginLine = int(sys.argv[2])
Nsteps = int(sys.argv[3])
direction = sys.argv[4]

if direction == 'heating':
	color = 'r'
elif direction == 'cooling':
	color = 'b'

if len(sys.argv) == 6:
	logfile2 = sys.argv[5]

	Time, T, Epair, Etot = LoadData(logfile,beginLine,Nsteps)
	Time2, T2, Epair2, Etot2 = LoadData(logfile2,beginLine,Nsteps)

	Plot4(T,T2,Epair,Etot,Epair2,Etot2,'temp',1,color,'Energy_temp')
	Plot4(Time,Time2,Epair,Etot,Epair2,Etot2,'time',1,color,'Energy_time')

else:
	# load data and then plot
	Time, T, Epair, Etot = LoadData(logfile,beginLine,Nsteps)

	Plot2(T,Epair,Etot,'temp',1,color,'Energy_temp')
	Plot2(Time,Epair,Etot,'time',1,color,'Energy_time')


print 'Time elapsed %.4f seconds' %(time.time() - t0)


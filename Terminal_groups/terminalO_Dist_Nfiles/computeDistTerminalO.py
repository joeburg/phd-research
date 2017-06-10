"""This class determines the network connectivity of the material given a LAMMPS .xyz file"""

import copy
import matplotlib.pyplot as plt
import numpy
import pylab
import scipy.interpolate

from scipy.stats import gaussian_kde

class TerminalODist:

	def __init__(self,inputfile,color):
		self.Natoms = 0
		self.NSi = 0
		self.NC = 0
		self.NO = 0

		self.Lx = 0
		self.Ly = 0
		self.Lz = 0

		self.O_stats = numpy.zeros(6) # freeO, non-bridging O, bridging O 
		self.Si_stats = numpy.zeros(6) # T0, T1, T2, T3
		self.Si_statsT = numpy.zeros(6) # only considers bridging O 

		self.p = 0 # Si-X-Si network connectivity
		self.q = 0 # condensation degree
		self.mSi = 0 # mean Si network connectivity
		self.rho = 0 # density 

		self.data = []
		self.CM = [] # Si-O bonding matrix
		self.CM_T = [] # only considers bridging O

		# analysis of network 
		self.LoadData(inputfile)
		print 'Data loaded for %s.' %inputfile

		print 'Computing Si-O bonds...'
		self.ComputeCM()

		print 'Finding terminal O...'
		terminalO = self.GetTerminalO()

		print 'Computing distance between terminal O...'
		terminalOdist = self.ComputeDist(terminalO)

		# print 'Plotting terminal O...'
		# self.PlotDensity(terminalOdist,color)

		print 'Writing data to file...'
		self.WriteData(terminalOdist,inputfile)

	def LoadData(self,inputfile):
		""" reads LAMMPS file, sorts data, and computes cell dimensions """ 
		f = open(inputfile)
		self.Natoms = int(f.readline())

		f.readline()

		while True:
			fields = f.readline().strip().split()
			if fields:

				atomtype = int(fields[0])
				xcoord = float(fields[1])
				ycoord = float(fields[2])
				zcoord = float(fields[3])

				# determine type of atom
				if atomtype == 1:
					self.NSi += 1
				elif atomtype == 2:
					self.NO += 1
				elif atomtype == 3:
					self.NC += 1
				else: 
					raise RuntimeError, "Incorrect atom type."

				# populate the data array
				self.data.append([atomtype,xcoord,ycoord,zcoord])

			else:
				break
		f.close()

		self.data = numpy.array(self.data)
		b = numpy.zeros((self.Natoms,4))

		# sort the atoms by Si, O, C
		Siidx = 0
		Oidx = self.NSi
		Cidx = self.NSi + self.NO
		for i in range(self.Natoms):
			if self.data[i][0] == 1:
				b[Siidx] = self.data[i]
				Siidx += 1
			elif self.data[i][0] == 2:
				b[Oidx] = self.data[i]
				Oidx += 1
			elif self.data[i][0] == 3:
				b[Cidx] = self.data[i]
				Cidx += 1

		self.data = b

		# find the dimensions of the simulation cell
		self.Lx = max(self.data[:,1] - min(self.data[:,1]))
		self.Ly = max(self.data[:,2] - min(self.data[:,2]))
		self.Lz = max(self.data[:,3] - min(self.data[:,3]))

	def Distance(self,dx,dy,dz):
		return (dx*dx + dy*dy + dz*dz)**0.5

	def ComputeCM(self):
		""" Computes all Si-O bonds and populates a CM matrix and an ID matrix""" 

		self.CM = numpy.zeros((self.NSi,self.NO))
	
		for i in range(self.NSi):
			for j in range(self.NSi,self.NSi+self.NO-1):

				dx = abs(self.data[j][1] - self.data[i][1])
				dy = abs(self.data[j][2] - self.data[i][2])
				dz = abs(self.data[j][3] - self.data[i][3])

				if self.Distance(dx,dy,dz) < 2.3:
					self.CM[i][j-self.NSi] = 1
				
				# consider PBCs
				if dx > self.Lx - 2.3:
					if self.Distance(dx + self.Lx, dy, dz) < 2.3:
						self.CM[i][j-self.NSi] = 1

					if self.Distance(dx - self.Lx, dy, dz) < 2.3:
						self.CM[i][j-self.NSi] = 1

				if dy > self.Ly - 2.3:
					if self.Distance(dx ,dy + self.Ly, dz) < 2.3:
						self.CM[i][j-self.NSi] = 1

					if self.Distance(dx, dy - self.Ly, dz) < 2.3:
						self.CM[i][j-self.NSi] = 1

				if dz > self.Lz - 2.3:
					if self.Distance(dx, dy, dz + self.Lz) < 2.3:
						self.CM[i][j-self.NSi] = 1

					if self.Distance(dx, dy, dz - self.Lz) < 2.3:
						self.CM[i][j-self.NSi] = 1

		self.CM_T = copy.deepcopy(self.CM)

	def GetTerminalO(self):
		""" Compute free and terminal O and write out""" 
		terminalO = []

		for i in range(self.NO):
			NO = numpy.sum(self.CM[:,i])
			if  NO == 1: # terminal O
				# index corresponding to data[] array
				idx = i + self.NSi
				terminalO.append(idx)

		return terminalO

	def ComputeDist(self, terminalO):
		terminalOdist = []
		for i in range(len(terminalO)):
			i_idx = terminalO[i]

			for j in range(i+1,len(terminalO)):				
				j_idx = terminalO[j]
 
				dx = abs(self.data[j_idx][1] - self.data[i_idx][1])
				dy = abs(self.data[j_idx][2] - self.data[i_idx][2])
				dz = abs(self.data[j_idx][3] - self.data[i_idx][3])

				d = self.Distance(dx,dy,dz)
				terminalOdist.append(d)

		terminalOdist = numpy.array(terminalOdist)
		return terminalOdist

	def GetDensity(self, data):
	    density = gaussian_kde(data)
	    density.covariance_factor = lambda : 0.075
	    density._compute_covariance()
	    return density     

	def PlotDensity(self,data,color):
		# plot the density function 
		density = self.GetDensity(data)
		xs = numpy.linspace(0,200,200)
		plt.plot(xs,density(xs),color)

	def WriteData(self,data,inputfile):
		# generate the bins for collecting the data frequency
		bins = numpy.linspace(0,200,201)

		# generate the histogram for the data set
		# note that len(hist) = len(data)-1
		histdata = numpy.histogram(data,bins=bins,normed=True)
		hist = histdata[0]

		# write out the density values
		fname = 'terminalODist_{}.csv'.format(inputfile[:-4])

		f = open(fname,'w')
		for i in range(len(hist)):
			f.write('%.2f,%.8f\n'%(bins[i],hist[i]))
		f.close()




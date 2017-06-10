"""This class determines the network connectivity of the material given a LAMMPS .xyz file"""

import copy
import matplotlib.pyplot as plt
import numpy
import pylab
import scipy.interpolate

from scipy.stats import gaussian_kde

class AtomDistance:

	def __init__(self,filenames,colors):
		self.Natoms = 0
		self.NSi = 0
		self.NC = 0
		self.NO = 0

		self.Lx = 0
		self.Ly = 0
		self.Lz = 0

		self.data = []

		self.ProcessFiles(filenames,colors)

	def ProcessFiles(self,filenames,colors):
		# set font type and axes thickness
		arialfont = {'fontname':'Arial'}
		plt.rcParams['axes.linewidth'] = 2
		plt.rcParams.update({'fontname':'Arial','font.size': 16})

		# process Si-Si distances 
		plt.figure(1)
		for i in range(len(filenames)):
			filename = filenames[i]
			color = colors[i]

			# analysis of network 
			self.LoadData(filename)
			print 'Data loaded for %s.' %filename

			print 'Computing Si-Si distance...'
			SiSi_dist = self.ComputeDist('Si','Si')

			print 'Plotting results...'
			self.PlotDensity(SiSi_dist,color)
			
			print 'Writing data to file...'
			fname = 'SiSi_dist_{}.csv'.format(filename[:-4])
			self.WriteData(SiSi_dist,fname)

			# clear the data array and simulation attributes
			self.Natoms = 0
			self.NSi = 0
			self.NC = 0
			self.NO = 0

			self.Lx = 0
			self.Ly = 0
			self.Lz = 0

			self.data = []

		plt.xlabel('Si-Si Distance, $\delta$ ($\AA$)',size=20,**arialfont)
		plt.ylabel('Frequency',size=20,**arialfont)
		plt.xlim((0,120))
		plt.axes().set_aspect(1./plt.axes().get_data_ratio())
		plt.savefig('SiSi_dist.pdf',bbox_inches='tight')
		plt.close(1)	


		# process C-C distances 
		plt.figure(2)
		for i in range(len(filenames)):
			filename = filenames[i]
			color = colors[i]

			# analysis of network 
			self.LoadData(filename)
			print 'Data loaded for %s.' %filename

			# analysis of network 
			print 'Computing C-C distance...'
			CC_dist = self.ComputeDist('C','C')

			print 'Plotting results...'
			self.PlotDensity(CC_dist,color)

			print 'Writing data to file...'
			fname = 'CC_dist_{}.csv'.format(filename[:-4])
			self.WriteData(CC_dist,fname)

			# clear the data array and simulation attributes
			self.Natoms = 0
			self.NSi = 0
			self.NC = 0
			self.NO = 0

			self.Lx = 0
			self.Ly = 0
			self.Lz = 0

			self.data = []

		plt.xlabel('C-C Distance, $\delta$ ($\AA$)',size=20,**arialfont)
		plt.ylabel('Frequency',size=20,**arialfont)
		plt.xlim((0,120))
		plt.axes().set_aspect(1./plt.axes().get_data_ratio())
		plt.savefig('CC_dist.pdf',bbox_inches='tight')
		plt.close(2)	


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

	def ComputeDist(self, type1, type2):
		# data ordered by Si, O, and C type atoms
		if type1 == 'Si' and type2 == 'Si':
			Ndist = (self.NSi*(self.NSi-1))/2
			N = self.NSi
			idx_shift = 0

		elif type1 == 'C' and type2 == 'C':
			Ndist = (self.NC*(self.NC-1))/2
			N = self.NC
			idx_shift = self.NSi + self.NO

		atomDist = numpy.zeros(Ndist)

		n = 0
		for i in range(N):
			i_idx = i + idx_shift

			for j in range(i+1,N):				
				j_idx = j + idx_shift
 	
				dx = abs(self.data[j_idx][1] - self.data[i_idx][1])
				dy = abs(self.data[j_idx][2] - self.data[i_idx][2])
				dz = abs(self.data[j_idx][3] - self.data[i_idx][3])

				d = self.Distance(dx,dy,dz)
				
				atomDist[n] = d
				n += 1

		return atomDist

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

	def WriteData(self,data,fname):
		# generate the bins for collecting the data frequency
		bins = numpy.linspace(0,200,201)

		# generate the histogram for the data set
		# note that len(hist) = len(data)-1
		histdata = numpy.histogram(data,bins=bins,normed=False)
		hist = histdata[0]

		# write out the histogram values
		f = open(fname,'w')
		for i in range(len(hist)):
			f.write('%.2f,%.2f\n'%(bins[i],hist[i]))
		f.close()




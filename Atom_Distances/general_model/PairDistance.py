import math
import numpy
import sys
import time


#-----------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------#
def Distance(dx, dy, dz):
	return (dx*dx + dy*dy + dz*dz)**0.5

#-----------------------------------------------------------------------------------------#

class PairDistance:

	def __init__(self, filename, atom_type1, atom_type2, L):
		if L:
			self.Lx = L
			self.Ly = L
			self.Lz = L
		else:
			self.Lx = 0
			self.Ly = 0
			self.Lz = 0

		data = self.LoadData(filename)

		data_type1, data_type2, same_types = self.SortData(data, atom_type1, atom_type2)

		distances = self.ComputeDistances(data_type1, data_type2, same_types)

		self.WriteData(distances, filename, atom_type1, atom_type2)

	#------------------------------------------------------------------------------------#
	def LoadData(self, filename):
		f = open(filename)

		# the first 2 lines are not relevant 
		f.readline()
		f.readline()

		# store only the relevant data 
		data = []

		while True:
			fields = f.readline().strip().split()
			if fields:
				atom_type = int(fields[0])
				x = float(fields[1])
				y = float(fields[2])
				z = float(fields[3])

				data.append((atom_type, x, y, z))

			else:
				break
		f.close()

		# get the cell dimensions
		data = numpy.array(data)

		# if the cell dimensions are not yet set, then compute them
		if not self.Lx: 
			x_coords = data[:,1]
			y_coords = data[:,2]
			z_coords = data[:,3]

			self.Lx = self.AverageMax(x_coords, 50) - self.AverageMin(x_coords, 50)
			self.Ly = self.AverageMax(y_coords, 50) - self.AverageMin(y_coords, 50)
			self.Lz = self.AverageMax(z_coords, 50) - self.AverageMin(z_coords, 50)

		return data


	def AverageMax(self, data, N):
		data = numpy.sort(data)
		return numpy.mean(data[-N:])

	def AverageMin(self, data, N):
		data = numpy.sort(data)
		return numpy.mean(data[:N])

	#------------------------------------------------------------------------------------#
	def SortData(self, data, atom_type1, atom_type2):
		data_type1 = []
		data_type2 = []

		if atom_type1 == atom_type2:
			for atom in data:
				atom_type = atom[0]
				if atom_type == atom_type1:
					data_type1.append(atom)

			data_type2 = data_type1
			same_types = True	

		else:
			for atom in data:
				atom_type = atom[0]
				if atom_type == atom_type1:
					data_type1.append(atom)
				elif atom_type == atom_type2:
					data_type2.append(atom)

			same_types = False

		return data_type1, data_type2, same_types

	#------------------------------------------------------------------------------------#
	def addDistance(self, dx, dy, dz, distances, n):
		d = Distance(dx, dy, dz)
		distances[n] = d
		n += 1	
		return n

	def ComputeDistances(self, data_type1, data_type2, same_types):
		if same_types:
			# N1 = N2 = N; N*(N-1)
			Ninteractions = len(data_type1)*(len(data_type1) - 1)
		else:
			# N1 * N2
			Ninteractions = len(data_type1)*len(data_type2)
		
		distances = numpy.zeros(Ninteractions)

		n = 0
		for i in range(len(data_type1)):
			atom_type1, x1, y1, z1 = data_type1[i]

			for j in range(len(data_type2)):
				atom_type2, x2, y2, z2 = data_type2[j]

				dx = abs(x1 - x2)
				dy = abs(y1 - y2)
				dz = abs(z1 - z2)

				if same_types:
					if not i == j:
						n = self.addDistance(dx, dy, dz, distances, n)
				else:
					n = self.addDistance(dx, dy, dz, distances, n)

		return distances

	#------------------------------------------------------------------------------------#
	def WriteData(self, data, filename, atom_type1, atom_type2):
		# write the data to a file
		filename = '%s_%s-%s_dist.csv' %(filename[:-4], atom_type1, atom_type2)

		f = open(filename, 'w')
		for d in data:
			f.write('%.6f\n' %d)
		f.close()


		# write the data as a histrogram 
		filename = '%s_%s-%s_hist.csv' %(filename[:-4], atom_type1, atom_type2)

		# generate the bins for collecting the data frequency
		bins = numpy.linspace(0, 1000, 1001)

		# generate the histogram for the data set
		# note that len(hist) = len(data)-1
		histdata = numpy.histogram(data, bins=bins, normed=False)
		hist = histdata[0]

		# write out the histogram values
		f = open(filename,'w')
		for i in range(len(hist)):
			# only write the non-zero terms
			if hist[i]:
				f.write('%.2f,%.2f\n'%(bins[i],hist[i]))
		f.close()

#-----------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------#
# main program 
if len(sys.argv) < 4:
	print 'Usage:'
	print '  python %s <filename> <atom type 1> <atom type 2> [cell length = False]' %sys.argv[0]
	exit()


t0 = time.time()

filename = sys.argv[1]
atom_type1 = int(sys.argv[2])
atom_type2 = int(sys.argv[3])

L = False
if len(sys.argv) == 5:
	L = float(sys.argv[4])

print '\nWorking with %s...' %filename

PairDistance(filename, atom_type1, atom_type2, L)

print '\nComputed distances in %.4f seconds.' %(time.time()-t0)

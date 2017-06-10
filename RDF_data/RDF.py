''' This program computes the RDF data of a LAMMPS simulation '''

import numpy

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
class rdf:

	def __init__(self,inputfile):
		self.inputfile = inputfile

		# initialize the data array for the .xyz file; sorted by atoms type: Si (1), O (2), then C (3)
		self.data = [] 

		# initialize the atom counts and simulation cell dimensions
		self.Natoms = 0
		self.NSi = 0
		self.NC = 0
		self.NO = 0

		self.Lx = 0 
		self.Ly = 0
		self.Lz = 0

		# load the .xyz file 
		self.LoadData(inputfile)
		print '\nData loaded...'

		# compute the RDF and write out the results
		print 'Computing the RDF...'
		rdf_data = self.ComputeRDF(inputfile)

	#------------------------------------------------------------------------------#

	def LoadData(self, inputfile):
		""" reads LAMMPS file, sorts data, and computes cell dimensions """ 
		f = open(inputfile)
		self.Natoms = int(f.readline())

		f.readline()

		self.data = numpy.zeros((self.Natoms,4))

		for i in range(self.Natoms):
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
				self.data[i] = [atomtype,xcoord,ycoord,zcoord]
		f.close()

		# find the dimensions of the simulation cell
		self.Lx = max(self.data[:,1] - min(self.data[:,1]))
		self.Ly = max(self.data[:,2] - min(self.data[:,2]))
		self.Lz = max(self.data[:,3] - min(self.data[:,3]))

		# sort the rows by atom type: Si (1), O (2), then C (3)
		self.data = numpy.array(sorted(self.data, key=lambda a_entry: a_entry[0]))


	#------------------------------------------------------------------------------#	
	
	def Distance(self,dx,dy,dz):
		return (dx*dx + dy*dy + dz*dz)**0.5


	def getBond(self,atom1,atom2,cutoff):
		''' computes if 2 atoms are bonded; uses periodic boundary conditions '''
		dx = abs(atom1[1] - atom2[1])
		dy = abs(atom1[2] - atom2[2])
		dz = abs(atom1[3] - atom2[3])

		d = self.Distance(dx,dy,dz)
		if d < cutoff:
			return d
		
		# consider PBCs in x and z directions; not in y-dir due to source and sink
		if dx > self.Lx - cutoff:
			d = self.Distance(dx + self.Lx, dy, dz)
			if d < cutoff:
				return d

			d = self.Distance(dx - self.Lx, dy, dz)
			if d < cutoff:
				return d

		if dy > self.Ly - cutoff:
			d = self.Distance(dx, dy + self.Ly, dz) 
			if d < cutoff:
				return d

			d = self.Distance(dx, dy - self.Ly, dz)
			if d < cutoff:
				return d

		if dz > self.Lz - cutoff:
			d = self.Distance(dx, dy, dz + self.Lz) 
			if d < cutoff:
				return d

			d = self.Distance(dx, dy, dz - self.Lz)
			if d < cutoff:
				return d

		# else return 0 to indicate that there's no bond
		return 0


	def ComputeRDF(self, inputfile):
		SiO_bonds = []
		CC_bonds = []
		SiC_bonds = []

		OO_dist = []
		SiSi_dist = []

		# compute the Si-O bonds
		for i in range(self.NSi):
			atom1 = self.data[i]
			Nbonds = 0
			for j in range(self.NSi,self.NSi+self.NO):
				atom2 = self.data[j]

				d = self.getBond(atom1,atom2,1.74)
				if d:
					SiO_bonds.append(d)
					Nbonds += 1

				# can have a maximum of 3 Si-O bonds
				if Nbonds == 3:
					break

		# compute the Si-C bonds
		for i in range(self.NSi):
			atom1 = self.data[i]
			for j in range(self.NSi+self.NO,self.Natoms):
				atom2 = self.data[j]

				d = self.getBond(atom1,atom2,2.01)
				if d:
					SiC_bonds.append(d)
					break # only 1 possible Si-C bond

		# compute the C-C bonds
		for i in range(self.NSi+self.NO,self.Natoms):
			atom1 = self.data[i]
			for j in range(i+1,self.Natoms):
				atom2 = self.data[j]

				d = self.getBond(atom1,atom2,1.64)
				if d:
					CC_bonds.append(d)
					break # only 1 possible C-C bond

		# compute the O-O distance 
		for i in range(self.NSi,self.NSi+self.NO):
			atom1 = self.data[i]
			for j in range(i+1,self.NSi+self.NO):
				atom2 = self.data[j]

				d = self.getBond(atom1,atom2,3.0)
				if d:
					OO_dist.append(d)

		# compute the Si-Si distance 
		for i in range(self.NSi):
			atom1 = self.data[i]
			for j in range(i+1,self.NSi):
				atom2 = self.data[j]

				d = self.getBond(atom1,atom2,3.5)
				if d:
					SiSi_dist.append(d)


		# print "Number of bonds = %d" %(len(SiO_bonds)+len(SiC_bonds)+len(CC_bonds))
		
		# compute the average bond lengths 
		SiO_bonds = numpy.array(SiO_bonds)
		SiC_bonds = numpy.array(SiC_bonds)
		CC_bonds = numpy.array(CC_bonds)

		OO_dist = numpy.array(OO_dist)
		SiSi_dist = numpy.array(SiSi_dist)
		
		SiO_avg = numpy.average(SiO_bonds)
		SiO_std = numpy.std(SiO_bonds)

		SiC_avg = numpy.average(SiC_bonds)
		SiC_std = numpy.std(SiC_bonds)

		CC_avg = numpy.average(CC_bonds)
		CC_std = numpy.std(CC_bonds)

		OO_avg = numpy.average(OO_dist)
		OO_std = numpy.std(OO_dist)

		SiSi_avg = numpy.average(SiSi_dist)
		SiSi_std = numpy.std(SiSi_dist)

		# compute the mode of the Si-O, O-O and Si-Si peaks
		# choose bin width to be 0.01 A
		SiO_mode = self.ComputeMode(SiO_bonds,174)
		OO_mode = self.ComputeMode(OO_dist,300)
		SiSi_mode = self.ComputeMode(SiSi_dist,350)

		# Write out the results
		self.WriteResults(SiO_avg, SiO_std, SiC_avg, SiC_std, CC_avg, CC_std, 
							OO_avg, OO_std, SiSi_avg, SiSi_std,
							SiO_mode, OO_mode, SiSi_mode, inputfile)

		# Write out the RDF data
		self.WriteRDFData(SiO_bonds, inputfile, 'SiO')
		self.WriteRDFData(SiC_bonds, inputfile, 'SiC')
		self.WriteRDFData(CC_bonds, inputfile, 'CC')
		self.WriteRDFData(OO_dist, inputfile, 'OO')
		self.WriteRDFData(SiSi_dist, inputfile, 'SiSi')

	#------------------------------------------------------------------------------#

	def ComputeMode(self, data, Nbins):
		histdata = numpy.histogram(data,bins=Nbins)

		# find the index of the most probable value
		idx = numpy.argmax(histdata[0])

		# take the distance to be the average of the bin limits 
		most_prob_dist = 0.5*(histdata[1][idx] + histdata[1][idx+1])

		return most_prob_dist

	def WriteResults(self, SiO_avg, SiO_std, SiC_avg, SiC_std, CC_avg, CC_std, 
						OO_avg, OO_std, SiSi_avg, SiSi_std,
						SiO_mode, OO_mode, SiSi_mode, inputfile):
		outputfile = '{}_bond_lengths.csv'.format(inputfile[:-4])
		f = open(outputfile,'w')
		f.write('%s,%.8f,%.8f\n' %('SiO',SiO_avg,SiO_std))
		f.write('%s,%.8f,%.8f\n' %('SiC',SiC_avg,SiC_std))
		f.write('%s,%.8f,%.8f\n' %('CC',CC_avg,CC_std))
		f.write('%s,%.8f,%.8f\n' %('OO',OO_avg,OO_std))
		f.write('%s,%.8f,%.8f\n' %('SiSi',SiSi_avg,SiSi_std))
		f.write('%s,%.8f\n' %('Si-O mode',SiO_mode))
		f.write('%s,%.8f\n' %('O-O mode',OO_mode))
		f.write('%s,%.8f\n' %('Si-Si mode',SiSi_mode))
		f.close()


	def WriteRDFData(self, data, inputfile, bond_type):
		outputfile = '{}_{}_rdf.csv'.format(inputfile[:-4], bond_type)
		f = open(outputfile,'w')

		bins = numpy.linspace(0,4,600)
		histdata = numpy.histogram(data,bins=bins,normed=False)
		density = histdata[0].astype(float) / numpy.sum(histdata[0])

		for i in range(len(density)):
			if density[i] > 0.0000000:
				f.write('%.8f,%.10f\n' % (bins[i],density[i]))

		f.close()





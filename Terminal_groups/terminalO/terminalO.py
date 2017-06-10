"""This class find the terminal OH groups in a material"""

import copy
import numpy


class TerminalO:

	def __init__(self,inputfile):
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
		print 'Data loaded.'

		print 'Computing Si-O bonds...'
		self.ComputeCM()

		print 'Finding terminal O...'
		self.GetTerminalO()
		# self.ComputeDensity() # must compute density after connecitivty
		# self.WriteSolution(inputfile)

	def LoadData(self,inputfile):
		""" reads LAMMPS file, sorts data, and computes cell dimensions """ 
		f = open(inputfile)
		self.Natoms = int(f.readline())

		f.readline()

		ID = 1
		Oid = 2 # O atom ID
		while True:
			fields = f.readline().strip().split()
			if fields:

				atomtype = int(fields[0])
				xcoord = float(fields[1])
				ycoord = float(fields[2])
				zcoord = float(fields[3])

				if atomtype == Oid:
					# jump to O atom IDs 
					ID = 9001
					Oid = 0

				atomID = ID 
				ID += 1

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
				self.data.append([atomtype,xcoord,ycoord,zcoord,atomID])
				# print [atomtype,xcoord,ycoord,zcoord,atomID]

			else:
				break
		f.close()

		self.data = numpy.array(self.data)
		b = numpy.zeros((self.Natoms,5))

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
		self.OID = numpy.zeros(self.NO)
		n = 0
		for i in range(self.NSi,self.NSi+self.NO):
			self.OID[n] = self.data[i][4]
			n+=1

		# print self.OID

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

		outfile = 'terminalO.txt'
		f = open(outfile, 'w')

		for i in range(self.NO):
			NO = numpy.sum(self.CM[:,i])
			if  NO == 0:
				# free O
				f.write('%d\n' %self.OID[i])
			elif NO == 1:
				# non-bridging O
				f.write('%d\n' %self.OID[i])
		f.close()

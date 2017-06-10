import copy
import glob
import numpy
import os

''' finds the terminal O atoms situated at pore surfaces '''

class PoreSurface:

	def __init__(self, inputfile, atomData, terminalOH, terminalGroups, dimensions, cutoff=5.0):
		self.Lx, self.Ly, self.Lz = dimensions

		print 'Loading porogen data...'
		trial = self.GetTrial(inputfile)
		porefile = self.GetPoreFile(trial)
		porogenAtoms = self.LoadData(porefile)

		print 'Finding terminal O on pore surface...'
		poresurface_OH = self.ComputeGroupsOnSurface(atomData, terminalOH, porogenAtoms, cutoff)
		poresurface = self.ComputeGroupsOnSurface(atomData, terminalGroups, porogenAtoms, cutoff)

		print 'Writing out results...'
		self.AnalyzeResults(poresurface, poresurface_OH, terminalOH, terminalGroups, trial, cutoff)

	def GetTrial(self, inputfile):
		# use the last '_' to get the trial 
		idx = inputfile.rfind('_')
		# remove the .xyz from the string
		trial = inputfile[idx+1 : -4]
		return trial

	def GetPoreFile(self, trial):
		porogenfiles = glob.glob('porogen_*.xyz')
		for porefile in porogenfiles:
			if trial in porefile:
				return porefile


	def LoadData(self,inputfile):
		""" reads LAMMPS file, sorts data, and computes cell dimensions """ 
		data = []

		f = open(inputfile)
		Natoms = int(f.readline())

		f.readline()

		while True:
			fields = f.readline().strip().split()
			if fields:
				atomtype = int(fields[0])
				xcoord = float(fields[1])
				ycoord = float(fields[2])
				zcoord = float(fields[3])
				# populate the data array
				data.append([atomtype,xcoord,ycoord,zcoord])
			else:
				break
		f.close()

		return data


	def Distance(self,dx,dy,dz):
		return (dx*dx + dy*dy + dz*dz)**0.5


	def isOnSurface(self, atomData, terminal_atom, porogenAtoms, cutoff):
		for porogen_atom in porogenAtoms:

			dx = abs(atomData[terminal_atom][1] - porogen_atom[1])
			dy = abs(atomData[terminal_atom][2] - porogen_atom[2])
			dz = abs(atomData[terminal_atom][3] - porogen_atom[3])

			d = self.Distance(dx,dy,dz)
			if d < cutoff:
				return True

			# consider PBCs
			if dx > self.Lx - cutoff:
				d = self.Distance(dx + self.Lx, dy, dz)
				if d < cutoff:
					return True

				d = self.Distance(dx - self.Lx, dy, dz)
				if d < cutoff:
					return True

			if dy > self.Ly - cutoff:
				d = self.Distance(dx ,dy + self.Ly, dz)
				if d < cutoff:
					return True

				d = self.Distance(dx, dy - self.Ly, dz)
				if d < cutoff:
					return True

			if dz > self.Lz - cutoff:
				d = self.Distance(dx, dy, dz + self.Lz)
				if d < cutoff:
					return True

				d = self.Distance(dx, dy, dz - self.Lz)
				if d < cutoff:
					return True

		# default is to return false
		return False


	def ComputeGroupsOnSurface(self, atomData, terminalGroups, porogenAtoms, cutoff):
		''' find terminal groups that are within a cutoff distance of porogen atoms; 
			we assume this to mean they are on the pore surface '''

		poresurface = []

		for terminal_group in terminalGroups:
			if type(terminal_group) == type(int()):
				terminal_atom = terminal_group
				if self.isOnSurface(atomData, terminal_atom, porogenAtoms, cutoff):
					poresurface.append(terminal_group)

			else:
				for terminal_atom in terminal_group:
					if self.isOnSurface(atomData, terminal_atom, porogenAtoms, cutoff):
						poresurface.append(terminal_group)
						# if one of atoms in the terminal group is found to 
						# be on the pore surface, don't need to check the others
						break

		return poresurface

	def AnalyzeResults(self, poresurface, poresurface_OH, terminalOH, terminalGroups, trial, cutoff):
		# write out the fraction of terminal O on pore surfaces
		fname = 'pore_surface_stats.txt'

		# self.WritePoreSurfaceStats(poresurface,terminalO,fname0,porecutoff)

		N_termGroups_on_surface = len(poresurface)
		N_termGroups = len(terminalGroups)
		N_termOH_on_surface = len(poresurface_OH)
		N_termOH = len(terminalOH)

		frac_surface = (float(N_termOH_on_surface) + float(N_termGroups_on_surface)) /\
						(N_termOH + N_termGroups)

		f = open(fname, 'a')
		f.write('\n\nTrial %s:\n' %trial)
		f.write('Cutoff = %.4f\n' %cutoff)
		f.write('Number of terminal OH = %d\n' %N_termOH)
		f.write('Number of terminal groups = %d\n' %N_termGroups)
		f.write('Number of terminal OH on pore surface = %d\n' %N_termOH_on_surface)
		f.write('Number of terminal groups on pore surface = %d\n' %N_termGroups_on_surface)
		f.write('Fraction terminal groups on pore surface = %.6f\n' %frac_surface)
		f.close()















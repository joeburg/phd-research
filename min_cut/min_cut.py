''' This class determines the min-cut of a material given a LAMMPS .xyz file '''

import graph_tool.all as gt
import numpy
import scipy.interpolate
import scipy.spatial
import time
import yaml

# from matplotlib.pyplot import plot_surface
# from mpl_toolkits.mplot3d import Axes3D

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

class MinCut:

	def __init__(self,inputfile,delta):
		self.inputfile = inputfile

		# initialize the graph; object with vertex and edge attributes 
		self.graph = gt.Graph()

		# initialize the number of verticies in the graph
		self.Nvertices = 0

		# the edges are weighted by the bond rupture force of the respective bond type
		self.E_SiC_bond = 2.8 # weight of Si-C bonds in pN
		self.E_SiO_bond = 3.3 # weight of Si-O bonds in pN
		self.E_CC_bond = 4.1 # weight of C-C bonds in pN
		self.E_interface = 10 # weight of bonds between source/sink and network

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
		self.x_min = 0
		self.x_max = 0
		self.z_min = 0
		self.z_max = 0
		self.y_avg = 0 


		# run the associated methods with computing the min-cut and fracture path
		# load the .xyz file 
		self.LoadData(inputfile)
		print 'Data loaded.'
		print 'Cell lengths = %.2f, %.2f nm' %(self.Lx*0.1, self.Lz*0.1)

		# decompose the simulation cell into the source, sink and network 
		print 'Decomposing the simulation cell into a flow network...'
		source, network, sink = self.formFlowNetwork(delta)
		decomposition = (len(source),len(network),len(sink))
		print '\nAtoms in source = %d' %len(source)
		print 'Atoms in network = %d' %len(network)
		print 'Atoms in sink = %d\n' %len(sink)

		# updates the graph with the atoms as vertices
		print 'Transforming the atoms into verticies...'
		source_vert, sink_vert = self.computeVerticies(source, network, sink)

		# updates the graph property, positions
		print 'Transforming the atom positions as a graph property...'
		positions, atomIDs = self.computePositions(source, network, sink) 

		# computes the bonds in the network and updates the graph with the bonds as edges
		print 'Computing the bonds in the network...'
		capacity = self.computeEdges(source, network, sink)
		print '\nNumber of edges = %d' %len(list(self.graph.edges()))
		print 'Number of vertices = %d\n' %len(list(self.graph.vertices()))

		# computes the min-cut of the graph between the source and sink
		print 'Computing the min-cut of the graph...'
		residual,max_flow,partition,mincut = self.computeMinCut(source_vert, sink_vert, capacity)

		# use this method to visualize the min-cut partition
		# self.visualizeGraphPartition(positions, partition, residual, capacity)

		# computes the set of bonds broken by the min-cut and the associated atoms
		print 'Analyzing the 3D crack path...'
		crack_path_results = self.AnalyzeCrackPath(partition, positions, atomIDs, delta)

		# # compute the Gc of the material
		Gc = 0
		# Gc = self.computeGc(bond_densities)
		# print 'Gc = %.2f J/m^2' %Gc
		

		# write out the results 
		self.WriteResults(delta, decomposition, crack_path_results, Gc)


	#------------------------------------------------------------------------------#

	def LoadData(self,inputfile):
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

		# sort the rows by atom type: Si (1), O (2), then C (3)
		self.data = numpy.array(sorted(self.data, key=lambda a_entry: a_entry[0]))

		# find the dimensions of the simulation cell 
		self.x_min = min(self.data[:,1])
		self.x_max = max(self.data[:,1])
		self.Lx = self.x_max - self.x_min

		self.z_min = min(self.data[:,3])
		self.z_max = max(self.data[:,3])
		self.Lz = self.z_max - self.z_min

		y_min = min(self.data[:,2])
		y_max = max(self.data[:,2])
		self.Ly = y_max - y_min
		self.y_avg = 0.5*(y_max + y_min)

	#------------------------------------------------------------------------------#
	''' methods associated with decomposing the simiulation cell into a flow network ''' 

	def formRegion(self, atom_idx_list):
		region = numpy.zeros((len(atom_idx_list),4))

		for i in range(len(atom_idx_list)):
			region[i] = self.data[atom_idx_list[i]]

		# sort the rows by atom type: Si (1), O (2), then C (3)
		region = numpy.array(sorted(region, key=lambda a_entry: a_entry[0]))

		return region

	def formFlowNetwork(self,delta):
		''' decomposes the simulation cell into the source, network and sink given 
		the min-cut height, delta; source/network/sink are arrays of type,x,y,z as in self.data'''

		source_idx = [] 
		network_idx = []
		sink_idx = []

		for i in range(self.Natoms):
			ycoord = self.data[i][2]

			# the source is all atoms above yavg + delta/2 
			if ycoord > self.y_avg + delta*0.5:
				source_idx.append(i)

			# the sink is all atoms below yavg - delta/2 
			elif ycoord < self.y_avg - delta*0.5:
				sink_idx.append(i)

			# all other atoms are within the network
			else:
				network_idx.append(i)

		source = self.formRegion(source_idx)
		network = self.formRegion(network_idx)
		sink = self.formRegion(sink_idx)

		return (source, network, sink)

	#------------------------------------------------------------------------------#
	''' methods associated with computing the graph verticies '''

	def computeVerticies(self, source, network, sink):
		# add the number of vertices in the network to the graph in addition to the 
		# source and sink verticies; Nvertices = Nnetwork_atoms + 2
		self.Nvertices = len(network) + 2

		self.graph.add_vertex(self.Nvertices)

		# we will define the source and sink verticies to be 0 and 1 respectively
		source_vert = self.graph.vertex(0)
		sink_vert = self.graph.vertex(1)

		return (source_vert, sink_vert)


	#------------------------------------------------------------------------------#
	''' methods associated with computing the graph positions '''

	def computeAvgPos(self, region):
		# region data rows are formmatted as type,x,y,z
		xavg = 0.5*(max(region[:,1]) + min(region[:,1]))
		yavg = 0.5*(max(region[:,2]) + min(region[:,2]))
		zavg = 0.5*(max(region[:,3]) + min(region[:,3]))

		return numpy.array([xavg, yavg, zavg])


	def computePositions(self, source, network, sink):
		# associate the x,y,z coordinates with each vertex; use vector<float> property map
		positions = self.graph.new_vertex_property('vector<float>')
		atomIDs = self.graph.new_vertex_property('int')

		# take the position of the source to be the average position of all the source atoms
		# give source atomID = 4
		positions[self.graph.vertex(0)] = self.computeAvgPos(source)
		atomIDs[self.graph.vertex(0)] = 4

		# take the position of the sink to be the average position of all the sink atoms
		# give sink atomID = 5
		positions[self.graph.vertex(1)] = self.computeAvgPos(sink)
		atomIDs[self.graph.vertex(1)] = 5

		# take the position of the network atoms from the data array
		for i in range(len(network)):
			# the network verticies begin at index 2
			positions[self.graph.vertex(i+2)] = network[i][1:]
			atomIDs[self.graph.vertex(i+2)] = network[i][0]

		return (positions, atomIDs)

	#------------------------------------------------------------------------------#
	''' methods associated with computing the graph edges '''

	def Distance(self,dx,dy,dz):
		return (dx*dx + dy*dy + dz*dz)**0.5


	def computeNatoms(self, region):
		NSi = 0
		NO = 0
		NC = 0

		for atom in region:
			atomtype = atom[0]
			# determine type of atom
			if atomtype == 1:
				NSi += 1
			elif atomtype == 2:
				NO += 1
			elif atomtype == 3:
				NC += 1

		Natoms = NSi + NO + NC
		return (Natoms, NSi, NO, NC)


	def isBond(self,atom1,atom2):
		''' computes if 2 atoms are bonded; uses periodic boundary conditions in the 
		x and z directions; the y direction is the min-cut axis so no PBCs needed '''
		dx = abs(atom1[1] - atom2[1])
		dy = abs(atom1[2] - atom2[2])
		dz = abs(atom1[3] - atom2[3])

		if self.Distance(dx,dy,dz) < 2.0:
			return True
		
		# consider PBCs in x and z directions; not in y-dir due to source and sink
		if dx > self.Lx - 2.0:
			if self.Distance(dx + self.Lx, dy, dz) < 2.0:
				return True

			if self.Distance(dx - self.Lx, dy, dz) < 2.0:
				return True

		if dz > self.Lz - 2.0:
			if self.Distance(dx, dy, dz + self.Lz) < 2.0:
				return True

			if self.Distance(dx, dy, dz - self.Lz) < 2.0:
				return True

		return False

	def addBond(self, vert1, vert2, E_bond, capacity):
		# add directed edge in both directions 
		e1 = self.graph.add_edge(self.graph.vertex(vert1), self.graph.vertex(vert2))
		e2 = self.graph.add_edge(self.graph.vertex(vert2), self.graph.vertex(vert1))

		# add the capacity of each edge 
		capacity[e1] = E_bond
		capacity[e2] = E_bond

		return capacity


	def computeNetworkBonds(self, capacity, network):
		''' computes the bonds within the network region of the simulation cell
		and correspondingly adds edges to the graph with the respective capacities '''
		Natoms, NSi, NO, NC = self.computeNatoms(network)

		# compute the Si-O bonds
		for i in range(NSi):
			atom1 = network[i]
			Nbonds = 0
			for j in range(NSi,NSi+NO):
				atom2 = network[j]

				if self.isBond(atom1,atom2):
					# network verts start at index 2
					capacity = self.addBond(i+2, j+2, self.E_SiO_bond, capacity)

					Nbonds += 1

				# can have a maximum of 3 Si-O bonds
				if Nbonds == 3:
					break

		# compute the Si-C bonds
		for i in range(NSi):
			atom1 = network[i]
			for j in range(NSi+NO,Natoms):
				atom2 = network[j]

				if self.isBond(atom1,atom2):
					capacity = self.addBond(i+2, j+2, self.E_SiC_bond, capacity)
					break # only 1 possible Si-C bond

		# compute the C-C bonds
		for i in range(NSi+NO,Natoms):
			atom1 = network[i]
			for j in range(i+1,Natoms):
				atom2 = network[j]

				if self.isBond(atom1,atom2):
					capacity = self.addBond(i+2, j+2, self.E_CC_bond, capacity)
					break # only 1 possible C-C bond

		return capacity


	def computeInterfaceBonds(self, capacity, network, region, region_vert_idx):
		''' computes the bonds that cross the slicing plane between the source/sink region
		and the network region; the source/sink region maps to a single vertex; the edges and 
		capacities are added accordingly '''

		Natoms_n, NSi_n, NO_n, NC_n = self.computeNatoms(network)
		Natoms_r, NSi_r, NO_r, NC_r = self.computeNatoms(region)

		# compute the Si-O bonds; Si in region; O in network
		for i in range(NSi_r):
			atom1 = region[i]
			Nbonds = 0
			for j in range(NSi_n,NSi_n+NO_n):
				atom2 = network[j]

				if self.isBond(atom1,atom2):
					# add directed edge in both directions (network verts start at index 2)
					capacity = self.addBond(region_vert_idx, j+2, self.E_interface, capacity)

					Nbonds += 1

				# can have a maximum of 3 Si-O bonds
				if Nbonds == 3:
					break

		# compute the Si-O bonds; O in region; Si in network
		for i in range(NSi_r,NSi_r+NO_r):
			atom1 = region[i]
			Nbonds = 0
			for j in range(NSi_n):
				atom2 = network[j]

				if self.isBond(atom1,atom2):
					# add directed edge in both directions (network verts start at index 2)
					capacity = self.addBond(region_vert_idx, j+2, self.E_interface, capacity)

					Nbonds += 1

				# can have a maximum of 3 Si-O bonds
				if Nbonds == 3:
					break	

		# compute the Si-C bonds; Si in region; C in network
		for i in range(NSi_r):
			atom1 = region[i]
			for j in range(NSi_n+NO_n,Natoms_n):
				atom2 = network[j]

				if self.isBond(atom1,atom2):
					# add directed edge in both directions (network verts start at index 2)
					capacity = self.addBond(region_vert_idx, j+2, self.E_interface, capacity)
					break # only 1 possible Si-C bond

		# compute the Si-C bonds; C in region; Si in network
		for i in range(NSi_r+NO_r,Natoms_r):
			atom1 = region[i]
			for j in range(NSi_n):
				atom2 = network[j]

				if self.isBond(atom1,atom2):
					# add directed edge in both directions (network verts start at index 2)
					capacity = self.addBond(region_vert_idx, j+2, self.E_interface, capacity)
					break # only 1 possible Si-C bond

		# compute the C-C bonds; C in region; C in network
		for i in range(NSi_r+NO_r,Natoms_r):
			atom1 = region[i]
			for j in range(NSi_n+NO_n,Natoms_n):
				atom2 = network[j]

				if self.isBond(atom1,atom2):
					# add directed edge in both directions (network verts start at index 2)
					capacity = self.addBond(region_vert_idx, j+2, self.E_interface, capacity)
					break # only 1 possible Si-C bond

		return capacity


	def computeEdges(self, source, network, sink):
		# each edge has a capacity according to its respective bond rupture force; use a property map
		capacity =  self.graph.new_edge_property('double')

		# compute the bonds in the network atoms
		capacity = self.computeNetworkBonds(capacity, network)

		# compute the bonds between the source and network atoms
		capacity = self.computeInterfaceBonds(capacity, network, source, 0)

		# compute the bonds between the sink and network atoms
		capacity = self.computeInterfaceBonds(capacity, network, sink, 1)

		return capacity	
		
	#------------------------------------------------------------------------------#	

	def computeMinCut(self, source_vert, sink_vert, capacity):
		''' computes the max flow, min-cut, and partition set of the graph '''

		residual = gt.push_relabel_max_flow(self.graph, source_vert, sink_vert, capacity)
		max_flow = sum(residual[e] for e in sink_vert.in_edges())

		partition = gt.min_st_cut(self.graph, source_vert, capacity, residual)
		mincut = sum([capacity[e] - residual[e] for e in self.graph.edges() if partition[e.source()] != partition[e.target()]]) 

		return (residual,max_flow,partition,mincut)

	#------------------------------------------------------------------------------#
	''' methods associated with analyzing the 3D crack path '''

	def getBrokenBonds(self, partition):
		# for each edge in the graph, if the source and target are in the same partition
		# then the bond was not broken; if they are in different partitions, 
		# then the bond was broken

		broken_bonds = set([])
		# atoms_with_bb = set([]) # atoms that participate in broken bonds

		for edge in self.graph.edges():
			source_vert = edge.source()
			target_vert = edge.target()

			if partition[source_vert] != partition[target_vert]:
				broken_bonds.add(edge)
				# atoms_with_bb.add(source_vert)
				# atoms_with_bb.add(target_vert)

		return broken_bonds


	def typesBrokenBonds(self, broken_bonds, atomIDs):

		NSiO_bb = 0
		NSiC_bb = 0
		NCC_bb = 0

		for edge in broken_bonds:
			atom_src = atomIDs[edge.source()]
			atom_tgt = atomIDs[edge.target()]

			if (atom_src == 1 and atom_tgt == 2) or (atom_src == 2 and atom_tgt == 1):
				NSiO_bb += 1

			elif (atom_src == 1 and atom_tgt == 3) or (atom_src == 3 and atom_tgt == 1):
				NSiC_bb += 1

			elif atom_src == 3 and atom_tgt == 3:
				NCC_bb += 1

		# the number of bonds is 1/2 of the computed edges since we used a directed 
		# graph with edges in both directions between atoms 
		NSiO_bb = 0.5*NSiO_bb
		NSiC_bb = 0.5*NSiC_bb
		NCC_bb = 0.5*NCC_bb
		return (NSiO_bb, NSiC_bb, NCC_bb)


	def getSmaller(self, y_atom, y_min):
		if y_atom < y_min:
			return y_atom
		else:
			return y_min

	def getLarger(self, y_atom, y_max):
		if y_atom > y_max:
			return y_atom
		else:
			return y_max

	def averagePoint(self, point1, point2):
		x1, y1, z1 = point1
		x2, y2, z2 = point2
		return numpy.array([0.5*(x1+x2), 0.5*(y1+y2), 0.5*(z1+z2)])


	def mincutHeight(self, broken_bonds, positions):
		y_min = 10000
		y_max = 0

		for edge in broken_bonds:
			y_atom1 = positions[edge.source()][1]
			y_atom2 = positions[edge.target()][1]

			y_min = self.getSmaller(y_atom1, y_min)
			y_min = self.getSmaller(y_atom2, y_min)

			y_max = self.getLarger(y_atom1, y_max)
			y_max = self.getLarger(y_atom2, y_max)

		return (y_max - y_min)


	def WriteSurface(self, x, y, z, filename):
		data = {'x' : x.tolist(), 'y' : y.tolist(), 'z' : z.tolist()}
		with open(filename, 'w') as f:
			yaml.dump(data,f)


	def sideLength(self, point1, point2):
		dx = abs(point1[0] - point2[0])
		dy = abs(point1[1] - point2[1])
		dz = abs(point1[2] - point2[2])
		return 	self.Distance(dx, dy, dz)


	def computeTriangleArea(self, point1, point2, point3):
		d1 = self.sideLength(point1, point2)
		d2 = self.sideLength(point1, point3)
		d3 = self.sideLength(point2, point3)

		s = 0.5*(d1 + d2 + d3)
		return (s*abs(s-d1)*abs(s-d2)*abs(s-d3))**0.5


	def generateFractureSurfMesh(self, points):
		# generate a mesh with the fracture points to create a surface (xz-plane)
		x = bb_points[:,0]
		y = bb_points[:,1]
		z = bb_points[:,2]

		points = numpy.column_stack((x,z))
		grid_x, grid_z = numpy.mgrid[self.x_min:self.x_max:1000j, self.z_min:self.z_max:1000j]

		Y1 = scipy.interpolate.griddata(points, y, (grid_x, grid_z), method='nearest')
		Y2 = scipy.interpolate.griddata(points, y, (grid_x, grid_z), method='linear')
		Y3 = scipy.interpolate.griddata(points, y, (grid_x, grid_z), method='cubic')

		self.WriteSurface(grid_x, Y1, grid_z, 'fracture_surface_nearest.yml')
		self.WriteSurface(grid_x, Y2, grid_z, 'fracture_surface_linear.yml')
		self.WriteSurface(grid_x, Y3, grid_z, 'fracture_surface_cubic.yml')


	def generateFractureSurface(self, broken_bonds, positions):
		''' generates a fracture plane based on the bonds broken '''

		# take the midpoints between the atoms that participate in a broken bond
		bb_points = numpy.zeros((len(broken_bonds), 3))

		i = 0
		for edge in broken_bonds:
			atom1 = positions[edge.source()]
			atom2 = positions[edge.target()]

			bb_points[i] = self.averagePoint(atom1, atom2)
			i += 1

		# use the average y coordinate in the fracture surface when adding boundary points
		y_avg = 0.5*(min(bb_points[:,1]) + max(bb_points[:,1]))

		# make boundaries around the xz plane to perform "contrained" Delaunay triangulation;
		# scipy.spatial.Delaunay does not support constrained triangulations so we add a dense arrays 
		# of points around the boundary; 
		# this ensures that the resulting area computed is not underestimated

		n = 20 # define number of points on each edge 
		# n points = 4*(n-2)+4 so as not to repeat the corner points
		BCs = numpy.zeros((4*(n-2)+4, 3))

		# left side 
		i = 0
		for z in numpy.linspace(self.z_min, self.z_max, num=n):
			BCs[i] = numpy.array([self.x_min, y_avg, z])
			i += 1

		# top side 
		for x in numpy.linspace(self.x_min+1, self.x_max, num=n-1): # dont repeat the corner point
			BCs[i] = numpy.array([x, y_avg, self.z_max])
			i += 1

		# right side 
		for z in numpy.linspace(self.z_max-1, self.z_min, num=n-1): # dont repeat the corner point
			BCs[i] = numpy.array([self.x_max, y_avg, z])
			i += 1

		# bottom side
		for x in numpy.linspace(self.x_max-1, self.x_min+1, num=n-2): # dont repeat the corner points
			BCs[i] = numpy.array([x, y_avg, self.z_min])
			i += 1


		# use the xz-plane projections to compute the contrained Delaunay triangulation with the BCs
		bb_points_BCs = numpy.concatenate((BCs, bb_points), axis=0)

		triangulation = scipy.spatial.Delaunay(bb_points_BCs[:,[0,2]])
		Nsimplices = len(triangulation.simplices)

		# given the trianuglation of the points, compute the area of each triangle
		surf_area = 0
		for simplex in triangulation.simplices:
			point1 = bb_points_BCs[simplex[0]]
			point2 = bb_points_BCs[simplex[1]]
			point3 = bb_points_BCs[simplex[2]]

			surf_area += self.computeTriangleArea(point1, point2, point3)

		# convert A^2 to nm^2
		surf_area = 0.01*surf_area
		
		# use this method to generate fracture surfaces meshes
		# self.generateFractureSurfMesh(bb_points_BCs)

		# use this method to visualize the delaunay triangulation 
		self.visualizeFractureSurfacePoints(bb_points_BCs[:,[0,2]])

		# use this method to write out positions on the surface
		# self.WriteFracSurftoVMD(bb_points)

		return (surf_area, Nsimplices)


	def computeBondDensity(self, Nbond_types, surf_area):
		NSiO_bb, NSiC_bb, NCC_bb = Nbond_types

		SiO_rho = float(NSiO_bb)/surf_area
		SiC_rho = float(NSiC_bb)/surf_area
		CC_rho = float(NCC_bb)/surf_area
		return (SiO_rho, SiC_rho, CC_rho)


	def AnalyzeCrackPath(self, partition, positions, atomIDs, delta):

		broken_bonds = self.getBrokenBonds(partition)

		Nbond_types = self.typesBrokenBonds(broken_bonds, atomIDs)
		print '\nNumber of Si-O bonds broken = %d' %Nbond_types[0]
		print 'Number of Si-C bonds broken = %d' %Nbond_types[1]
		print 'Number of C-C bonds broken = %d' %Nbond_types[2]

		surf_area, Nsimplices = self.generateFractureSurface(broken_bonds, positions)
		print '\nNumber of simplicies = %d' %Nsimplices
		print 'Fracture surface area = %.2f nm^2' %surf_area

		bond_densities = self.computeBondDensity(Nbond_types, surf_area)
		print '\nTotal fracture bond density = %.2f nm^-2' %(sum(bond_densities))
		print 'Si-O fracture bond density = %.2f nm^-2' %bond_densities[0]
		print 'Si-C fracture bond density = %.2f nm^-2' %bond_densities[1]
		print 'C-C fracture bond density = %.2f nm^-2' %bond_densities[2]

		h_mincut = self.mincutHeight(broken_bonds, positions)
		print '\nmin-cut height, h = %.2f A\n' %h_mincut

		self.WriteBBtoVMD(broken_bonds, atomIDs, positions, delta)

		return (Nbond_types, surf_area, Nsimplices, bond_densities, h_mincut)


	#------------------------------------------------------------------------------#
	def computeGc(self, bond_densities):
		rho = sum(bond_densities) # in nm^-2






	#------------------------------------------------------------------------------#
	''' methods to write out and visualize data '''

	def visualizeFractureSurfacePoints(self, points):
		graph_del, pos_del = gt.triangulation(points, type='delaunay')

		cap = graph_del.new_edge_property("double")
		for edge in graph_del.edges():
			cap[edge] = 10

		gt.graph_draw(graph_del, pos=pos_del,output="delaunay_frac_surf.pdf")


	def getAtomLabel(self, atomID):
		if atomID == 1:
			return 'Si'
		elif atomID == 2:
			return 'O'
		elif atomID == 3:
			return 'C'

	def WriteResults(self, delta, decomposition, crack_path_results, Gc):
		Nsource, Nnetwork, Nsink = decomposition
		Nbond_types, surf_area, Nsimplices, bond_densities, h_mincut = crack_path_results
		NSiO_bb, NSiC_bb, NCC_bb = Nbond_types
		SiO_rho, SiC_rho, CC_rho = bond_densities

		filename = '{}_mincut_results_{}.txt'.format(self.inputfile[:-4],delta)
		f = open(filename, 'w')
		f.write('Results for min-cut analysis:\n\n')
		f.write('min-cut height, delta, = %.2f A\n' %delta)
		f.write('Natoms in source = %d\n' %Nsource)
		f.write('Natoms in network = %d\n' %Nnetwork)
		f.write('Natoms in sink = %d\n' %Nsink)
		f.write('Number of Si-O bonds broken = %d\n' %NSiO_bb)
		f.write('Number of Si-C bonds broken = %d\n' %NSiC_bb)
		f.write('Number of C-C bonds broken = %d\n' %NCC_bb)
		f.write('Number of simplices = %d\n' %Nsimplices)
		f.write('Simulation cell area = %.2f nm^2\n' %(self.Lx*self.Lz*0.01))
		f.write('Fracture surface area = %.2f nm^2\n' %surf_area)
		f.write('Total fracture bond density = %.2f nm^-2\n' %(sum(bond_densities)))
		f.write('Si-O fracture bond density = %.2f nm^-2\n' %SiO_rho)
		f.write('Si-C fracture bond density = %.2f nm^-2\n' %SiC_rho)
		f.write('C-C fracture bond density = %.2f nm^-2\n' %CC_rho)
		f.write('fracture (min-cut) height, h, = %.2f nm\n' %(h_mincut*0.1))
		f.write('Gc = %.2f J/m^2' %Gc)
		f.close()


	def WriteFracSurftoVMD(self, points):
		''' write the midpoints (in the fracture plane) of the bonds to VMD format '''

		filename = '{}_frac_plane_VMD.xyz'.format(self.inputfile[:-4])
		f = open(filename, 'w')
		f.write('%d\n' %(len(points)))
		f.write('Atoms\n')
		for point in points:
			x, y, z = point
			f.write('%s  %.4f  %.4f  %.4f\n' % ('Fr',x,y,z))
		f.close()		


	def WriteBBtoVMD(self, broken_bonds, atomIDs, positions, delta):
		''' write the atoms comprising the broken bonds to VMD format '''

		filename = '{}_crack_path_VMD_{}.xyz'.format(self.inputfile[:-4], delta)
		f = open(filename, 'w')
		f.write('%d\n' %(2*len(broken_bonds)))
		f.write('Atoms\n')
		for edge in broken_bonds:
			atom_src = self.getAtomLabel(atomIDs[edge.source()])
			atom_tgt = self.getAtomLabel(atomIDs[edge.target()])

			atom_src_pos = positions[edge.source()]
			atom_tgt_pot = positions[edge.target()]

			f.write('%s  %.4f  %.4f  %.4f\n' % (atom_src, atom_src_pos[0], 
				atom_src_pos[1],atom_src_pos[2]))

			f.write('%s  %.4f  %.4f  %.4f\n' % (atom_tgt, atom_tgt_pot[0], 
				atom_tgt_pot[1],atom_tgt_pot[2]))

		f.close()


	def visualizeGraphPartition(self, positions, partition, residual, capacity):
		ofile = '{}_graph_partition.png'.format(self.inputfile[:-4])
		gt.graph_draw(self.graph, pos=positions, edge_pen_width=gt.prop_to_size(capacity, mi=0.5, ma=1, power=1),
			edge_text=residual, vertex_fill_color=partition, vertex_text=self.graph.vertex_index,
			vertex_font_size=2, edge_font_size=0, output=ofile)




















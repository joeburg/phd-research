import graph_tool.all as gt

import numpy

# from numpy.random import seed, random
# from scipy.linalg import norm

# gt.seed_rng(42)
# seed(42)
# points = random((400, 2))

# points[0] = [0, 0]
# points[1] = [1, 1]
# g, pos = gt.triangulation(points, type="delaunay")

# print g
# print pos

# g.set_directed(True)
# edges = list(g.edges())
# # reciprocate edges
# for e in edges:
#    g.add_edge(e.target(), e.source())
# # The capacity will be defined as the inverse euclidean distance
# cap = g.new_edge_property("double")
# for e in g.edges():
#     cap[e] = min(1.0 / norm(pos[e.target()].a - pos[e.source()].a), 10)
# # g.edge_properties["cap"] = cap
# # g.vertex_properties["pos"] = pos
# # g.save("flow-example.xml.gz")
# gt.graph_draw(g, pos=pos, edge_pen_width=gt.prop_to_size(cap, mi=0, ma=3, power=1),output="flow-example.pdf")

class MinCut:

	def __init__(self,inputfile,molType):

		if molType not in set(['EtOCS', 'EtOCSMe', '135Benz','SiO2']):
			raise RuntimeError, 'Invalid precursor type.'
		else:
			# precursor type: Et-OCS, Et-OCS(Me), 135-Benzene, SiO2
			self.molType = molType 

		 # initialize the graph; object with vertex and edge attributes 
		self.graph = gt.Graph()

		# associate the x,y,z coordinates with each vertex; use vector<float> property map
		self.positions = self.graph.new_vertex_property('vector<float>')

		# each edge has a capacity according to its respective bond rupture force; use a property map
		self.capacity =  self.graph.new_edge_property('double')

		self.E_SiC_bond = 5000001 # weight of Si-C bonds
		self.E_SiO_bond = 6000100 # weight of Si-O bonds
		self.E_CC_bond = 7000001 # weight of C-C bonds
		self.E_bd = 100000000 # weight of network bonds and atoms in source and sink

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
		self.y_avg = 0 


		self.LoadData(inputfile)
		print 'Data loaded.'

		self.TestGraph()

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
		self.Lx = max(self.data[:,1]) - min(self.data[:,1])
		self.Lz = max(self.data[:,3]) - min(self.data[:,3])

		y_min = min(self.data[:,2])
		y_max = max(self.data[:,2])
		self.Ly = y_max - y_min
		self.y_avg = 0.5*(y_max + y_min)

	def TestGraph(self):
		# add the atoms as verticies in the graph
		self.graph.add_vertex(self.Natoms)

		# associate the x,y,z coordinates with each vertex
		for i in range(self.Natoms):
			self.positions[self.graph.vertex(i)] = self.data[i][1:]

		for i in range(self.Natoms):
			print self.positions[i]

		# gt.graph_draw(self.graph, vertex_text=self.graph.vertex_index, vertex_font_size=18,
		# 	output_size=(200, 200), output="test.png")



#------------------------------------------------------------------------------#
# ifile = 'test.xyz'
# mol = 'EtOCS'
# MinCut(ifile,mol)


x = numpy.arange(10)

g = gt.Graph()
g.add_vertex(len(x))
# print g
e_info = g.new_edge_property('int')
i = 0
for v1 in g.vertices():
	for v2 in g.vertices():
		if v2 > v1:
			e1 = g.add_edge(g.vertex(v1), g.vertex(v2))
			e2 = g.add_edge(g.vertex(v2), g.vertex(v1))

			e_info[e1] = i
			e_info[e2] = i

			i+=1

for e in g.edges():
	print e_info[e]


# gt.graph_draw(g, vertex_text=g.vertex_index, vertex_font_size=18,output_size=(200, 200), output="test.png")

# vprop_vfloat = g.new_vertex_property('vector<float>')
# vprop_vfloat[g.vertex(8)] = [0.5,0.3,0.4]

# for v in g.vertices():
# 	print v

# print vprop_vfloat
# for v in vprop_vfloat:
# 	print v 

# print g.vertex(1)
# gt.graph_draw(g, pos=pos, edge_pen_width=gt.prop_to_size(cap, mi=0, ma=3, power=1),output="flow-example.pdf")






import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D

import numpy
import yaml


def surfacePlot(X,Y,Z,filename):
	x_min = numpy.min(X)
	x_max = numpy.max(X)

	y_min = numpy.min(Y)
	y_max = numpy.max(Y)

	plt.figure()
	plt.subplot(111)
	plt.imshow(Z.T, extent=(x_min, x_max, y_min, y_max), origin='lower')
	plt.savefig(filename)
	plt.close()

	# ax = fig.add_subplot(111,projection='3d')
	# ax.plot_surface(X,Y,Z)
	# plt.savefig(filename)
	# plt.close()


def LoadData(filename):
	with open(filename) as f:
		data = yaml.load(f)

	grid_x = numpy.array(data['x'])
	Y = numpy.array(data['y'])
	grid_z = numpy.array(data['z'])

	return (grid_x, Y, grid_z)



f1 = 'fracture_surface_nearest.yml'
f2 = 'fracture_surface_linear.yml'
f3 = 'fracture_surface_cubic.yml'

files = [f1, f2, f3]

for f in files:
	ofile = '{}.png'.format(f[:-4])
	grid_x, Y, grid_z = LoadData(f)
	surfacePlot(grid_x, grid_z, Y, ofile)
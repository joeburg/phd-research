##import numpy as np
##from mpl_toolkits.mplot3d import Axes3D
##import matplotlib.pyplot as plt
##
##def randrange(n, vmin, vmax):
##    return (vmax-vmin)*np.random.rand(n) + vmin
##
##fig = plt.figure()
##ax = fig.add_subplot(111, projection='3d')
##n = 100
##for c, m, zl, zh in [('r', 'o', -50, -25), ('b', '^', -30, -5)]:
##    xs = randrange(n, 23, 32)
##    ys = randrange(n, 0, 100)
##    zs = randrange(n, zl, zh)
##    ax.scatter(xs, ys, zs, c=c, marker=m)
##
##ax.set_xlabel('X Label')
##ax.set_ylabel('Y Label')
##ax.set_zlabel('Z Label')
##
##plt.show()



import mpl_toolkits.mplot3d as a3
import matplotlib.colors as colors
import pylab as pl
import scipy as sp

ax = a3.Axes3D(pl.figure())
for i in range(10000):
    vtx = sp.rand(3,3)
    tri = a3.art3d.Poly3DCollection([vtx])
    tri.set_color(colors.rgb2hex(sp.rand(3)))
    tri.set_edgecolor('k')
    ax.add_collection3d(tri)
pl.show()

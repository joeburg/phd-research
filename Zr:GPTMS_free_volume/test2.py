import numpy as np
from scipy.spatial import ConvexHull
from scipy.spatial import Delaunay


def convex_hull_volume(pts):
    ch = ConvexHull(pts)
    dt = Delaunay(pts[ch.vertices])
    tets = dt.points[dt.simplices]
    return np.sum(tetrahedron_volume(tets[:, 0], tets[:, 1],
                                     tets[:, 2], tets[:, 3]))


def tetrahedron_volume(a, b, c, d):
    return np.abs(np.einsum('ij,ij->i', a-d, np.cross(b-d, c-d))) / 6

pts = np.array([[0,0,0],[1.1311,1,1.21412],[1.432324,0,0],[0,1.34423,0],[0,0,1.41324213],[1.3232,1.42321,0],[1,0,1],[0,1,1]])

print convex_hull_volume(pts)

"""
Purpose: finds the circumcenter of an array of points
"""

import numpy
import numpy.linalg

def is_wellcentered(pts, tol=1e-8):
    """Determines whether a set of points defines a well-centered simplex.
    """
    barycentric_coordinates = circumcenter_barycentric(pts)    
    return min(barycentric_coordinates) > tol

def circumcenter_barycentric(pts):
    """Barycentric coordinates of the circumcenter of a set of points.
    
    Parameters
    ----------
    pts : array-like
        An N-by-K array of points which define an (N-1)-simplex in K dimensional space.
        N and K must satisfy 1 <= N <= K + 1 and K >= 1.

    Returns
    -------
    coords : ndarray
        Barycentric coordinates of the circumcenter of the simplex defined by pts.
        Stored in an array with shape (K,)
        
    Examples
    --------
    >>> from pydec.math.circumcenter import *
    >>> circumcenter_barycentric([[0],[4]])           # edge in 1D
    array([ 0.5,  0.5])
    >>> circumcenter_barycentric([[0,0],[4,0]])       # edge in 2D
    array([ 0.5,  0.5])
    >>> circumcenter_barycentric([[0,0],[4,0],[0,4]]) # triangle in 2D
    array([ 0. ,  0.5,  0.5])
   
    References
    ----------
    Uses an extension of the method described here:
    http://www.ics.uci.edu/~eppstein/junkyard/circumcenter.html
    """
    pts = numpy.asarray(pts)
    rows,cols = pts.shape

    assert(rows <= cols + 1)    

    A = numpy.bmat( [[ 2*numpy.dot(pts,pts.T), numpy.ones((rows,1)) ],
               [  numpy.ones((1,rows)) ,  numpy.zeros((1,1))  ]] )

    b = numpy.hstack((numpy.sum(pts * pts, axis=1),numpy.ones((1))))
    x = numpy.linalg.solve(A,b)
    bary_coords = x[:-1]  

    return bary_coords
    
def circumcenter(pts):
    """Circumcenter and circumradius of a set of points.
    
    Parameters
    ----------
    pts : array-like
        An N-by-K array of points which define an (N-1)-simplex in K dimensional space.
        N and K must satisfy 1 <= N <= K + 1 and K >= 1.

    Returns
    -------
    center : ndarray
        Circumcenter of the simplex defined by pts.  Stored in an array with shape (K,)
    radius : float
        Circumradius of the circumsphere that circumscribes the points defined by pts.
        
    Examples
    --------
    >>> circumcenter([[0],[1]])             # edge in 1D
    (array([ 0.5]), 0.5)
    >>> circumcenter([[0,0],[1,0]])         # edge in 2D
    (array([ 0.5,  0. ]), 0.5)
    >>> circumcenter([[0,0],[1,0],[0,1]])   # triangle in 2D
    (array([ 0.5,  0.5]), 0.70710678118654757)
   
    References
    ----------
    Uses an extension of the method described here:
    http://www.ics.uci.edu/~eppstein/junkyard/circumcenter.html
    """
    pts = numpy.asarray(pts)      
    bary_coords = circumcenter_barycentric(pts)
    center = numpy.dot(bary_coords,pts)
    radius = numpy.linalg.norm(pts[0,:] - center)
    return (center,radius)

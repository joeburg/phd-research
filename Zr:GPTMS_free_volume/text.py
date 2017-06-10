from random import random
from chull import Vector,Hull

##sphere = []
##for i in range(10):  
##    x,y,z = 2*random()-1,2*random()-1,2*random()-1
##    if x*x+y*y+z*z < 1.0:
##        sphere.append(Vector(x,y,z))  

pts = [Vector(1.1,2.3,3.232321),Vector(3,4,5),Vector(3,4,2),Vector(45,42,3),Vector(4,5,3),Vector(45,3,2)]



h=Hull(pts)  

h.Print()

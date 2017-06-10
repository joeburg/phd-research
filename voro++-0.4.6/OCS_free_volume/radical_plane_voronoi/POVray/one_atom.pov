#version 3.6;

// Right-handed coordinate system where the z-axis points upwards
camera {
	location <30,-50,25>
	sky z
	right -0.15*x*image_width/image_height
	up 0.15*z
	look_at <0,0,2.8>
}

// White background
background{rgb 1}

// Two lights with slightly different colors
light_source{<-8,-20,30> color rgb <0.77,0.75,0.75>}
light_source{<25,-12> color rgb <0.43,0.45,0.45>}

// The radius of the cylinders to be used when drawing the Voronoi cells
#declare r=0.08;

// Different colors for the particles
#declare f1=finish{reflection 0.17 specular 0.3 ambient 0.42}
#declare tSi=texture{pigment{rgb <1,1,0>} finish{f1}}
#declare tC=texture{pigment{rgb <0.25,0.25,0.25>} finish{f1}}
#declare tO=texture{pigment{rgb <1,0,0>} finish{f1}}

// The polydisperse particle packing
union{
#include "one_atom_p.pov"
}

// The Voronoi cells for the packing, computed using the radical Voronoi
// tessellation
union{
#include "one_atom_v.pov"
	pigment{rgb <0.75,0.75,0.75>} finish{specular 0.5 ambient 0.42}
}

#version 3.6;
#include "colors.inc"
global_settings { assumed_gamma 1.0 }
camera {
    location <60, 60, 60>
    look_at <50,50,50>
    angle 36
  }
background   { color rgb <1,1,1> }
//light_source { <70, 70, 70> color rgb <1, 1, 1> translate <-5, 5, -5> }
                
sphere{<49.4418,49.7341,46.0169>, 2 rgb <1,0,0>}
sphere{<50.7794,48.4176,44.8693>, 1.5 rgb <0,0,1>}

// sphere{<51.8612,49.5895,46.7834>,texture{ pigment{color Yellow} }}
// sphere{<49.5718,48.2845,47.377>,texture{ pigment{color Gray25} }}
// cylinder{<49.5718,48.2845,47.377>,<50.7794,48.4176,44.8693>,texture{ pigment{color Turquoise} }}
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <stdio.h>

#include "voro++.hh"
using namespace voro;

int main(int argc, char *argv[])
{
    if (argc < 3)
    {
        std::cout << "Usage:" << std::endl;
        std::cout << "  " << argv[0] << " <input file> <boundary file>" << std::endl;
        return 0;
    }
    const char *inputfile = argv[1];
    const char *BCfile = argv[2];
    
    /* make name for output file */
    const char *suffix = ".vol";
    char buf[100];
    strcpy(buf,inputfile);
    strcat(buf,suffix);
    const char *outputfile = buf;
    
    /* read x,y,z limits from boundary file */
    double xmin,xmax,ymin,ymax,zmin,zmax;
    FILE * pFile;
    pFile = fopen(BCfile,"r");
    fscanf(pFile,"%lg %lg %lg %lg %lg %lg",&xmin,&xmax,&ymin,&ymax,&zmin,&zmax);
    fclose(pFile);
     
    /* setup the number of blocks the container is divided into;
     set spacing to 2 angstroms in each direction */
    double L = xmax - xmin;
    const int Nx = int(L/8), Ny=Nx, Nz=Nx;
    
    /* create a container with the geometry and make it periodic. Allocate
     space for 8 atoms within each computational block. */
    container con(xmin,xmax,ymin,ymax,zmin,zmax,Nx,Ny,Nz,true,true,true,8);
    
    /* add the atoms to the container */
    con.import(inputfile);
    
    /* compute Voronoi cell network and write out volumes of 
     Voronoi cells to output file */
    const char *format = "%i  %x  %y  %z  %v";
    con.print_custom(format,outputfile);

    std::cout << "Successfully computed the Voronoi cell network!" << std::endl;
    
    return 0;
}

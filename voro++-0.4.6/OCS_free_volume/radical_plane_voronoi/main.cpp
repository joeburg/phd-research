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
        std::cout << "  " << argv[0] << " <input file> <boundary file> [POV output] [wall files]" << std::endl;
        return 0;
    }
    const char *inputfile = argv[1];
    const char *BCfile = argv[2];
    
    bool POV = false;
    if (argc == 4) POV = true;
    
    bool walls = false;
    const char *wallFile = NULL;
    if (argc == 5)
    {
        wallFile = argv[4];
        walls = true;
    }
    
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
    if (pFile == NULL)
    {
        std::cerr << "ERROR: Could not open file: " << *BCfile << std::endl;
        exit(1);
    }
    
    fscanf(pFile,"%lg %lg %lg %lg %lg %lg",&xmin,&xmax,&ymin,&ymax,&zmin,&zmax);
    fclose(pFile);
     
    /* setup the number of blocks the container is divided into;
     set spacing to 2 angstroms in each direction */
    double L = xmax - xmin;
    const int Nx = int(L/8), Ny=Nx, Nz=Nx;
    
    /* create a container with the geometry and make it periodic. Allocate
     space for 8 atoms within each computational block. */
    container_poly con(xmin,xmax,ymin,ymax,zmin,zmax,Nx,Ny,Nz,true,true,true,8);
    
    /* add the atoms to the container */
    con.import(inputfile);
    
    /* add walls at the boundaries container */
    if (walls)
    {
        FILE *fp;
        fp = fopen(wallFile,"r");
        
        if (fp == NULL)
        {
            std::cerr << "ERROR: Could not open file: " << *wallFile << std::endl;
            exit(1);
        }
        
        double a,b,c,d;
        unsigned int ID=0;
        while (fscanf(fp,"%lg %lg %lg %lg",&a,&b,&c,&d)!=EOF)
        {
            wall_plane plane(a,b,c,d,ID);
            con.add_wall(plane);
            ID++;
        }
        fclose(fp);
    }
        
    /* compute Voronoi cell network and write out volumes of 
     Voronoi cells to output file */
    const char *format = "%i  %x  %y  %z  %v  %r  %w  %g  %s";
    con.print_custom(format,outputfile);
    
    /* output the tessellation in POV-Ray format */
    if (POV) {
        /* make name for output pov file */
        const char *suffix = "_particles.pov";
        char buf[100];
        strcpy(buf,inputfile);
        strcat(buf,suffix);
        const char *POVparticlefile = buf;
        
         con.draw_particles_pov(POVparticlefile);
        
        const char *suffix2 = "_cells.pov";
        char buf2[100];
        strcpy(buf2,inputfile);
        strcat(buf2,suffix2);
        const char *POVcellfile = buf2;
        
        con.draw_cells_pov(POVcellfile);
    }
    
    std::cout << "Successfully computed the Voronoi cell network!" << std::endl;
    
    return 0;
}

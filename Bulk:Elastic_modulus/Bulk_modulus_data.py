#Purpose: obtain average pressure and dilitation from LAMMPS log file
# to obtain the bulk modulus; output file in mathematica format for plotting;
# also use linear regression to calculate modulus
#Joe Burg

import numpy
import sys

from scipy import stats
from Tkinter import Tk
from tkFileDialog import askopenfilename
Tk().withdraw()

#################################################################################

def load_data(filename):
    data = []
    f = open(filename)
    for line in f:
        line = line.strip().split()
        data.append(line)
    f.close()
    return data

def analyze_press_data(p_data,V0,cutoff,P_avg,Dil_avg):
    p = []
    v = []
    dilitation = []
    
    for i in range(len(p_data)):
        if i > cutoff:
            if V0:
                #print p_data[i]
                p.append(float(p_data[i][4])*-0.0001)

                dil = (float(p_data[i][8])-V0)/V0
                if abs(dil) < 0.00001:
                    dilitation.append(0.0)
                else:
                    dilitation.append(dil)
            else:
                #need to find V0 using press0 data
                v.append(float(p_data[i][8]))
                     
    if V0:
        p = numpy.array(p)
        dilitation = numpy.array(dilitation)

        pavg = numpy.mean(p)
        perror = numpy.std(p)
        dilavg = numpy.mean(dilitation)
        dilerror = numpy.std(dilitation)

        P_avg.append(pavg)
        Dil_avg.append(dilavg)
        
        f1.write('{{%.10f,%.10f}, ErrorBar[%.10f,%.10f]},'%(dilavg,pavg,dilerror,perror))
        f2.write('{%.10f,%.10f},'%(dilavg,pavg))

        return P_avg,Dil_avg
    else:
        #return V0 
        return numpy.mean(v)

def get_elastic_mod(bulk_mod):
    return 3*bulk_mod*(1-2*0.25)

################################################################################
# filename = askopenfilename()
# data = load_data(filename)

# start = input("What line does the data for 5000bar start? ")
# start -= 1

# N_press = input("How many pressure steps are there? ")
# outputfile = raw_input('Provide the output file name: ')
# cutoff = input("At what timestep do you want to cutoff the press data? ")

if len(sys.argv) < 6:
    print 'Usage:'
    print '  python %s <inputfile> <start line> <number pressure steps> <outputfile> <cutoff>' %sys.argv[0]
    exit()

inputfile = sys.argv[1]
start = int(sys.argv[2])
start -= 1
N_press = int(sys.argv[3])
outputfile = sys.argv[4]
cutoff = int(sys.argv[5])

# load input data
data = load_data(inputfile)

# order of variables -> press5000, press2500, press1500, press1000, press500
#press0, pull500, pull1000, pull1500, pull2500, pull5000
g = globals()
p_vars = []
for i in range(N_press):
    varname = 'p'+str(i)
    p_vars.append(varname)
    g[varname] = data[start:start+202]
    start += 201+28

#find V0 using press0 data
V0 = analyze_press_data(g[p_vars[N_press/2]],0,cutoff,0,0)

#run data for each pressure step
P_avg_press = [] 
Dil_avg_press = []
P_avg0 = []
Dil_avg0 = []
P_avg_pull = []
Dil_avg_pull = []
f1 = open(outputfile, 'w')
f2 = open(outputfile[:-4]+'_2.txt','w')
for i in range(N_press):
    if i < N_press/2:
        P_avg_press,Dil_avg_press = analyze_press_data(g[p_vars[i]],V0,cutoff,\
                                                      P_avg_press,Dil_avg_press)
    elif i == N_press/2:
        P_avg0,Dil_avg0 = analyze_press_data(g[p_vars[i]],V0,cutoff,\
                                                      P_avg0,Dil_avg0)
    else:
        P_avg_pull,Dil_avg_pull = analyze_press_data(g[p_vars[i]],V0,cutoff,\
                                                      P_avg_pull,Dil_avg_pull)
f1.close()
f2.close()

#add 0 press case to both compressive and tensile data
P_avg_press += P_avg0
Dil_avg_press += Dil_avg0
P_avg_pull += P_avg0
Dil_avg_pull += Dil_avg0


#use a linear regression to calculate the bulk modulus; compute elastic
#modulus from bulk modulus and write to file
#only consider volumetric strains of 4.5% or less
Dil_press = []
Dil_pull = []
P_pull = []
P_press = []
for i in range(len(Dil_avg_press)):
    if abs(Dil_avg_press[i]) < 0.045:
        Dil_press.append(Dil_avg_press[i])
        P_press.append(P_avg_press[i])

    if abs(Dil_avg_pull[i]) < 0.045:
        Dil_pull.append(Dil_avg_pull[i])
        P_pull.append(P_avg_pull[i])

compression = numpy.polyfit(numpy.array(Dil_press),numpy.array(P_press),1)
tension = numpy.polyfit(numpy.array(Dil_pull),numpy.array(P_pull),1)

bulk_mod_c = compression[0]
elastic_mod_c = get_elastic_mod(bulk_mod_c)

bulk_mod_t = tension[0]
elastic_mod_t = get_elastic_mod(bulk_mod_t)


f = open(outputfile[:-4]+'_modulus.txt','w')
f.write('Bulk modulus (compression) = %.4f\n' % bulk_mod_c)
f.write('Bulk modulus (tension) = %.4f\n' % bulk_mod_t)
f.write('Elastic modulus (compression) = %.4f\n' % elastic_mod_c)
f.write('Elastic modulus (tension) = %.4f' % elastic_mod_t)
f.close()

print 'All done!'

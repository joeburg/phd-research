import glob
import numpy
import sys
import time

#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#

def LoadData(inputfile):
	""" intelligently reads log.lammps file without user input """ 

	# end of document flag if looped simulation
	stop_flag = 'jump'

	# if not looped simulation, count if 10 consecutive blank lines and then break
	Nblanks = 0

	# flag to start looking for pressure steps
	start_flag = 'measure modulus'
	start_reading = False

	# flags used to start/stop reading pressure step
	pressure_step_start = 'Step Temp'
	pressure_step_stop = 'Loop time'
	pressure_step_active = False

	# data structure to store pressure steps
	pressure_steps = []
	pressure_step = []

	f = open(inputfile)

	while True:
		line = f.readline().strip()

		# look for end of doc flag
		if stop_flag in line:
			break

		# for log files without a stop flag, count the number of empty lines
		if not line:
			Nblanks += 1
			if Nblanks > 10:
				break
		else:
			Nblanks = 0

		# search for start flag before activating data collection
		if start_flag in line:
			start_reading = True
			continue

		# once start flag has been reached, begin collecting pressure steps
		if start_reading:
			if pressure_step_start in line:
				pressure_step_active = True
				pressure_step = [] # reinitialize data structure 
				continue

			elif pressure_step_stop in line:
				pressure_step_active = False
				pressure_steps.append(numpy.array(pressure_step))
				continue

			if pressure_step_active:
				fields = line.split()
				try:
					pressure = float(fields[4])
					volume = float(fields[8])
					pressure_step.append((pressure, volume))
				except:
					# if log data not present, it means that there was an error during the simulation
					print "\nError with %s. The simulation did not finish property. "%inputfile+\
							"Possible 'Lost Bonds, Angles or Dihedrals' error in LAMPMS.\n"
					break

	f.close()
	return pressure_steps


def getEquilibVolume(pressure_steps, cutoff):
	for step in pressure_steps:
		P_avg = numpy.mean(step[:,0])
		V_avg = numpy.mean(step[:,1])
		# look for the pressure near 0 bar
		if abs(P_avg) < cutoff: 
			return V_avg
	# if no volume close to zero, raise an exception
	raise RuntimeError, "Error: Could not find equilibrium volume!"


def shiftDataToOrigin(data, data_eq):
	''' this function shifts the data to the origin '''
	# data_eq has the form (P_eq, Dil_eq)
	P_eq, Dil_eq = data_eq

	# data has the form [(P_avg, P_std, Dil_avg, Dil_std), ]
	if P_eq > 0:
		data[:,0] = data[:,0] - P_eq
	else:
		data[:,0] = data[:,0] + abs(P_eq)

	if Dil_eq > 0:
		data[:,2] = data[:,2] - Dil_eq
	else:
		data[:,2] = data[:,2] + abs(Dil_eq)
	
	return data


def ComputeAvgPressureAvgDilitation(pressure_steps, Veq, eq_cutoff):
	# store average compression/tension values in an array 
	# [(P_avg, P_std, Dil_avg, Dil_std)]
	data_C = []
	data_T = []
	data_eq = []
	
	for step in pressure_steps:
		# convert pressure from bar to GPa and 
		# define -P as compression and +P as tension
		step[:,0] = step[:,0]*(-0.0001)

		# convert the volume to a dilitation 
		step[:,1] = (step[:,1] - Veq)/Veq

		# compute the averages and std of the pressure and dilitation
		P_avg = numpy.mean(step[:,0])
		P_std = numpy.std(step[:,0])

		Dil_avg = numpy.mean(step[:,1])
		Dil_std = numpy.std(step[:,1])

		# add the equilibrium pressure/dilitation to both the 
		# compression and tension data
		if abs(P_avg) < eq_cutoff*0.0001:
			data_C.append((P_avg, P_std, Dil_avg, Dil_std))
			data_T.append((P_avg, P_std, Dil_avg, Dil_std))
			data_eq = (P_avg, Dil_avg)

		# if the pressure is negative, add to the compression data
		elif P_avg < 0: 
			data_C.append((P_avg, P_std, Dil_avg, Dil_std))

		# if the pressure is positive, add to the tension data
		else:
			data_T.append((P_avg, P_std, Dil_avg, Dil_std))

	# convert data lists to numpy arrays
	data_C = numpy.array(data_C)
	data_T = numpy.array(data_T)

	# shift the data to the origin
	data_C = shiftDataToOrigin(data_C, data_eq)
	data_T = shiftDataToOrigin(data_T, data_eq)

	return (data_C, data_T)


def ComputeModulus(data, poisson):
	# assumes data is a numpy array
	# structure of compression/tension data [(P_avg, P_std, Dil_avg, Dil_std)]
	# only consider volumetric strains of 4.5% or less
	data_clean = []
	for step in data:
		if step[2] < 0.045:
			data_clean.append(step)
	data_clean = numpy.array(data_clean)

	# linear fit of dilitation vs pressure data to compute the bulk modulus
	K = numpy.polyfit(data[:,2],data[:,0],1)[0]
	# convert the bulk modulus to the elastic modulus for an isotropic material
	E = 3*K*(1 - 2*poisson)
	return (K, E)


def WriteModulusResults(inputfile, data_C, data_T, poisson):
	# compute the compressive bulk and elastic moduli
	K_C, E_C = ComputeModulus(data_C, poisson)

	# compute the tensile bulk and elastic moduli
	K_T, E_T = ComputeModulus(data_T, poisson)

	# write the results
	# inputfiles have format .lammps
	filename = '%s_modulus_results.txt' %inputfile[:-7]

	f = open(filename, 'w')
	f.write('# modulus results\n')
	f.write('# Poisson Ratio: %.4f\n\n' %poisson)
	f.write('Bulk modulus (compression) = %.4f GPa\n' % K_C)
	f.write('Bulk modulus (tension) = %.4f GPa\n' % K_T)
	f.write('Elastic modulus (compression) = %.4f GPa\n' % E_C)
	f.write('Elastic modulus (tension) = %.4f GPa' % E_T)
	f.close()

def WriteDataMathematicaFormat(inputfile, data_C, data_T):
	''' this function writes the data files in mathematica format 
		for both the ListPlot format and ListErrorPlot format '''

	# inputfiles have format .lammps
	ofileListPLot = '%s_data_ListPlot.txt' %inputfile[:-7]
	ofileListErrorPlot = '%s_data_ListErrorPlot.txt' %inputfile[:-7]

	fPlt = open(ofileListPLot, 'w')
	fErrPlt = open(ofileListErrorPlot, 'w')

	fPlt.write('# avg dilidation vs avg pressure in Mathematica ListPlot format\n')
	fPlt.write('# the dilidation and pressure are averaged over 20,000 timesteps\n\n')

	fErrPlt.write('# avg dilidation vs avg pressure in Mathematica ListErrorPlot format\n')
	fErrPlt.write('# the dilidation and pressure are averaged over 20,000 timesteps\n\n')

	fPlt.write('compression data:\n')
	fErrPlt.write('compression data:\n')
	
	for step in data_C:
		P_avg, P_std, Dil_avg, Dil_std = step

		fPlt.write('{%.10f,%.10f},' %(Dil_avg, P_avg))
		fErrPlt.write('{{%.10f,%.10f}, ErrorBar[%.10f,%.10f]},' %(Dil_avg, P_avg,\
																	Dil_std, P_std))

	fPlt.write('\n\ntension data:\n')
	fErrPlt.write('\n\ntension data:\n')
	
	for step in data_T:
		P_avg, P_std, Dil_avg, Dil_std = step

		fPlt.write('{%.10f,%.10f},' %(Dil_avg, P_avg))
		fErrPlt.write('{{%.10f,%.10f}, ErrorBar[%.10f,%.10f]},' %(Dil_avg, P_avg,\
																	Dil_std, P_std))
	fPlt.close()
	fErrPlt.close()

#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------#
# main program

if len(sys.argv) < 2:
	print 'Usage:'
	print '  python %s <cutoff (N of initial steps to exlude from each pressure step)> [Change Poisson Ratio or Eq Pressure (1 == yes)]' %sys.argv[0]
	exit()

t0 = time.time()

cutoff = int(sys.argv[1])
if len(sys.argv) == 3:
	change = input('\nWhat do you want to change?\n\n'+\
								'\t1) Poisson Ratio (default = 0.25)\n'+\
								'\t2) Equilibrium Pressure Range (default = 75bar)\n\n'+\
								'\tchoice: ')

	if change == 1:
		poisson = input('Give the new Poisson Ratio: ')
		eq_cutoff = 75
	elif change == 2:
		eq_cutoff = input('Give the new equilibrium pressure range: ')
		poisson = 0.25
	else: 
		raise RuntimeError, "Error: '%d' is not a choice!" %change
else:
	poisson = 0.25
	eq_cutoff = 75

# grab all the log files in the current directory
inputfiles = glob.glob("log*.lammps")

# iterate over each log file in the directory
for inputfile in inputfiles:
	print 'Working with %s...' %inputfile
	pressure_steps = LoadData(inputfile)
	Veq = getEquilibVolume(pressure_steps, eq_cutoff)
	data_C, data_T = ComputeAvgPressureAvgDilitation(pressure_steps, Veq, eq_cutoff)

	# Write out results
	WriteModulusResults(inputfile, data_C, data_T, poisson)
	WriteDataMathematicaFormat(inputfile, data_C, data_T)

print 'Computed moduli in %.4f seconds.' %(time.time()-t0)

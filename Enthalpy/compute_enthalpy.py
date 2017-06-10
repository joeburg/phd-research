import sys

#------------------------------------------------------------------------#
def load_log_file(filename, start_step):
	f = open(filename)

	# read lines in log file until the thermo data is reached
	line = 0
	for i in range(start_step-1):
		f.readline()

	# load in data (1000 steps is about 1 K)
	data = []

	N = 0
	temp = 0
	press = 0 
	vol = 0 
	energy = 0
	pressvol = 0
	for i in range(250001):
		fields = f.readline().strip().split()

		# average T, P, V, E every 1000 timesteps (~1 K)
		T = float(fields[2])
		P = float(fields[4])
		V = float(fields[8])
		Etot = float(fields[10])
		PV = P*V
		N += 1

		if N < 1000:
			temp += T 
			press += P
			vol += V
			energy += Etot
			pressvol += PV

		else:
			# recored the ensemble average 
			# convert P from bar to Pa, Vol from A^3 to m^3, E from eV to J
			data.append((temp/N, press/N*100000, vol/N*(10**-30), energy/N*(1.60218*10**-19),pressvol/N*100000*(10**-30)))

			# reset the T, P, V, E, N
			temp = T
			press = P
			vol = V
			energy = Etot
			pressvol = PV 
			N = 1

	f.close()
	return data

def compute_enthalpy(data):
	# set the number of atoms in the system
	NSi = 3000
	NO = 5340
	NC = 3000
	N = NSi+NO+NC

	Na = 6.022e23 # avagadros number 

	# # mSi = 4.6637066e-26 # in kg
	# # mO = 2.6567626e-26 # in kg
	# # mC = 1.9944325e-26 # in kg
	# M = NSi*mSi + NO*mO + NC*mC

	#values at 50 K 
	T0, P0, V0, E0, PV0 = data[0]

	# values at 300 K
	T2, P2, V2, E2, PV2 = data[-1]

	# compute the change in enthalpy between 300K and 50K 
	delta_U = (E2 - E0)*Na/N/1000 # in kJ/mol
	delta_PV = (PV2 - PV0)*Na/N/1000 # in kJ/mol
	delta_H = delta_U + delta_PV # in kJ/mol

	return (delta_U, delta_PV, delta_H)


def write_enthalpy(delta_H, load):
	f = open('enthalpy.txt', 'a')
	f.write('{%s, %.8f}, ' %(load,delta_H))

def write_internal_energy(delta_U, load):
	f = open('internal_energy.txt', 'a')
	f.write('{%s, %.8f}, ' %(load,delta_U))

def write_work_done(delta_PV, load):
	f = open('work_done.txt', 'a')
	f.write('{%s, %.8f}, ' %(load,delta_PV))

#------------------------------------------------------------------------#

if len(sys.argv) < 4:
	print 'Usage:'
	print '  python %s <filename> <start step> <applied load>' %sys.argv[0]
	exit()

filename = sys.argv[1]
start_step = int(sys.argv[2])
load = sys.argv[3]


# compute the enthalpy 
data = load_log_file(filename, start_step)
delta_U, delta_PV, delta_H = compute_enthalpy(data)

# write out data
write_enthalpy(delta_H, load)
write_internal_energy(delta_U, load)
write_work_done(delta_PV, load)




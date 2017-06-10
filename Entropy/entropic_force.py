import sys

#------------------------------------------------------------------------#
def load_log_file(filename, start_step):
	f = open(filename)

	# read lines in log file until the thermo data is reached
	line = 0
	for i in range(start_step-1):
		f.readline()

	# load in data (1000 steps is about 1 K)
	data_P = []

	N = 0
	temp = 0
	press = 0 
	vol = 0 
	energy = 0
	for i in range(250001):
		fields = f.readline().strip().split()
		if i == 0:
			Teq = float(fields[2])
			Peq = float(fields[4])*100000
			Veq = float(fields[8])*(10**-30)
			Eeq = float(fields[10])*(1.6*10**-19)

			data_eq = (Teq, Peq, Veq, Eeq)

		else: 
			Tave = float(fields[2])
			Pave = float(fields[4])
			V = float(fields[8])
			Etot = float(fields[10])
			N += 1

			if N < 1000:
				temp += Tave 
				press += Pave
				vol += V
				energy += Etot

			else:
				# recored the ensemble average 
				# convert P from bar to Pa, Vol from A^3 to m^3, E from eV to J
				data_P.append((temp/N, press/N*100000, vol/N*(10**-30), energy/N*(1.6*10**-19)))

				# reset the T, P, V, E, N
				temp = Tave
				press = Pave
				vol = V
				energy = Etot 
				N = 1

	f.close()
	return data_eq, data_P

def compute_force(data_eq, data_P):
	# compute the force between thermodynamics at 300K and 50K 
	# Teq = data_eq[0]
	# Peq = data_eq[1]
	# Veq = data_eq[2]
	# Eeq = data_eq[3]

	Teq = data_P[0][0]
	Peq = data_P[0][1]
	Veq = data_P[0][2]
	Eeq = data_P[0][3]

	T2 = data_P[-1][0]
	P2 = data_P[-1][1]
	V2 = data_P[-1][2]
	E2 = data_P[-1][3]

	delta_S = (Eeq + Peq*Veq)/Teq - (E2 + P2*V2)/T2 # in J/K
	delta_U = Eeq - E2
	delta_V = Veq - V2

	# only consider entropic term in force
	N = 11340*2 # number of bonds 
	Na = 6.02*10**23

	# f = -Teq*delta_S/delta_V
	# f = -f/Na*(10**12) #in pN/mol and define resistive force as positive 
	# fin = delta_U/delta_V/Na*(10**12) # in pN/mol 

	delta_L = abs(delta_V)**(1.0/3)

	f = abs(Teq*delta_S/delta_L)
	f = f/N*1e9 # in nN/bond

	fin = abs(delta_U/delta_L/N*1e9) # in nN/bond
	return f, fin, delta_S

def write_force(force, load):
	f = open('entropic_forces.txt', 'a')
	f.write('{%s, %.20f}, ' %(load,force))
	f.close()

def write_force_internal(force, load):
	f = open('internal_forces.txt', 'a')
	f.write('{%s, %.6f}, ' %(load,force))
	f.close()

def write_entropy(delta_S, load):
	f = open('entropy.txt', 'a')
	f.write('{%s, %.2f}, ' %(load,delta_S*(10**17)))

#------------------------------------------------------------------------#

if len(sys.argv) < 4:
	print 'Usage:'
	print '  python %s <filename> <start step> <applied load>' %sys.argv[0]
	exit()

filename = sys.argv[1]
start_step = int(sys.argv[2])
load = sys.argv[3]

data_eq, data_P = load_log_file(filename, start_step)
force, force_in, delta_S = compute_force(data_eq, data_P)
write_force(force,load)
write_entropy(delta_S,load)
write_force_internal(force_in,load)




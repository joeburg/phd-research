# this file allows to easily add new models without 
# directly changing the source code
# to add a new precursor, update
	# 1) precursor type
	# 2) condensation degree
	# 3) network connectivity / relative network connectivity 
	# 4) mean Si network connectivity
	# 5) density

#-----------------------------------------------------------------------------------#
# 1) precursor type
ValidMolTypes = ['EtOCS',\
				'OCSEt',\
				'EtOCSMe',\
				'135Benz',\
				'135Benzene',\
				'13Benz',\
				'13Benzene',\
				'14Benz',\
				'14Benzene',\
				'SiO2',\
				'EtOCSMethyl',\
				'EtOCSVinyl',\
				'EtOCSPhenyl',\
				'TVS',\
				'TVSR1',\
				'TVSR2',\
				'TVSR3',\
				'TVSR4',\
				'TVSR6',\
				'EtOCS_EtOCSPhenyl']

#-----------------------------------------------------------------------------------#
# 2) condensation degree

def q_SiO3(NT0, NT1, NT2, NT3, NT4):
	return (1.0/3)*(0*NT0 + 1*NT1 + 2*NT2 + 3*NT3)

def q_SiO4(NT0, NT1, NT2, NT3, NT4):
	return (1.0/4)*(0*NT0 + 1*NT1 + 2*NT2 + 3*NT3 + 4*NT4)

CondensationDegree = {'EtOCS': q_SiO3,
						'OCSEt': q_SiO3,
						'EtOCSMe': q_SiO3,
						'135Benz': q_SiO3,
						'135Benzene': q_SiO3,
						'13Benz': q_SiO3,
						'13Benzene': q_SiO3,
						'14Benz': q_SiO3,
						'14Benzene': q_SiO3,
						'SiO2': q_SiO4,
						'EtOCSMethyl': q_SiO3,
						'EtOCSVinyl': q_SiO3,
						'EtOCSPhenyl': q_SiO3,
						'TVS': q_SiO3,
						'TVSR1': q_SiO3,
						'TVSR2': q_SiO3,
						'TVSR3': q_SiO3,
						'TVSR4': q_SiO3,
						'TVSR6': q_SiO3,
						'EtOCS_EtOCSPhenyl': q_SiO3}

#-----------------------------------------------------------------------------------#
# 2) network connectivity

def bridged_silane(q):
	return (1.0/4)*(1.0 + 3.0*q)	

def bridged_silane_w_terminal_group(q):
	return (1.0/4)*(1.0 + 2.5*q)

def hyperconnected_silane(q):
	return (1.0/5)*(2.0 + 3.0*q)

def silicon(q):
	return 1.0*q

def mixed(q, Nprec1, Nprec2, moltype1, moltype2):
	a = float(Nprec1)/(Nprec1+Nprec2)
	b = float(Nprec2)/(Nprec1+Nprec2)
	return a*Connectivity[moltype1](q) + b*Connectivity[moltype2](q)

Connectivity = {'EtOCS': bridged_silane,
				'OCSEt': bridged_silane,
				'EtOCSMe': bridged_silane_w_terminal_group,
				'135Benz': hyperconnected_silane,
				'135Benzene': hyperconnected_silane,
				'13Benz': bridged_silane,
				'13Benzene': bridged_silane,
				'14Benz': bridged_silane,
				'14Benzene': bridged_silane,
				'SiO2': q_SiO4,
				'EtOCSMethyl': bridged_silane_w_terminal_group,
				'EtOCSVinyl': bridged_silane_w_terminal_group,
				'EtOCSPhenyl': bridged_silane_w_terminal_group,
				'TVS': hyperconnected_silane,
				'TVSR1': hyperconnected_silane,
				'TVSR2': hyperconnected_silane,
				'TVSR3': hyperconnected_silane,
				'TVSR4': hyperconnected_silane,
				'TVSR6': hyperconnected_silane,
				'EtOCS_EtOCSPhenyl': mixed}

#-----------------------------------------------------------------------------------#
# 4) mean Si network connectivity

def mSi_3_5(p):
	return 3.5*p

def mSi_4(p):
	return 4.0*p

def mSi_5(p):
	return 5.0*p

def mSi_mixed(p, Nprec1, Nprec2, moltype1, moltype2):
	a = float(Nprec1)/(Nprec1+Nprec2)
	b = float(Nprec2)/(Nprec1+Nprec2)
	return a*MeanSiConnectivity[moltype1](p) + b*MeanSiConnectivity[moltype2](p)

MeanSiConnectivity = {'EtOCS': mSi_4,
					'OCSEt': mSi_4,
					'EtOCSMe': mSi_3_5,
					'135Benz': mSi_5,
					'135Benzene': mSi_5,
					'13Benz': mSi_4,
					'13Benzene': mSi_4,
					'14Benz':mSi_4,
					'14Benzene':mSi_4,
					'SiO2': q_SiO4,
					'EtOCSMethyl': mSi_3_5,
					'EtOCSVinyl': mSi_3_5,
					'EtOCSPhenyl': mSi_3_5,
					'TVS': mSi_5,
					'TVSR1': mSi_5,
					'TVSR2': mSi_5,
					'TVSR3': mSi_5,
					'TVSR4': mSi_5,
					'TVSR6': mSi_5,
					'EtOCS_EtOCSPhenyl': mSi_mixed}

#-----------------------------------------------------------------------------------#
# 5) density parameters 

def silane_3Si(NSi):
	# the silane precursors are related to the number of Si atoms
	return (1.0/3.0)*NSi

def silane_2Si(NSi):
	return 0.5*NSi

def silane_1Si(NSi):
	return 1.0*NSi

def silane_mixed(NSi, Nprec1, Nprec2, moltype1, moltype2):
	a = float(Nprec1)/(Nprec1+Nprec2)
	b = float(Nprec2)/(Nprec1+Nprec2)
	return a*Nprecursors[moltype1](NSi) + b*Nprecursors[moltype2](NSi)

def massPrecursor_mixed(Nprec1, Nprec2, moltype1, moltype2):
	a = float(Nprec1)/(Nprec1+Nprec2)
	b = float(Nprec2)/(Nprec1+Nprec2)
	return a*MassPrecursor[moltype1] + b*MassPrecursor[moltype2]	


Mass = {'Si' : 28.09,
		'C' : 12.01,
		'CH' : 13.02,
		'CH2' : 14.03,
		'CH3' : 15.04,
		'O' : 15.999,
		'OH' : 17.009}

MassPrecursor = {'EtOCS' : 2*Mass['Si'] + 2*Mass['CH2'],
				'OCSEt' : 2*Mass['Si'] + 2*Mass['CH2'],
				'EtOCSMethyl' : 2*Mass['Si'] + 2*Mass['CH2'] + Mass['CH3'],
				'EtOCSMe' : 2*Mass['Si'] + 2*Mass['CH2'] + Mass['CH3'],
				'EtOCSVinyl' : 2*Mass['Si'] + 3*Mass['CH2'] + Mass['CH'],
				'EtOCSPhenyl' : 2*Mass['Si'] + 2*Mass['CH2'] + Mass['C'] + 5*Mass['CH'],
				'135Benz' : 3*Mass['Si'] + 3*Mass['C'] + 3*Mass['CH'],
				'135Benzene' : 3*Mass['Si'] + 3*Mass['C'] + 3*Mass['CH'],
				'13Benz' : 2*Mass['Si'] + 2*Mass['C'] + 4*Mass['CH'],
				'13Benzene' : 2*Mass['Si'] + 2*Mass['C'] + 4*Mass['CH'],
				'14Benz' : 2*Mass['Si'] + 2*Mass['C'] + 4*Mass['CH'],
				'14Benzene' : 2*Mass['Si'] + 2*Mass['C'] + 4*Mass['CH'],
				'SiO2' : Mass['Si'],
				'TVS' : 3*Mass['Si'] + Mass['C'] + 6*Mass['CH2'] + Mass['CH3'],
				'TVSR1' : 3*Mass['Si'] + Mass['C'] + 3*Mass['CH2'] + Mass['CH3'],
				'TVSR2' : 3*Mass['Si'] + Mass['C'] + 6*Mass['CH2'] + Mass['CH3'],
				'TVSR3' : 3*Mass['Si'] + Mass['C'] + 9*Mass['CH2'] + Mass['CH3'],
				'TVSR4' : 3*Mass['Si'] + Mass['C'] + 12*Mass['CH2'] + Mass['CH3'],
				'TVSR6' : 3*Mass['Si'] + Mass['C'] + 18*Mass['CH2'] + Mass['CH3'],
				'EtOCS_EtOCSPhenyl': massPrecursor_mixed}

Nprecursors = {'EtOCS' : silane_2Si,
				'OCSEt' : silane_2Si,
				'EtOCSMethyl' : silane_2Si,
				'EtOCSMe' : silane_2Si,
				'EtOCSVinyl' : silane_2Si,
				'EtOCSPhenyl' : silane_2Si,
				'135Benz' : silane_3Si,
				'135Benzene' : silane_3Si,
				'13Benz' : silane_2Si,
				'13Benzene' : silane_2Si,
				'14Benz' : silane_2Si,
				'14Benzene' : silane_2Si,
				'SiO2' : silane_1Si,
				'TVS' : silane_3Si,
				'TVSR1' : silane_3Si,
				'TVSR2' : silane_3Si,
				'TVSR3' : silane_3Si,
				'TVSR4' : silane_3Si,
				'TVSR6' : silane_3Si,
				'EtOCS_EtOCSPhenyl': silane_mixed}



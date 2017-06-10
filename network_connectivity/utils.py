
def namePrecursors():
	molType1 = raw_input('\nGive the name of precursor 1: ')
	molType2 = raw_input('Give the name of precursor 2: ')
	return (molType1, molType2)

def getNPrecursors(molType1, molType2):
	Nprecursors1 = input('\nGive the number of %s precursors: ' %molType1)
	Nprecursors2 = input('Give the number of %s precursors: ' %molType2)
	return (Nprecursors1, Nprecursors2)

def getPrecursorParams(molType, Nprecursors1=None, Nprecursors2=None):
	molType1 = None
	molType2 = None

	# if the number of precursors are given, assume the molType is in
	# stanford form and extract molType1 and molType2
	if Nprecursors1:
		idx = molType.find('_')
		molType1 = molType[:idx]
		molType2 = molType[idx+1:]

	else:
		# try to find the moltypes 
		idx = molType.find('_')
		if idx >= 0:
			molType1 = molType[:idx]
			molType2 = molType[idx+1:]

			mixed = input('Is this a mixed precuror model with "%s" and "%s"?\n'
							%(molType1,molType2)+\
							'\t1) Yes\n'+\
							'\t2) Yes, but the precursor names are incorrect.\n'+\
							'\t2) No\n\n'+\
							'\tchoice: ')

			dispatcher = {1: (molType1, molType2),
							2: namePrecursors,
							3: False}

			try:
				response = dispatcher[mixed]()
			except:
				response = dispatcher[mixed]

			if response:
				molType1, molType2 = response
				Nprecursors1, Nprecursors2 = getNPrecursors(molType1, molType2)

	# setup the parameter dictionary
	PrecursorParams = {'molType': molType,
						'molType1': molType1,
						'molType2': molType2,
						'Nprecursors1': Nprecursors1,
						'Nprecursors2': Nprecursors2}
	return PrecursorParams

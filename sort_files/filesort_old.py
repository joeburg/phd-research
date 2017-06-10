trials = ['2.1','2.2','2.3','3.1','3.2','3.3','4.1','4.2','4.3','5.1','5.2','5.3','6.1','6.2','6.3']
#------------------------------------------------------------------------------#

file_base_xyz = '135Benz_{}_{}.xyz'

for i in range(170000,395000,5000):
	for trial in trials:

		# make the trial directories and sub directories
		top_dir = 'trial{}/'.format(trial[ : trial.find('.')])
		sub_dir = '{}trial{}/'.format(top_dir, trial)
		data_sub_dir = '{}data/'.format(sub_dir)

		# make the top directory if it dosent exit
		if not os.path.exists(top_dir):
			os.makedirs(top_dir)

		# make the sub-directory 
		if not os.path.exists(sub_dir):
			os.makedirs(sub_dir)

		# make the data sub-directory
		if not os.path.exists(data_sub_dir):
			os.makedirs(data_sub_dir)

		# rename the file
		file_name = file_base_xyz.format(i,trial)
		file_name_new = '{}{}'.format(data_sub_dir,file_name)
		os.rename(file_name, file_name_new)

#------------------------------------------------------------------------------#

file_base_rdf = '135Benz_rdf_{}.txt'
file_base_por = 'porogen_170000_{}.xyz'
file_base_benz_por = '135BenzPorogen_170000_{}.xyz'
file_base_log = 'log_{}.lammps'

for trial in trials:
	# make the trial directories and sub directories
	top_dir = 'trial{}/'.format(trial[ : trial.find('.')])
	sub_dir = '{}trial{}/'.format(top_dir, trial)
	data_sub_dir = '{}data/'.format(sub_dir)

	# make the top directory if it dosent exit
	if not os.path.exists(top_dir):
		os.makedirs(top_dir)

	# make the sub-directory 
	if not os.path.exists(sub_dir):
		os.makedirs(sub_dir)

	# make the data sub-directory
	if not os.path.exists(data_sub_dir):
		os.makedirs(data_sub_dir)

	# rename the file
	file_name_rdf = file_base_rdf.format(trial)
	file_name_rdf_new = '{}{}'.format(sub_dir,file_name_rdf)
	os.rename(file_name_rdf, file_name_rdf_new)

	# rename the file
	file_name_por = file_base_por.format(trial)
	file_name_por_new = '{}{}'.format(data_sub_dir,file_name_por)
	os.rename(file_name_por, file_name_por_new)

	# rename the file
	file_name_benz_por = file_base_benz_por.format(trial)
	file_name_benz_por_new = '{}{}'.format(data_sub_dir,file_name_benz_por)
	os.rename(file_name_benz_por, file_name_benz_por_new)

	# rename the file
	file_name_log = file_base_log.format(trial)
	file_name_log_new = '{}{}'.format(sub_dir,file_name_log)
	os.rename(file_name_log, file_name_log_new)	
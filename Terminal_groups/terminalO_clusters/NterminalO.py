''' sums a histogram to get the number of terminal O '''

import glob
import time

#--------------------------------------------------------------------------#

def LoadDataHist(inputfile):
		""" reads histogram file """ 
		bins = []
		counts = []

		f = open(inputfile)

		while True:
			fields = f.readline().strip().split(',')
			if fields[0]:
				bin = float(fields[0])
				count = float(fields[1])

				bins.append(bin)
				counts.append(count)
			else:
				break
		f.close()

		return bins, counts

#--------------------------------------------------------------------------#
t0 = time.time()

filenames = glob.glob('data/clusters_hist_cutoff_3.6_OCSEtpore_140000_*.csv')


fout = 'NterminalO_cutoff_3.6.txt'
f = open(fout, 'w')

for fname in filenames:
	NtermO = 0
	bins, counts = LoadDataHist(fname)

	# add up terminal O
	for i in range(len(counts)):
		NtermO += bins[i]*counts[i]

	# write out results
	f.write('%s\n' % fname)
	f.write('N terminal O = %d\n\n' % NtermO)

f.close()

print 'Found the number of terminal O in %.4f seconds.' %(time.time()-t0) 



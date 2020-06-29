import numpy as np
from obspy import read, Stream
from glob import glob
import os, sys

####
# tally up all the unique station pairs
# for each pair:
#	find all the monthly cross-correlation stacks
#	read and stack
#	write overall stack as sac file for aftan, with ndays properly set as user0 etc.
####

corfiles = np.array(glob('20??.???/COR/*/*SAC'))  # all monthly crosscorr files
pairs = np.array(['_'.join(e.split('/')[-1].split('.')[0].split('_')[1:]) for e in corfiles])
u_pairs = np.unique(pairs)

odir = 'COR/'

for p in u_pairs:
	print(p)
	if p.split('_')[0] == p.split('_')[1]:
		continue  # skip autocorrelations

	flist = corfiles[np.where(pairs == p)]
	for i in range(len(flist)):
		if i == 0:
			st = read(flist[i]) # start the stream with one trace
		else:
			st1 = read(flist[i])
			st[0].data = st[0].data + st1[0].data
			st[0].stats.sac.user0 = st[0].stats.sac.user0 + st1[0].stats.sac.user0

	st[0].data = st[0].data/len(flist)

	ofile = os.path.join(odir,'COR_%s.SAC' % p)
	st.write(ofile,'SAC')

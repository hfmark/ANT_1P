import numpy as np
from obspy import read, Stream
from glob import glob
from datetime import datetime
import os, sys

####
# tally up all the unique station pairs
# for each pair:
#	find all the monthly cross-correlation stacks
#	read and stack
#	write overall stack as sac file for aftan, with ndays properly set as user0 etc.
# ALSO:
#	group monthly cross-correlations into threes, stack
#	write with some kind of useful nomenclature
#	np.array(list(zip(a,a[1:],a[2:])))
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
	for i in range(len(flist)):  # FULL stack
		if i == 0:
			st = read(flist[i]) # start the stream with one trace
			st[0].data = st[0].data*st[0].stats.sac.user0  # weight by ndays
		else:
			st1 = read(flist[i])
			st[0].data = st[0].data + st1[0].data*st1[0].stats.sac.user0
			st[0].stats.sac.user0 = st[0].stats.sac.user0 + st1[0].stats.sac.user0

	st[0].data = st[0].data/st[0].stats.sac.user0

	ofile = os.path.join(odir,'COR_%s.SAC' % p)
	st.write(ofile,'SAC')

	# now sort in order of time, get overlapping 3-month segments, and stack those
	if len(flist) < 4:  # if we can't get two 3-month stacks, skip this
		continue

	# month/year sorting
	dates = np.array([datetime.strptime(e.split('/')[0], '%Y.%b') for e in flist])
	inds = np.argsort(dates)
	flist = flist[inds]

	# ordered overlapping triples
	tri_month = np.array(list(zip(flist,flist[1:],flist[2:])))

	# loop triples, sum, write
	for i in range(len(tri_month)):
		mons = np.array([datetime.strptime(e.split('/')[0].split('.')[1],'%b').strftime('%m') \
			for e in tri_month[i]])
		ofile = os.path.join(odir,'COR_%s_%s.%s.%s.SAC' % (p,mons[0],mons[1],mons[2]))
		for j in range(3):
			if j == 0:
				st = read(tri_month[i][j])
				st[0].data = st[0].data*st[0].stats.sac.user0
			else:
				st1 = read(tri_month[i][j])
				st[0].data = st[0].data + st1[0].data*st1[0].stats.sac.user0
				st[0].stats.sac.user0 = st[0].stats.sac.user0 + st1[0].stats.sac.user0

		st[0].data = st[0].data/st[0].stats.sac.user0

		st.write(ofile,'SAC')

import numpy as np
from glob import glob
from datetime import datetime
import os, sys

####
# list all dates that we have any data for, then:
#	find all miniseed files for a given day (all networks)
#	cat mseed files
#	call rdseed to make a seed dayfile with combined dataless
#	rename dayfile by date
####

# get a list of ALL dates with data (check 1P, YJ, EN)
P_dirs = np.array(glob('seed/1P/*/'))
P_days = np.array([datetime.strptime('.'.join(a.split('_')[-1].split('.')[:3]),'%Y.%b.%d') \
				for a in P_dirs])

Y_dirs = np.array(glob('seed/YJ/*/'))
Y_days = np.array([datetime.strptime('.'.join(a.split('_')[-1].split('.')[:3]),'%Y.%b.%d')\
				for a in Y_dirs])

E_file = np.array(glob('seed/ENAP/*.mseed'))
E_days = np.array([datetime.strptime('.'.join(a.split('-')[-1].split('.')[:3]),'%Y.%m.%d') \
				for a in E_file])

days = np.unique(np.hstack((P_days,Y_days,E_days)))
dataless_file = 'seed/dataless/combined.dataless'

trouble = []
print(len(days),'days')
for d in days:
	if d.year < 2018 or d.year >= 2019:  # for now, just do 2018 (testing for MG* and AY*)
		continue
	mseed_list = []
	if d in P_days:
		fdir = P_dirs[P_days == d][0]
		flist = glob(fdir+'*.mseed')
		for ffilei in flist:
			mseed_list.append(ffile)
	if d in Y_days:
		fdir = Y_dirs[Y_days == d][0]
		flist = glob(fdir+'*.mseed')
		for i in flist:
			mseed_list.append(ffile)
	if d in E_days:
		ffile = E_file[E_days==d][0]
		mseed_list.append(ffile)

	cat_str = ' '.join(mseed_list)  # string of miniseed files to be catted
	mseed_ofile = 'seed/dayfiles/PAT_' + d.strftime('%Y.%m.%d') + '.mseed'
	seed_ofile = 'seed/dayfiles/PAT_' + d.strftime('%Y.%m.%d') + '.seed'

	iq = os.system('cat %s > %s' % (cat_str, mseed_ofile)
	if iq != 0:
		trouble.append(d)
		continue
	iq = os.system('rdseed -d -o 5 -f %s -g %s > junk' % (mseed_ofile, dataless_file))
	if iq != 0:
		trouble.append(d)
		continue
	os.rename('seed.rdseed',seed_ofile)

if len(trouble) > 0:
	f = open('trouble.dat','w')
	for d in trouble:
		f.write(d.strftime('%Y.%m.%d'))
		f.write('\n')
	f.close()

import numpy as np
from obspy.core.stream import Stream, read
from glob import glob
from datetime import datetime
import os, sys

####
# get daily ENAP files from /DB, downsample to 1Hz, save in day dirs as miniseed
####

mseed_dir_in = '/DB/PATAGONIA/ENAP/data_files/'
mseed_dir_out = '/P/hmark/ANT_1P/seed/ENAP/'

# get list of all relevant files
dayfiles = np.array(glob(mseed_dir_in + '20??/*/*.HH?.*'))

# (messy) qc: clean ANMA files for bad timing days
bad_files_1 = np.array([mseed_dir_in + '2019/ANMA/ANMA.EN..HH%s.2019.%03d' % (j,k) for j in ['E','N','Z'] for k in np.append(np.arange(314,347),365)])
bad_files_2 = np.array([mseed_dir_in + '2020/ANMA/ANMA.EN..HH%s.2020.%03d' % (j,k) for j in ['E','N','Z'] for k in np.arange(0,78)])
bad_files = np.append(bad_files_1,bad_files_2)

# get list of all days and of unique days
days_all = np.array([datetime.strptime(a[-8:],'%Y.%j').strftime('%Y.%m.%d') for a in dayfiles])
days = np.unique(days_all)

for dy in days:
	file_list = dayfiles[np.where(days_all == dy)]  # files for this day, all stations

	sys.exit()

	st = Stream()
	for fl in file_list:
		if fl not in bad_files and os.stat(fl).st_size != 0:
			st += read(fl)
	if len(st) > 0:
		# decimate or resample to 1Hz
		if np.all([tr.stats.sampling_rate == 100 for tr in st]):
			st.decimate(100, no_filter=True)
		else:
			st.resample(1.0, no_filter=True)

		# write miniseed
		ofile = os.path.join(mseed_dir_out,'EN-%s.%s.%s.mseed' % \
				(dy.split('.')[0],dy.split('.')[1],dy.split('.')[2]))
		st.write(ofile, format='MSEED')

from obspy import read
from glob import glob
import os, sys


####
# read in cross-correlation from sac file (seed2cor), symmetrize, write symmetric sac file
####

#cor_list = glob('COR/*/*LHZ*LHZ*SAC')  # skip other components and derivative files
cor_list = glob('COR/*SAC')

for ifile in cor_list:
	xc = read(ifile,'SAC')

	n = xc[0].count()
	mid = int((n-1)/2)
	if n%2 != 1:  # if length is wrong, skip
		print('cross corr %s cannot be symmetrized' % (ifile.split('/')[-1].split('.')[0]))
		continue

	for tr in xc:
		# symmetrize trace data
		h1 = tr.data[mid:]
		h2 = tr.data[mid::-1]
		tr.data = (h1 + h2)/2.

		tr.stats.starttime = tr.stats.endtime  # hopefully this will auto-adjust?

	# write symmetrized sac file
	ofile = ifile+'_s'
	xc.write(ofile,'SAC')

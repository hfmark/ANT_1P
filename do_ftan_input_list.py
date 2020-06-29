import numpy as np
from glob import glob
from scipy.interpolate import interp1d
import os, sys

####
# write list input file for aftan
#
# each line, for processing one symmetrized cross correlation, looks like:
#	-1 1.5 5 7 100 20 1 0.5 0.2 2 COR_109C_BAK.SAC
# where the parameters are, left to right:
#	piover4, vmin, vmax, tmin, tmax, thresh, ffact, taperl, snr, fmatch, filename
#
# piover4: sign for pi/4 term in phase velocity calculation
# vmin, vmax: min and max velocity for ftan
# tmin, tmax: min and max period for ftan
# thresh: threshold for catching dispersion curve jumps, usually 10?
# ffact: multiplicative factor for ftan filter width [alpha = ffact*20.*sqrt(dist/1000.)]
#	NOTE that this is hard-coded to 1 in the driver
# taperl: factor for left-end seismogram tapering: taper = taperl*tmax
# snr: phase-matching filter parameter, spectra ratio for cutting point
# fmatch: factor for length of phase matching window
# filename: name of sac file with cross correlation stack
####

cor_list = np.sort(glob('COR/*SAC'))  # list files that we want to process
sta1_lst = np.array([e.split('/')[-1].split('_')[1] for e in cor_list])
sta2_lst = np.array([e.split('/')[-1].split('_')[2].split('.')[0] for e in cor_list])
twosta = sta1_lst != sta2_lst  # don't use autocorrelations
cor_list = cor_list[twosta]

phfile = 'ak135_phvel.dat'

minp = 3; maxp = 45
pp = np.arange(minp,maxp+1)
minv = 1.7; maxv = 5.2;
minSNR = 5;

thresh = 10 # 20; 
ffact = 1; taperl = 0.5; fsnr = 0.2; fmatch = 2;

f = open('aftan.lst','w')  # output file for input to aftan

for cf in cor_list[62:92]:
	# find corresponding snr file, get period range for dispersion curve calc
	sf = cf+'_s_snr.cv.p.txt'
	per,snr = np.loadtxt(sf,usecols=(0,1),unpack=True)
	resamp = interp1d(per,snr,fill_value='extrapolate')
	snr_rs = resamp(pp)
	usable = pp[np.where(snr_rs > minSNR)]

	f.write('-1 %.1f %.1f %.1f %.1f %.1f %.1f %.1f %.1f %.1f %s %s\n' % \
		(minv, maxv, min(usable), max(usable), thresh, ffact, taperl, fsnr, fmatch, cf, phfile))

f.close()

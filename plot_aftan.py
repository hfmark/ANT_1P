import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.backends.backend_pdf import PdfPages
import shapefile
from scipy.interpolate import interp1d, griddata
from copy import copy
from glob import glob
from obspy import read, read_inventory
import os, sys

####
# read dispersion curves and amplitude maps from aftan output, plot
# (maybe with station info/map/xcors as in pysismo?)
####

def read_AMP(sta1,sta2,iper,nf=2,root='COR/',return_dist=False):
	"""
	read amplitudes from aftan output file
	return array where each row is a velocity and each column is a period
	"""
	# find amplitude file
	fname = os.path.join(root,'COR_%s_%s_full.SAC_%i_AMP' % (sta1,sta2,nf))
	if not os.path.isfile(fname):
		sta1, sta2 = sta2, sta1  # swap station order and try again
		fname = os.path.join(root,'COR_%s_%s_full.SAC_%i_AMP' % (sta1,sta2,nf))
		if not os.path.isfile(fname):
			print('no AMP file for %s, %s' % (sta1,sta2))
			return

	# read grid parameters
	f=open(fname,'r')
	line = f.readline()
	f.close()
	nper,nvel,dt,dist = [float(e) for e in line.split()]

	# get amplitudes from file, resample on a regular velocity grid
	per,time,ampT = np.loadtxt(fname,usecols=(0,1,2),unpack=True,skiprows=1)
	vel = dist/np.unique(time)  # times are even, vels are uneven!
	# periods are also uneven, fix that too
	iper_long = iper[per.astype(int)-1]
	grid_x,grid_y = np.mgrid[min(iper):max(iper):np.complex(nper), min(vel):max(vel):np.complex(nvel)]
	amp = griddata(np.vstack((iper_long,dist/time)).T,ampT,(grid_x,grid_y),method='linear')
	amp = np.flipud(amp.T)

	# normalize each column (better for plotting)
	for i in range(amp.shape[1]):
		amp[:,i] = amp[:,i]/max(abs(amp[:,i]))
	
	minT = min(np.unique(time)); maxT = max(np.unique(time))
	vmin = dist/maxT; vmax = dist/minT
	
	if return_dist:
		return amp, vmin, vmax, dist
	return amp, vmin, vmax

def read_DISP(sta1,sta2,nf=2,itr=1,root='COR/'):
	"""
	read dispersion curves from aftan output file
	return nominal and instantaneous periods, group and phase velocities
	"""
	# find disp file:
	# nf=1 for pre-phase matched filtering, nf=2 for filtered
	# itr=0 for initial curves without jump correction, itr=1 for final curves
	fname = os.path.join(root,'COR_%s_%s_full.SAC_%i_DISP.%i' % (sta1,sta2,nf,itr))
	if not os.path.isfile(fname):
		sta1, sta2 = sta2, sta1
		fname = os.path.join(root,'COR_%s_%s_full.SAC_%i_DISP.%i' % (sta1,sta2,nf,itr))
		if not os.path.isfile(fname):
			print('no DISP file for %s, %s' % (sta1, sta2))
			return

	# read outputs
	cper, iper, gvel, pvel = np.loadtxt(fname,usecols=(1,2,3,4),unpack=True)

	return cper, iper, gvel, pvel

def calc_cutoffperiod(sta1,sta2,dist,wvln=2.5,nf=2,itr=1,root='COR/',gvel=None,iper=None):
	"""
	find the max period we want to use for a given station pair based on 
	station spacing and a minimum number of wavelengths and the estimated group velocities
	"""
	# TODO: iper or cper???? probably iper
	if gvel is None or iper is None:  # get group vel curve if not given
		_,iper,gvel,_ = read_DISP(sta1,sta2,nf=nf,itr=itr,root=root)

	try:
		cutoffperiod = iper[dist/(wvln*gvel) >= iper].max()
	except ValueError: # no max, basically
		cutoffperiod = min(iper) - 1  # no good periods here

	return cutoffperiod

def read_snr_precalc(sta1,sta2,root='COR/'):
	"""
	read snr calculated prior to ftan (referenced to center periods)
	"""
	# get snr filename
	fname = os.path.join(root,'COR_%s_%s_full.SAC_s_snr.cv.p.txt' % (sta1,sta2))
	if not os.path.exists(fname):
		sta1, sta2 = sta2, sta1
		fname = os.path.join(root,'COR_%s_%s_full.SAC_s_snr.cv.p.txt' % (sta1,sta2))
		if not os.path.exists(fname):
			print('no snr file found for %s, %s' % (sta1, sta2))
			return

	# load SNR (use post, not precursor)
	cper,snr = np.loadtxt(fname,usecols=(0,1),unpack=True)
	return cper, snr

def find_sta_network(sta,inv):
	"""
	find the network code for a given station name, since that info doesn't make it
	into the seed2cor workflow for both stations of a pair
	"""
	return inv.select(station=sta).networks[0].code

def plot_basemap(ax=None,coast_shp='shapefiles/SouthAmericaCoasts.shp',\
		   bbox=[-76,-65,-56,-42]):
	"""
	plot basemap with coastlines
	"""
	assert ax, 'need to provide an axis for plotting'

	sf = shapefile.Reader(coast_shp)
	for shape in sf.shapes():
		parts = list(shape.parts) + [len(shape.points)]
		partlims = zip(parts[:-1],parts[1:])
		for i1, i2 in partlims:
			points = shape.points[i1:i2]
			x, y = zip(*points)
			ax.plot(x,y,'-',lw=0.75,color='k')

	ax.set_aspect('equal')
	ax.set_xlim(bbox[:2])
	ax.set_ylim(bbox[2:])
	ax.grid()

	return

def get_sta_coords(sta,listfile='stations.lst'):
	"""
	retrieve station coords from list used as seed2cor input
	"""
	assert os.path.isfile(listfile),'%s not found' % listfile

	stnm = np.loadtxt(listfile,usecols=(0,),dtype=str)
	lat,lon = np.loadtxt(listfile,usecols=(1,2),unpack=True)

	return lat[stnm==sta][0],lon[stnm==sta][0]

def get_ndays(sta1,sta2,root='COR/'):
	"""
	get number of stacked days in xcor for a pair of stations from sac header
	"""
	fname = os.path.join(root,'COR_%s_%s_full.SAC_s' % (sta1,sta2))
	if not os.path.isfile(fname):
		sta1, sta2 = sta2, sta1
		fname = os.path.join(root,'COR_%s_%s_full.SAC_s' % (sta1,sta2))
		if not os.path.isfile(fname):
			print('no symmetric xcor file found for %s, %s' % (sta1, sta2))
			return

	xc = read(fname,'SAC')

	return xc[0].stats.sac.user0
	
if __name__ == '__main__':

	done_list = np.sort(glob('COR/*full*2_AMP'))  # only plot sets that *have* a final result
	sta1_list = [e.split('/')[-1].split('_')[1] for e in done_list]
	sta2_list = [e.split('/')[-1].split('_')[2].split('.')[0] for e in done_list]

	pdfpath = 'output/FTAN_plots.pdf'
	if os.path.isfile(pdfpath):
		iq = input('%s exists. Overwrite? y/[n]' % pdfpath) or 'n'
		if iq == 'n': sys.exit()

	pdf = PdfPages(pdfpath)

	inv = read_inventory('seed/dataless/combined.xml')

	for i in range(len(sta1_list)):
		sta1 = sta1_list[i]; sta2 = sta2_list[i]
		if sta1 == sta2:
			continue
		# check for all necessary files:
		disp = glob('COR/COR_%s_%s*full*DISP*' % (sta1, sta2))
		amp = glob('COR/COR_%s_%s*full*AMP*' % (sta1, sta2))
		if len(disp) != 4 or len(amp) != 2:  # file is probably missing
			print('skipping %s-%s, some files missing' % (sta1,sta2))
			continue
		print(sta1,sta2)

		# start plotting
		fig = plt.figure(figsize=(12,5))
		xlm = [4,45]

		#### un-phase-match-filtered
		# amplitudes, periods, velocities
		cper,iper,gvel,pvel = read_DISP(sta1,sta2,nf=1,itr=1)
		_,iper_amp,_,_ = read_DISP(sta1,sta2,nf=1,itr=0)
		amp,vmin,vmax,dist = read_AMP(sta1,sta2,iper_amp,nf=1,return_dist=True)
		if min(iper) == max(iper):  # singularity!! oops
			plt.close()
			continue
		cfp = calc_cutoffperiod(sta1,sta2,dist,nf=1,itr=1,gvel=gvel,iper=iper)
		cper_0 = copy(cper); iper_0 = copy(iper)

		gs2 = gridspec.GridSpec(1,1,wspace=0.2,hspace=0)
		ax2 = fig.add_subplot(gs2[0,0])

		ax2.imshow(amp,extent=[min(iper_amp),max(iper_amp),vmin,vmax],aspect='auto')
		ax2.plot(iper,gvel,color='k',label='group')
		ax2.plot(iper,pvel,color='r',label='phase')
		ax2.plot([cfp,cfp],[vmin,vmax],color='grey')

		ax2.set_xlim(xlm)
		ax2.set_ylim(vmin,vmax)
		ax2.set_xlabel('Instantaneous period [s]')
		ax2.set_ylabel('Velocity [km/s]')
		ax2.set_title('Raw FTAN')

		#### clean FTAN and clean curves
		# amplitudes, periods, velocities
		cper,iper,gvel,pvel = read_DISP(sta1,sta2,nf=2,itr=1)
		_,iper_amp,_,_ = read_DISP(sta1,sta2,nf=2,itr=0)
		amp,vmin,vmax = read_AMP(sta1,sta2,iper_amp,nf=2)
		cfp = calc_cutoffperiod(sta1,sta2,dist,nf=2,itr=1,gvel=gvel,iper=iper)

		gs3 = gridspec.GridSpec(1,1,wspace=0.2,hspace=0)
		ax3 = fig.add_subplot(gs3[0,0])

		ax3.imshow(amp,extent=[min(iper_amp),max(iper_amp),vmin,vmax],aspect='auto')
		ax3.plot(iper,gvel,color='k',label='group')
		ax3.plot(iper,pvel,color='r',label='phase')
		ax3.plot([cfp,cfp],[vmin,vmax],color='grey')

		# add SNR (calculated before aftan)
		psnr,snr = read_snr_precalc(sta1,sta2)
		nom_to_inst = interp1d(cper,iper,fill_value='extrapolate')
		ipsnr = nom_to_inst(psnr)
		more_pts = interp1d(ipsnr,snr,fill_value='extrapolate')
		snr_fine = more_pts(iper)

		ax33 = ax3.twinx()
		ax33.plot(iper,snr_fine,color='m')
		for tl in ax33.get_yticklabels():
			tl.set_color('m')

		ax3.set_xlabel('Instantaneous period [s]')
		ax3.legend(fontsize=9,loc='upper right')
		ax3.set_title('Clean FTAN')
		ax3.set_xlim(xlm)
		ax3.set_ylim(vmin,vmax)

		#### station map
		gs4 = gridspec.GridSpec(1,1,wspace=0.2,hspace=0.0)
		ax4 = fig.add_subplot(gs4[0,0])
		plot_basemap(ax4)
		ll1 = get_sta_coords(sta1); ll2 = get_sta_coords(sta2); nm = (sta1,sta2)
		x = (ll1[1],ll2[1]); y = (ll1[0],ll2[0])
		ax4.plot(x,y,'^-',color='k',ms=10,mfc='w',mew=1,zorder=100)
		for lon,lat,label in zip(x, y, nm):
			ax4.text(lon, lat, label, ha='center', va='bottom', fontsize=7, weight='bold')

		#### nominal vs observed periods
		gs5 = gridspec.GridSpec(1,1,wspace=0.2,hspace=0.0)
		ax5 = fig.add_subplot(gs5[0,0])
		ax5.plot(cper_0,iper_0,label='raw')
		ax5.plot(cper,iper,label='clean')
		ax5.grid()
		ax5.set_xlabel('Nominal period [s]')
		ax5.set_ylabel('Instantaneous period [s]')
		ax5.legend(fontsize=7,loc='lower right')

		# make things not overlap
		gs2.update(left=0.1,right=0.4)
		gs3.update(left=0.45,right=0.75)
		gs4.update(left=0.84,right=0.98,bottom=0.50)
		gs5.update(left=0.84,right=0.98,top=0.44)

		# plot title
		ndays = get_ndays(sta1,sta2)
		fig.suptitle('%s.%s - %s.%s: %.2f km, %i days' % (find_sta_network(sta1,inv), sta1, \
					find_sta_network(sta2,inv), sta2, dist, ndays), fontsize=14)

		pdf.savefig(fig)
		plt.close()

	pdf.close()

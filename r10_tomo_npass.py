import numpy as np
import pickle
import mods.PAT.utils as ut
import mods.PAT.files as vr
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from matplotlib.backends.backend_pdf import PdfPages
import ScientificColourMaps6 as SCM6
from copy import copy
import os, sys

####
# tomography-ify dispersion curves with n passes/cleaning
#   input files get made and remade iteratively
####

# dirs
idir = 'rt_inputs'
odir = 'rt_outputs'
rt_exe = os.path.expanduser('~/Patagonia/bin/RayTomo/bin/tomo_sp_cu_s')

# for plotting
plt.ioff()
ext = [-77,-68,-54.5,-42]
min_res = 0.02  # may need adjusting? but seems ok for now
min_dens = 10

# tomo parameters:
ia = 'ani'
vtypes = ['phase'] #,'group']
periods = [8.,14.,20.,26.,32.,40.]
npass = 4
grid_steps = [0.3,0.3,0.3,0.3]
minSNR = 5
corr_lengths = [50,50,50,50]  # sigma?
alphas = [600,400,250,150]  # smoothing
betas = [25,25,25,25]  # damping?
fancy_names = ['1st','2nd','3rd','4th']
std_factor = 3
# for anisotropy smoothing:
alpha_ani = 600; corr_ani = 300


assert np.all([len(grid_steps)==npass,len(corr_lengths)==npass,\
               len(alphas)==npass,len(betas)==npass,len(fancy_names)>=npass]),\
               'not enough parameters for %i passes' % npass

skip_stations = [] # ['MG01','VTDF','DSPA','VOH01','DGER','RGND']
skip_pairs = []  # ['sta1-sta2']

# set up pdfs for plots
pdfs = {}
for vtype in vtypes:
    pdfs[vtype] = PdfPages('../Plots/%i-pass-rttomo-%s_%s.pdf' % (npass,ia,vtype))

# load dispersion curves
pickle_file = 'COR/disp_curves.pickle'
f = open(pickle_file,'rb')
curves = pickle.load(f)
f.close()
sta1 = np.array([c.station1.name for c in curves])
sta2 = np.array([c.station2.name for c in curves])


# loop periods/vtypes, do all passes
for period in periods:
    for vtype in vtypes:
        print('Doing {} s {} velocities'.format(period,vtype))
###### TODO: recalculate alphas per period/vtype?
        # start with clean skip_pairs for this period/vtype
        skippairs = copy(skip_pairs)
        # loop passes
        for i in range(npass):
            print('\t{} pass, skipping {} pairs'.format(fancy_names[i],len(skippairs)))
            # write input files for phase and group at all periods, excluding selected stations and pairs
            ofile = open('%s/data_%s_%s_%i_%i.rt' % (idir,ia,vtype,i,period),'w') # data
            tfile = open('%s/inds_%s_%s_%i_%i.rt' % (idir,ia,vtype,i,period),'w') # pair tracking
            count = 0  # curve index tracking
            # loop curves
            for j in range(len(curves)):
                c = curves[j]
                if sta1[j] in skip_stations or sta2[j] in skip_stations or \
                    '%s-%s' % (sta1[j],sta2[j]) in skippairs:  # skip bad stations and pairs
                    continue
                v,_,snr = c.filtered_vel_sdev_SNR(period,vtype=vtype)
                if not np.isnan(v) and snr > minSNR:
                    ofile.write('%6i %8.3f %8.3f %8.3f %8.3f %8.4f 1 1\n' % \
                        (count,c.station1.coord[1],c.station1.coord[0],\
                        c.station2.coord[1],c.station2.coord[0],v))
                    tfile.write('%6i %8s %8s\n' % (count,c.station1.name,c.station2.name))
                    count += 1
            ofile.close()
            tfile.close()

            # run raytomo
            dfile = '%s/data_%s_%s_%i_%i.rt' % (idir,ia,vtype,i,period)
            ofile_root = '%s/tomo_%s_%s_%i' % (odir,ia,vtype,i)

            if ia == 'iso':
                to_run = 'me\n4\n5\n-55.2 -42.9 0.3\n6\n-75.9 -67.5 0.3\n10\n0.25\n2.0\nR\n%s\n0.3\n1.0\n12\n%.1f %.1f %.1f 999\n19\n25\n26\nx\ngo\nEOF' % \
            (vtype.capitalize()[0],alphas[i],betas[i],corr_lengths[i])
            elif ia == 'ani':
                to_run = 'me\n4\n5\n-55.2 -42.9 0.3\n6\n-75.9 -67.5 0.3\n10\n0.25\n2.0\nR\n%s\n0.5\n2.0\n11\n1\n12\n%.1f %.1f %.1f 999\n13\n%.1f %.1f 999\n19\n25\n26\nx\ngo\nEOF' % \
            (vtype.capitalize()[0],alphas[i],betas[i],corr_lengths[i],alpha_ani,corr_ani)

            os.system('%s %s %s %i << EOF 1>junk 2>rt_err \n%s' % \
                        (rt_exe,dfile,ofile_root,period,to_run))

            if i == npass - 1: # last pass, plot and move on to next period
                fig = plt.figure(figsize=(11.3,7.2))  # 9.2
                ax1 = plt.subplot2grid((1,2),(0,0))  # 3
                ax2 = plt.subplot2grid((1,2),(0,1))  # 3

                # read isotropic velocities
                lon,lat,vph = np.loadtxt('%s_%i.1' % (ofile_root,period),usecols=(0,1,2),unpack=True)

                # set up grid
                xi = np.arange(min(lon),max(lon),0.3)
                yi = np.arange(min(lat),max(lat),0.3)
                Xi,Yi = np.meshgrid(xi,yi)

                # triangulate, interpolate velocities
                triang = tri.Triangulation(lon,lat)
                interp_vph = tri.LinearTriInterpolator(triang,vph)
                vph_grid = interp_vph(Xi,Yi)

                # read and interpolate response amplitudes from resolution analysis
                lon,lat,rea = np.loadtxt('%s_%i.rea' % (ofile_root,period),usecols=(0,1,4),unpack=True)
                triang = tri.Triangulation(lon,lat)
                interp_rea = tri.LinearTriInterpolator(triang,rea)
                rea_grid = interp_rea(Xi,Yi)

                lon,lat,dens = np.loadtxt('%s_%i.res' % (ofile_root,period),usecols=(0,1,2),unpack=True)
                triang = tri.Triangulation(lon,lat)
                interp_dens = tri.LinearTriInterpolator(triang,dens)
                dens_grid = interp_dens(Xi,Yi)

                # mask velocities where no resolution
                vph_grid.mask[rea_grid < min_res] = True
                vph_grid.mask[dens_grid < min_dens] = True

                # plot masked array
                con1 = ax1.imshow(vph_grid,extent=[min(xi),max(xi),min(yi),max(yi)],\
                        cmap=SCM6.roma,origin='lower')
                ut.plot_basemap(ax1,bbox=ext)
                cb1 = plt.colorbar(con1,ax=ax1,fraction=0.04,orientation='horizontal')
                con1.set_clim(vr.vlims[period])

                if ia == 'ani':
                    lon,lat,a_amp,a_psi = np.loadtxt('%s_%i.1' % (ofile_root,period),\
                            usecols=(0,1,5,6),unpack=True)
                    #u = a_amp*np.cos(np.degrees(a_psi))*1e4
                    #v = a_amp*np.sin(np.degrees(a_psi))*1e4
                    #for il in range(len(lon)):
                    #    ax1.plot([lon[il] - u[il]/2,lon[il] + u[il]/2],[lat[il] - v[il]/2,lat[il] + v[il]/2],color='k',lw=0.4)
                    ax1.quiver(lon,lat,a_amp*np.cos(np.degrees(a_psi)),a_amp*np.sin(np.degrees(a_psi)))

                con2 = ax2.imshow(dens_grid,extent=[min(xi),max(xi),min(yi),max(yi)],\
                        cmap=SCM6.hawaii,origin='lower')
                ut.plot_basemap(ax2,bbox=ext)
                cb2 = plt.colorbar(con2,ax=ax2,fraction=0.04,orientation='horizontal')

                plt.suptitle('%i sec, %i used/%i skipped, alpha %i' % \
                            (period,count,len(skippairs),alphas[i]))

                pdfs[vtype].savefig(fig)
                plt.close()

            # otherwise, read resids, track outlier paths, and add to skip_pairs
            rid,resids = np.loadtxt('%s_%i.resid' % (ofile_root,period),usecols=(0,8),unpack=True)
            rstd = np.std(resids)
            to_skip = np.where(resids > rstd*std_factor)[0]

            s1,s2 = np.loadtxt('%s/inds_%s_%s_%i_%i.rt' % (idir,ia,vtype,i,period),\
                                usecols=(1,2),unpack=True,dtype=str)
            for j in to_skip:
                skippairs.append('%s-%s' % (s1[j],s2[j]))

os.system('rm -f junk')
os.system('rm -f rt_err')

for vtype in vtypes:
    pdfs[vtype].close()

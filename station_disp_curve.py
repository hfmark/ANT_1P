import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
import mods.PAT.anttomo as ant
import mods.PAT.files as vr
from matplotlib.backends.backend_pdf import PdfPages
import pickle
import os, sys

####
# read maps from EQ and ANT tomo, read stations, get disp curves for station locations
####

plot_group = False

# read in EQ tomography outputs
aphv = sio.loadmat(os.path.expanduser('~/Patagonia/EQ_tomo/helmholtz_stack_all.mat'),variable_names='avgphv')  # current last version, with station exclusions set
#aphv = sio.loadmat(os.path.expanduser('~/Patagonia/EQ_tomo/helmholtz_stack_ENAPnoamp_noMG.mat'),variable_names='avgphv')
aphv = aphv['avgphv']
emaps = {aphv[:,i]['period'][0][0][0]:ant.EQVelocityMap(aphv[:,i]) for i in range(len(aphv[0]))}

# read in ANT outputs
f = open('output/4-pass-tomography_phase.pickle','rb')
vmaps = pickle.load(f)
f.close()

if plot_group:
    f = open('output/4-pass-tomography_group.pickle','rb')
    vmaps_group = pickle.load(f)
    f.close()

opdf = '../Plots/disp_curve_stations.pdf'
pdf = PdfPages(opdf)

# get station info (coords/names for all stations)
stnm = np.loadtxt(vr.sta_list, usecols=(0,), dtype=str)
slon,slat = np.loadtxt(vr.sta_list, usecols=(1,2), unpack=True)

# periods for both
eper = [emaps[k].period for k in emaps.keys()]
aper = [vmaps[k].period for k in vmaps.keys()]
if plot_group:
    gper = [vmaps_group[k].period for k in vmaps_group.keys()]

plt.ioff()
# for a station, get velocities from ANT and EQ and plot both? to see how they compare?
for i in range(len(stnm)):
    try:
        evel = [emaps[k].interp_velocity(slon[i],slat[i])[0] for k in emaps.keys()]
        estd = [emaps[k].interp_velocity(slon[i],slat[i])[1] for k in emaps.keys()]
        avel = [vmaps[k].interp_velocity(slon[i],slat[i]) for k in vmaps.keys()]
        astd = [max(2./vmaps[k].density_interp(slon[i],slat[i])[0], 0.05) for k in vmaps.keys()]
        astd = [min(0.2,e) for e in astd]
        if plot_group:
            gvel = [vmaps_group[k].interp_velocity(slon[i],slat[i]) for k in vmaps_group.keys()]
            # NOTE this assumes set of ant phase periods includes all of the group periods
            gstd = [max(2./vmaps[k].density_interp(slon[i],slat[i])[0], 0.05) for k in vmaps_group.keys()]
            gstd = [min(0.2,e) for e in gstd]
            gstd = 2.5*np.array(gstd)
    except:
        print(stnm[i])
        continue  # probably a bounds error for interpolation

    fig = plt.figure()
    #plt.plot(aper, avel, '.-', label='ambient, phase')
    plt.errorbar(aper, avel, yerr=astd, fmt='.-', label='ambient, phase')
    if plot_group:
        #plt.plot(gper, gvel, '.-', label='ambient, group')
        plt.errorbar(gper, gvel, yerr=gstd, fmt='.-', label='ambient, group')
    #plt.plot(eper, evel, '.-', label='EQ')
    plt.errorbar(eper, evel, yerr=estd, fmt='.-', label='EQ, phase')
    plt.legend(fontsize=9)
    plt.title(stnm[i])
    plt.xlabel('Period [s]')
    plt.ylabel('Phase velocity [km/s]')
    plt.xlim(5,105)
    #plt.ylim(2.0,4.7)
    plt.ylim(0.5,4.5)
    pdf.savefig(fig)
    plt.close()

pdf.close()

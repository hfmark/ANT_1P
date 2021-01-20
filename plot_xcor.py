import numpy as np
import matplotlib.pyplot as plt
from obspy import read
import os, sys

####
# plot xcor filtered at a few different bands
# compare linear stack and pws/various wu values
####

plt.ion()

period_bands = [[4, 7], [7, 15], [10, 22], [15, 30], [20, 35]]  # from pysismo
signal_vmin = 1.7; signal_vmax = 5.2
trail = 200.
nplot = len(period_bands) + 1

pws = read('test_pws/COR/CURI_MG02/tf-pws_0.5.sac')
lin = read('test_pws/COR/CURI_MG02/tl.sac')

lin[0].data = lin[0].data - np.mean(lin[0].data)

tmax = pws[0].stats.sac.dist/signal_vmin + trail*2
tmax = min(1.05*tmax, pws[0].times().max())
xlim = (0,tmax)

fig = plt.figure()
axlist = [fig.add_subplot(nplot,1,i) for i in range(1, nplot+1)]
axlist[0].get_figure().subplots_adjust(hspace=0)

# limits of yaxis from min/max of the window we actually care about
mask = (pws[0].times() >= min(pws[0].stats.sac.dist/signal_vmax,tmax)) & (pws[0].times() <= tmax)

#ylim = (lin[0].data[mask].min(), lin[0].data[mask].max())
ylim = (-1,1)
# plot full xcor
axlist[0].plot(lin[0].times(),lin[0].data/max(abs(lin[0].data)))

#ylim = (pws[0].data[mask].min(), pws[0].data[mask].max())
# plot full xcor
axlist[0].plot(pws[0].times(),pws[0].data/max(abs(pws[0].data)))

# signal window
tsig_min = pws[0].stats.sac.dist/signal_vmax
tsig_max = pws[0].stats.sac.dist/signal_vmin

axlist[0].plot(2*[tsig_min], ylim, color='k')
axlist[0].plot(2*[tsig_max], ylim, color='k')

# axis formatting
axlist[0].set_xlim(xlim)
axlist[0].set_ylim(ylim)
axlist[0].grid(True)
axlist[0].set_xticklabels([])

# plot for different bands:
for ax, (tmin,tmax) in zip(axlist[1:], period_bands):
    lastplot = ax is axlist[-1]

    fltp = pws.copy().filter('bandpass',freqmin=1/tmax,freqmax=1/tmin)
    fltl = lin.copy().filter('bandpass',freqmin=1/tmax,freqmax=1/tmin)

    # limits of yaxis from min/max of the window we actually care about
    #ylim = (fltl[0].data[mask].min(), fltl[0].data[mask].max())

    # plot full xcor
    ax.plot(fltl[0].times(),fltl[0].data/max(abs(fltl[0].data)))
    ax.plot(fltp[0].times(),fltp[0].data/max(abs(fltp[0].data)))

    ax.plot(2*[tsig_min], ylim, color='k')
    ax.plot(2*[tsig_max], ylim, color='k')

    # axis formatting
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.grid(True)
    if not lastplot:
        axlist[0].set_xticklabels([])

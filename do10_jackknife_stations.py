import numpy as np
import mods.PAT.anttomo as ant
import pickle
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from copy import copy
import itertools as it
import os, sys

####
# test tomography with one set of params, jackknifing excluded station choice
####

# load dispersion curves
pickle_file = 'COR/disp_curves.pickle'
f = open(pickle_file, 'rb')
curves = pickle.load(f)
f.close()

# get list of all stations involved
stations = np.unique(np.hstack(([c.station1.name for c in curves],[c.station2.name for c in curves])))

# choose tomo parameters
period = 26.0
vtype = 'phase'
grid_step = 0.25
minspectsnr = 5.0
corr_length = 50
alpha = 250
beta = 50
lam = 0.3
minresheight = 0.02

# set up output pdf
opdf = 'output/station-jk-%.1f_%s.pdf' % (period, vtype)
if os.path.exists(opdf):
    iq = input('outfile already present. Replace? [y]/n') or 'y'
    if iq == 'n' or iq == 'N':
        sys.exit()
pdf = PdfPages(opdf)

skippairs = []

# loop stations for skip_stations, do tomo, plot with a title that says what was skipped
for sta in stations:
    skip_stations = [sta]
    print('testing: %s, %.1f s' % (sta, period))

    v = ant.VelocityMap(dispersion_curves=curves,
                                    period=period,
                                    skipstations=skip_stations,
                                    skippairs=skippairs,
                                    resolution_fit='gaussian',
                                    min_resolution_height=minresheight,
                                    verbose=False,
                                    lonstep=grid_step,
                                    latstep=grid_step,
                                    minspectSNR=minspectsnr,
                                    correlation_length=corr_length,
                                    alpha=alpha,
                                    beta=beta,
                                    lambda_=lam,
                                    vtype=vtype)

    title = ("SKIPPING {0}: {1} s, {2} x {2} deg, min SNR = {3}, corr. length "
             "= {4} km, alpha = {5}, beta = {6}, lambda = {7} ({8} paths)")
    title = title.format(sta, period, grid_step, minspectsnr, corr_length,
                         alpha, beta, lam, len(v.paths))
    fig = v.plot(title=title, showplot=False)

    # add to pdf
    pdf.savefig(fig)
    plt.close()

# close pdf
pdf.close()

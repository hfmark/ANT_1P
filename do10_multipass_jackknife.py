import numpy as np
import mods.PAT.anttomo as ant
import pickle
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from copy import copy
import os, sys

####
# tomography-ify dispersion curves with n passes/cleaning
####

# load dispersion curves
pickle_file = 'COR/disp_curves.pickle'
f = open(pickle_file, 'rb')
curves = pickle.load(f)
f.close()

# get list of all stations involved
stations = np.unique(np.hstack(([c.station1.name for c in curves],[c.station2.name for c in curves])))

plt.ioff()

# set tomography parameters to loop over for n passes
period = 15.0

_vtype = 'phase'
_npass = 4
_grid_steps = 0.3*np.ones(_npass)
_minspectsnrs = 5.0*np.ones(_npass)
_corr_lengths = 50*np.ones(_npass) # 50 seems better than 100
_alphas = (600,400,250,150)
_betas = 25*np.ones(_npass)
_lambdas = 0.3*np.ones(_npass)
_fancy_names = ('1st','2nd','3rd','4th')
minresheight = 0.02
lonmin = -75.9  # override automatic gridding so that things match EQ tomo
latmin = -55.2
nlon = 29
nlat = 42

assert np.all([len(_grid_steps)==_npass,len(_minspectsnrs)==_npass,len(_corr_lengths)==_npass,\
               len(_alphas)==_npass,len(_betas)==_npass,len(_lambdas)==_npass,\
               len(_fancy_names)>=_npass]), 'not enough parameters for %i passes' % _npass
# set up output pdf
opdf = '../Plots/ANT_station-jk-%.1f_%s.pdf' % (period, _vtype)
if os.path.exists(opdf):
    iq = input('outfile already present. Replace? [y]/n') or 'y'
    if iq == 'n' or iq == 'N':
        sys.exit()
pdf = PdfPages(opdf)

_skip_pairs = []

#_skip_stations = []
_skip_stations = ['MG01','VTDF','DSPA','VOH01','DGER','RGND']

for sta in stations:
    if sta in _skip_stations:
        continue  # skip the ones that are always skipped
    skip_stations = copy(_skip_stations)
    skip_stations.append(sta)
    print('testing: %s, %.1f s' % (sta, period))

    skippairs = copy(_skip_pairs)
    for passnb in range(_npass):
        s = ("{} pass (rejecting {} pairs): grid step = {}, min SNR = {}, "
            "corr. length = {} km, alpha = {}, beta = {}, lambda = {}")
        print(s.format(_fancy_names[passnb], len(skippairs),
                        _grid_steps[passnb], _minspectsnrs[passnb],
                        _corr_lengths[passnb], _alphas[passnb],
                        _betas[passnb], _lambdas[passnb]))

        # Performing the tomographic inversion to produce a velocity map
        # at period = *period* , with parameters given above:
        # - *lonstep*, *latstep* control the internode distance of the grid
        # - *minnbtrimester*, *maxsdev*, *minspectSNR*, *minspectSNR_nosdev*
        #   correspond to the selection criteria
        # - *alpha*, *corr_length* control the spatial smoothing term
        # - *beta*, *lambda_* control the weighted norm penalization term
        #
        # Note that if no value is given for some parameter, then the
        # inversion will use the default value defined in the configuration
        # file.
        #
        # (See doc of VelocityMap for a complete description of the input
        # arguments.)

        try:
            v = ant.VelocityMap(dispersion_curves=curves,
                                    period=period,
                                    skipstations=skip_stations, # not _skip_stations!
                                    skippairs=skippairs,
                                    resolution_fit='gaussian',
                                    min_resolution_height=minresheight,
                                    verbose=False,
                                    lonstep=_grid_steps[passnb],
                                    latstep=_grid_steps[passnb],
                                    minspectSNR=_minspectsnrs[passnb],
                                    correlation_length=_corr_lengths[passnb],
                                    alpha=_alphas[passnb],
                                    beta=_betas[passnb],
                                    lambda_=_lambdas[passnb],
                                    vtype=_vtype,
                                    lonmin=lonmin,
                                    latmin=latmin,
                                    nlon=nlon,
                                    nlat=nlat)
        except ant.CannotPerformTomoInversion as err:
            print("Cannot perform tomo inversion: {}".format(err))
            for fig in periodfigs:
                plt.close(fig)
            # next period
            break

        if passnb < _npass-1:
            # pairs whose residual is > 3 times the std dev of
            # the residuals are rejected from the next pass
            residuals = v.traveltime_residuals()
            maxresidual = 3 * residuals.std()
            badpairs = [(c.station1.name, c.station2.name)
                        for c, r in zip(v.disp_curves, residuals)
                        if abs(float(r)) > maxresidual]
            for ib in range(len(badpairs)):
                skippairs.append(badpairs[ib])

    title = ("SKIPPING {0}: {1} s, {2} x {2} deg, min SNR = {3}, corr. length "
             "= {4} km, alpha = {5}, beta = {6}, lambda = {7} ({8} paths)")
    title = title.format(sta, period, _grid_steps[-1], _minspectsnrs[-1], _corr_lengths[-1],
                         _alphas[-1], _betas[-1], _lambdas[-1], len(v.paths))
    fig = v.plot(title=title, showplot=False)

    # add to pdf
    pdf.savefig(fig)
    plt.close()


# closing pdf file
#pagenbs = range(pdf.get_pagecount())
pdf.close()

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

plt.ioff()

# set tomography parameters to loop over for n passes
#_vtype = 'phase'
_vtype = 'group'

if _vtype == 'phase':
    _periods = [8., 9., 10., 11., 12., 13., 14., 15., 16., 17., 18., 19., 20., 22., 24., 26., 28., 30., 32., 34.]
elif _vtype == 'group':
    _periods = [8., 9., 10., 11., 12., 13., 14., 15., 16., 17., 18., 19., 20., 22., 24., 26., 28., 30.]

#_periods = [8.0, 14.0, 20.0, 25.0, 32.0]  # tests only

_npass = 4
_grid_steps = 0.3*np.ones(_npass)
_minspectsnrs = 5.0*np.ones(_npass)
_corr_lengths = 50*np.ones(_npass) # 50 seems better than 100
_alphas = (600,400,300,200)
_betas = 25*np.ones(_npass)
_lambdas = 0.3*np.ones(_npass)
_fancy_names = ('1st','2nd','3rd','4th')
minresheight = 0.02
lonmin = -76.  # override automatic gridding so that things match EQ tomo
latmin = -55.
nlon = 31
nlat = 42

assert np.all([len(_grid_steps)==_npass,len(_minspectsnrs)==_npass,len(_corr_lengths)==_npass,\
               len(_alphas)==_npass,len(_betas)==_npass,len(_lambdas)==_npass,\
               len(_fancy_names)>=_npass]), 'not enough parameters for %i passes' % _npass

#_skip_pairs = [('AY01','GUMN')]
#_skip_stations = []
#_skip_stations = ['CHN01','VOH01'] # 'RRS01' 'COC01'  # from previous version with linear stack
_skip_stations = ['ANMA','BQON','DGER','LAVG','LOSC','MJTA','RGND','TROP']
if _vtype == 'phase':
#    _skip_pairs = [('ANMA','DGER'),('COYC','TAPA'),('CURI','RPTE'),('DGER','GO08'),('MG04','RGND'),\
#                    ('RRS01','VOH01'),('RMG01','VCC01'),('AMG01','COC01'),('CHN01','VOH01'),\
#                    ('AY01','LSMN'),('LSR01','VOH01'),('GO08','GRAF')]
                    # 8s(x3), 14s(x4), 20s(x4), 26s(x1)
#    _skip_pairs = [('ANMA','DGER'),('COYC','TAPA'),('CURI','RPTE'),('DGER','GO08'),('MG04','RGND'),\
#                    ('COYC','MG04'),('GO08','GRAF')]
    _skip_pairs = [('RMG01','RPR01')]  # this was from pws but maybe still holds?
elif _vtype == 'group':
    _skip_pairs = []

# set up output pdf file and pickle file
opdf = '../Plots/%i-pass-tomography_%s.pdf' % (_npass, _vtype)
opickle = 'output/%i-pass-tomography_%s.pickle' % (_npass, _vtype)
if os.path.exists(opdf) or os.path.exists(opickle):
    iq = input('outfile(s) already present. replace? [y]/n') or 'y'
    if iq == 'n' or iq == 'N':
        sys.exit()
pdf = PdfPages(opdf)

vmaps = {} # dict for final maps at each period

for period in _periods:
    print("\nDoing period {} s".format(period))

    periodfigs = []

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
                                    skipstations=_skip_stations,
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
        else:
            # adding velocity map to the dict of final maps
            vmaps[period] = v
            maxresidual = None

        # creating a figure summing up the results of the inversion:
        # - 1st panel = map of velocities or velocity anomalies
        # - 2nd panel = map of interstation paths and path densities
        #               (highlighting paths with large diff
        #                between obs/predicted travel-time)
        # - 3rd panel = resolution map
        #
        # See doc of VelocityMap.plot(), VelocityMap.plot_velocity(),
        # VelocityMap.plot_pathdensity(), VelocityMap.plot_resolution()
        # for a detailed description of the input arguments.

        title = ("Period = {0} s, {1} pass, grid {2} x {2} deg, "
                 "min SNR = {3}, corr. length = {4} km, alpha = {5}, "
                 "beta = {6}, lambda = {7} ({8} paths, {9} rejected)")
        title = title.format(period, _fancy_names[passnb],
                             _grid_steps[passnb], _minspectsnrs[passnb],
                             _corr_lengths[passnb], _alphas[passnb],
                             _betas[passnb], _lambdas[passnb], len(v.paths),
                             len(skippairs))

        if passnb == _npass-1:
            # we highlight paths that will be rejected
            fig = v.plot(title=title, showplot=False,
                         highlight_residuals_gt=maxresidual)
            # appending fig to figures of period
            periodfigs.append(fig)

    else:
#       # if we did not break from loop:
#       # let's compare the 2-pass tomography with a one-pass tomography
#       s = ("One-pass tomography: grid step = {}, min SNR = {}, "
#            "corr. length = {} km, alpha = {}, beta = {}, lambda = {}")
#       print(s.format(_grid_steps[-1], _minspectsnrs[-1], _corr_lengths[-1],
#                       _alphas[-1], _betas[-1], _lambdas[-1]))
#
#       # tomographic inversion
#       try:
#           v = ant.VelocityMap(dispersion_curves=curves,
#                                   period=period,
#                                   verbose=False,
#                                   resolution_fit='gaussian',
#                                   min_resolution_height=minresheight,
#                                   lonstep=_grid_steps[-1],
#                                   latstep=_grid_steps[-1],
#                                   minspectSNR=_minspectsnrs[-1],
#                                   correlation_length=_corr_lengths[-1],
#                                   alpha=_alphas[-1],
#                                   beta=_betas[-1],
#                                   lambda_=_lambdas[-1],
#                                   vtype=_vtype,
#                                   lonmin=lonmin,
#                                   latmin=latmin,
#                                   nlon=nlon,
#                                   nlat=nlat)
#       except CannotPerformTomoInversion as err:
#           print("Cannot perform tomo inversion: {}".format(err))
#           for fig in periodfigs:
#               plt.close(fig)
#           # next period
#           continue
#
#       # figure
#       title = ("Period = {0} s, one pass, grid {1} x {1} deg, "
#                "min SNR = {2}, corr. length = {3} km, alpha = {4}, "
#                "beta = {5}, lambda = {6} ({7} paths)")
#       title = title.format(period, _grid_steps[-1], _minspectsnrs[-1],
#                            _corr_lengths[-1], _alphas[-1],
#                            _betas[-1], _lambdas[-1], len(v.paths))
#       fig = v.plot(title=title, showplot=False)
#
#       # appending fig to figures of period
#       periodfigs.append(fig)

        # exporting figs to pdf
        for fig in periodfigs:
            pdf.savefig(fig)
            plt.close(fig)

# closing pdf file
#pagenbs = range(pdf.get_pagecount())
pdf.close()

# exporting final maps (using pickle) as a dict:
# {period: instance of ant.VelocityMap}

print("\nExporting final velocity maps to file: " + opickle)
f = open(opickle,'wb')
pickle.dump(vmaps, f, protocol=2)
f.close()

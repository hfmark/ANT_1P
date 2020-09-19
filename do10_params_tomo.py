import numpy as np
import mods.PAT.anttomo as ant
import pickle
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from copy import copy
import itertools as it
import os, sys

####
# test tomography on dispersion curves with a variety of parameters
####

# load dispersion curves
pickle_file = 'COR/disp_curves.pickle'
f = open(pickle_file, 'rb')
curves = pickle.load(f)
f.close()

# set tomography parameters to loop over for n passes
_periods = [8.0, 20.0, 25.0, 32.0]
_vtype = 'phase'
_grid_steps = [0.25, 0.5]
_minspectsnrs = [5.0]
_corr_lengths = [25,50,100,200]
_alphas = [100,200,400,800]
_betas = [20,50,100]
_lambdas = [0.3]
minresheight = 0.02

_skip_stations = []
#_skip_pairs = [] #[('AY01','GUMN')]
_skip_pairs = [('ANMA','DGER'),('COYC','TAPA'),('CURI','RPTE'),('DGER','GO08'),('MG04','RGND'),\
		 ('RRS01','VOH01'),('RMG01','VCC01'),('AMG01','COC01'),('CHN01','VOH01'),\
		 ('AY01','LSMN'),('LSR01','VOH01'),('GO08','GRAF')]

# set up output pdf file and pickle file
opdf = '../Plots/test-tomography_%s.pdf' % (_vtype)
oLcurve = 'output/Lcurve_data_%s.dat' % (_vtype)
if os.path.exists(opdf) or os.path.exists(oLcurve):
    iq = input('outfile(s) already present. replace? [y]/n') or 'y'
    if iq == 'n' or iq == 'N':
        sys.exit()
plt.ioff()
pdf = PdfPages(opdf)
fout = open(oLcurve, 'w')

param_lists = it.product(_periods, _grid_steps, _minspectsnrs, _corr_lengths, \
				_alphas, _betas, _lambdas)
param_lists = list(param_lists)

for period, grid_step, minspectSNR, corr_length, alpha, beta, lambda_ in param_lists:
    s = ("Period = {} s, grid step = {}, min SNR = {}, corr. length "
             "= {} km, alpha = {}, beta = {}, lambda = {}")
    print(s.format(period, grid_step, minspectSNR, corr_length, alpha, beta, lambda_))

    skippairs = copy(_skip_pairs)
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

    v = ant.VelocityMap(dispersion_curves=curves,
                                    period=period,
                                    skipstations=_skip_stations,
                                    skippairs=skippairs,
                                    resolution_fit='gaussian',
                                    min_resolution_height=minresheight,
                                    verbose=False,
                                    lonstep=grid_step,
                                    latstep=grid_step,
                                    minspectSNR=minspectSNR,
                                    correlation_length=corr_length,
                                    alpha=alpha,
                                    beta=beta,
                                    lambda_=lambda_,
                                    vtype=_vtype)

    # get misfit and norm for L curves
    misfit = v.velocity_residuals()
    norm = v.model_norm()

    fout.write('%5.1f%5.1f%5.1f%8.1f%8.1f%8.1f%8.2f%7.3f%10.3f%10.3f\t%i\n' % \
            (period, grid_step, minspectSNR, corr_length, alpha, beta, lambda_, \
             minresheight, abs(misfit).sum(), norm, len(misfit)))

    # creating a figure summing up the results of the inversion:
    # - 1st panel = map of velocities or velocity anomalies
    # - 2nd panel = map of interstation paths and path densities
    # - 3rd panel = resolution map
    #
    # See doc of VelocityMap.plot(), VelocityMap.plot_velocity(),
    # VelocityMap.plot_pathdensity(), VelocityMap.plot_resolution()
    # for a detailed description of the input arguments.

    title = ("Period = {0} s, grid {1} x {1} deg, min SNR = {2}, corr. length "
             "= {3} km, alpha = {4}, beta = {5}, lambda = {6} ({7} paths)")
    title = title.format(period, grid_step, minspectSNR, corr_length,
                         alpha, beta, lambda_, len(v.paths))
    fig = v.plot(title=title, showplot=False)

    # exporting plot to pdf
    pdf.savefig(fig)
    plt.close()

# closing pdf file and Lcurve file
pdf.close()
fout.close()

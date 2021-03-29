import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os, sys

####
# read in misfit and norm for different tomo parameters, plot L curves and find minima
####

A = np.loadtxt('output/Lcurve_data_phase.dat')
dat = pd.DataFrame(A,columns=['period', 'grid_step', 'minsnr', 'corr_length', 'alpha', 'beta', 'lam', 'minresh', 'misfit', 'norm', 'npaths'])

plt.ion()
p1 = dat.period == 8.0
plt.figure()
plt.scatter(dat[p1].misfit/dat[p1].npaths, dat[p1].norm/dat[p1].npaths, c=dat[p1].corr_length)
plt.colorbar(fraction=0.04)

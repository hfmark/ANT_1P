import numpy as np
import mods.PAT.anttomo as ant
from obspy import read_inventory
from glob import glob
import pickle
import os, sys

####
# read aftan outputs, construct DispersionCurve() objects, pickle
####

curves = []

# find all station pairs with stacks
corfiles = np.array(glob('COR/*_full.SAC_2_DISP.1'))  # full crosscorr only
pairs = np.array(['_'.join(e.split('/')[-1].split('_')[1:3]) for e in corfiles])
inv = read_inventory('seed/dataless/combined.xml')

# loop pairs
for p in pairs:
    # check: do we have the final aftan outputs? *_full.SAC_2_DISP.1
    disp2file = 'COR/COR_%s_full.SAC_2_DISP.1' % p
    if not os.path.isfile(disp2file):
        print('%s is unfinished' % p)
        continue  # skip this pair

    print(p)

    # if so, get that filename, snr filename, and trifile list
    snrfile = 'COR/COR_%s_full.SAC_s_snr.cv.p.txt' % p
    trifiles = glob('COR/COR_%s_??.??.??.SAC_2_DISP.1' % p)
    sta1 = p.split('_')[0]; sta2 = p.split('_')[1]

    # make a dispersion curve and add to list of curves
    cv = ant.DispersionCurve(sta1, sta2, inv, disp2file, snrfile, trifiles)
    curves.append(cv)
    
# pickle the curves
ofile = 'COR/disp_curves.pickle'
f = open(ofile, 'wb')
pickle.dump(curves, f, protocol=2)
f.close()

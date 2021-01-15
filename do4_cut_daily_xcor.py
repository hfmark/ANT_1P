import numpy as np
from obspy import read
from glob import glob
from datetime import datetime
import os, sys

####
# tally up all unique station pairs
# for each pair:
    # check that there are at least [N] daily xcor files
    # if so, loop day files:
        # read each, slice in half and write positive and reversed-negative as separate files
        # put those in a directory named for the pair
####

# NOTE that the list of corfiles is LONG
# this could perhaps be made slightly less insane by listing month.year dirs and looping those,
# then sub-looping pairs for each one. os.makedirs would be ok with that, pairs would
# get added as the code progressed.

corfiles = np.array(glob('????.???/COR_D/*/*.SAC'))  # all daily crosscorr files
pairs = np.array(['_'.join(e.split('/')[-1].split('.')[0].split('_')[1:3]) for e in corfiles])
u_pairs = np.unique(pairs)

odir = 'COR/'

for p in u_pairs:
    if p.split('_')[0] == p.split('_')[1]:
        continue  # skip autocorrelations

    print(p)
    flist = corfiles[np.where(pairs == p)]  # all the full sym/asym daily crosscorrelation files
    p_odir = os.path.join(odir,p)
    os.makedirs(p_odir,exist_ok=True)
    for i in range(len(flist)):
        st = read(flist[i])  # read daily xcor file

        st1 = st.copy()  # copy and keep causaul half
        st1[0].data = st1[0].data[3000:]  # hard-coded to 1/2 assumed length of trace (6001)
        st2 = st.copy()  # copy and keep acausal half
        st2[0].data = st2[0].data[::-1]   # reverse
        st2[0].data = st2[0].data[3000:]  # and trim

        # rewrite halves
        day_str = '.'.join([flist[i].split('/')[0],flist[i].split('/')[-1].split('.')[0].split('_')[-1]])
        f_pos = '_'.join([p,day_str]) + '.SACpos'
        o_pos = os.path.join(p_odir,f_pos)
        st1.write(o_pos,'SAC')
        f_neg = '_'.join([p,day_str]) + '.SACneg'
        o_neg = os.path.join(p_odir,f_neg)
        st2.write(o_neg,'SAC')

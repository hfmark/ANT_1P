import numpy as np
from glob import glob
from datetime import datetime
import os, sys

####
# list all the station pairs
# for each pair:
    # glob all daily cross-correlation files
    # get dates for each, sort by date
    # write ALL to a file
    # get list of unique months
    # for each overlapping set of 3 months, write a list of just those daily files
# write some lists of lists for running things?
####

pairs = np.array(glob('COR/*/'))

for p in pairs:
    corfiles = np.array(glob(p+'*.SACpos'))  # this is only half the files, note
    daystrs = np.array([datetime.strptime('.'.join(e.split('/')[-1].split('_')[-1].split('.')[:-1]), '%Y.%b.%d') for e in corfiles])
    order = np.argsort(daystrs)
    corfiles = corfiles[order]
    daystrs = daystrs[order]

    yrmonth = np.array([e.strftime('%Y.%m') for e in daystrs])  # all *months* with data
    u_month = np.sort(np.unique(yrmonth))
    if len(u_month) < 4:  # not enough for two 3-month stacks
        continue
    
    tri_month = np.array(list(zip(u_month, u_month[1:], u_month[2:])))

    for i in range(len(tri_month)):  # loop trios of months, write list for each
        tri_inds = np.where(np.logical_or(yrmonth==tri_month[i][0],np.logical_or(yrmonth==tri_month[i][1],yrmonth==tri_month[i][2])))[0]
        fout = open('tfpws_in/%s_%s.%s.%s.in' % (p.split('/')[-2],tri_month[i][0],tri_month[i][1],tri_month[i][2]),'w')
        for j in range(len(tri_inds)):
            fpos = corfiles[tri_inds[j]]
            fneg = fpos[:-3]+'neg'
            fout.write('/master-ssd/hmark/%s\n/master-ssd/hmark/%s\n' % (fpos, fneg))

        fout.close()

    # write a list of everything
    fout = open('tfpws_in/%s_full.in' % p.split('/')[-2],'w')
    for j in range(len(corfiles)):
        fpos = corfiles[j]
        fneg = fpos[:-3]+'neg'
        fout.write('/master-ssd/hmark/%s\n/master-ssd/hmark/%s\n' % (fpos, fneg))

    fout.close()

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


# now that we have a list of daily xcor files for each pair of stations, make *lists of lists*
# that are roughly equal in total size for running on the cluster

# UPDATE 01.04.2021: write many more lists of [per_list] lines each because that seems to be
# the limit for what xargs will run before just stopping?

# write out list of #lines-per-file for all the inputs
os.system('wc -l tfpws_in/*in > lines_per_list.dat')

# read in that list
nln = np.loadtxt('lines_per_list.dat',usecols=(0,))
fnames = np.loadtxt('lines_per_list.dat',usecols=(1,),dtype=str)

total = nln[-1]; nln = nln[:-1]; fnames = fnames[:-1]  # get rid of "total" line

total = len(nln)  # instead of lines of lines, just lines

#nlist  = 12; per_list = int(total/nlist)  # approximate # lines per list
per_list = 600; nlist = int(total/per_list) + 1  # assuming not a multiple of 1200 exactly


# temporary addition for AY03/ANMA, to speed things up a bit
i = 0; lsum = 0
fout = open('tellus_lists/list_AA_%i.dat' % i, 'w')
for j in range(len(fnames)):
    sta1 = fnames[j].split('/')[1].split('_')[0]
    sta2 = fnames[j].split('/')[1].split('_')[1]
    middle = fnames[j].split('/')[-1][:-3]
    if sta1 in ['AY03','ANMA'] or sta2 in ['AY03','ANMA']:
        fg.write('%s rm wu=0.5 tls=/master-ssd/hmark/tl_%s.sacn tfpws=/master-ssd/hmark/tf_%s.sacn\n' % (fnames[j],middle,middle))
        lsum += 1
    if lsum == per_list:
        lsum = 0
        i += 1
        fout.close()
        fout = open('tellus_lists/list_AA_%i.dat' % i, 'w')
fout.close()

sys.exit()

# temporary addition to make stacking lists just for specific stations that need restacking
fg = open('tellus_lists/list_GUMN.dat', 'w')
fo = open('tellus_lists/list_OHRS.dat', 'w')
for j in range(len(fnames)):
    sta1 = fnames[j].split('/')[1].split('_')[0]
    sta2 = fnames[j].split('/')[1].split('_')[1]
    middle = fnames[j].split('/')[-1][:-3]
    if sta1 == 'GUMN' or sta2 == 'GUMN':
        fg.write('%s rm wu=0.5 tls=/master-ssd/hmark/tl_%s.sacn tfpws=/master-ssd/hmark/tf_%s.sacn\n' % (fnames[j],middle,middle))
    elif sta1 == 'OHRS' or sta2 == 'OHRS':
        fo.write('%s rm wu=0.5 tls=/master-ssd/hmark/tl_%s.sacn tfpws=/master-ssd/hmark/tf_%s.sacn\n' % (fnames[j],middle,middle))

fg.close()
fo.close()



sys.exit()


j = 0
for i in range(nlist):
    fout = open('tellus_lists/list_%i.dat' % i,'w')

    lsum = 0
    while lsum < per_list and j < len(fnames):
        middle = fnames[j].split('/')[-1][:-3]
        firststa = middle.split('_')[0]
        #if len(firststa) == 4:  # this is NOT a YJ station
        fout.write('%s rm wu=0.5 tls=/master-ssd/hmark/tl_%s.sacn tfpws=/master-ssd/hmark/tf_%s.sacn\n' % (fnames[j],middle,middle))

        lsum += 1
        j += 1

    fout.close()

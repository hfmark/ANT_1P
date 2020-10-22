import numpy as np
from glob import glob
from obspy import read
import os, sys, shutil

####
# move miniseed files for permanent stations downloaded separately into established 1P dirs for
# the same days
####

# list new dirs and get day bits of names
new_dirs = np.array(glob('../seed/new/*'))
new_days = np.array(['.'.join(e.split('/')[-1].split('_')[-1].split('.')[:-1]) for e in new_dirs])

# list old dirs and get day bits of names
old_dirs = np.array(glob('../seed/1P/*'))
old_days = np.array(['.'.join(e.split('/')[-1].split('_')[-1].split('.')[:-1]) for e in old_dirs])

# loop new dirs/DSPA miniseed files
for i,ndir in enumerate(new_dirs):
    # find old dir with same day as new one
    nday = new_days[i]
    if nday in old_days:
        oind = np.where(old_days == nday)[0]
        # list the filename
        fMH = glob(ndir+'/*.mseed')[0]  # For DSPA, this is always just one file long
        print(fMH)
        try:
            st = read(fMH)
            print('read')
        except:
            pass
        print(oind, old_dirs[oind][0])
        print(os.path.join(old_dirs[oind][0],fMH.split('/')[-1]))

#        try:
#            # read in the file, resample to 1Hz, change channel to LHZ and rewrite in 1P dir
#            st = read(fMH)
#            st.merge(method=1,fill_value=0)
#            st.resample(1.)
#            for tr in st:
#                tr.stats.channel = 'L' + tr.stats.channel[1:]
#            ofile = os.path.join(old_dirs[oind],fMH.split('/')[-1])
#            st.write(ofile,encoding=5)
#        except:
#            print('error?')
#            continue
#    else:
#        print('no dir for', nday)

import numpy as np
from glob import glob
import os, sys

####
# move miniseed files for permanent stations downloaded separately into established 1P dirs for
# the same days
####

# list new dirs and get day bits of names
new_dirs = np.array(glob('../seed/new/*'))
new_days = np.array(['.'.join(e.split('/')[-1].split('.')[:-1]) for e in new_dirs])

# list old dirs and get day bits of names
old_dirs = np.array(glob('../seed/1P/*'))
old_days = np.array(['.'.join(e.split('/')[-1].split('.')[:-1]) for e in old_dirs])

# loop new dirs/MG01.mseed files
for i,ndir in enumerate(new_dirs):
    print(i,ndir)
# for each file, find old dir with same day as new one
# move the file there

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
print(new_days)

# list old dirs and get day bits of names

# loop new dirs/MG01.mseed files
# for each file, find old dir with same day as new one
# move the file there

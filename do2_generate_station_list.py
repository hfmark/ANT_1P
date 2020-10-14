import numpy as np
from glob import glob
from obspy import read_inventory
import os, sys

####
# write a station list (name, lat, lon) for seed2cor
# NOTE that this is for the old version of seed2cor; I think the new version wants (name, lon, lat)??
####

# get all the station info
inv = read_inventory('seed/dataless/combined.xml')

# find unique stations (channel doesn't matter; they're all co-located anyway)
ok_sid = np.array(inv.get_contents()['channels'])
_,inds = np.unique(np.array([e.split('.')[1] for e in ok_sid]),return_index=True)
sid = ok_sid[inds]

# open station list file for writing
fs = open('stations.lst','w')

# loop stations and write
for s in sid:
    crd = inv.get_coordinates(s)
    sta = s.split('.')[1]
    fs.write('%s  %.4f  %.4f\n' % (sta,crd['longitude'],crd['latitude']))  # reversed for new s2c

fs.close()

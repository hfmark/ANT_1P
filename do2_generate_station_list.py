import numpy as np
from glob import glob
import os, sys

####
# write a station list (name, lat, lon) for seed2cor
####

# get list of unique stations for 1P/etc and ENAP, along with a representative seed file for each
# rdseed -S -f [filename] > junk
# from rdseed.stations, get name and coords and write to another file

# open station list file for writing
fs = open('stations.lst','w')

# find a single seed file for each unique station
sfiles = glob('/P/hmark/ANT_1P/seed/1P/*/*.seed')
stations = [a.split('/')[-1].split('.')[0] for a in sfiles]
stations,inds = np.unique(stations,return_index=True)

for i in range(len(stations)):
	os.system('rdseed -S -f %s > junk' % sfiles[inds[i]])

	with open('rdseed.stations','r') as f: line = f.readline()

	line = line.split()
	fs.write('%s  %s  %s\n' % (line[0],line[2],line[3]))

fs.close()
sys.exit()






sfiles = glob('/P/hmark/ANT_1P/seed/ENAP/*.seed')
stations = [a.split('/')[-1].split('.')[1] for a in sfiles]
stations,inds = np.unique(stations,return_index=True)

for i in range(len(stations)):
	os.system('rdseed -S -f %s > junk' % sfiles[inds[i]])

	with open('rdseed.stations','r') as f: line = f.readline()

	line = line.split()
	fs.write('%s  %s  %s\n' % (line[0],line[2],line[3]))

fs.close()

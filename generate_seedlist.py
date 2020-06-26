import numpy as np
from datetime import datetime
from glob import glob
import os, sys

####
# write a list of [path yyyy mm dd] for seedfiles to input to seed2cor
# include both 1P/C/C1/G station files and ENAP seed files (once those exist)
####

# open seedlist file for writing
fs = open('seedfiles.lst','w')

# 1P etc
flist = glob('/P/hmark/ANT_1P/seed/1P/*/*.seed')

# for each, get y/m/d, glob seed files, and write
for i in range(len(flist)):
	date = '.'.join(flist[i].split('/')[-2].split('_')[-1].split('.')[:3])
	d = datetime.strptime(date,'%Y.%b.%d')
	if d.year == 2018 and d.month < 3:
		fs.write('%s  %s  %02d  %02d\n' % (flist[i],str(d.year),d.month,d.day))

fs.close()
sys.exit()





# ENAP:
flist = glob('/P/hmark/ANT_1P/seed/ENAP/*.seed')
for i in range(len(flist)):
	date = flist[i].split('/')[-1].split('_')[-1][:10]
	d = datetime.strptime(date,'%Y.%m.%d')
	fs.write('%s  %s  %02d  %02d\n' % (flist[i],str(d.year),d.month,d.day))

fs.close()

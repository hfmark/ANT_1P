import numpy as np
from datetime import datetime
from glob import glob
import os, sys

####
# write a list of [path yyyy mm dd] for seedfiles to input to seed2cor
# dayfiles include 1P/C/C1/G/ENAP all catted together
####

# open seedlist file for writing
fs = open('seedfiles.lst','w')

# Everything!
flist = glob('/P/hmark/ANT_1P/seed/dayfiles/*.seed')

# for each, get y/m/d, glob seed files, and write
for i in range(len(flist)):
	date = '.'.join(flist[i].split('/')[-1].split('_')[-1].split('.')[:3])
	d = datetime.strptime(date,'%Y.%m.%d')
	fs.write('%s  %s  %02d  %02d\n' % (flist[i],str(d.year),d.month,d.day))

fs.close()

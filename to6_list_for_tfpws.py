import numpy as np
import glob as glob
import os, sys

####
# count # of days to be stacked in all the list files, and roughly balance them into lists of 
# list files with ~equal lengths
# write lists of list files along with args for tf_pws0mp input
####

# write out list of #lines-per-file for all the inputs
os.system('wc -l tfpws_in/*in > lines_per_list.dat')

# read in that list
nln = np.loadtxt('lines_per_list.dat',usecols=(0,))
fnames = np.loadtxt('lines_per_list.dat',usecols=(1,),dtype=str)

total = nln[-1]; nln = nln[:-1]; fnames = fnames[:-1]  # get rid of "total" line

nlist  = 8; per_list = int(total/nlist)  # approximate # lines per list

j = 0
for i in range(nlist):
    fout = open('list_%i.dat' % i,'w')

    lsum = 0
    while lsum <= per_list and j < len(fnames):
        middle = fnames[j].split('/')[-1][:-3]
        fout.write('%s rm tls=/master-ssd/hmark/tl_%s.sacn tfpws=/master-ssd/hmark/tf_%s.sacn\n' % (fnames[j],middle,middle))

        lsum += nln[j]
        j += 1

    fout.close()

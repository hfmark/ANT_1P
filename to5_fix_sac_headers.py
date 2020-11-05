import numpy as np
from glob import glob
from obspy import read
import os, sys

####
# the sac files written by tf-pws can't be properly read by the snr code because of some issue with
# the headers
# stupid workaround: read each into python and rewrite (solves the distance problem) after
# putting b=0 and e=3000 into stats.sac (solves the time problem). The rest seems ok.
####

corfiles = np.sort(glob('COR/stacks/tf*sacn'))

for fname in corfiles:
    st = read(fname)
    fname_out = fname[:-5]+'_hdr.sacn'
    st[0].stats.sac['b'] = 0.
    st[0].stats.sac['e'] = 3000.
    st.write(fname_out,'SAC')

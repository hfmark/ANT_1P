import numpy as np
from glob import glob
import os, sys

####
# run aftan on all the separate list files
####

aftan_lists = glob('aftan_*.lst')

for lst in aftan_lists:
    print(lst)
    os.system('../bin/aftan/src/aftan_c_test %s > junk.ftan' % lst)


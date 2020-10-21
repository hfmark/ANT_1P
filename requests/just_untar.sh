#!/bin/bash
#

## untar, and get rid of archives

seeddir=./seed/YJ
#seeddir=./seed/temp

cd $seeddir

for tarfile in  *'.tar.mseed'; do
  tar -xf ${tarfile}
  rm -f ${tarfile}

done

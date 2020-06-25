#!/bin/bash
#

## convert miniseed files from IRIS back into full seed with headers for seed2cor

seeddir=./seed/1P

cd $seeddir


for d in */; do
  cd $d

  for mini in *'.mseed'; do
    rdseed -d -o 5 -f $mini -g '../../dataless/IRISDMC-'${mini%.mseed}'.dataless' >> junk

    mv seed.rdseed ${mini%.mseed}'.seed'
#    rm ${mini}
#    rm ${mini%.mseed}'.dataless'
#    rm ${mini%.mseed}'.meta'
#    rm ${mini%.mseed}'.xml'
  done

  cd ../
done




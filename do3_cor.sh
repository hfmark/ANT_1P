#!/bin/bash

export OMP_NUM_THREADS=12
echo 'remember to enter y to start correlation calculation'
/P/hmark/bin/Seed2Cor/src/Seed2Cor parameters_new.txt 12 > junk_cor

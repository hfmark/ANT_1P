#!/bin/bash

#ls COR/*/*SAC_s > sym_cor.lst
ls COR/*SAC_s > sym_cor.lst

/P/hmark/bin/xcor_snr/spectral_snr_s2c sym_cor.lst


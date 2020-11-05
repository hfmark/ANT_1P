#!/bin/bash

ls COR/stacks/tf*sacn > pws_cor.lst

/P/hmark/bin/xcor_snr/spectral_snr_s2c pws_cor.lst

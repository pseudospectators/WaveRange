#!/bin/bash
#
# Download compressed FluSI regular output data, reconstruct and compress with 1e-3 tolerance using the command files.
#
wget https://osf.io/5kcuq/download --output-document=ux_hit.enc.h5
./wrdec < outmeta
./wrenc < inmeta 

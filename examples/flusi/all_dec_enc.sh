#!/bin/bash
#
# Download compressed FluSI regular output data, reconstruct and compress with 1e-3 tolerance using the command files.
# The contents of the command file outmeta substitute for the following interactive input:
#   Enter compressed data file name []: ux_hit.enc.h5
#   Enter reconstructed file name []: ux_hit.h5
#   Enter file type (0: regular output; 1: backup) [0]: 0
#   Enter output data type (1: float; 2: double) [2]: 2
# The contents of the command file inmeta substitute for the following interactive input:
#   Enter input file name []: ux_hit.h5
#   Enter output file name []: ux_hit_lr.enc.h5
#   Enter file type (0: regular output; 1: backup) [0]: 0
#   Enter base cutoff relative tolerance [1e-16]: 1e-3
#
wget https://osf.io/5kcuq/download --output-document=ux_hit.enc.h5
./wrdec < outmeta
./wrenc < inmeta

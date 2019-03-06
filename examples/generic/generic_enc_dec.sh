#!/bin/bash
#
# Generate a sample Fortran unformatted sequential access binary file 'data.bin' 
# that contains a 32*32*32 double precision array, a 64*64*64 double precision array 
# and 1 single precision variable. Modify 'Makefile' if necessary and type 'make' to build an 
# executable 'create_in_field'. Compress the first array with 1e-6 relative tolerance,
# the second array with 1e-3 relative tolerance, and leave the last single precision variable 
# uncompressed. Finally, reconstruct the data from the compressed format. 
# This example assumes that the record length is 4 bit. No endian conversion is performed.
# The contents of the command file inmeta substitute for the following interactive input:
#   Enter input data file name [data.bin]: data.bin
#   Enter encoded data file name [data.wrb]: data.wrb
#   Enter encoding header file name [data.wrh]: data.wrh
#   Enter file type (0: Fortran sequential w 4-byte recl; 1: Fortran sequential w 8-byte recl; 2: C/C++) [0]: 0
#   Enter endian conversion (0: do not perform; 1: inversion) [0]: 0
#   Enter the number of fields in the file, nf [1]: 3
#   Field number 0
#   Enter input data type (1: float; 2: double) [2]: 2
#   Enter the number of data points in the first dimension, nx [16]: 32
#   Enter the number of data points in the second dimension, ny [16]: 32
#   Enter the number of data points in the third dimension, nz [16]: 32
#   Enter the number of data points in the higher (slowest) dimensions, nh [1]: 1
#   Invert the order of the dimensions? (0: no; 1: yes) [0]: 0
#   Enter compression flag (0: do not compress; 1: compress) [1]: 1
#   Enter base cutoff relative tolerance [1e-16]: 1e-6
#   Field number 1
#   Enter input data type (1: float; 2: double) [2]: 2
#   Enter the number of data points in the first dimension, nx [16]: 64
#   Enter the number of data points in the second dimension, ny [16]: 64
#   Enter the number of data points in the third dimension, nz [16]: 64
#   Enter the number of data points in the higher (slowest) dimensions, nh [1]: 1
#   Invert the order of the dimensions? (0: no; 1: yes) [0]: 0
#   Enter compression flag (0: do not compress; 1: compress) [1]: 1
#   Enter base cutoff relative tolerance [1e-16]: 1e-3
#   Field number 2
#   Enter input data type (1: float; 2: double) [2]: 1
#   Enter the number of data points in the first dimension, nx [16]: 1
#   Enter the number of data points in the second dimension, ny [16]: 1
#   Enter the number of data points in the third dimension, nz [16]: 1
#   Enter the number of data points in the higher (slowest) dimensions, nh [1]: 1
#   Invert the order of the dimensions? (0: no; 1: yes) [0]: 0
#   Enter compression flag (0: do not compress; 1: compress) [1]: 0
# The contents of the command file outmeta substitute for the following interactive input:
#   Enter encoded data file name [data.wrb]: data.wrb
#   Enter encoding header file name [data.wrh]: data.wrh
#   Enter extracted (output) data file name [datarec.bin]: datarec.bin
#   Enter file type (0: Fortran sequential w 4-byte recl; 1: Fortran sequential w 8-byte recl; 2: C/C++) [0]: 0
#   Enter endian conversion (0: do not perform; 1: inversion) [0]: 0
#
make
./create_in_field
./wrenc < inmeta
./wrdec < outmeta


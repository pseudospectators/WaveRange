# WaveRange
WaveRange: CFD data compression utility

I. GENERAL INFORMATION
  
Reference:
doc/cfdproc2017.pdf
Dmitry Kolomenskiy, Ryo Onishi and Hitoshi Uehara "Wavelet-Based Compression of CFD Big Data"
Proceedings of the 31st Computational Fluid Dynamics Symposium, Kyoto, December 12-14, 2017
Paper No. C08-1

This work is supported by the FLAGSHIP2020, MEXT within the priority study4 
(Advancement of meteorological and global environmental predictions utilizing 
observational “Big Data”).
WaveRange is a utility for compression of CFD data, inspired by image processing. It uses wavelet decomposition and subsequent range coding with quantization suitable for floating-point data.

Copyright (C) 2017  Dmitry Kolomenskiy

Copyright (C) 2017  Ryo Onishi

Copyright (C) 2017  JAMSTEC

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

II. COMPILING AND BUILDING

1) Modify the 'config.mk' file. Set the environment variables 'CC' and 'CXX' to point to the desired compilers, set the compiler flags in 'CFLAGS' and 'CXXFLAGS', and the archiver in 'AR'. 
* HDF5 (<https://www.hdfgroup.org/downloads/hdf5/>) and MPI are required for building the FluSI interface. If compiling with HDF5 support, modify the paths in 'HDF_INC' and 'HDF_LIB' to point to the valid library files. It may be necessary to use mpicxx or h5c++. 
* If not using HDF5, select a serial compiler and empty 'HDF_INC' and 'HDF_LIB'.
2) Type 'make' to build the executable files. To only build one of the interfaces, type 'make flusi' or 'make mssg'.
3) Executables will appear in 'bin/' directory. 'flusi/' will contain compression and reconstruction utilities for FluSI output data, 'mssg/' will contain similar utilities for MSSG data. The encoder executable file names contain 'enc', the decoder executable file names contain 'dec'.

III. USAGE

1) Copy the executables into the same directory with the data files. Sample compressed data files can be downloaded from <https://osf.io/pz4n8/>. 
2) If necessary, copy and modify the sample 'inmeta' and 'outmeta' control files that can be found in 'tests/'.
3) Run the utility with the parameters read from the control files or use parameter string.

   Example: ./wrenc < inmeta; ./wrdec < outmeta

IV. BRIEF DESCRIPTION OF THE FILES

* config.mk : build configuration file
* LICENSE.txt : license file
* Makefile : make file
* README.md : this readme file
* doc/cfdproc2017.pdf : brief theoretical introduction
* src/core/defs.h : constant parameter definitions 
* src/core/wrappers.cpp : encoding/decoding subroutines including wavelet transforn and range coding
* src/core/wrappers.h : header for wrappers.cpp
* src/flusi/hdf5_interfaces.cpp : subroutine for handling FluSI HDF5 files
* src/flusi/hdf5_interfaces.h : header for hdf5_interfaces.cpp
* src/flusi/main_dec.cpp : main FluSI decoder program
* src/flusi/main_enc.cpp : main FluSI encoder program
* src/mssg/ctrl_aux.cpp : subroutines for handling MSSG control files
* src/mssg/ctrl_aux.h : header for ctrl_aux.cpp
* src/mssg/mssg_dec.cpp : main MSSG decoder program
* src/mssg/mssg_enc.cpp : main MSSG encoder program
* src/rangecod/Makefile : range coder make file
* src/rangecod/port.h : range coder constant parameters
* src/rangecod/rangecod.c : range coder subroutines
* src/rangecod/rangecod.h : header for rangecod.c
* src/waveletcdf97_3d/Makefile : wavelet transform make file
* src/waveletcdf97_3d/waveletcdf97_3d.c : wavelet transform subroutines
* src/waveletcdf97_3d/waveletcdf97_3d.h : header for waveletcdf97_3d.c
* tests/flusi/all_enc_dec.sh : encode/decode sample bash script for FluSI
* tests/flusi/inmeta : example FluSI encoder parameter
* tests/flusi/outmeta : example FluSI decoder parameters
* tests/flusi/wrdec : symbolic link to WaveRange FluSI decoder executable
* tests/flusi/wrenc : symbolic link to WaveRange FluSI encoder executable
* tests/mssg/regout/all_enc_dec.sh : encode/decode sample bash script for MSSG regular output
* tests/mssg/regout/inmeta : example MSSG regular output encoder parameters
* tests/mssg/regout/outmeta : example MSSG regular output decoder parameters
* tests/mssg/regout/wrmssgdec : symbolic link to WaveRange MSSG decoder executable
* tests/mssg/regout/wrmssgenc : symbolic link to WaveRange MSSG encoder executable
* tests/mssg/separated/all_enc_dec.sh : encode/decode sample bash script for MSSG divided restart files
* tests/mssg/separated/wrmssgdec : symbolic link to WaveRange MSSG decoder executable
* tests/mssg/separated/wrmssgenc : symbolic link to WaveRange MSSG encoder executable
* tests/mssg/unified/all_enc_dec.sh : encode/decode sample bash script for MSSG united restart files
* tests/mssg/unified/wrmssgdec : symbolic link to WaveRange MSSG decoder executable
* tests/mssg/unified/wrmssgenc :symbolic link to WaveRange MSSG encoder executable

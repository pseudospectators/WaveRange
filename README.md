# WaveRange
WaveRange: CFD data compression utility

Copyright (C) 2017  Dmitry Kolomenskiy
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
  
Reference:
doc/cfdproc2017.pdf
Dmitry Kolomenskiy, Ryo Onishi and Hitoshi Uehara "Wavelet-Based Compression of CFD Big Data"
Proceedings of the 31st Computational Fluid Dynamics Symposium, Kyoto, December 12-14, 2017
Paper No. C08-1

This work is supported by the FLAGSHIP2020, MEXT within the priority study4 
(Advancement of meteorological and global environmental predictions utilizing 
observational “Big Data”).

I. GENERAL INFORMATION

WaveRange is a utility for compression of CFD data, inspired by image processing. It uses wavelet decomposition and subsequent range coding with quantization suitable for floating-point data.

II. COMPILING AND BUILDING

1) Modify the 'config.mk' file. MPI and HDF libraries only necessary for FluSI interface. If not using it, select a serial compiler and empty HDF path.
2) Type 'make' to build the executable files. To only build one of the interfaces, type 'make flusi' or 'make mssg'.
3) Executables will appear in 'bin/' directory. 'flusi/' will contain compression and reconstruction utilities for FluSI output data, 'mssg/' will contain similar utilities for MSSG data. The encoder executable file names contain 'enc', the decoder executable file names contain 'dec'.

III. USAGE

1) Copy the executables into the same directory with the data files.
2) If necessary, copy and modify the sample 'inmeta' and 'outmeta' control files that can be found in 'tests/'.
3) Run the utility with the parameters read from the control files or use parameter string.
examples: 
wrenc < inmeta
wrdec < outmeta

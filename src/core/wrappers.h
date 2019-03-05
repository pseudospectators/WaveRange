/*
    wrappers.h : This file is part of WaveRange CFD data compression utility

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
  
    Reference:
    doc/cfdproc2017.pdf
    Dmitry Kolomenskiy, Ryo Onishi and Hitoshi Uehara "Wavelet-Based Compression of CFD Big Data"
    Proceedings of the 31st Computational Fluid Dynamics Symposium, Kyoto, December 12-14, 2017
    Paper No. C08-1

    This work is supported by the FLAGSHIP2020, MEXT within the priority study4 
    (Advancement of meteorological and global environmental predictions utilizing 
    observational “Big Data”).
*/

/* C/C++ interface */
/* Encoding subroutine with wavelet transform and range coding 
    nx : (INPUT) number of elements of the input 3D field in the first (fastest) direction
    ny : (INPUT) number of elements of the input 3D field in the second direction
    nz : (INPUT) number of elements of the input 3D field in the third (slowest) direction
    fld_1d : (INPUT) input 3D field reshaped in a 1D array with the direction x being contiguous (as in Fortran)
    wtflag : (INPUT) wavelet transform flag: 0 if not transforming, 1 if transforming
    mx : (INPUT) number local cutoff subdomains in x, recommended mx=1
    my : (INPUT) number local cutoff subdomains in y, recommended my=1
    mz : (INPUT) number local cutoff subdomains in z, recommended mz=1
    cutoffvec : (INPUT) local cutoff vector, e.g. cutoffvec = new double[1]; cutoffvec[0] = relative_global_tolerance;
    tolabs : (OUTPUT) absolute global tolerance
    midval : (OUTPUT) mid-value, midval = minval+(maxval-minval)/2, where maxval is the maximum and minval is the minimum of fld_1d
    halfspanval : (OUTPUT) half-span, halfspanval = (maxval-minval)/2
    wlev : (OUTPUT) number of wavelet transform levels
    nlay : (OUTPUT) number of bit planes
    ntot_enc : (OUTPUT) total number of elemens of the encoded array data_enc
    deps_vec : (OUTPUT) quantization step size vector, double deps_vec[nlaymax]; where nlaymax is an output of setup_wr
    minval_vec : (OUTPUT) bit plane offset vector, double minval_vec[nlaymax];
    len_enc_vec : (OUTPUT) number of elements in the encoded bit planes, unsigned long int len_enc_vec[nlaymax];
    data_enc : (OUTPUT) range-encoded output data array, defined as, e.g., unsigned char *data_enc = new unsigned char[ntot_enc_max]; where ntot_enc_max is an output of setup_wr */ 
extern "C" void encoding_wrap(int nx, int ny, int nz, double *fld_1d, int wtflag, int mx, int my, int mz, double *cutoffvec, double& tolabs, double& midval, double& halfspanval, unsigned char& wlev, unsigned char& nlay, unsigned long int& ntot_enc, double *deps_vec, double *minval_vec, unsigned long int *len_enc_vec, unsigned char *data_enc);

/* Decoding subroutine with range decoding and inverse wavelet transform 
    nx : (INPUT) number of elements of the input 3D field in the first (fastest) direction
    ny : (INPUT) number of elements of the input 3D field in the second direction
    nz : (INPUT) number of elements of the input 3D field in the third (slowest) direction
    fld_1d : (OUTPUT) output reconstructed 3D field reshaped in a 1D array with the direction x being contiguous (as in Fortran)
    tolabs : (NOT USED)
    midval : (INPUT) mid-value, midval = minval+(maxval-minval)/2, where maxval is the maximum and minval is the minimum of fld_1d
    halfspanval : (INPUT) half-span, halfspanval = (maxval-minval)/2
    wlev : (INPUT) number of wavelet transform levels
    nlay : (INPUT) number of bit planes
    ntot_enc : (INPUT) total number of elemens of the encoded array data_enc
    deps_vec : (INPUT) quantization step size vector, double deps_vec[nlay];
    minval_vec : (INPUT) bit plane offset vector, double minval_vec[nlay];
    len_enc_vec : (INPUT) number of elements in the encoded bit planes, unsigned long int len_enc_vec[nlay];
    data_enc : (INPUT) range-encoded data array, defined as, e.g., unsigned char *data_enc = new unsigned char[ntot_enc]; */ 
extern "C" void decoding_wrap(int nx, int ny, int nz, double *fld_1d, double& tolabs, double& midval, double& halfspanval, unsigned char& wlev, unsigned char& nlay, unsigned long int& ntot_enc, double *deps_vec, double *minval_vec, unsigned long int *len_enc_vec, unsigned char *data_enc);

/* Return the number of bit planes and the required encoded data array size, as needed for memory allocation
    nlaymax : maximum allowed number of bit planes
    ntot_enc_max : maximum allowed total number of elemens of the encoded array data_enc */ 
extern "C" void setup_wr(int nx, int ny, int nz, unsigned char& nlaymax, unsigned long int& ntot_enc_max);

/* Fortran interface */
/* Encoding subroutine with wavelet transform and range coding 
    nx : (INPUT) number of elements of the input 3D field in the first (fastest) direction
    ny : (INPUT) number of elements of the input 3D field in the second direction
    nz : (INPUT) number of elements of the input 3D field in the third (slowest) direction
    fld : (INPUT) output reconstructed 3D field reshaped in a 1D array with the direction x being contiguous (as in Fortran)
    wtflag : (INPUT) wavelet transform flag: 0 if not transforming, 1 if transforming
    tolrel : (INPUT) relative global tolerance
    tolabs : (OUTPUT) absolute global tolerance
    midval : (OUTPUT) mid-value, midval = minval+(maxval-minval)/2, where maxval is the maximum and minval is the minimum of fld
    halfspanval : (OUTPUT) half-span, halfspanval = (maxval-minval)/2
    wlev : (OUTPUT) number of wavelet transform levels
    nlay : (OUTPUT) number of bit planes
    ntot_enc : (OUTPUT) total number of elemens of the encoded array data_enc
    deps_vec : (OUTPUT) quantization step size vector
    minval_vec : (OUTPUT) bit plane offset vector
    len_enc_vec : (OUTPUT) number of elements in the encoded bit planes
    data_enc : (OUTPUT) range-encoded output data array */ 
extern "C" void encoding_wrap_f(int *nx, int *ny, int *nz, double *fld, int *wtflag, double *tolrel, double& tolabs, double& midval, double& halfspanval, unsigned char& wlev, unsigned char& nlay, long int& ntot_enc, double *deps_vec, double *minval_vec, long int *len_enc_vec, unsigned char *data_enc);

/* Decoding subroutine with range decoding and inverse wavelet transform 
    nx : (INPUT) number of elements of the input 3D field in the first (fastest) direction
    ny : (INPUT) number of elements of the input 3D field in the second direction
    nz : (INPUT) number of elements of the input 3D field in the third (slowest) direction
    fld : (OUTPUT) output reconstructed 3D field reshaped in a 1D array with the direction x being contiguous (as in Fortran)
    midval : (INPUT) mid-value, midval = minval+(maxval-minval)/2, where maxval is the maximum and minval is the minimum of fld
    halfspanval : (INPUT) half-span, halfspanval = (maxval-minval)/2
    wlev : (INPUT) number of wavelet transform levels
    nlay : (INPUT) number of bit planes
    ntot_enc : (INPUT) total number of elemens of the encoded array data_enc
    deps_vec : (INPUT) quantization step size vector
    minval_vec : (INPUT) bit plane offset vector
    len_enc_vec : (INPUT) number of elements in the encoded bit planes
    data_enc : (INPUT) range-encoded data array */ 
extern "C" void decoding_wrap_f(int *nx, int *ny, int *nz, double *fld, double& midval, double& halfspanval, unsigned char& wlev, unsigned char& nlay, long int& ntot_enc, double *deps_vec, double *minval_vec, long int *len_enc_vec, unsigned char *data_enc);

/* Return the number of bit planes and the required encoded data array size, as needed for memory allocation
    nx : (INPUT) number of elements of the input 3D field in the first (fastest) direction
    ny : (INPUT) number of elements of the input 3D field in the second direction
    nz : (INPUT) number of elements of the input 3D field in the third (slowest) direction
    nlaymax : (OUTPUT) maximum allowed number of bit planes
    ntot_enc_max : (OUTPUT) maximum allowed total number of elemens of the encoded array data_enc */ 
extern "C" void setup_wr_f(int *nx, int *ny, int *nz, int& nlaymax, long int& ntot_enc_max);

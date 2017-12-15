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

/* Encoding subroutine with wavelet transform and range coding */ 
void encoding_wrap(int nx, int ny, int nz, double *fld_1d, int wtflag, int mx, int my, int mz, double *cutoffvec, double& tolabs, double& midval, double& halfspanval, unsigned char& wlev, unsigned char& nlay, unsigned long int& ntot_enc, double *deps_vec, double *minval_vec, unsigned long int *len_enc_vec, unsigned char *data_enc);
/* Decoding subroutine with range decoding and inverse wavelet transform*/ 
void decoding_wrap(int nx, int ny, int nz, double *fld_1d, double& tolabs, double& midval, double& halfspanval, unsigned char& wlev, unsigned char& nlay, unsigned long int& ntot_enc, double *deps_vec, double *minval_vec, unsigned long int *len_enc_vec, unsigned char *data_enc);
/* Use the range encoder for an array fld_q */ 
void range_encode(unsigned char *fld_q, unsigned long int ntot, unsigned char *enc_q, unsigned long int& len_out_q);
/* Decode the range-encoded data */
void range_decode(unsigned char *enc_q, unsigned long int len_out_q, unsigned char *dec_q, unsigned long int ntot);



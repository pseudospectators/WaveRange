/*
    gen_aux.h : This file is part of WaveRange CFD data compression utility

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

/* Write a field to an unformatted fortran/C/C++ binary file */
void write_field_gen( const char *filename, int idset, int ifiletype, int flag_convertendian, int nbytes, unsigned char *recl, int nx, int ny, int nz, int nh, int idinv, double *fld );
/* Read a field from an unformatted fortran/C/C++ binary file */
void read_field_gen( const char *filename, int idset, int ifiletype, int flag_convertendian, int nbytes, unsigned char *recl, int nx, int ny, int nz, int nh, int idinv, long *btpos, double *fld );
/* Write unsigned char type data set */
void write_field_gen_enc( const char *filename, unsigned char *fld, unsigned long int ntot_enc );
/* Read unsigned char type data set */
void read_field_gen_enc( std::ifstream &inputfile, unsigned char *fld, unsigned long int ntot_enc );
/* Write double type data set */
void write_field_gen_raw( const char *filename, int nbytes, double *fld, unsigned long int ntot );
/* Read double type data set */
void read_field_gen_raw( std::ifstream &inputfile, int nbytes, double *fld, unsigned long int ntot );
/* Write encoding header file */
void write_header_gen_enc( const char *filename, int idset, int *nbytes, unsigned char *recl, int *nx, int *ny, int *nz, int *nh, int *idinv, int *icomp, double *tol_base, double *tolabs, double *midval, double *halfspanval, unsigned char *wlev, unsigned char *nlay, unsigned long int *ntot_enc, double *deps_vec, double *minval_vec, unsigned long int *len_enc_vec );
/* Read encoding header file */
void read_header_gen_enc( std::ifstream &fs, int idset, int *nbytes, unsigned char *recl, int *nx, int *ny, int *nz, int *nh, int *idinv, int *icomp, double *tol_base, double *tolabs, double *midval, double *halfspanval, unsigned char *wlev, unsigned char *nlay, unsigned long int *ntot_enc, double *deps_vec, double *minval_vec, unsigned long int *len_enc_vec );

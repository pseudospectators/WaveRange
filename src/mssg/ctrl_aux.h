/*
    ctrl_aux.h : This file is part of WaveRange CFD data compression utility

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

/* Read restart control file */
void read_control_file( const char *filename, int &nx, int &ny, int &nz, int &nprocx, int &nprocy, char dsettab[NDSMAX][256], int &ndset );
/* Read regular output control file */
void read_control_file_grads( const char *filename, int &nx, int &ny, int &nz, int &nt, double &undef, char *dsetname );
/* Write a field to an unformatted fortran binary file */
void write_field_mssg( const char *filename, int flag_convertendian, int nbytes, int idset, int nx, int ny, int nz, int nxloc, int nyloc, int ixst, int iyst, double *fld );
/* Read a field from an unformatted fortran binary file */
void read_field_mssg( const char *filename, int flag_convertendian, int nbytes, int idset, int nx, int ny, int nz, int nxloc, int nyloc, int ixst, int iyst, double *fld );
/* Write unsigned char type data set */
void write_field_mssg_enc( const char *filename, unsigned char *fld, unsigned long int ntot_enc );
/* Read unsigned char type data set */
void read_field_mssg_enc( std::ifstream &inputfile, unsigned char *fld, unsigned long int ntot_enc );
/* Write encoding header file */
void write_header_mssg_enc( const char *filename, int idset, const char *dsetname, double *tolabs, double *midval, double *halfspanval, unsigned char *wlev, unsigned char *nlay, unsigned long int *ntot_enc, double *deps_vec, double *minval_vec, unsigned long int *len_enc_vec );
/* Read encoding header file */
void read_header_mssg_enc( std::ifstream &fs, int idset, char *dsetname, double *tolabs, double *midval, double *halfspanval, unsigned char *wlev, unsigned char *nlay, unsigned long int *ntot_enc, double *deps_vec, double *minval_vec, unsigned long int *len_enc_vec );

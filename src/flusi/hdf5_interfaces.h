/*
    hdf5_interfaces.h : This file is part of WaveRange CFD data compression utility

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
*/

/* Operator function for finding dataset names */
herr_t op_func ( hid_t loc_id, const char *name, const H5O_info_t *info, void *operator_data );
/* Write double type attribute string */
int write_attrib_dble( const char *filename, const char *dsetname, const char *aname, double *attribute, int nattr );
/* Read double type attribute string */
int read_attrib_dble( const char *filename, const char *dsetname, const char *aname, double *attribute );
/* Write int type attribute string */
int write_attrib_int( const char *filename, const char *dsetname, const char *aname, int *attribute, int nattr );
/* Read int type attribute string */
int read_attrib_int( const char *filename, const char *dsetname, const char *aname, int *attribute );
/* Write encoding attribute string */
int write_attrib_enc( const char *filename, const char *dsetname, double *tolabs, double *midval, double *halfspanval, unsigned char *wlev, unsigned char *nlay, unsigned long int *ntot_enc, double *deps_vec, double *minval_vec, unsigned long int *len_enc_vec );
/* Read encoding attribute string */
int read_attrib_enc( const char *filename, const char *dsetname, double *tolabs, double *midval, double *halfspanval, unsigned char *wlev, unsigned char *nlay, unsigned long int *ntot_enc, double *deps_vec, double *minval_vec, unsigned long int *len_enc_vec );
/* Write double type data set */
int write_field_hdf5( const char *filename, const char *dsetname, double *fld, int nx, int ny, int nz, hid_t outtype );
/* Read double type data set */
int read_field_hdf5( const char *filename, const char *dsetname, double *fld );
/* Write unsigned char type data set */
int write_field_hdf5_enc( const char *filename, const char *dsetname, unsigned char *fld, unsigned long int ntot_enc );
/* Read unsigned char type data set */
int read_field_hdf5_enc( const char *filename, const char *dsetname, unsigned char *fld );




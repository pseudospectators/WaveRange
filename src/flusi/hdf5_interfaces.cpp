/*
    hdf5_interfaces.cpp : This file is part of WaveRange CFD data compression utility

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

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>

#include <iostream>
#include <exception>
#include <memory>

#include "hdf5.h"
#include "../core/defs.h"
#include "hdf5_interfaces.h"

using namespace std;


/* Operator function for finding dataset names */
herr_t op_func ( hid_t loc_id, const char *name, const H5O_info_t *info, void *operator_data )
{
    // File name passed as operator data using type conversion
    char *od = (char *) operator_data;

    // Check if the current object is the root group, and if not print
    // the full path name and type. Modified from HDF5 tutorial.
    if (name[0] == '.')      
      printf ("  File contents:\n");
    else
        switch (info->type) {
            case H5O_TYPE_GROUP:
                printf ("%s  (Group)\n", name);
                break;
            case H5O_TYPE_DATASET:
                printf ("%s  (Dataset)\n", name);
                while ((*od++ = *name++) != '\0');
                break;
            case H5O_TYPE_NAMED_DATATYPE:
                printf ("%s  (Datatype)\n", name);
                break;
            default:
                printf ("%s  (Unknown)\n", name);
        }

    return 0;
}


/* Write double type attribute string */
int write_attrib_dble( const char *filename, const char *dsetname, const char *aname, double *attribute, int nattr )
{

   hid_t file_id, dset_id, attr_id, aspace_id;
   hsize_t dims[1];
   htri_t exists;
   herr_t status;
   int err = 0;

   // Open file
   file_id = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT);

   // Open dataset
   dset_id = H5Dopen2(file_id, dsetname, H5P_DEFAULT);

   // check if attribute exists
   exists = H5Aexists(dset_id, aname);
   if (exists) {
     // Open attribute (it exists already)
     attr_id = H5Aopen(dset_id, aname, H5P_DEFAULT);

     // Get dataspace
     aspace_id = H5Aget_space(attr_id);

   } else {

     // Determine the dataspace identifier aspace_id
     dims[0] = nattr;
     aspace_id = H5Screate_simple(1, dims, NULL);

     // set attr_id, ie create an attribute attached to the object dset_id
     attr_id = H5Acreate2(dset_id,aname,H5T_NATIVE_DOUBLE,aspace_id,H5P_DEFAULT,H5P_DEFAULT);
   }

   // Write the attribute data attribute to the attribute identifier attr_id.
   status = H5Awrite(attr_id, H5T_NATIVE_DOUBLE, attribute);
   if (status<0) err = 1;

   // Close the attribute
   status = H5Aclose(attr_id);
   if (status<0) err = 1;

   // Close the dataset
   status = H5Dclose(dset_id);
   if (status<0) err = 1;

   // Close the file
   status = H5Fclose(file_id);
   if (status<0) err = 1;

   // Exit
   return err;
}


/* Read double type attribute string */
int read_attrib_dble( const char *filename, const char *dsetname, const char *aname, double *attribute )
{
   hid_t faplist_id, file_id, dset_id, attr_id;
   htri_t exists;
   herr_t status;
   int err = 0;

   // Open file
   faplist_id = H5Pcreate (H5P_FILE_ACCESS);
   status = H5Pset_fapl_stdio (faplist_id);
   if (status<0) err = 1;
   file_id = H5Fopen (filename, H5F_ACC_RDONLY, faplist_id);

   // Open dataset
   dset_id = H5Dopen2(file_id, dsetname, H5P_DEFAULT);

   // check if attribute exists
   exists = H5Aexists(dset_id, aname);
   if (exists) {
     // open attribute
     attr_id = H5Aopen(dset_id, aname, H5P_DEFAULT);
     // read attribute data
     status = H5Aread(attr_id, H5T_NATIVE_DOUBLE, attribute);
     if (status<0) err = 1;
     // Close attribute
     status = H5Aclose(attr_id);
     if (status<0) err = 1;
     // return no error
     err = 0;
   } else {
     attribute = 0;
     // Return error
     err = 1;
   }

   // Close the dataset
   status = H5Dclose(dset_id);
   if (status<0) err = 1;

   // Close the file
   status = H5Fclose(file_id);
   if (status<0) err = 1;

   // Exit
   return err;
}


/* Write int type attribute string */
int write_attrib_int( const char *filename, const char *dsetname, const char *aname, int *attribute, int nattr )
{
   hid_t file_id, dset_id, attr_id, aspace_id;
   hsize_t dims[1];
   htri_t exists;
   herr_t status;
   int err = 0;

   // Open file
   file_id = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT);

   // Open dataset
   dset_id = H5Dopen2(file_id, dsetname, H5P_DEFAULT);

   // check if attribute exists
   exists = H5Aexists(dset_id, aname);
   if (exists) {
     // Open attribute (it exists already)
     attr_id = H5Aopen(dset_id, aname, H5P_DEFAULT);
     // Get dataspace
     aspace_id = H5Aget_space(attr_id);
   } else {
     // Determine the dataspace identifier aspace_id
     dims[0] = nattr;
     aspace_id = H5Screate_simple(1, dims, NULL);
     // set attr_id, ie create an attribute attached to the object dset_id
     attr_id = H5Acreate2(dset_id,aname,H5T_NATIVE_INT,aspace_id,H5P_DEFAULT,H5P_DEFAULT);
   }

   // Write the attribute data attribute to the attribute identifier attr_id.
   status = H5Awrite(attr_id, H5T_NATIVE_INT, attribute);
   if (status<0) err = 1;

   // Close the attribute
   status = H5Aclose(attr_id);
   if (status<0) err = 1;

   // Close the dataset
   status = H5Dclose(dset_id);
   if (status<0) err = 1;

   // Close the file
   status = H5Fclose(file_id);
   if (status<0) err = 1;

   // Exit
   return err;
}

/* Read int type attribute string */
int read_attrib_int( const char *filename, const char *dsetname, const char *aname, int *attribute )
{
   hid_t faplist_id, file_id, dset_id, attr_id;
   htri_t exists;
   herr_t status;
   int err = 0;

   // Open file
   faplist_id = H5Pcreate (H5P_FILE_ACCESS);
   status = H5Pset_fapl_stdio (faplist_id);
   if (status<0) err = 1;
   file_id = H5Fopen (filename, H5F_ACC_RDONLY, faplist_id);

   // Open dataset
   dset_id = H5Dopen2(file_id, dsetname, H5P_DEFAULT);

   // Check if attribute exists
   exists = H5Aexists(dset_id, aname);
   if (exists) {
     // Open attribute
     attr_id = H5Aopen(dset_id, aname, H5P_DEFAULT);
     // Read attribute data
     status = H5Aread(attr_id, H5T_NATIVE_INT, attribute);
     if (status<0) err = 1;
     // Close attribute
     status = H5Aclose(attr_id);
     if (status<0) err = 1;
     // Return no error
     err = 0;
   } else {
     attribute = 0;
     // Return error
     err = 1;
   }

   // Close the dataset
   status = H5Dclose(dset_id);
   if (status<0) err = 1;

   // Close the file
   status = H5Fclose(file_id);
   if (status<0) err = 1;

   // Exit
   return err;
}


/* Write encoding attribute string */
int write_attrib_enc( const char *filename, const char *dsetname, double *tolabs, double *midval, double *halfspanval, unsigned char *wlev, unsigned char *nlay, unsigned long int *ntot_enc, double *deps_vec, double *minval_vec, unsigned long int *len_enc_vec )
{
   hid_t file_id, dset_id, attr_id, aspace_id;
   hsize_t dims[1];
   herr_t status;
   int err = 0;
   int cv = CODER_VERSION;

   // Open file
   file_id = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT);

   // Open dataset
   dset_id = H5Dopen2(file_id, dsetname, H5P_DEFAULT);

   // -- Coder version --
   // Determine the dataspace identifier aspace_id
   dims[0] = 1;
   aspace_id = H5Screate_simple(1, dims, NULL);
   // set attr_id, ie create an attribute attached to the object dset_id
   attr_id = H5Acreate2(dset_id,"coder_version",H5T_NATIVE_INT,aspace_id,H5P_DEFAULT,H5P_DEFAULT);
   // Write the attribute data attribute to the attribute identifier attr_id.
   status = H5Awrite(attr_id, H5T_NATIVE_INT, &cv);
   if (status<0) err = 1;
   // Close the attribute
   status = H5Aclose(attr_id);
   if (status<0) err = 1;

   // -- tolabs --
   // Determine the dataspace identifier aspace_id
   dims[0] = 1;
   aspace_id = H5Screate_simple(1, dims, NULL);
   // set attr_id, ie create an attribute attached to the object dset_id
   attr_id = H5Acreate2(dset_id,"tolabs",H5T_NATIVE_DOUBLE,aspace_id,H5P_DEFAULT,H5P_DEFAULT);
   // Write the attribute data attribute to the attribute identifier attr_id.
   status = H5Awrite(attr_id, H5T_NATIVE_DOUBLE, tolabs);
   if (status<0) err = 1;
   // Close the attribute
   status = H5Aclose(attr_id);
   if (status<0) err = 1;

   // -- midval --
   // Determine the dataspace identifier aspace_id
   dims[0] = 1;
   aspace_id = H5Screate_simple(1, dims, NULL);
   // set attr_id, ie create an attribute attached to the object dset_id
   attr_id = H5Acreate2(dset_id,"midval",H5T_NATIVE_DOUBLE,aspace_id,H5P_DEFAULT,H5P_DEFAULT);
   // Write the attribute data attribute to the attribute identifier attr_id.
   status = H5Awrite(attr_id, H5T_NATIVE_DOUBLE, midval);
   if (status<0) err = 1;
   // Close the attribute
   status = H5Aclose(attr_id);
   if (status<0) err = 1;

   // -- halfspanval --
   // Determine the dataspace identifier aspace_id
   dims[0] = 1;
   aspace_id = H5Screate_simple(1, dims, NULL);
   // set attr_id, ie create an attribute attached to the object dset_id
   attr_id = H5Acreate2(dset_id,"halfspanval",H5T_NATIVE_DOUBLE,aspace_id,H5P_DEFAULT,H5P_DEFAULT);
   // Write the attribute data attribute to the attribute identifier attr_id.
   status = H5Awrite(attr_id, H5T_NATIVE_DOUBLE, halfspanval);
   if (status<0) err = 1;
   // Close the attribute
   status = H5Aclose(attr_id);
   if (status<0) err = 1;

   // -- wlev --
   // Determine the dataspace identifier aspace_id
   dims[0] = 1;
   aspace_id = H5Screate_simple(1, dims, NULL);
   // set attr_id, ie create an attribute attached to the object dset_id
   attr_id = H5Acreate2(dset_id,"wlev",H5T_NATIVE_UCHAR,aspace_id,H5P_DEFAULT,H5P_DEFAULT);
   // Write the attribute data attribute to the attribute identifier attr_id.
   status = H5Awrite(attr_id, H5T_NATIVE_UCHAR, wlev);
   if (status<0) err = 1;
   // Close the attribute
   status = H5Aclose(attr_id);
   if (status<0) err = 1;

   // -- nlay --
   // Determine the dataspace identifier aspace_id
   dims[0] = 1;
   aspace_id = H5Screate_simple(1, dims, NULL);
   // set attr_id, ie create an attribute attached to the object dset_id
   attr_id = H5Acreate2(dset_id,"nlay",H5T_NATIVE_UCHAR,aspace_id,H5P_DEFAULT,H5P_DEFAULT);
   // Write the attribute data attribute to the attribute identifier attr_id.
   status = H5Awrite(attr_id, H5T_NATIVE_UCHAR, nlay);
   if (status<0) err = 1;
   // Close the attribute
   status = H5Aclose(attr_id);
   if (status<0) err = 1;

   // -- ntot_enc --
   // Determine the dataspace identifier aspace_id
   dims[0] = 1;
   aspace_id = H5Screate_simple(1, dims, NULL);
   // set attr_id, ie create an attribute attached to the object dset_id
   attr_id = H5Acreate2(dset_id,"ntot_enc",H5T_NATIVE_ULONG,aspace_id,H5P_DEFAULT,H5P_DEFAULT);
   // Write the attribute data attribute to the attribute identifier attr_id.
   status = H5Awrite(attr_id, H5T_NATIVE_ULONG, ntot_enc);
   if (status<0) err = 1;
   // Close the attribute
   status = H5Aclose(attr_id);
   if (status<0) err = 1;

   // Only save the following attributes if the encoded data is non-trivial
   if (ntot_enc > 0)
     { 
       // -- deps_vec --
       // Determine the dataspace identifier aspace_id
       dims[0] = int(*nlay);
       aspace_id = H5Screate_simple(1, dims, NULL);
       // set attr_id, ie create an attribute attached to the object dset_id
       attr_id = H5Acreate2(dset_id,"deps_vec",H5T_NATIVE_DOUBLE,aspace_id,H5P_DEFAULT,H5P_DEFAULT);
       // Write the attribute data attribute to the attribute identifier attr_id.
       status = H5Awrite(attr_id, H5T_NATIVE_DOUBLE, deps_vec);
       if (status<0) err = 1;
       // Close the attribute
       status = H5Aclose(attr_id);
       if (status<0) err = 1;

       // -- minval_vec --
       // Determine the dataspace identifier aspace_id
       dims[0] = int(*nlay);
       aspace_id = H5Screate_simple(1, dims, NULL);
       // set attr_id, ie create an attribute attached to the object dset_id
       attr_id = H5Acreate2(dset_id,"minval_vec",H5T_NATIVE_DOUBLE,aspace_id,H5P_DEFAULT,H5P_DEFAULT);
       // Write the attribute data attribute to the attribute identifier attr_id.
       status = H5Awrite(attr_id, H5T_NATIVE_DOUBLE, minval_vec);
       if (status<0) err = 1;
       // Close the attribute
       status = H5Aclose(attr_id);
       if (status<0) err = 1;

       // -- len_enc_vec --
       // Determine the dataspace identifier aspace_id
       dims[0] = int(*nlay);
       aspace_id = H5Screate_simple(1, dims, NULL);
       // set attr_id, ie create an attribute attached to the object dset_id
       attr_id = H5Acreate2(dset_id,"len_enc_vec",H5T_NATIVE_ULONG,aspace_id,H5P_DEFAULT,H5P_DEFAULT);
       // Write the attribute data attribute to the attribute identifier attr_id.
       status = H5Awrite(attr_id, H5T_NATIVE_ULONG, len_enc_vec);
       if (status<0) err = 1;
       // Close the attribute
       status = H5Aclose(attr_id);
       if (status<0) err = 1;
     }

   // Close the dataset
   status = H5Dclose(dset_id);
   if (status<0) err = 1;

   // Close the file
   status = H5Fclose(file_id);
   if (status<0) err = 1;

   // Exit
   return err;
}


/* Read encoding attribute string */
int read_attrib_enc( const char *filename, const char *dsetname, double *tolabs, double *midval, double *halfspanval, unsigned char *wlev, unsigned char *nlay, unsigned long int *ntot_enc, double *deps_vec, double *minval_vec, unsigned long int *len_enc_vec )
{
   hid_t faplist_id, file_id, dset_id, attr_id;
   htri_t exists;
   herr_t status;
   int err = 0;

   // Open file
   faplist_id = H5Pcreate (H5P_FILE_ACCESS);
   status = H5Pset_fapl_stdio (faplist_id);
   file_id = H5Fopen (filename, H5F_ACC_RDONLY, faplist_id);

   // Open dataset
   dset_id = H5Dopen2(file_id, dsetname, H5P_DEFAULT);

   // -- tolabs --
   // check if attribute exists
   exists = H5Aexists(dset_id, "tolabs");
   if (exists) {
     // open attribute
     attr_id = H5Aopen(dset_id, "tolabs", H5P_DEFAULT);
     // read attribute data
     status = H5Aread(attr_id, H5T_NATIVE_DOUBLE, tolabs);
     if (status<0) err = 1;
     // Close attribute
     status = H5Aclose(attr_id);
     if (status<0) err = 1;
     // return no error
     err = 0;
   } else {
     tolabs = 0;
     // Return error
     err = 1;
   }

   // -- midval --
   // check if attribute exists
   exists = H5Aexists(dset_id, "midval");
   if (exists) {
     // open attribute
     attr_id = H5Aopen(dset_id, "midval", H5P_DEFAULT);
     // read attribute data
     status = H5Aread(attr_id, H5T_NATIVE_DOUBLE, midval);
     if (status<0) err = 1;
     // Close attribute
     status = H5Aclose(attr_id);
     if (status<0) err = 1;
     // return no error
     err = 0;
   } else {
     midval = 0;
     // Return error
     err = 1;
   }

   // -- halfspanval --
   // check if attribute exists
   exists = H5Aexists(dset_id, "halfspanval");
   if (exists) {
     // open attribute
     attr_id = H5Aopen(dset_id, "halfspanval", H5P_DEFAULT);
     // read attribute data
     status = H5Aread(attr_id, H5T_NATIVE_DOUBLE, halfspanval);
     if (status<0) err = 1;
     // Close attribute
     status = H5Aclose(attr_id);
     if (status<0) err = 1;
     // return no error
     err = 0;
   } else {
     halfspanval = 0;
     // Return error
     err = 1;
   }

   // -- wlev --
   // check if attribute exists
   exists = H5Aexists(dset_id, "wlev");
   if (exists) {
     // open attribute
     attr_id = H5Aopen(dset_id, "wlev", H5P_DEFAULT);
     // read attribute data
     status = H5Aread(attr_id, H5T_NATIVE_UCHAR, wlev);
     if (status<0) err = 1;
     // Close attribute
     status = H5Aclose(attr_id);
     if (status<0) err = 1;
     // return no error
     err = 0;
   } else {
     wlev = 0;
     // Return error
     err = 1;
   }

   // -- nlay --
   // check if attribute exists
   exists = H5Aexists(dset_id, "nlay");
   if (exists) {
     // open attribute
     attr_id = H5Aopen(dset_id, "nlay", H5P_DEFAULT);
     // read attribute data
     status = H5Aread(attr_id, H5T_NATIVE_UCHAR, nlay);
     if (status<0) err = 1;
     // Close attribute
     status = H5Aclose(attr_id);
     if (status<0) err = 1;
     // return no error
     err = 0;
   } else {
     nlay = 0;
     // Return error
     err = 1;
   }

   // -- ntot_enc --
   // check if attribute exists
   exists = H5Aexists(dset_id, "ntot_enc");
   if (exists) {
     // open attribute
     attr_id = H5Aopen(dset_id, "ntot_enc", H5P_DEFAULT);
     // read attribute data
     status = H5Aread(attr_id, H5T_NATIVE_ULONG, ntot_enc);
     if (status<0) err = 1;
     // Close attribute
     status = H5Aclose(attr_id);
     if (status<0) err = 1;
     // return no error
     err = 0;
   } else {
     ntot_enc = 0;
     // Return error
     err = 1;
   }

   // Only read the following attributes if the encoded data is non-trivial
   if (ntot_enc > 0)
     { 

       // -- deps_vec --
       // check if attribute exists
       exists = H5Aexists(dset_id, "deps_vec");
       if (exists) {
         // open attribute
         attr_id = H5Aopen(dset_id, "deps_vec", H5P_DEFAULT);
         // read attribute data
         status = H5Aread(attr_id, H5T_NATIVE_DOUBLE, deps_vec);
         if (status<0) err = 1;
         // Close attribute
         status = H5Aclose(attr_id); 
         if (status<0) err = 1;
         // return no error
         err = 0;
       } else {
         deps_vec = 0;
         // Return error
         err = 1;
       }

       // -- minval_vec --
       // check if attribute exists
       exists = H5Aexists(dset_id, "minval_vec");
       if (exists) {
         // open attribute
         attr_id = H5Aopen(dset_id, "minval_vec", H5P_DEFAULT);
         // read attribute data
         status = H5Aread(attr_id, H5T_NATIVE_DOUBLE, minval_vec);
         if (status<0) err = 1;
         // Close attribute
         status = H5Aclose(attr_id); 
         if (status<0) err = 1;
         // return no error
         err = 0;
       } else {
         minval_vec = 0;
         // Return error
         err = 1;
       }

       // -- len_enc_vec --
       // check if attribute exists
       exists = H5Aexists(dset_id, "len_enc_vec");
       if (exists) {
         // open attribute
         attr_id = H5Aopen(dset_id, "len_enc_vec", H5P_DEFAULT);
         // read attribute data
         status = H5Aread(attr_id, H5T_NATIVE_ULONG, len_enc_vec);
         if (status<0) err = 1;
         // Close attribute
         status = H5Aclose(attr_id);
         if (status<0) err = 1;
         // return no error
         err = 0;
       } else {
         len_enc_vec = 0;
         // Return error
         err = 1;
       }
     }

   // Close the dataset
   status = H5Dclose(dset_id);
   if (status<0) err = 1;

   // Close the file
   status = H5Fclose(file_id);
   if (status<0) err = 1;

   // Exit
   return err;
}


/* Write float or double type data set */
int write_field_hdf5( const char *filename, const char *dsetname, double *fld, int nx, int ny, int nz, hid_t outtype )
{
   hid_t file_id, dset_id, dataspace_id, plist_id;
   hsize_t dims[3];
   //hsize_t chdims[3]; // Chunking is not necessary, removed for better portability
   herr_t status;
   int err = 0;

   // Open file
   file_id = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT);

   // Create the data space for the dataset
   dims[0] = hsize_t(nz); 
   dims[1] = hsize_t(ny);
   dims[2] = hsize_t(nx);  
   dataspace_id = H5Screate_simple(3, dims, NULL);

   // Create dataset.
   plist_id = H5Pcreate(H5P_DATASET_CREATE);
   //chdims[0] = hsize_t(min(nz,4096)); 
   //chdims[1] = hsize_t(min(ny,4096));
   //chdims[2] = hsize_t(min(nx,4096));  
   //status = H5Pset_chunk(plist_id, 3, chdims);
   status = H5Pset_fill_value(plist_id, outtype, 0);
   status = H5Pset_alloc_time(plist_id, H5D_ALLOC_TIME_EARLY);

   // Create the dataset
   dset_id = H5Dcreate2(file_id, dsetname, outtype, dataspace_id, 
                          H5P_DEFAULT, plist_id, H5P_DEFAULT);

   // Close the property list
   status = H5Pclose(plist_id);
   if (status<0) err = 1;

   // Write dataset
   status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, dataspace_id, H5P_DEFAULT, fld);
   if (status<0) err = 1;

   // Close the dataset
   status = H5Dclose(dset_id);
   if (status<0) err = 1;

   // Terminate access to the data space
   status = H5Sclose(dataspace_id);
   if (status<0) err = 1;

   // Close the file
   status = H5Fclose(file_id);
   if (status<0) err = 1;

   // Exit
   return err;
}


/* Read double type data set */
int read_field_hdf5( const char *filename, const char *dsetname, double *fld )
{
   hid_t faplist_id, file_id, dset_id;
   herr_t status;
   int err = 0;

   // Open file
   faplist_id = H5Pcreate (H5P_FILE_ACCESS);
   status = H5Pset_fapl_stdio (faplist_id);
   file_id = H5Fopen (filename, H5F_ACC_RDONLY, faplist_id);

   // Open dataset
   dset_id = H5Dopen2(file_id, dsetname, H5P_DEFAULT);

   // Read dataset
   status = H5Dread(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, fld);
   if (status<0) err = 1;

   // Close the dataset
   status = H5Dclose(dset_id);
   if (status<0) err = 1;

   // Close the file
   status = H5Fclose(file_id);
   if (status<0) err = 1;

   // Exit
   return err;
}


/* Write unsigned char type data set */
int write_field_hdf5_enc( const char *filename, const char *dsetname, unsigned char *fld, unsigned long int ntot_enc )
{
   hid_t file_id, dset_id, dataspace_id;
   hsize_t dims[1];
   herr_t status;
   int err = 0;

   // Open file
   file_id = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT);

   // Create the data space for the dataset
   dims[0] = hsize_t(ntot_enc); 
   dataspace_id = H5Screate_simple(1, dims, NULL);

   // Create the dataset
   dset_id = H5Dcreate2(file_id, dsetname, H5T_NATIVE_UCHAR, dataspace_id, 
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

   // Write dataset
   status = H5Dwrite(dset_id, H5T_NATIVE_UCHAR, dataspace_id, dataspace_id, H5P_DEFAULT, fld);
   if (status<0) err = 1;

   // Close the dataset
   status = H5Dclose(dset_id);
   if (status<0) err = 1;

   // Terminate access to the data space
   status = H5Sclose(dataspace_id);
   if (status<0) err = 1;

   // Close the file
   status = H5Fclose(file_id);
   if (status<0) err = 1;

   // Exit
   return err;
}


/* Read unsigned char type data set */
int read_field_hdf5_enc( const char *filename, const char *dsetname, unsigned char *fld )
{
   hid_t faplist_id, file_id, dset_id;
   herr_t status;
   int err = 0;

   // Open file
   faplist_id = H5Pcreate (H5P_FILE_ACCESS);
   status = H5Pset_fapl_stdio (faplist_id);
   file_id = H5Fopen (filename, H5F_ACC_RDONLY, faplist_id);

   // Open dataset
   dset_id = H5Dopen2(file_id, dsetname, H5P_DEFAULT);

   // Read dataset
   status = H5Dread(dset_id, H5T_NATIVE_UCHAR, H5S_ALL, H5S_ALL, H5P_DEFAULT, fld);
   if (status<0) err = 1;

   // Close the dataset
   status = H5Dclose(dset_id);
   if (status<0) err = 1;

   // Close the file
   status = H5Fclose(file_id);
   if (status<0) err = 1;

   // Exit
   return err;
}

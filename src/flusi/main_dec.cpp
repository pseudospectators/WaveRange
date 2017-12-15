/*
    main_dec.cpp : This file is part of WaveRange CFD data compression utility

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
#include <fstream>
#include <exception>
#include <memory>
#include <string>
#include <sstream>

#include "hdf5.h"
#include "../core/defs.h"
#include "../core/wrappers.h"
#include "hdf5_interfaces.h"

using namespace std;


// Main code for decoding
int main( int argc, char *argv[] )
{
    // Dataset dimensions
    int nx;
    int ny;
    int nz;

    // More variable declarations
    double *fld_1d_rec;
    double tolabs;
    double midval, halfspanval;
    unsigned char wlev, nlay;
    unsigned long int ntot_enc;
    double deps_vec[NLAYMAX];
    double minval_vec[NLAYMAX];
    unsigned long int len_enc_vec[NLAYMAX];

    // I/O variable declarations
    int ifiletype = 0, iouttype = 1;
    hid_t faplist_id, file_id, outtype;
    herr_t status;
    htri_t exists;
    int err = 0;

    // parameters variable definitons
    string in_name;
    string out_name;
    string bar, bar2;

    cout << "usage: ./wrdec compressed_000.h5 decompressed_000.h5 TYPE PRECISION\n";
    cout << "where TYPE=(0: regular output; 1: backup) and PRECISION=(1:single 2:double)\n";
    cout << "interactive mode if not enough arguments are passed.\n";

    if ( argc == 5 )
    {
      cout << "automatic mode.";
      in_name = argv[1];
      out_name = argv[2];
      bar = argv[3];
      bar2 = argv[4];
    }
    else
    {
      /* Prepare for decoding */
      // Read file name
      cout << "Enter compressed data file name []: ";
      getline (cin,in_name);
      cout << "Enter reconstructed file name []: ";
      getline (cin,out_name);
      cout << "Enter file type (0: regular output; 1: backup) [0]: ";
      getline (cin,bar);
      cout << "Enter output data type (1: float; 2: double) [2]: ";
      getline (cin,bar2);
    }

    stringstream(bar) >> ifiletype;
    stringstream(bar2) >> iouttype;


    switch (iouttype) {
      case 1: {  outtype = H5T_NATIVE_FLOAT; break; }
      case 2: {  outtype = H5T_NATIVE_DOUBLE; break; }
      default: ;
    }

    // Print out metadata
    cout << endl << "=== Decoding parameters ===" << endl;
    cout << "Input file name: " << in_name << endl;
    cout << "Output file name: " << out_name << endl;
    cout << "File type (0: regular output; 1: backup): " << ifiletype << endl;
    cout << "Output data type (1: float; 2: double): " << iouttype << endl;

    // Start HDF5
    status = H5open();

    // Create a new output file using default properties
    file_id = H5Fcreate(out_name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    // Close the file
    status = H5Fclose(file_id);

    /* Encoding */
    switch (ifiletype) {

      // For regular output files
      case 0:
        {
          /* Start reading from file */
          // Dataset name
	  char dsetname[13];

          // Attributes
          double time, nu, epsi, domain_size[3];
          int nx, ny, nz, nxyz[3];

          // Open file
          faplist_id = H5Pcreate (H5P_FILE_ACCESS);
          status = H5Pset_fapl_stdio (faplist_id);
          file_id = H5Fopen (in_name.c_str(), H5F_ACC_RDONLY, faplist_id);

          // Find dataset name
	  status = H5Ovisit (file_id,H5_INDEX_NAME,H5_ITER_NATIVE,op_func,(void *) &dsetname);

          // Close file
          status = H5Fclose(file_id);


          /* Read dataset */
          // Read attributes
          err = read_attrib_dble( in_name.c_str(), dsetname, "time", &time );
          err = read_attrib_dble( in_name.c_str(), dsetname, "viscosity", &nu );
          err = read_attrib_dble( in_name.c_str(), dsetname, "epsi", &epsi );
          err = read_attrib_dble( in_name.c_str(), dsetname, "domain_size", domain_size );
          err = read_attrib_int( in_name.c_str(), dsetname, "nxyz", nxyz );

          // Size of the dataset
          nx = nxyz[0];
          ny = nxyz[1];
          nz = nxyz[2];
          unsigned long int ntot = (unsigned long int)(nx)*(unsigned long int)(ny)*(unsigned long int)(nz);

          // Diagnostics
          cout << " dset=" << dsetname << " nx=" << nx  << " ny=" << ny  << " nz=" << nz << endl;

          // Read coding attributes
          err = read_attrib_enc( in_name.c_str(), dsetname, &tolabs, &midval, &halfspanval, &wlev, &nlay, &ntot_enc, deps_vec, minval_vec, len_enc_vec );

          // Allocate encoded data array
          unsigned char *data_enc = new unsigned char[ntot_enc];

          // Read data if the compressed data set is non-trivial
          if (ntot_enc)
            err = read_field_hdf5_enc( in_name.c_str(), dsetname, data_enc );

          // Allocate data for wavelet reconstruction
          fld_1d_rec = new double[ntot];

          /* Do decoding */
          // Apply decoding routine
          decoding_wrap(nx,ny,nz,fld_1d_rec,tolabs,midval,halfspanval,wlev,nlay,ntot_enc,deps_vec,minval_vec,len_enc_vec,data_enc);
          cout << "  decode: fld_1d_rec[0]=" << fld_1d_rec[0] << " fld_1d_rec[last]=" << fld_1d_rec[ntot-1UL] << endl;

          // Calculate min and max
          double maxval = fld_1d_rec[0];
          double minval = fld_1d_rec[0];
          for(unsigned long int j1 = 0; j1 < ntot; j1++)
            {
              if (fld_1d_rec[j1] < minval) minval = fld_1d_rec[j1];
              if (fld_1d_rec[j1] > maxval) maxval = fld_1d_rec[j1];
            }

          // Print min and max
          cout << "        min=" << minval << " max=" << maxval << endl;

          // Deallocate memory
          delete [] data_enc;

          /* Write reconstructed data to a file */
          // Write data
          err = write_field_hdf5( out_name.c_str(), dsetname, fld_1d_rec, nx, ny, nz, outtype );

          // Write attributes
          err = write_attrib_dble( out_name.c_str(), dsetname, "time", &time, 1 );
          err = write_attrib_dble( out_name.c_str(), dsetname, "viscosity", &nu, 1 );
          err = write_attrib_dble( out_name.c_str(), dsetname, "epsi", &epsi, 1 );
          err = write_attrib_dble( out_name.c_str(), dsetname, "domain_size", domain_size, 3 );
          err = write_attrib_int( out_name.c_str(), dsetname, "nxyz", nxyz, 3 );

          // Deallocate memory
          delete [] fld_1d_rec;

          break;
        }

      // For backup files
      case 1:
        {
          // Dataset parameters
          double attributes[8];
          const int ndset = 50;
          int dset_exist[ndset];
          char dsettab[ndset][13] = {"ux","uy","uz","nlkx0","nlky0","nlkz0","nlkx1","nlky1","nlkz1",
          "bx","by","bz","bnlkx0","bnlky0","bnlkz0","bnlkx1","bnlky1","bnlkz1",
          "scalar1","scalar1_nlk0","scalar1_nlk1",
          "scalar2","scalar2_nlk0","scalar2_nlk1",
          "scalar3","scalar3_nlk0","scalar3_nlk1",
          "scalar4","scalar4_nlk0","scalar4_nlk1",
          "scalar5","scalar5_nlk0","scalar5_nlk1",
          "scalar6","scalar6_nlk0","scalar6_nlk1",
          "scalar7","scalar7_nlk0","scalar7_nlk1",
          "scalar8","scalar8_nlk0","scalar8_nlk1",
          "scalar9","scalar9_nlk0","scalar9_nlk1",
          "uavgx","uavgy","uavgz","ekinavg","Z_avg"};

           // Loop for all datasets
          for (int j=0; j<ndset; j++)
          {
            /* Read compressed data from a file */
            // Check if dataset exists
            faplist_id = H5Pcreate (H5P_FILE_ACCESS);
            status = H5Pset_fapl_stdio (faplist_id);
            file_id = H5Fopen (in_name.c_str(), H5F_ACC_RDONLY, faplist_id);
            exists = H5Lexists( file_id, dsettab[j], H5P_DEFAULT );
            status = H5Fclose(file_id);
            if (exists) dset_exist[j] = 1; else dset_exist[j] = 0;

            // Proceed only with those datasets that exist
            if (dset_exist[j] == 1)
            {
              /* Read dataset */
              // Read attributes
              // (/time,dt1,dt0,dble(n1),dble(it),dble(nx),dble(ny),dble(nz)/)
              err = read_attrib_dble( in_name.c_str(), dsettab[j], "bckp", attributes );

              // Size of the dataset
              nx = int(attributes[5]);
              ny = int(attributes[6]);
              nz = int(attributes[7]);
              unsigned long int ntot = (unsigned long int)(nx)*(unsigned long int)(ny)*(unsigned long int)(nz);

              // Diagnostics
              cout << " dset=" << dsettab[j] << " nx=" << nx  << " ny=" << ny  << " nz=" << nz << endl;

              // Read coding attributes
              err = read_attrib_enc( in_name.c_str(), dsettab[j], &tolabs, &midval, &halfspanval, &wlev, &nlay, &ntot_enc, deps_vec, minval_vec, len_enc_vec );

              // Allocate encoded data array
              unsigned char *data_enc = new unsigned char[ntot_enc];

              // Read data if the compressed data set is non-trivial
              if (ntot_enc)
                err = read_field_hdf5_enc( in_name.c_str(), dsettab[j], data_enc );

              // Allocate data for wavelet reconstruction
              fld_1d_rec = new double[ntot];

              /* Do decoding */
              // Apply decoding routine
              decoding_wrap(nx,ny,nz,fld_1d_rec,tolabs,midval,halfspanval,wlev,nlay,ntot_enc,deps_vec,minval_vec,len_enc_vec,data_enc);
              cout << "  decode: fld_1d_rec[0]=" << fld_1d_rec[0] << " fld_1d_rec[last]=" << fld_1d_rec[ntot-1UL] << endl;

              // Calculate min and max
              double maxval = fld_1d_rec[0];
              double minval = fld_1d_rec[0];
              for(unsigned long int j1 = 0; j1 < ntot; j1++)
                {
                  if (fld_1d_rec[j1] < minval) minval = fld_1d_rec[j1];
                  if (fld_1d_rec[j1] > maxval) maxval = fld_1d_rec[j1];
                }

              // Echo min and max
              cout << "        min=" << minval << " max=" << maxval << endl;

              // Deallocate memory
              delete [] data_enc;

              /* Write reconstructed data to a file */
              // Write data
              err = write_field_hdf5( out_name.c_str(), dsettab[j], fld_1d_rec, nx, ny, nz, outtype );

              // Write attributes
              err = write_attrib_dble( out_name.c_str(), dsettab[j], "bckp", attributes, 8 );

              // Deallocate memory
              delete [] fld_1d_rec;
            }
          }
          break;
        }

      default:
        // Display error message
        cout << "Error: unknown file type" << endl;

    }

    // Stop HDF5
    status = H5close();

    return 0;
}

/*
    main_enc.cpp : This file is part of WaveRange CFD data compression utility

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


// Main code for encoding
int main( int argc, char *argv[] )
{
    // Dataset dimensions
    int nx = 0;
    int ny = 0;
    int nz = 0;
    unsigned long int ntot = 0;

    // Base tolerance, applied as relative to max(fabs(fld_1d))
    double tol_base = 1e-16;

    // Local cutoff tolerance parameters
    int mx, my, mz;
    unsigned int mtot;
    double *cutoffvec;

    // Data variable declarations
    double *fld_1d;
    double tolabs;
    double midval, halfspanval;
    unsigned char wlev, nlay;
    unsigned long int ntot_enc;
    double deps_vec[NLAYMAX];
    double minval_vec[NLAYMAX];
    unsigned long int len_enc_vec[NLAYMAX];

    // I/O variable declarations
    int ifiletype = 0;
    hid_t faplist_id, file_id;
    herr_t status;
    htri_t exists;
    int err = 0;

    string in_name;
    string out_name;
    string bar;
    string bar2;

    cout << "usage: ./wrenc original_000.h5 compressed_000.h5 TYPE TOLERANCE\n";
    cout << "where TYPE=(0: regular output; 1: backup) and TOLERANCE=(e.g. 1.0e-5)\n";
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
      /* Prepare for encoding */
      // Read metadata
      cout << "Enter input file name []: ";
      getline (cin,in_name);
      cout << "Enter output file name []: ";
      getline (cin,out_name);
      cout << "Enter file type (0: regular output; 1: backup) [0]: ";
      getline (cin,bar);
      cout << "Enter base cutoff relative tolerance [1e-16]: ";
      getline (cin,bar2);
    }

    stringstream(bar) >> ifiletype;
    stringstream(bar2) >> tol_base;

    // Print out metadata
    cout << endl << "=== Compression parameters ===" << endl;
    cout << "Input file name: " << in_name << endl;
    cout << "Output file name: " << out_name << endl;
    cout << "File type (0: regular output; 1: backup): " << ifiletype << endl;
    cout << "Base cutoff relative tolerance: " << tol_base << endl;

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
          /* Define uniform cutoff */
          mx = 1;
          my = 1;
          mz = 1;
          mtot = 1;
          cutoffvec = new double[1];
          cutoffvec[0] = tol_base;

          /* Start reading from file */
          // Dataset name
	  char dsetname[13];

          // Attributes
          double time, nu, epsi, domain_size[3];
          int nxyz[3];

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
          ntot = (unsigned long int)(nx)*(unsigned long int)(ny)*(unsigned long int)(nz);

          // Diagnostics
          cout << " dset=" << dsetname << " nx=" << nx  << " ny=" << ny  << " nz=" << nz << endl;

          // Allocate array
          fld_1d = new double[ntot];

          // Read data
          err = read_field_hdf5( in_name.c_str(), dsetname, fld_1d );
          cout << "  read: fld_1d[0]=" << fld_1d[0] << " fld_1d[last]=" << fld_1d[ntot-1UL] << endl;

          // Calculate min and max
          double minval = fld_1d[0];
          double maxval = fld_1d[0];
          for(unsigned long int j1 = 0; j1 < ntot; j1++)
            {
              minval = fmin(minval,fld_1d[j1]);
              maxval = fmax(maxval,fld_1d[j1]);
            }

	  // Print min and max
          cout << "        min=" << minval << " max=" << maxval << endl;

          // Allocate encoded data array (will be stored in a file)
          // Encoded array may be longer than the original
          unsigned char *data_enc = new unsigned char[SAFETY_BUFFER_FACTOR*NLAYMAX*(ntot<1024UL?1024UL:ntot)];

          /* Do encoding */
          // Apply encoding routine
          encoding_wrap(nx,ny,nz,fld_1d,1,mx,my,mz,cutoffvec,tolabs,midval,halfspanval,wlev,nlay,ntot_enc,deps_vec,minval_vec,len_enc_vec,data_enc);

          // Print efficient global cutoff
          cout << "        tolabs=" << tolabs << endl;

          // Deallocate memory
          delete [] fld_1d;

          /* Write compressed data to a file */
          // Write data if the compressed data set is non-trivial
          err = write_field_hdf5_enc( out_name.c_str(), dsetname, data_enc, ntot_enc );

          // Write attributes
          err = write_attrib_dble( out_name.c_str(), dsetname, "time", &time, 1 );
          err = write_attrib_dble( out_name.c_str(), dsetname, "viscosity", &nu, 1 );
          err = write_attrib_dble( out_name.c_str(), dsetname, "epsi", &epsi, 1 );
          err = write_attrib_dble( out_name.c_str(), dsetname, "domain_size", domain_size, 3 );
          err = write_attrib_int( out_name.c_str(), dsetname, "nxyz", nxyz, 3 );

          // Write coding attributes
          err = write_attrib_enc( out_name.c_str(), dsetname, &tolabs, &midval, &halfspanval, &wlev, &nlay, &ntot_enc, deps_vec, minval_vec, len_enc_vec );

          // Deallocate memory
          delete [] data_enc;

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

          /* Determine local cutoff */
          // Separate branches for uniform and non-uniform cutoff
          if (UNIFORM_CUTOFF)
          {
            // Define uniform cutoff
            mx = 1;
            my = 1;
            mz = 1;
            mtot = 1;
            cutoffvec = new double[1];
            cutoffvec[0] = tol_base;
          }
          else
          {
            // Read velocity attributes
            err = read_attrib_dble( in_name.c_str(), "ux", "bckp", attributes );

            // Size of the velocity dataset
            nx = int(attributes[5]);
            ny = int(attributes[6]);
            nz = int(attributes[7]);
            ntot = (unsigned long int)(nx)*(unsigned long int)(ny)*(unsigned long int)(nz);

            // Allocate full size array
            fld_1d = new double[ntot];

            // Desired subdomain block size
            const int bxyz = DS_BLOCK;

            // Number of blocks
            mx = nx/bxyz;
            my = ny/bxyz;
            mz = nz/bxyz;
            mtot = mx*my*mz;

            // Allocate cutoff array
            cutoffvec = new double[mtot];

            // Define and allocate downscaled data structures
            double *um[3];
            double *wabs;
            int *nm;
            for (int j=0; j<3; j++) um[j] = new double[mtot];
            nm = new int[mtot];

            // Loop for the three velocity components
            for (int j=0; j<3; j++)
            {
              // Read velocity component data
              err = read_field_hdf5( in_name.c_str(), dsettab[j], fld_1d );

              // Initialize downscaled dataset
              for (unsigned int k=0; k<mtot; k++)
              {
                um[j][k] = 0;
                nm[k] = 0;
              }

              // Fill in the downscaled data
              for (unsigned long int i=0; i<ntot; i++)
              {
                // Cartesian coordinates of the local precision cutoff block
                int kx = int(double(i%nx)/double(nx)*double(mx));
                int ky = int(double((i/nx)%ny)/double(ny)*double(my));
                int kz = int(double(i/nx/ny)/double(nz)*double(mz));

                // 1D index
                int k = kx+mx*ky+mx*my*kz;

                // Sum up the values
                um[j][k] += fld_1d[i];
                nm[k] += 1;
              }

              // Calculate the averages
              for (unsigned int k=0; k<mtot; k++) um[j][k] /= nm[k];
            }

            // Deallocate memory
            delete [] fld_1d;
            delete [] nm;

            // Allocate scaled vorticity data
            wabs = new double[mtot];

            // Compute scaled vorticity and its maximum
            double wabsmax = 0;
            for (int kz=0; kz<mz; kz++)
              for (int ky=0; ky<my; ky++)
                for (int kx=0; kx<mx; kx++)
                  {
                    double wx = (um[2][kx+mx*((ky+1)%my)+mx*my*kz]-um[2][kx+mx*((ky-1)%my)+mx*my*kz]) -
                                (um[1][kx+mx*ky+mx*my*(kz+1)%mz]-um[1][kx+mx*ky+mx*my*((kz-1)%mz)]);
                    double wy = (um[0][kx+mx*ky+mx*my*(kz+1)%mz]-um[0][kx+mx*ky+mx*my*((kz-1)%mz)]) -
                                (um[2][((kx+1)%mx)+mx*ky+mx*my*kz]-um[2][((kx-1)%mx)+mx*ky+mx*my*kz]);
                    double wz = (um[1][((kx+1)%mx)+mx*ky+mx*my*kz]-um[1][((kx-1)%mx)+mx*ky+mx*my*kz]) -
                                (um[0][kx+mx*((ky+1)%my)+mx*my*kz]-um[0][kx+mx*((ky-1)%my)+mx*my*kz]);
                    wabs[kx+mx*ky+mx*my*kz] = sqrt(wx*wx+wy*wy+wz*wz);
                    if (wabs[kx+mx*ky+mx*my*kz] > wabsmax) wabsmax = wabs[kx+mx*ky+mx*my*kz];
                  }

            // Print out local vorticity values
            // cout << "Downscaled velocity values" << endl;
            // for (unsigned int k=0; k<mtot; k++) cout << " k=" << k << " wabs_norm=" << wabs[k]/wabsmax << " um1=" << um[0][k] << " um2=" << um[1][k] << " um3=" << um[2][k] << endl;

            // Calculate local cutoff for each block
            for (unsigned int k=0; k<mtot; k++)
            //  if (wabs[k] > 1.0/256.0*wabsmax) cutoffvec[k] = tol_base; else cutoffvec[k] = tol_base*128.0;
              if (wabs[k] > 1.0/128.0*wabsmax) cutoffvec[k] = tol_base; else cutoffvec[k] = tol_base*16.0;

            // Deallocate
            for (int j=0; j<3; j++) delete [] um[j];
            delete [] wabs;

            // Print out local cutoff values
            // cout << "Local cutoff values" << endl;
            // for (unsigned int k=0; k<mtot; k++) cout << " k=" << k << " tol=" << cutoffvec[k] << endl;
          }

          /* Encoding: loop for all datasets */
          for (int j=0; j<ndset; j++)
          {
            /* Read the input data from file */
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
              ntot = (unsigned long int)(nx)*(unsigned long int)(ny)*(unsigned long int)(nz);

              // Diagnostics
              cout << " dset=" << dsettab[j] << " nx=" << nx  << " ny=" << ny  << " nz=" << nz << endl;

              // Allocate array
              fld_1d = new double[ntot];

              // Read data
              err = read_field_hdf5( in_name.c_str(), dsettab[j], fld_1d );
              cout << "  read: fld_1d[0]=" << fld_1d[0] << " fld_1d[last]=" << fld_1d[ntot-1UL] << endl;

              // Calculate min and max
              double minval = fld_1d[0];
              double maxval = fld_1d[0];
              for(unsigned long int j1 = 0; j1 < ntot; j1++)
                {
                  minval = fmin(minval,fld_1d[j1]);
                  maxval = fmax(maxval,fld_1d[j1]);
                }

              // Print min and max
              cout << "        min=" << minval << " max=" << maxval << endl;

              // Allocate encoded data array (will be stored in a file)
              // Encoded array may be longer than the original
              unsigned char *data_enc = new unsigned char[SAFETY_BUFFER_FACTOR*NLAYMAX*(ntot<1024UL?1024UL:ntot)];

              /* Do encoding */
              // Apply encoding routine
              encoding_wrap(nx,ny,nz,fld_1d,1,mx,my,mz,cutoffvec,tolabs,midval,halfspanval,wlev,nlay,ntot_enc,deps_vec,minval_vec,len_enc_vec,data_enc);

              // Print efficient global cutoff
              cout << "        tolabs=" << tolabs << endl;

              // Deallocate memory
              delete [] fld_1d;

              /* Write compressed data to a file */
              // Write data if the compressed data set is non-trivial
              if (ntot_enc > 0)
                err = write_field_hdf5_enc( out_name.c_str(), dsettab[j], data_enc, ntot_enc );

              // Write attributes
              err = write_attrib_dble( out_name.c_str(), dsettab[j], "bckp", attributes, 8 );

              // Write coding attributes
              err = write_attrib_enc( out_name.c_str(), dsettab[j], &tolabs, &midval, &halfspanval, &wlev, &nlay, &ntot_enc, deps_vec, minval_vec, len_enc_vec );

              // Deallocate memory
              delete [] data_enc;
            }
          }
          break;
        }

      default:
        // Display error message
        cout << "Error: unknown file type" << endl;
    }

    // Deallocate memory
    delete [] cutoffvec;

    // Stop HDF5
    status = H5close();

    // Display a message on exit
    cout << "=== End of compression ===\n";

    return 0;
}

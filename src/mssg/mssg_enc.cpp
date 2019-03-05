/*
    mssg_enc.cpp : This file is part of WaveRange CFD data compression utility

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
#include <sstream>
#include <iomanip>
#include <limits>
#include <cassert>

#include "../core/defs.h"
#include "../core/wrappers.h"
#include "ctrl_aux.h"

using namespace std;


// Main code for encoding
int main( int argc, char *argv[] )
{
    // Dataset dimensions
    int nx = 0;
    int ny = 0;
    int nz = 0;
    unsigned long int ntot = 0;

    // Number of subdomains
    int nprocx, nprocy;

    // This proc id
    int thisproc = 0;

    // Local field size
    int nxloc, nyloc;

    // Base tolerance, applied as relative to max(abs(fld_1d))
    double tol_base = 1e-16;

    // Floating point input file endian conversion flag (0: do not convert; 1: convert)
    int flag_convertendian = 1;

    // Floating point input file precision (4: single; 8: double)
    int nbytes;

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
    int ifiletype = 0, iintype = 2;
    int err = 0;
    string prefix_name, ext_name = ".enc", in_name, out_name, header_name, control_name;
    string bar, bar2, bar3, bar4, bar5;

    // Interactive mode help string
    cout << "usage: ./wrmssgenc FILE_NAME_PREFIX ENCODED_NAME_EXT TYPE PRECISION ENDIANFLIP TOLERANCE PROCID\n";
    cout << "where TYPE=(0: regular output; 1: backup united; 2: backup divided), PRECISION=(1:single; 2:double), ENDIANFLIP=(0:no; 1:yes), TOLERANCE=(e.g. 1.0e-16) and PROCID=(this proc id)\n";
    cout << "interactive mode if not enough arguments are passed.\n";

    /* Prepare for encoding */
    if ( argc == 8 )
    {
      // Read metadata from parameter string
      cout << "automatic mode.";
      prefix_name = argv[1];
      ext_name = argv[2];
      bar = argv[3];
      bar2 = argv[4];
      bar3 = argv[5];
      bar4 = argv[6];
      bar5 = argv[7];
    }
    else
    {
      // Read metadata
      cout << "Enter data file name prefix []: ";
      getline (cin,prefix_name);
      cout << "Enter encoded file extension name [.enc]: ";
      getline (cin,ext_name);
      cout << "Enter file type (0: regular output; 1: backup merged; 2: backup separated) [0]: ";
      getline (cin,bar);
      cout << "Enter input data type (1: float; 2: double) [2]: ";
      getline (cin,bar2);
      cout << "Enter endian conversion (0: do not perform; 1: inversion) [1]: ";
      getline (cin,bar3);
      cout << "Enter base cutoff relative tolerance [1e-16]: ";
      getline (cin,bar4);
      cout << "Enter id of this proc [0]: ";
      getline (cin,bar5);

    }
    stringstream(bar) >> ifiletype;
    stringstream(bar2) >> iintype;
    if (iintype==1) nbytes = 4; else nbytes = 8;
    stringstream(bar3) >> flag_convertendian;
    stringstream(bar4) >> tol_base;
    stringstream(bar5) >> thisproc;

    // Print out metadata
    cout << endl << "=== Compression parameters ===" << endl;
    cout << "Data file name prefix: " << prefix_name << endl;
    cout << "Encoded file extension name: " << ext_name << endl;
    cout << "File type (0: regular output; 1: backup merged; 2: backup separated): " << ifiletype << endl;
    cout << "Input files contain " << nbytes << "-byte floating point data" << endl;
    if (flag_convertendian) cout << "Convert big endian to little endian or vice versa" << endl;
    cout << "Base cutoff relative tolerance: " << tol_base << endl;
    cout << "This proc id: " << thisproc << endl;

    /* Encoding */
    switch (ifiletype) {

      // For regular output files
      case 0:
        {
          // Data set name and file path
          char dsetname[256];

          // Number of fields in a file
          int nt;

          // Read control file
          control_name = prefix_name + ".ctl";
          read_control_file_grads(control_name.c_str(),nx,ny,nz,nt,dsetname);

          // Size of the dataset
          ntot = (unsigned long int)(nx)*(unsigned long int)(ny)*(unsigned long int)(nz);

          // Define uniform cutoff
          mx = 1;
          my = 1;
          mz = 1;
          mtot = 1;
          cutoffvec = new double[1];
          cutoffvec[0] = tol_base;

          // Diagnostics
          cout << " dset=" << dsetname << " nx=" << nx  << " ny=" << ny << " nz=" << nz << " nt=" << nt << endl;

          // Output file names
          header_name = prefix_name + "_h"+ ext_name;
          out_name = prefix_name + "_f"+ ext_name;

          // Create header file
          fstream fheader;
          fheader.open(header_name.c_str(), fstream::out | fstream::trunc);
          assert(fheader.is_open());
          fheader << " ===== Header file for compressed MSSG regular output data =====" << endl;
          fheader << " Coder version: " << CODER_VERSION << endl;
          fheader << " File name prefix: " << prefix_name.c_str() << endl;
          fheader << " Encoded file extension name: " << ext_name << endl;
          fheader << " File type (0: regular output; 1: backup merged; 2: backup separated): " << ifiletype << endl;
          fheader << " Input files contained " << nbytes << "-byte floating point data" << endl;
          if (flag_convertendian) fheader << " Converted big endian to little endian or vice versa" << endl;
          fheader << " Base cutoff relative tolerance: " << tol_base << endl;
          fheader.close();

          // Create a new encoded data file. Overwrite if exists
          ofstream foutput;
          foutput.open(out_name.c_str(), ios::binary|ios::out|ios::trunc);
          assert(foutput.is_open());
          foutput.close();

          // Loop for all time instants in the dataset
          for (int it=0; it<nt; it++)
            {
              // Print field number on the screen
              cout << "Field number it=" << it << endl;

              // Allocate array
              fld_1d = new double[ntot];

              /* Read data */
              read_field_mssg(dsetname,flag_convertendian,nbytes,it,nx,ny,nz,nx,ny,0,0,fld_1d);
              cout << "  read: fld_1d[0]=" << fld_1d[0] << " fld_1d[last]=" << fld_1d[ntot-1UL] << endl;

              // Calculate min and max
              double maxval = fld_1d[0];
              double minval = fld_1d[0];
              for(unsigned long int j1 = 0; j1 < ntot; j1++)
                {
                  if (fld_1d[j1] < minval) minval = fld_1d[j1];
                  if (fld_1d[j1] > maxval) maxval = fld_1d[j1];
                }

              // Print min and max
              cout << "        min=" << minval << " max=" << maxval << endl;

              // Detect masking and encode it
              if (minval < MSSG_MASK_THRESHOLD) 
                {
                  // Compute the mean value used for padding
                  double fld_pad = 0;
                  int jmean = 0;
                  for(unsigned long int j1 = 0; j1 < ntot; j1++)
                    if (fld_1d[j1] >= MSSG_MASK_THRESHOLD) 
                      {
                        fld_pad += fld_1d[j1];
                        jmean++;
                      }
                  fld_pad /= jmean;

                  // Allocate array
                  double *mask_1d = new double[ntot];

                  // Separate the mask and the field
                  for(unsigned long int j1 = 0; j1 < ntot; j1++)
                    if (fld_1d[j1] < MSSG_MASK_THRESHOLD) 
                      {
                        fld_1d[j1] = fld_pad;
                        mask_1d[j1] = minval;
                      }
                    else mask_1d[j1] = 0;

                  // Print masking info
                  cout << " Masking detected, padding with fld_pad=" << fld_pad << ", mask min=" << minval << endl;

                  // Allocate encoded data array (will be stored in a file)
                  // Encoded array may be longer than the original
                  unsigned char *data_enc = new unsigned char[SAFETY_BUFFER_FACTOR*NLAYMAX*(ntot<1024UL?1024UL:ntot)];

                  // Tolerance of the mask function encoding
                  double *cutoffvecmask = new double[1];
                  cutoffvecmask[0] = MSSG_MASK_TOLREL; 

                  // Apply encoding routine to the mask
                  encoding_wrap(nx,ny,nz,mask_1d,0,mx,my,mz,cutoffvecmask,tolabs,midval,halfspanval,wlev,nlay,ntot_enc,deps_vec,minval_vec,len_enc_vec,data_enc);

                  // Deallocate memory
                  delete [] mask_1d;
                  delete [] cutoffvecmask;

                  // Write compressed mask to a file
                  // Append the header file with coding attributes
                  write_header_mssg_enc(header_name.c_str(),it,"mask",&tolabs,&midval,&halfspanval,&wlev,&nlay,&ntot_enc,deps_vec,minval_vec,len_enc_vec);
  
                  // Write data if the compressed data set is non-trivial
                  if (ntot_enc > 0)
                      write_field_mssg_enc(out_name.c_str(),data_enc,ntot_enc);

                  // Text output
                  cout << " Mask done, encoding the main field..." << endl;

                  // Deallocate memory
                  delete [] data_enc;
                }

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
              // Append the header file with coding attributes
              write_header_mssg_enc(header_name.c_str(),it,dsetname,&tolabs,&midval,&halfspanval,&wlev,&nlay,&ntot_enc,deps_vec,minval_vec,len_enc_vec);
  
              // Write data if the compressed data set is non-trivial
              if (ntot_enc > 0)
                  write_field_mssg_enc(out_name.c_str(),data_enc,ntot_enc);

              // Deallocate memory
              delete [] data_enc;
            }

          break;
        }

      // For backup files
      case 1: // Merge binary output in one file
      case 2: // Separate output for each subdomain 
        {
          // Dataset parameters
          int ndset;
          char dsettab[NDSMAX][256];

          // Read control file
          control_name = prefix_name + ".nmlst";
          read_control_file(control_name.c_str(),nx,ny,nz,nprocx,nprocy,dsettab,ndset);

          // Calculate domain decomposition parameters
          nxloc = nx / nprocx;
          nyloc = ny / nprocy;

          // Print out the control file data
          cout << endl << "=== Parameters read from control file ===" << endl;
          cout << "nx(=nlg+i_over*2) = " << nx << "; ny(=npg+j_over*2) = " << ny << "; nr(=nz) = " << nz << "; dim_size(=nprocx,nprocy) = " << nprocx << ", " << nprocy << "; ndset = " << ndset << endl;
          for (int j = 0; j<ndset; j++) cout << "record number = " << j+1 << "; field = " << dsettab[j] << endl;

          // Size of the velocity dataset
          if (ifiletype == 1) 
            ntot = (unsigned long int)(nx)*(unsigned long int)(ny)*(unsigned long int)(nz);
          else
            ntot = (unsigned long int)(nxloc)*(unsigned long int)(nyloc)*(unsigned long int)(nz);

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
            // Display error message and stop
            cout << "Local cutoff for MSSG restart not implemented" << endl;
            throw std::exception();
          }

          // Proc label
          stringstream lbl;
          lbl << setw(MSSG_FILE_DIG) << setfill('0') << thisproc;

          // Output file names
          if (ifiletype == 1) 
            {
              header_name = prefix_name + "_h"+ ext_name;
              out_name = prefix_name + "_f"+ ext_name;
            }
          else
            {
              header_name = prefix_name + "_h"+ lbl.str() + ext_name;
              out_name = prefix_name + "_f"+ lbl.str() + ext_name;
            }
          //cout << "Output file names: " << out_name << " " << header_name << endl;

          // Create header file
          fstream fheader;
          fheader.open(header_name.c_str(), fstream::out | fstream::trunc);
          assert(fheader.is_open());
          fheader << " ===== Header file for compressed MSSG restart data =====" << endl;
          fheader << " Coder version: " << CODER_VERSION << endl;
          fheader << " File name prefix: " << prefix_name.c_str() << endl;
          fheader << " Encoded file extension name: " << ext_name << endl;
          fheader << " File type (0: regular output; 1: backup merged; 2: backup separated): " << ifiletype << endl;
          fheader << " Input files contained " << nbytes << "-byte floating point data" << endl;
          if (flag_convertendian) fheader << " Converted big endian to little endian or vice versa" << endl;
                             else fheader << " Did not perform endian conversion" << endl;
          fheader << " Base cutoff relative tolerance: " << tol_base << endl;

          // Write the first ('time') dataset in the header file
          fheader << " -----" << endl;
          fheader << "1" << endl;
          fheader << " Data set name = " << dsettab[0] << endl;

          /* Read time from the first dataset */
          fld_1d = new double[nxloc*nyloc*nz];
          in_name = prefix_name + ".p_" + lbl.str();
          read_field_mssg(in_name.c_str(),flag_convertendian,nbytes,0,nxloc,nyloc,nz,nxloc,nyloc,0,0,fld_1d);

          // Write first MSSG_TIME_REC_LEN elements of time record
          const int time_rec_len = MSSG_TIME_REC_LEN;
          fheader << " first " << time_rec_len << " elements of time record" << endl;
          for (int j=0; j<time_rec_len; j++) 
            fheader << setprecision(numeric_limits<long double>::digits10 + 1) << fld_1d[j] << " ";
          fheader << endl;

          // Close header file and deallocate memory
          fheader.close();
          delete [] fld_1d;

          // Create a new encoded data file. Overwrite if exists
          ofstream foutput;
          foutput.open(out_name.c_str(), ios::binary|ios::out|ios::trunc);
          assert(foutput.is_open());
          foutput.close();

          /* Encoding: loop for all datasets */
          for (int idset=1; idset<ndset; idset++)
          {
            // Allocate array
            fld_1d = new double[ntot];

            /* Read the input data */
            // Read from all files or just one file, depending on ifiletype
            if (ifiletype == 1)
              {   
                // Loop for all subdomains
                for (int iprocy=0; iprocy<nprocy; iprocy++)
                  for (int iprocx=0; iprocx<nprocx; iprocx++)
                    {
                      // 1D array index
                      int iproc = iprocx + nprocx*iprocy;

                      // Calculate local start indexes
                      int ixst = iprocx*nxloc;
                      int iyst = iprocy*nyloc;

                      // Open file
                      stringstream lbliproc;
                      lbliproc << setw(MSSG_FILE_DIG) << setfill('0') << iproc;
                      in_name = prefix_name + ".p_" + lbliproc.str();
                      //cout << "Reading from " << in_name << endl;

                      // Read data in the subdomain
                      read_field_mssg(in_name.c_str(),flag_convertendian,nbytes,idset,nx,ny,nz,nxloc,nyloc,ixst,iyst,fld_1d);
                    }

                // Diagnostics
                cout << " dset=" << dsettab[idset] << " nx=" << nx << " ny=" << ny << " nz=" << nz << endl;
                cout << "  read: fld_1d[0]=" << fld_1d[0] << " fld_1d[last]=" << fld_1d[ntot-1UL] << endl;
              }
            else
              {
                // Read data in the local domain
                read_field_mssg(in_name.c_str(),flag_convertendian,nbytes,idset,nxloc,nyloc,nz,nxloc,nyloc,0,0,fld_1d);

                // Diagnostics
                cout << " dset=" << dsettab[idset] << " nxloc=" << nxloc << " nyloc=" << nyloc << " nz=" << nz << endl;
                cout << "  read: fld_1d[0]=" << fld_1d[0] << " fld_1d[last]=" << fld_1d[ntot-1UL] << endl;
              }

            // Calculate min and max
            double maxval = fld_1d[0];
            double minval = fld_1d[0];
            for(unsigned long int j1 = 0; j1 < ntot; j1++)
              {
                if (fld_1d[j1] < minval) minval = fld_1d[j1];
                if (fld_1d[j1] > maxval) maxval = fld_1d[j1];
              }

            // Print min and max
            cout << "        min=" << minval << " max=" << maxval << endl;

            // Allocate encoded data array (will be stored in a file)
            // Encoded array may be longer than the original
            unsigned char *data_enc = new unsigned char[SAFETY_BUFFER_FACTOR*NLAYMAX*(ntot<1024UL?1024UL:ntot)];

            /* Do encoding */
            // Apply encoding routine
            if (ifiletype == 1)
                encoding_wrap(nx,ny,nz,fld_1d,1,mx,my,mz,cutoffvec,tolabs,midval,halfspanval,wlev,nlay,ntot_enc,deps_vec,minval_vec,len_enc_vec,data_enc);
            else
                encoding_wrap(nxloc,nyloc,nz,fld_1d,1,mx,my,mz,cutoffvec,tolabs,midval,halfspanval,wlev,nlay,ntot_enc,deps_vec,minval_vec,len_enc_vec,data_enc);

            // Print efficient global cutoff
            cout << "        tolabs=" << tolabs << endl;

            // Deallocate memory
            delete [] fld_1d;

            /* Write compressed data to a file */
            // Append the header file with coding attributes
            write_header_mssg_enc(header_name.c_str(),idset,dsettab[idset],&tolabs,&midval,&halfspanval,&wlev,&nlay,&ntot_enc,deps_vec,minval_vec,len_enc_vec);

            // Write data if the compressed data set is non-trivial
            if (ntot_enc > 0)
                write_field_mssg_enc(out_name.c_str(),data_enc,ntot_enc);

            // Deallocate memory
            delete [] data_enc;
          }
          break;
        }

      default:
        // Display error message
        cout << "Error: unknown file type" << endl;
    }

    // Deallocate memory
    delete [] cutoffvec;

    return 0;
}

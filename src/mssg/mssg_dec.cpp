/*
    mssg_dec.cpp : This file is part of WaveRange CFD data compression utility

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
#include <string.h>
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


// Main code for decoding
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

    // Floating point input file endian conversion flag (0: do not convert; 1: convert)
    int flag_convertendian;

    // Floating point input file precision (4: single; 8: double)
    int nbytes;

    // Masking parameter value
    double undef;

    // More variable declarations
    double *fld_1d_rec;
    double *mask_1d_rec;
    double tolabs;
    double midval, halfspanval;
    unsigned char wlev, nlay;
    unsigned long int ntot_enc;
    double deps_vec[NLAYMAX];
    double minval_vec[NLAYMAX];
    unsigned long int len_enc_vec[NLAYMAX];

    // I/O variable declarations
    int ifiletype = 0, iouttype = 1;
    int err = 0;

    // Parameters variable definitons
    string in_prefix_name, ext_name, out_prefix_name, in_name, out_name, header_name, control_name;
    string bar, bar2, bar3, bar4;

    cout << "usage: ./wrmssgdec ENCODED_NAME_PREFIX ENCODED_NAME_EXT EXTRACTED_NAME_PREFIX TYPE PRECISION ENDIANFLIP PROCID\n";
    cout << "where TYPE=(0: regular output; 1: backup united; 2: backup divided), PRECISION=(1:single; 2:double), ENDIANFLIP=(0:no; 1:yes) and PROCID=(this proc id)\n";
    cout << "interactive mode if not enough arguments are passed.\n";

    /* Prepare for decoding */
    if ( argc == 8 )
    {
      // Read metadata from parameter string
      cout << "automatic mode.";
      in_prefix_name = argv[1];
      ext_name = argv[2];
      out_prefix_name = argv[3];
      bar = argv[4];
      bar2 = argv[5];
      bar3 = argv[6];
      bar4 = argv[7];
    }
    else
    {
      // Read metadata
      cout << "Enter encoded data file name prefix []: ";
      getline (cin,in_prefix_name);
      cout << "Enter encoded data file extension name [.enc]: ";
      getline (cin,ext_name);
      cout << "Enter extracted data file name prefix []: ";
      getline (cin,out_prefix_name);
      cout << "Enter file type (0: regular output; 1: backup merged; 2: backup separated) [0]: ";
      getline (cin,bar);
      cout << "Enter extracted data type (1: float; 2: double) [2]: ";
      getline (cin,bar2);
      cout << "Enter endian conversion (0: do not perform; 1: inversion) [1]: ";
      getline (cin,bar3);
      cout << "Enter id of this proc [0]: ";
      getline (cin,bar4);
    }
    stringstream(bar) >> ifiletype;
    stringstream(bar2) >> iouttype;
    if (iouttype==1) nbytes = 4; else nbytes = 8;
    stringstream(bar3) >> flag_convertendian;
    stringstream(bar4) >> thisproc;

    // Print out metadata
    cout << endl << "=== Decoding parameters ===" << endl;
    cout << "Encoded file name prefix: " << in_prefix_name << endl;
    cout << "Encoded file extension name: " << ext_name << endl;
    cout << "Extracted file name prefix: " << out_prefix_name << endl;
    cout << "File type (0: regular output; 1: backup merged; 2: backup separated): " << ifiletype << endl;
    cout << "Output files contain " << nbytes << "-byte floating point data" << endl;
    if (flag_convertendian) cout << "Convert big endian to little endian or vice versa" << endl;
    cout << "This proc id: " << thisproc << endl;

    /* Decoding */
    switch (ifiletype) {

      // For regular output files
      case 0:
        {
          /* Start reading from file */
          // Dataset name
          char dsetname[256];

          // Number of fields in a file
          int nt;

          // Read control file
          control_name = in_prefix_name + ".ctl";
          read_control_file_grads(control_name.c_str(),nx,ny,nz,nt,undef,dsetname);

          // Print out the control file data
          cout << endl << "=== Parameters read from control file ===" << endl;
          cout << " dset=" << dsetname << " nx=" << nx  << " ny=" << ny  << " nz=" << nz << " nt=" << nt << " undef=" << undef << endl;

          // Size of the dataset
          ntot = (unsigned long int)(nx)*(unsigned long int)(ny)*(unsigned long int)(nz);

          // Copy the control file with a new name
          if (strcmp(in_prefix_name.c_str(),out_prefix_name.c_str()) != 0 )
            {
              string copy_control_name = out_prefix_name + ".ctl";
              ifstream ictrl;
              ofstream octrl;
              ictrl.open(control_name.c_str(), fstream::in);
              assert(ictrl.is_open());
              octrl.open(copy_control_name.c_str(), fstream::out|fstream::trunc);
              assert(octrl.is_open());
              char c;
              while (ictrl.get(c)) octrl << c;
              ictrl.close();
              octrl.close();  
            }

          // Output file name
          out_name = out_prefix_name + ".grd";

          // Proc label
          stringstream lbl;
          lbl << setw(MSSG_FILE_DIG) << setfill('0') << thisproc;

          // Input file names
          header_name = in_prefix_name + "_h"+ ext_name;
          in_name = in_prefix_name + "_f"+ ext_name;
          //cout << "Input file names: " << in_name << " " << header_name << endl;

          // Open header file
          ifstream fheader;
          fheader.open(header_name.c_str(), fstream::in);
          assert(fheader.is_open());

          // Skip first 8 lines from the header file
          string str; 
          for (int j=0; j<9; j++) getline(fheader, str);

          // Open encoded data file name
          ifstream finput;
          finput.open(in_name.c_str(), ios::binary|ios::in);
          assert(finput.is_open());

          // Loop for all time instants in the dataset
          for (int it=0; it<nt; it++)
            {
              // Allocate array
              fld_1d_rec = new double[ntot];


              // Read from the header file with coding attributes
              char dsetnamehdr[256];
              read_header_mssg_enc(fheader,it,dsetnamehdr,&tolabs,&midval,&halfspanval,&wlev,&nlay,&ntot_enc,deps_vec,minval_vec,len_enc_vec);
              //cout << " dsetnamehdr=" << dsetnamehdr << " it=" << it << endl;

              // Initialize array
              for (unsigned long int j=0; j<ntot; j++) fld_1d_rec[j] = midval;
 
              // Mask presence flag and middle value
              int mask_flag = 0;
              double mask_midval = 0;

              // Check if this field is a mask
              if (strcmp(dsetnamehdr,"mask") == 0)
                {
                  // Set mask presence flag to true
                  mask_flag = 1;

                  // Allocate array
                  mask_1d_rec = new double[ntot];

                  // Initialize array
                  for (unsigned long int j=0; j<ntot; j++) mask_1d_rec[j] = midval;

                  // Reconstruct mask
                  if (ntot_enc > 0)
                    {
                      // Allocate encoded data array
                      unsigned char *data_enc = new unsigned char[ntot_enc];
 
                      // Read from file
                      read_field_mssg_enc(finput,data_enc,ntot_enc);

                      // Apply decoding routine
                      cout << "  decoding mask_1d_rec, it=" << it << endl;
                      decoding_wrap(nx,ny,nz,mask_1d_rec,tolabs,midval,halfspanval,wlev,nlay,ntot_enc,deps_vec,minval_vec,len_enc_vec,data_enc);
                      cout << "  decode: mask_1d_rec[0]=" << mask_1d_rec[0] << " mask_1d_rec[last]=" << mask_1d_rec[ntot-1UL] << endl;

                      // Reconstruct binary mask
                      mask_midval = midval;
                      for(unsigned long int j1 = 0; j1 < ntot; j1++)
                          if (mask_1d_rec[j1] < midval) mask_1d_rec[j1] = undef; else mask_1d_rec[j1] = 0;

                      // Echo min and max
                      cout << "        min=" << undef << " max=" << 0 << endl;

                      // Deallocate memory
                      delete [] data_enc;

                      // Read header once again
                      read_header_mssg_enc(fheader,it,dsetnamehdr,&tolabs,&midval,&halfspanval,&wlev,&nlay,&ntot_enc,deps_vec,minval_vec,len_enc_vec);
                      //cout << " dsetnamehdr=" << dsetnamehdr << " it=" << it << endl;
                    }
                }

              // Reconstruct field
              if (ntot_enc > 0)
                {
                  // Allocate encoded data array
                  unsigned char *data_enc = new unsigned char[ntot_enc];
 
                  // Read from file
                  read_field_mssg_enc(finput,data_enc,ntot_enc);

                  // Apply decoding routine
                  cout << "  decoding fld_1d_rec, it=" << it << endl;
                  decoding_wrap(nx,ny,nz,fld_1d_rec,tolabs,midval,halfspanval,wlev,nlay,ntot_enc,deps_vec,minval_vec,len_enc_vec,data_enc);
                  cout << "  decode: fld_1d_rec[0]=" << fld_1d_rec[0] << " fld_1d_rec[last]=" << fld_1d_rec[ntot-1UL] << endl;

                  // Calculate min and max
                  double minval = fld_1d_rec[0];
                  double maxval = fld_1d_rec[0];
                  for(unsigned long int j1 = 0; j1 < ntot; j1++)
                    {
                      minval = fmin(minval,fld_1d_rec[j1]);
                      maxval = fmax(maxval,fld_1d_rec[j1]);
                    }

                  // Echo min and max
                  cout << "        min=" << minval << " max=" << maxval << endl;

                  // Deallocate memory
                  delete [] data_enc;
                }

              // Combine the field with the mask
              if (mask_flag) for (unsigned long int j=0; j<ntot; j++) 
                {
                  if (mask_1d_rec[j]<mask_midval) fld_1d_rec[j] = mask_1d_rec[j];
                }

              // Write data in the local domain
              write_field_mssg(out_name.c_str(),flag_convertendian,nbytes,it,nx,ny,nz,nx,ny,0,0,fld_1d_rec);

              // Diagnostics
              cout << "  wrote: fld_1d_rec[0]=" << fld_1d_rec[0] << " fld_1d_rec[last]=" << fld_1d_rec[ntot-1UL] << endl;

              // Deallocate memory
              if (mask_flag) delete [] mask_1d_rec;
              delete [] fld_1d_rec;
            }

          // Close encoded data file
          finput.close();

          // Close headerfile
          fheader.close();

          break;
        }

      // For backup files
      case 1:
      case 2:
        {
          // Dataset parameters
          int ndset;
          char dsettab[NDSMAX][256];

          // Read control file
          control_name = in_prefix_name + ".nmlst";
          read_control_file(control_name.c_str(),nx,ny,nz,nprocx,nprocy,dsettab,ndset);

          // Calculate domain decomposition parameters
          nxloc = nx / nprocx;
          nyloc = ny / nprocy;

          // Print out the control file data
          cout << endl << "=== Parameters read from control file ===" << endl;
          cout << "nx = " << nx << "; ny = " << ny << "; nr(=nz) = " << nz << "; dim_size(=nprocx,nprocy) = " << nprocx << ", " << nprocy << "; ndset = " << ndset << endl;
          for (int j = 0; j<ndset; j++) cout << "record number = " << j+1 << "; field = " << dsettab[j] << endl;

          // Size of the velocity dataset
          if (ifiletype == 1) 
            ntot = (unsigned long int)(nx)*(unsigned long int)(ny)*(unsigned long int)(nz);
          else
            ntot = (unsigned long int)(nxloc)*(unsigned long int)(nyloc)*(unsigned long int)(nz);

          // Copy the control file with a new name
          if (strcmp(in_prefix_name.c_str(),out_prefix_name.c_str()) != 0 )
            {
              string copy_control_name = out_prefix_name + ".nmlst";
              ifstream ictrl;
              ofstream octrl;
              ictrl.open(control_name.c_str(), fstream::in);
              assert(ictrl.is_open());
              octrl.open(copy_control_name.c_str(), fstream::out|fstream::trunc);
              assert(octrl.is_open());
              char c;
              while (ictrl.get(c)) octrl << c;
              ictrl.close();
              octrl.close();  
            }

          // Proc label
          stringstream lbl;
          lbl << setw(MSSG_FILE_DIG) << setfill('0') << thisproc;

          // Input file names
          if (ifiletype == 1) 
            {
              header_name = in_prefix_name + "_h"+ ext_name;
              in_name = in_prefix_name + "_f"+ ext_name;
            }
          else
            {
              header_name = in_prefix_name + "_h"+ lbl.str() + ext_name;
              in_name = in_prefix_name + "_f"+ lbl.str() + ext_name;
            }
          //cout << "Input file names: " << in_name << " " << header_name << endl;

          // Open header file
          ifstream fheader;
          fheader.open(header_name.c_str(), fstream::in);
          assert(fheader.is_open());

          // Open encoded data file name
          ifstream finput;
          finput.open(in_name.c_str(), ios::binary|ios::in);
          assert(finput.is_open());

          // Loop for all datasets
          for (int idset=0; idset<ndset; idset++)
            {
              // Allocate array
              fld_1d_rec = new double[ntot];

              // Initialize array
              for (unsigned long int j=0; j<ntot; j++) fld_1d_rec[j] = 0;

              // Diagnostics
              if (ifiletype == 1) 
                cout << " dset=" << dsettab[idset] << " nx=" << nx << " ny=" << ny << " nz=" << nz << endl;
              else
                cout << " dset=" << dsettab[idset] << " nxloc=" << nxloc << " nyloc=" << nyloc << " nz=" << nz << endl;

              /* Read data from files */
              if (idset == 0)
                {
                  /* Read time record from the header file */
                  // Read time data
                  string str;
                  const int time_rec_len = MSSG_TIME_REC_LEN;
                  for (int j=0; j<12; j++) getline(fheader, str);
                  for (int j=0; j<time_rec_len; j++) fheader >> fld_1d_rec[j];
                  getline(fheader, str);

                  // Print time data on the standard output
                  cout << "  'time' record = ";
                  for (int j=0; j<time_rec_len; j++) cout << " " << setprecision(numeric_limits<long double>::digits10 + 1) << fld_1d_rec[j];
                  cout << endl;

                  // Copy time record to each subdomain, if read from unified encoded file
                  if (ifiletype == 1) 
                    for (int iprocy=0; iprocy<nprocy; iprocy++)
                      for (int iprocx=0; iprocx<nprocx; iprocx++)
                        {
                          // For all elements that contain time record
                          for (int ix = 0; ix < time_rec_len; ix++)
                          {
                            // 1D index
                              unsigned long int j = (unsigned long int)(ix+iprocx*nxloc) + 
                                                    (unsigned long int)(nx)*(unsigned long int)(iprocy*nyloc);
                            // Copy data element 
                            if (iprocx+iprocy>0) fld_1d_rec[j] = fld_1d_rec[ix];
                          }
                        }
                }
              else 
                {
                  /* Read all subsequent records */
                  // Read from the header file with coding attributes
                  char dsetnamehdr[256];
                  read_header_mssg_enc(fheader,idset,dsetnamehdr,&tolabs,&midval,&halfspanval,&wlev,&nlay,&ntot_enc,deps_vec,minval_vec,len_enc_vec);
                }

              /* Read and decode compressed data */
              // Reconstruct data if the compressed data set is non-trivial
              if (idset > 0)
                { 
                if (ntot_enc > 0)
                  {
                    // Allocate encoded data array
                    unsigned char *data_enc = new unsigned char[ntot_enc];
 
                    // Read from file
                    read_field_mssg_enc(finput,data_enc,ntot_enc);

                    /* Do decoding */
                    // Apply decoding routine
                    if (ifiletype == 1) 
                      decoding_wrap(nx,ny,nz,fld_1d_rec,tolabs,midval,halfspanval,wlev,nlay,ntot_enc,deps_vec,minval_vec,len_enc_vec,data_enc);
                    else
                      decoding_wrap(nxloc,nyloc,nz,fld_1d_rec,tolabs,midval,halfspanval,wlev,nlay,ntot_enc,deps_vec,minval_vec,len_enc_vec,data_enc);
                    cout << "  decode: fld_1d_rec[0]=" << fld_1d_rec[0] << " fld_1d_rec[last]=" << fld_1d_rec[ntot-1UL] << endl;

                    // Calculate min and max
                    double minval = fld_1d_rec[0];
                    double maxval = fld_1d_rec[0];
                    for(unsigned long int j1 = 0; j1 < ntot; j1++)
                      {
                        minval = fmin(minval,fld_1d_rec[j1]);
                        maxval = fmax(maxval,fld_1d_rec[j1]);
                      }

                    // Echo min and max
                    cout << "        min=" << minval << " max=" << maxval << endl;

                    // Deallocate memory
                    delete [] data_enc;
                  }
                else // ntot_enc == 0
                  {
                    // All elements are equal
                    for (unsigned long int j=0; j<ntot; j++) fld_1d_rec[j] = midval;
                  }
                }

              /* Append binary floating-point files */
              // Write nproc files or just one file, depending on ifiletype
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

                      // Output file name
                      stringstream lbliproc;
                      lbliproc << setw(MSSG_FILE_DIG) << setfill('0') << iproc;
                      out_name = out_prefix_name + ".p_" + lbliproc.str();
                      //cout << "Writing into " << out_name << endl;

                      // Write data in the subdomain
                      write_field_mssg(out_name.c_str(),flag_convertendian,nbytes,idset,nx,ny,nz,nxloc,nyloc,ixst,iyst,fld_1d_rec);
                    }
                }
              else
                {
                  // Output file name
                  out_name = out_prefix_name + ".p_" + lbl.str();

                  // Write data in the local domain
                  write_field_mssg(out_name.c_str(),flag_convertendian,nbytes,idset,nxloc,nyloc,nz,nxloc,nyloc,0,0,fld_1d_rec);
                }

              // Diagnostics
              cout << "  wrote: fld_1d_rec[0]=" << fld_1d_rec[0] << " fld_1d_rec[last]=" << fld_1d_rec[ntot-1UL] << endl;

              // Deallocate memory
              delete [] fld_1d_rec;
            }

          // Close encoded data file
          finput.close();

          // Close headerfile
          fheader.close();

          break;
        }

      default:
        // Display error message
        cout << "Error: unknown file type" << endl;

    }

    // Display a message on exit
    cout << "=== End of decompression ===\n";

    return 0;
}

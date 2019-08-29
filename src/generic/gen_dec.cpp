/*
    gen_dec.cpp : This file is part of WaveRange CFD data compression utility

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
#include "gen_aux.h"

using namespace std;


// Main code for decoding
int main( int argc, char *argv[] )
{
    // Number of fields in a file
    int nf = 0;

    // Dataset dimensions
    int nx = 0;
    int ny = 0;
    int nz = 0;
    unsigned long int ntot = 0;

    // Higher-dimensional datasets are treated as 3d with the size in the third dimension 
    // equal to the higher-dimensional size nh times the z-size nz
    int nh;
    int nzh;

    // Base tolerance, applied as relative to max(fabs(fld_1d))
    double tol_base;

    // Floating point input file endian conversion flag (0: do not convert; 1: convert)
    int flag_convertendian = 0;

    // Floating point input file precision (4: single; 8: double)
    int nbytes;

    // Data variable declarations
    double *fld_1d_rec;
    double tolabs;
    double midval, halfspanval;
    unsigned char wlev, nlay;
    unsigned long int ntot_enc;
    double deps_vec[NLAYMAX];
    double minval_vec[NLAYMAX];
    unsigned long int len_enc_vec[NLAYMAX];

    // I/O variable declarations
    int ifiletype = 0, idinv, icomp;

    // I/O file names
    string in_name = "data.wrb", header_name = "data.wrh", out_name = "datarec.bin";

    // I/O read buffer string, Fortran record length
    string bar;
    unsigned char recl[8];
    for (int j = 0; j < 8; j++) recl[j] = 0;

    cout << "usage: ./wrdec ENCODED_FILE HEADER_FILE EXTRACTED_FILE TYPE ENDIANFLIP\n";
    cout << "where TYPE=(0: Fortran sequential w 4-byte recl; 1: Fortran sequential w 8-byte recl; 2: C/C++) and ENDIANFLIP=(0:no; 1:yes)\n";
    cout << "interactive mode if not enough arguments are passed.\n";

    /* Prepare for decoding */
    if ( argc == 6 )
    {
      // Read metadata from parameter string
      cout << "automatic mode.";
      in_name = argv[1];
      header_name = argv[2];
      out_name = argv[3];
      bar = argv[4];
      stringstream(bar) >> ifiletype;
      bar = argv[5];
      stringstream(bar) >> flag_convertendian;
    }
    else
    {
      // Read metadata
      cout << "Enter encoded data file name [data.wrb]: ";
      getline (cin,in_name);
      if (in_name.empty()) in_name = "data.wrb";
      cout << "Enter encoding header file name [data.wrh]: ";
      getline (cin,header_name);
      if (header_name.empty()) header_name = "data.wrh";
      cout << "Enter extracted (output) data file name [datarec.bin]: ";
      getline (cin,out_name);
      if (out_name.empty()) out_name = "datarec.bin";
      cout << "Enter file type (0: Fortran sequential w 4-byte recl; 1: Fortran sequential w 8-byte recl; 2: C/C++) [0]: ";
      getline (cin,bar);
      if (!bar.empty()) stringstream(bar) >> ifiletype;
      cout << "Enter endian conversion (0: do not perform; 1: inversion) [0]: ";
      getline (cin,bar);
      if (!bar.empty()) stringstream(bar) >> flag_convertendian;
    }

    // Print out metadata
    cout << endl << "=== Decoding parameters ===" << endl;
    cout << "Encoded data file name " << in_name << endl;
    cout << "Encoding header file name " << header_name << endl;
    cout << "Extracted (output) data file name: " << out_name << endl;
    cout << "File type (0: Fortran sequential w 4-byte recl; 1: Fortran sequential w 8-byte recl; 2: C/C++): " << ifiletype << endl;
    if (flag_convertendian) cout << "Convert big endian to little endian or vice versa" << endl;

    /* Encoding */
    switch (ifiletype) {

      // Fortran sequential access file with 4-byte space for record length
      case 0:
      // Fortran sequential access file with 8-byte space for record length
      case 1:
      // C/C++
      case 2:
        {
          /* Start reading from file */
          // Open header file
          ifstream fheader;
          fheader.open(header_name.c_str(), fstream::in);
          assert(fheader.is_open());

          // Skip first 5 lines from the header file
          string str; 
          for (int j=0; j<5; j++) getline(fheader, str);

          // Read the number of fields
          getline(fheader, str);
          str.erase(0,34);
          stringstream(str) >> nf;

          // Open encoded data file name
          ifstream finput;
          finput.open(in_name.c_str(), ios::binary|ios::in);
          assert(finput.is_open());

          // Loop for all fields in the dataset
          for (int it=0; it<nf; it++)
            {
              // Read from the header file with coding attributes
              read_header_gen_enc(fheader,it,&nbytes,recl,&nx,&ny,&nz,&nh,&idinv,&icomp,&tol_base,&tolabs,&midval,&halfspanval,&wlev,&nlay,&ntot_enc,deps_vec,minval_vec,len_enc_vec);

              // Print number of data points
              cout << "  contains " << nbytes << "-byte floating point data" << endl;
              cout << "  nx=" << nx << "  ny=" << ny << "  nz=" << nz << "  nh=" << nh;
              if (idinv) cout << " and reordering" << endl; else cout << endl;

              // The third and all higher dimensions are concatenated
              nzh = nz*nh;

              // Size of the floating-point array
              ntot = (unsigned long int)(nx)*(unsigned long int)(ny)*(unsigned long int)(nz)*(unsigned long int)(nh);

              // Allocate array
              fld_1d_rec = new double[ntot];

              // If compression flag is true for this field, read and reconstruct
              // Otherwise, read the original field from the file
              if (icomp) 
                {

                  // Initialize array
                  for (unsigned long int j=0; j<ntot; j++) fld_1d_rec[j] = midval;

                  // Reconstruct field
                  if (ntot_enc > 0)
                    {
                      // Allocate encoded data array
                      unsigned char *data_enc = new unsigned char[ntot_enc];
 
                      // Read from file
                      read_field_gen_enc(finput,data_enc,ntot_enc);
    
                      // Apply decoding routine
                      cout << "  decoding fld_1d_rec, field number " << it << endl;
                      decoding_wrap(nx,ny,nzh,fld_1d_rec,tolabs,midval,halfspanval,wlev,nlay,ntot_enc,deps_vec,minval_vec,len_enc_vec,data_enc);
                      cout << "  decode: fld_1d_rec[0]=" << fld_1d_rec[0] << " fld_1d_rec[last]=" << fld_1d_rec[ntot-1UL] << endl;

                      // Deallocate memory
                      delete [] data_enc;
                    }
                }
              else
                {
                  // Read an uncompressed field
                  read_field_gen_raw(finput,nbytes,fld_1d_rec,ntot);
                }

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

              // Write data in the local domain
              write_field_gen(out_name.c_str(),it,ifiletype,flag_convertendian,nbytes,recl,nx,ny,nz,nh,idinv,fld_1d_rec);

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

    return 0;
}

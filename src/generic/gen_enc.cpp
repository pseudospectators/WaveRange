/*
    gen_enc.cpp : This file is part of WaveRange CFD data compression utility

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
//NECs 2019/10/02
#include "../core/trim_split.h"
#include <algorithm> // for std::transform
//NECe 2019/10/02
#include "gen_aux.h"

using namespace std;



// Main code for encoding
int main( int argc, char *argv[] )
{
    // Number of fields in a file
    int nf = 1;

    // Dataset dimensions
    int nx = 16;
    int ny = 16;
    int nz = 16;
    unsigned long int ntot = 16*16*16;

    // Higher-dimensional datasets are treated as 3d with the size in the third dimension 
    // equal to the higher-dimensional size nh times the z-size nz
    int nh = 1;
    int nzh = 16;

    // Base tolerance, applied as relative to max(fabs(fld_1d))
    double tol_base = 1e-16;

    // Floating point input file endian conversion flag (0: do not convert; 1: convert)
    int flag_convertendian = 0;

    // Floating point input file precision (4: single; 8: double)
    int nbytes;

    // Local cutoff tolerance parameters
    int mx, my, mz;
    double *cutoffvec;

    // Field parameter arrays
    int *nbytes_vec, *nx_vec, *ny_vec, *nz_vec, *nh_vec, *idinv_vec, *icomp_vec;
    double *tol_base_vec;

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
    int ifiletype = 0, iintype = 2, idinv = 0, icomp = 1;
    string in_name = "data.bin", out_name = "data.wrb", header_name = "data.wrh";

    // I/O read buffer string, current position in the input file, Fortran record length
    string bar;
    long btpos = 0L;
    unsigned char recl[8];
    for (int j = 0; j < 8; j++) recl[j] = 0;

    //NECs 2019/10/02
    char file_name[] = "inmeta";
    int read_chk_flag = 0;    
    std::ifstream ifs(file_name);

    //check format of inmeta
    if (!ifs.fail()) //check file existence
    {
       std::string str;
       std::string sbuf[11]; // buffer for read params
       // read parameters from file(inmeta)
       std::cout << "==== " << file_name << " exists. ====" << std::endl;
       while (getline(ifs, str))
       {
          std::string str_trim = trim(str, " \t\v\r\n");
          if ( (int)(str_trim.length()) != 0 )
          {
             if ( str_trim.at(0) == '&' ) 
             {
                char delim[] = "=";
                std::vector<std::string> str_split = split(str_trim, delim);
                std::size_t size = str_split.size();
                if ((int)(size) == 2) 
                {
                    read_chk_flag = 1;
                    //std::string var_name = trim(str_split[0].substr(1,str_split[0].length()), " \t\v\r\n"); // 先頭の&を除く
                    std::string var_name = trim(str_split[0], " \t\v\r\n"); 
                    std::string var_value = trim(str_split[1], " \t\v\r\n");
                    std::transform(var_name.begin(), var_name.end(), var_name.begin(), ::tolower);
                    if (var_name == "&in_name") { in_name = var_value;}
                    if (var_name == "&out_name") { out_name = var_value;}
                    if (var_name == "&header_name") { header_name = var_value;}
                    if (var_name == "&file_type") { sbuf[0] = var_value;}
                    if (var_name == "&endian_conversion") { sbuf[1] = var_value;}
                    if (var_name == "&number_of_field") { sbuf[2] = var_value;}

                    //std::cout << var_name << " = " << var_value << std::endl;
                }
                else 
                {
                   /* eeror handling */
                   if ((int)(size) > 1) {
                      std::cout << "==== Error : '=' exists twice in a sentence :" << str_trim << " ====" << std::endl;
                   } else {
                      std::cout << "==== Error : 'value' is missing in a sentence :" << str_trim << " ====" << std::endl;
                   }
                   return -1;
                }
             }
          }
       }
       //
       if (read_chk_flag == 1) 
       {
          std::cout << "in_name = " << in_name << std::endl;
          std::cout << "out_name = " << out_name << std::endl;
          std::cout << "header_name = " << header_name << std::endl;
          std::cout << "file_type = " << sbuf[0] << std::endl;
          std::cout << "endian_conversion = " << sbuf[1] << std::endl;
          std::cout << "number_of_field = " << sbuf[2] << std::endl;
          std::cout << "" << std::endl;
          if (in_name.empty()) in_name = "data.bin";
          if (out_name.empty()) out_name = "data.wrb";
          if (header_name.empty()) header_name = "data.wrh";
          if (!sbuf[0].empty()) stringstream(sbuf[0]) >> ifiletype;
          if (!sbuf[1].empty()) stringstream(sbuf[1]) >> flag_convertendian;
          if (!sbuf[2].empty()) stringstream(sbuf[2]) >> nf;

       
          // Allocate arrays for the parameters of each field
          nbytes_vec = new int[nf];
          nx_vec = new int[nf];
          ny_vec = new int[nf];
          nz_vec = new int[nf];
          nh_vec = new int[nf];
          idinv_vec = new int[nf];
          icomp_vec = new int[nf];
          tol_base_vec = new double[nf];
       
          // read parameters of each field from file(inmeta)
          int chk_nf=0;
          int field_id = -1;
          std::ifstream ifs(file_name);
          while (getline(ifs, str))
          {
             std::string str_trim = trim(str, " \t\v\r\n");
             if ( (int)(str_trim.length()) != 0 )
             {
                if ( str_trim.at(0) == '%' ) 
                {
                   char delim[] = "=";
                   std::vector<std::string> str_split = split(str_trim, delim);
                   std::size_t size = str_split.size();
                   if ((int)(size) == 2) 
                   {
                      std::string field_name = trim(str_split[0], " \t\v\r\n"); 
                      std::string field_value = trim(str_split[1], " \t\v\r\n");
                      std::transform(field_name.begin(), field_name.end(), field_name.begin(), ::tolower);
                      if (field_name == "%field") {
                         if (!field_value.empty()) { 
                            stringstream(field_value) >> field_id;
                            std::cout << "==== read parameters for field " << field_id << " ====" << std::endl;
                            chk_nf = chk_nf + 1;
                         }
                      }
                   }
                }
                if ( str_trim.at(0) == '&' ) 
                {
                   char delim[] = "=";
                   std::vector<std::string> str_split = split(str_trim, delim);
                   std::size_t size = str_split.size();
                   if ((int)(size) == 2) 
                   {
                      std::string var_name = trim(str_split[0], " \t\v\r\n"); 
                      std::string var_value = trim(str_split[1], " \t\v\r\n");
                      std::transform(var_name.begin(), var_name.end(), var_name.begin(), ::tolower);
                      if (var_name == "&input_data_type") { sbuf[3] = var_value;}
                      if (var_name == "&nx") { sbuf[4] = var_value;}
                      if (var_name == "&ny") { sbuf[5] = var_value;}
                      if (var_name == "&nz") { sbuf[6] = var_value;}
                      if (var_name == "&nh") { sbuf[7] = var_value;}
                      if (var_name == "&order") { sbuf[8] = var_value;}
                      if (var_name == "&compress") { sbuf[9] = var_value;}
                      if (var_name == "&tolerance") { sbuf[10] = var_value;}
                   }
                }
                if ( str_trim.at(0) == '/' )  // end of field parameters
                {
                   std::cout << "input_data_type = " << sbuf[3] << std::endl;
                   std::cout << "nx = " << sbuf[4] << std::endl;
                   std::cout << "ny = " << sbuf[5] << std::endl;
                   std::cout << "nz = " << sbuf[6] << std::endl;
                   std::cout << "nh = " << sbuf[7] << std::endl;
                   std::cout << "order = " << sbuf[8] << std::endl;
                   std::cout << "compress = " << sbuf[9] << std::endl;
                   std::cout << "tolerance = " << sbuf[10] << std::endl;
                   std::cout << "" << std::endl;
                   if (!sbuf[3].empty()) stringstream(sbuf[3]) >> iintype;
                   if (!sbuf[4].empty()) stringstream(sbuf[4]) >> nx;
                   if (!sbuf[5].empty()) stringstream(sbuf[5]) >> ny;
                   if (!sbuf[6].empty()) stringstream(sbuf[6]) >> nz;
                   if (!sbuf[7].empty()) stringstream(sbuf[7]) >> nh;
                   if (!sbuf[8].empty()) stringstream(sbuf[8]) >> idinv;
                   if (!sbuf[9].empty()) stringstream(sbuf[9]) >> icomp;
                   if (!sbuf[10].empty()) stringstream(sbuf[10]) >> tol_base;
                   if (iintype==1) nbytes = 4; else nbytes = 8;
                   // Fill the arrays of parameters for field
                   nbytes_vec[field_id] = nbytes;
                   nx_vec[field_id] = nx;
                   ny_vec[field_id] = ny;
                   nz_vec[field_id] = nz;
                   nh_vec[field_id] = nh;
                   idinv_vec[field_id] = idinv;
                   icomp_vec[field_id] = icomp;
                   tol_base_vec[field_id] = tol_base;
                }
             }
          }
          if (chk_nf != nf) {
             std::cout << "==== Number of fields is " << nf << " ====" << std::endl;
             std::cout << "==== Number of blocks for field parameters are not sufficient. ====" << std::endl;
             std::cout << "==== Numbef of block = " << chk_nf << " ====" << std::endl;
             return -1;
          }
       }
       else //read_chk_flag
       {
           std::cout << "==== read parameters from " << file_name << " as old format. ====" << std::endl;
           // read parameters from old format file.
           std::ifstream ifs(file_name);
           getline(ifs, in_name);
           getline(ifs, out_name);
           getline(ifs, header_name);
           getline(ifs, sbuf[0]); 
           getline(ifs, sbuf[1]); 
           getline(ifs, sbuf[2]); 
           std::cout << "in_name = " << in_name << std::endl;
           std::cout << "out_name = " << out_name << std::endl;
           std::cout << "header_name = " << header_name << std::endl;
           std::cout << "file_type = " << sbuf[0] << std::endl;
           std::cout << "endian_conversion = " << sbuf[1] << std::endl;
           std::cout << "number_of_field = " << sbuf[2] << std::endl;
           std::cout << "" << std::endl;
           if (in_name.empty()) in_name = "data.bin";
           if (out_name.empty()) out_name = "data.wrb";
           if (header_name.empty()) header_name = "data.wrh";
           if (!sbuf[0].empty()) stringstream(sbuf[0]) >> ifiletype;
           if (!sbuf[1].empty()) stringstream(sbuf[1]) >> flag_convertendian;
           if (!sbuf[2].empty()) stringstream(sbuf[2]) >> nf;
       
           // Allocate arrays for the parameters of each field
           nbytes_vec = new int[nf];
           nx_vec = new int[nf];
           ny_vec = new int[nf];
           nz_vec = new int[nf];
           nh_vec = new int[nf];
           idinv_vec = new int[nf];
           icomp_vec = new int[nf];
           tol_base_vec = new double[nf];
       
           for (int it=0; it<nf; it++)
           {
              getline(ifs, sbuf[3]); 
              getline(ifs, sbuf[4]); 
              getline(ifs, sbuf[5]); 
              getline(ifs, sbuf[6]); 
              getline(ifs, sbuf[7]); 
              getline(ifs, sbuf[8]); 
              getline(ifs, sbuf[9]); 
              getline(ifs, sbuf[10]); 
              std::cout << "input_data_type = " << sbuf[3] << std::endl;
              std::cout << "nx = " << sbuf[4] << std::endl;
              std::cout << "ny = " << sbuf[5] << std::endl;
              std::cout << "nz = " << sbuf[6] << std::endl;
              std::cout << "nh = " << sbuf[7] << std::endl;
              std::cout << "order = " << sbuf[8] << std::endl;
              std::cout << "compress = " << sbuf[9] << std::endl;
              std::cout << "tolerance = " << sbuf[10] << std::endl;
              std::cout << "" << std::endl;
              if (!sbuf[3].empty()) stringstream(sbuf[3]) >> iintype;
              if (!sbuf[4].empty()) stringstream(sbuf[4]) >> nx;
              if (!sbuf[5].empty()) stringstream(sbuf[5]) >> ny;
              if (!sbuf[6].empty()) stringstream(sbuf[6]) >> nz;
              if (!sbuf[7].empty()) stringstream(sbuf[7]) >> nh;
              if (!sbuf[8].empty()) stringstream(sbuf[8]) >> idinv;
              if (!sbuf[9].empty()) stringstream(sbuf[9]) >> icomp;
              if (!sbuf[10].empty()) stringstream(sbuf[10]) >> tol_base;
              if (iintype==1) nbytes = 4; else nbytes = 8;
              // Fill the arrays of parameters for field
              nbytes_vec[it] = nbytes;
              nx_vec[it] = nx;
              ny_vec[it] = ny;
              nz_vec[it] = nz;
              nh_vec[it] = nh;
              idinv_vec[it] = idinv;
              icomp_vec[it] = icomp;
              tol_base_vec[it] = tol_base;
           }
       }
    }
    else
    {
       std::cout << "==== " << file_name << " doesn't exists. ====" << std::endl;
       //NECe 2019/10/02
        
       // Interactive mode help string
       cout << "usage: ./wrenc INPUT_FILE ENCODED_FILE HEADER_FILE TYPE ENDIANFLIP NF PRECISION NX NY NZ TOLERANCE\n";
       cout << "where TYPE=(0: Fortran sequential w 4-byte recl; 1: Fortran sequential w 8-byte recl; 2: C/C++),\n";
       cout << "      ENDIANFLIP=(0:no; 1:yes), NF=(how many fields, e.g. 1), PRECISION=(1:single; 2:double),\n";
       cout << "      NX=(e.g. 16), NY=(e.g. 16), NZ=(e.g. 16) and TOLERANCE=(e.g. 1.0e-16)\n";
       cout << "interactive mode if not enough arguments are passed.\n";
         
       /* Prepare for encoding */
       if ( argc == 12 )
       {
           // Read metadata from parameter string
           // It is assumed that all fields are three-dimensional, have the same size,
           // are stored using the same floating-point format and will be compressed using the same relative tolerance
           cout << "automatic mode.";
           in_name = argv[1];
           out_name = argv[2];
           header_name = argv[3];
           bar = argv[4];
           stringstream(bar) >> ifiletype;
           bar = argv[5];
           stringstream(bar) >> flag_convertendian;
           bar = argv[6];
           stringstream(bar) >> nf;
           bar = argv[7];
           stringstream(bar) >> iintype;
           if (iintype==1) nbytes = 4; else nbytes = 8;
           bar = argv[8];
           stringstream(bar) >> nx;
           bar = argv[9];
           stringstream(bar) >> ny;
           bar = argv[10];
           stringstream(bar) >> nz;
           bar = argv[11];
           stringstream(bar) >> tol_base;
           // Allocate arrays for the parameters of each field
           nbytes_vec = new int[nf];
           tol_base_vec = new double[nf];
           nx_vec = new int[nf];
           ny_vec = new int[nf];
           nz_vec = new int[nf];
           nh_vec = new int[nf];
           idinv_vec = new int[nf];
           icomp_vec = new int[nf];
           // Fill the arrays of parameters in a loop
           for (int it=0; it<nf; it++)
           {
             nbytes_vec[it] = nbytes;
             tol_base_vec[it] = tol_base;
             nx_vec[it] = nx;
             ny_vec[it] = ny;
             nz_vec[it] = nz;
             nh_vec[it] = nh;
             idinv_vec[it] = idinv;
             icomp_vec[it] = icomp;
          }
       }
       else
       {
          // Read metadata that are in common for all fields in the file
          // In this case, different fields may have different shape, etc.
          cout << "Enter input data file name [data.bin]: ";
          getline (cin,in_name);
          if (in_name.empty()) in_name = "data.bin";
          cout << "Enter encoded data file name [data.wrb]: ";
          getline (cin,out_name);
          if (out_name.empty()) out_name = "data.wrb";
          cout << "Enter encoding header file name [data.wrh]: ";
          getline (cin,header_name);
          if (header_name.empty()) header_name = "data.wrh";
          cout << "Enter file type (0: Fortran sequential w 4-byte recl; 1: Fortran sequential w 8-byte recl; 2: C/C++) [0]: ";
          getline (cin,bar);
          if (!bar.empty()) stringstream(bar) >> ifiletype;
          cout << "Enter endian conversion (0: do not perform; 1: inversion) [0]: ";
          getline (cin,bar);
          if (!bar.empty()) stringstream(bar) >> flag_convertendian;
          // Read number of fields
          cout << "Enter the number of fields in the file, nf [1]: ";
          getline (cin,bar);
          if (!bar.empty()) stringstream(bar) >> nf;
          // Allocate arrays for the parameters of each field
          nbytes_vec = new int[nf];
          tol_base_vec = new double[nf];
          nx_vec = new int[nf];
          ny_vec = new int[nf];
          nz_vec = new int[nf];
          nh_vec = new int[nf];
          idinv_vec = new int[nf];
          icomp_vec = new int[nf];
          // Read parameters of each field in a loop
          for (int it=0; it<nf; it++)
          {
             cout << "Field number " << it << endl;
             cout << "Enter input data type (1: float; 2: double) [2]: ";
             getline (cin,bar);
             if (!bar.empty()) stringstream(bar) >> iintype;
             if (iintype==1) nbytes_vec[it] = 4; else nbytes_vec[it] = 8;
             cout << "Enter the number of data points in the first dimension, nx [16]: ";
             getline (cin,bar);
             if (!bar.empty()) stringstream(bar) >> nx;
             nx_vec[it] = nx;
             cout << "Enter the number of data points in the second dimension, ny [16]: ";
             getline (cin,bar);
             if (!bar.empty()) stringstream(bar) >> ny;
             ny_vec[it] = ny;
             cout << "Enter the number of data points in the third dimension, nz [16]: ";
             getline (cin,bar);
             if (!bar.empty()) stringstream(bar) >> nz;
             nz_vec[it] = nz;
             cout << "Enter the number of data points in the higher (slowest) dimensions, nh [1]: ";
             getline (cin,bar);
             if (!bar.empty()) stringstream(bar) >> nh;
             nh_vec[it] = nh;
             cout << "Invert the order of the dimensions? (0: no; 1: yes) [0]: ";
             getline (cin,bar);
             if (!bar.empty()) stringstream(bar) >> idinv;
             idinv_vec[it] = idinv;
             cout << "Enter compression flag (0: do not compress; 1: compress) [1]: ";
             getline (cin,bar);
             if (!bar.empty()) stringstream(bar) >> icomp;
             icomp_vec[it] = icomp;
             if (icomp)
             {
                cout << "Enter base cutoff relative tolerance [1e-16]: ";
                getline (cin,bar);
                if (!bar.empty()) stringstream(bar) >> tol_base;
                tol_base_vec[it] = tol_base;
             }
             else tol_base_vec[it] = 0;
          }
       }
    } 
    
    // Print out metadata
    cout << endl << "=== Compression parameters ===" << endl;
    cout << "Input data file name: " << in_name << endl;
    cout << "Encoded data file name: " << out_name << endl;
    cout << "Encoding header file name: " << header_name << endl;
    cout << "File type (0: Fortran sequential w 4-byte recl; 1: Fortran sequential w 8-byte recl; 2: C/C++): " << ifiletype << endl;
    if (flag_convertendian) cout << "Convert big endian to little endian or vice versa" << endl;
    cout << "Number of fields in the file, nf: " << nf << endl;

    // Define uniform cutoff
    mx = 1;
    my = 1;
    mz = 1;
    cutoffvec = new double[1];
    cutoffvec[0] = tol_base;

    // Diagnostics
//    cout << " nx=" << nx  << " ny=" << ny << " nz=" << nz << " nf=" << nf << endl;

    // Create header file
    fstream fheader;
    fheader.open(header_name.c_str(), fstream::out | fstream::trunc);
    assert(fheader.is_open());
    fheader << " ===== Header file for compressed data =====" << endl;
    fheader << " Coder version: " << CODER_VERSION << endl;
    fheader << " Encoded data file name: " << out_name.c_str() << endl;
    fheader << " File type (0: Fortran sequential w 4-byte recl; 1: Fortran sequential w 8-byte recl; 2: C/C++): " << ifiletype << endl;
    if (flag_convertendian) fheader << " Converted big endian to little endian or vice versa" << endl; 
      else fheader << " No endian conversion" << endl; 
    fheader << " Number of fields in the file, nf: " << nf << endl;
    fheader.close();

    // Create a new encoded data file. Overwrite if exists
    ofstream foutput;
    foutput.open(out_name.c_str(), ios::binary|ios::out|ios::trunc);
    assert(foutput.is_open());
    foutput.close();

    /* Encoding */
    switch (ifiletype) {

      // Fortran sequential access file with 4-byte space for record length
      case 0:
      // Fortran sequential access file with 8-byte space for record length
      case 1:
      // C/C++
      case 2:
        {
          // Loop for all fields in the dataset
          for (int it=0; it<nf; it++)
            {
              // Print field number on the screen
              cout << "Field number " << it << endl;

              // Set up parameters of the current field
              nbytes = nbytes_vec[it];
              tol_base = tol_base_vec[it];
              idinv = idinv_vec[it];
              icomp = icomp_vec[it];
              nx = nx_vec[it];
              ny = ny_vec[it];
              nz = nz_vec[it];
              nh = nh_vec[it];

              // Print number of data points
              cout << "  contains " << nbytes << "-byte floating point data" << endl;
              cout << "  nx=" << nx << "  ny=" << ny << "  nz=" << nz << "  nh=" << nh;
              if (idinv) cout << " and reordering" << endl; else cout << endl;

              // The third and all higher dimensions are concatenated
              nzh = nz*nh;

              // Size of the floating-point array
              ntot = (unsigned long int)(nx)*(unsigned long int)(ny)*(unsigned long int)(nz)*(unsigned long int)(nh);

              // Allocate array
              fld_1d = new double[ntot];

              /* Read data */
              read_field_gen(in_name.c_str(),it,ifiletype,flag_convertendian,nbytes,recl,nx,ny,nz,nh,idinv,&btpos,fld_1d);
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

              // If compression flag is true for this field, compress and write the compressed data in the file
              // Otherwise, write the original field in the file
              if (icomp) 
                {
                  // Print compression status
                  cout << "  Compression enabled with base relative tolerance " << tol_base << endl;

                  // Allocate encoded data array (will be stored in a file)
                  // Encoded array may be longer than the original
                  unsigned char *data_enc = new unsigned char[SAFETY_BUFFER_FACTOR*NLAYMAX*(ntot<1024UL?1024UL:ntot)];

                  /* Do encoding */
                  // Apply encoding routine
                  encoding_wrap(nx,ny,nzh,fld_1d,1,mx,my,mz,cutoffvec,tolabs,midval,halfspanval,wlev,nlay,ntot_enc,deps_vec,minval_vec,len_enc_vec,data_enc);

                  // Deallocate memory
                  delete [] fld_1d;

                  // Print efficient global cutoff
                  cout << "        tolabs=" << tolabs << endl;

                  /* Write compressed data to a file */
                  // Append the header file with coding attributes
                  write_header_gen_enc(header_name.c_str(),it,&nbytes,recl,&nx,&ny,&nz,&nh,&idinv,&icomp,&tol_base,&tolabs,&midval,&halfspanval,&wlev,&nlay,&ntot_enc,deps_vec,minval_vec,len_enc_vec);

                  // Write data if the compressed data set is non-trivial
                  if (ntot_enc > 0)
                      write_field_gen_enc(out_name.c_str(),data_enc,ntot_enc);

                  // Deallocate memory
                  delete [] data_enc;
                }
              else
                {
                  // Print compression status
                  cout << "  Compression disabled" << endl;

                  /* Write uncompressed data to a file */
                  // Append the header file with coding attributes
                  write_header_gen_enc(header_name.c_str(),it,&nbytes,recl,&nx,&ny,&nz,&nh,&idinv,&icomp,&tol_base,&tolabs,&midval,&halfspanval,&wlev,&nlay,&ntot_enc,deps_vec,minval_vec,len_enc_vec);

                  // Write data in the local domain in the uncompressed C format
                  write_field_gen_raw(out_name.c_str(),nbytes,fld_1d,ntot);

                  // Deallocate memory
                  delete [] fld_1d;
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
    delete [] nbytes_vec;
    delete [] tol_base_vec;
    delete [] nx_vec;
    delete [] ny_vec;
    delete [] nz_vec;
    delete [] nh_vec;
    delete [] idinv_vec;
    delete [] icomp_vec;

    // Display a message on exit
    cout << "=== End of compression ===\n";

    return 0;
}


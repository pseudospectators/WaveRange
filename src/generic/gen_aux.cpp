/*
    gen_aux.cpp : This file is part of WaveRange CFD data compression utility

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

#include <stdlib.h>
#include <string>
#include <cstring>
#include <iostream>
#include <fstream>
#include <exception>
#include <iomanip>
#include <limits>

#include "../core/defs.h"
#include "gen_aux.h"

using namespace std;


/* Write a field to an unformatted fortran binary file */
void write_field_gen( const char *filename, int idset, int ifiletype, int flag_convertendian, int nbytes, unsigned char *recl, int nx, int ny, int nz, int nh, int idinv, double *fld )
{
    // I/O variable declarations
    ofstream outputfile;
    unsigned char foo[ifiletype==0?4:8], temp1[nbytes];
    double buf;

    // Check nbytes value
    if ( (nbytes != 4) && (nbytes != 8) ) 
      {
        // Display error message
        cout << "Generic input nbytes must be equal to 4 or 8" << endl;
        throw std::exception();
      }

    // Print out file name and dataset id 
    //cout << "Output data file name: " << filename << endl; 
    //cout << "Dataset id: " << idset << endl;  

    // Open input file
    if (idset == 0)
      outputfile.open(filename, ios::out | ios::binary | ios::trunc);
    else
      outputfile.open(filename, ios::out | ios::binary | ios::app);

    // Temporary array for conversion
    unsigned char *temp2;

    // If fortran sequential access file, read the record length.
    // See http://gcc.gnu.org/onlinedocs/gcc-3.4.4/g77/Portable-Unformatted-Files.html#fn-1
    // Unformatted sequential records consist of
    //   1. A number giving the length of the record contents;
    //   2. the length of record contents again (for backspace). 
    // The record length is of C type long; this means that it is 8 bytes on 64-bit systems such as Alpha GNU/Linux and 
    // 4 bytes on other systems, such as x86 GNU/Linux. Consequently such files cannot be exchanged between 64-bit and 32-bit systems, 
    // even with the same basic number format. 
    if (ifiletype == 0) // Fortran sequential with 4-byte record length
      {
        if (flag_convertendian)
          for (int j1 = 0; j1 < 4; j1++) foo[j1] = recl[4-1-j1];
        else
          for (int j1 = 0; j1 < 4; j1++) foo[j1] = recl[j1];
        outputfile.write(reinterpret_cast<char*>(foo), 4); 
      }
    else if (ifiletype == 1) // Fortran sequential with 8-byte record length
      { 
        if (flag_convertendian)
          for (int j1 = 0; j1 < 8; j1++) foo[j1] = recl[8-1-j1];
        else
          for (int j1 = 0; j1 < 8; j1++) foo[j1] = recl[j1];
        outputfile.write(reinterpret_cast<char*>(foo), 8); 
      }

    // Write the field in a file
    if (!idinv) // idinv == 0
      {
        // Write in the direct order
        for (int ih = 0; ih < nh; ih++)
          {
          for (int iz = 0; iz < nz; iz++)
            {
            for (int iy = 0; iy < ny; iy++)
              {
              for (int ix = 0; ix < nx; ix++)
                {
                  // 1D index 
                  unsigned long int j = (unsigned long int)(ix) + 
                                        (unsigned long int)(nx)*(unsigned long int)(iy) +
                                        (unsigned long int)(nx)*(unsigned long int)(ny)*(unsigned long int)(iz) +
                                        (unsigned long int)(nx)*(unsigned long int)(ny)*(unsigned long int)(nz)*(unsigned long int)(ih);

                  // Load one element from the data array
                  buf = fld[j];

                  // Reinterpret as float or double
                  if (nbytes == 4) 
                    {
                      float foo = float(buf);
                      float *buf4 = &foo;
                      temp2 = reinterpret_cast<unsigned char*>(buf4);
                    }
                  else
                    {
                      double *buf8 = &buf;
                      temp2 = reinterpret_cast<unsigned char*>(buf8);
                    }

                  // Endian conversion
                  if (flag_convertendian)
                    for (int j1 = 0; j1 < nbytes; j1++) temp1[j1] = temp2[nbytes-1-j1];
                  else
                    for (int j1 = 0; j1 < nbytes; j1++) temp1[j1] = temp2[j1];

                  // Write data element in the file
                  outputfile.write(reinterpret_cast<char*>(temp1), nbytes);
                }
              }
            }
          }
      }
    else  // idinv != 0
      {
        for (int ix = 0; ix < nx; ix++)
          {
          for (int iy = 0; iy < ny; iy++)
            {
            for (int iz = 0; iz < nz; iz++)
              {
              for (int ih = 0; ih < nh; ih++)
                {
                  // 1D index 
                  unsigned long int j = (unsigned long int)(ix) + 
                                        (unsigned long int)(nx)*(unsigned long int)(iy) +
                                        (unsigned long int)(nx)*(unsigned long int)(ny)*(unsigned long int)(iz) +
                                        (unsigned long int)(nx)*(unsigned long int)(ny)*(unsigned long int)(nz)*(unsigned long int)(ih);

                  // Load one element from the data array
                  buf = fld[j];

                  // Reinterpret as float or double
                  if (nbytes == 4) 
                    {
                      float foo = float(buf);
                      float *buf4 = &foo;
                      temp2 = reinterpret_cast<unsigned char*>(buf4);
                    }
                  else
                    {
                      double *buf8 = &buf;
                      temp2 = reinterpret_cast<unsigned char*>(buf8);
                    }

                  // Endian conversion
                  if (flag_convertendian)
                    for (int j1 = 0; j1 < nbytes; j1++) temp1[j1] = temp2[nbytes-1-j1];
                  else
                    for (int j1 = 0; j1 < nbytes; j1++) temp1[j1] = temp2[j1];

                  // Write data element in the file
                  outputfile.write(reinterpret_cast<char*>(temp1), nbytes);
                }
              }
            }
          }
      }

    // If fortran sequential access file, write the record length again.
    if (ifiletype == 0) // Fortran sequential with 4-byte record length
      {
        if (flag_convertendian)
          for (int j1 = 0; j1 < 4; j1++) foo[j1] = recl[4-1-j1];
        else
          for (int j1 = 0; j1 < 4; j1++) foo[j1] = recl[j1];
        outputfile.write(reinterpret_cast<char*>(foo), 4); 
      }
    else if (ifiletype == 1) // Fortran sequential with 8-byte record length
      { 
        if (flag_convertendian)
          for (int j1 = 0; j1 < 8; j1++) foo[j1] = recl[8-1-j1];
        else
          for (int j1 = 0; j1 < 8; j1++) foo[j1] = recl[j1];
        outputfile.write(reinterpret_cast<char*>(foo), 8); 
      }

    // Close output file
    outputfile.close();
}


/* Read a field from an unformatted fortran binary file */
void read_field_gen( const char *filename, int idset, int ifiletype, int flag_convertendian, int nbytes, unsigned char *recl, int nx, int ny, int nz, int nh, int idinv, unsigned long int *btpos, double *fld )
{
    // I/O variable declarations
    ifstream inputfile;
    unsigned char foo[ifiletype==0?4:8], temp1[nbytes], temp2[nbytes];
    double buf;

    // Check nbytes value
    if ( (nbytes != 4) && (nbytes != 8) ) 
      {
        // Display error message
        cout << "Generic input nbytes must be equal to 4 or 8" << endl;
        throw std::exception();
      }

    // Print out file name and dataset id 
    //cout << "Input data file name: " << filename << endl; 
    //cout << "Dataset id: " << idset << endl;  

    // Open input file
    inputfile.open(filename, ios::in|ios::binary);

    // Skip preceding datasets, read by blocks of 1 byte
    for (unsigned long int j = 0UL; j < *btpos/4; j++ )
      inputfile.read(reinterpret_cast<char*>(foo), 4);

    // If fortran sequential access file, read the record length.
    // See http://gcc.gnu.org/onlinedocs/gcc-3.4.4/g77/Portable-Unformatted-Files.html#fn-1
    // Unformatted sequential records consist of
    //   1. A number giving the length of the record contents;
    //   2. the length of record contents again (for backspace). 
    // The record length is of C type long; this means that it is 8 bytes on 64-bit systems such as Alpha GNU/Linux and 
    // 4 bytes on other systems, such as x86 GNU/Linux. Consequently such files cannot be exchanged between 64-bit and 32-bit systems, 
    // even with the same basic number format. 
    if (ifiletype == 0) // Fortran sequential with 4-byte record length
      {
        *btpos += 4UL;
        inputfile.read(reinterpret_cast<char*>(foo), 4); 
        if (flag_convertendian)
          for (int j1 = 0; j1 < 4; j1++) recl[j1] = foo[4-1-j1];
        else
          for (int j1 = 0; j1 < 4; j1++) recl[j1] = foo[j1];
      }
    else if (ifiletype == 1) // Fortran sequential with 8-byte record length
      { 
        *btpos += 8UL;
        inputfile.read(reinterpret_cast<char*>(foo), 8); 
        if (flag_convertendian)
          for (int j1 = 0; j1 < 8; j1++) recl[j1] = foo[8-1-j1];
        else
          for (int j1 = 0; j1 < 8; j1++) recl[j1] = foo[j1];
      }
 
    // Read from file
    if (!idinv) // idinv == 0
      {
        // Read and store in the direct order
        for (int ih = 0; ih < nh; ih++)
          {
          for (int iz = 0; iz < nz; iz++)
            {
            for (int iy = 0; iy < ny; iy++)
              {
              for (int ix = 0; ix < nx; ix++)
                {
                  // Read data element from file
                  inputfile.read(reinterpret_cast<char*>(temp1), nbytes);

                  // Endian conversion
                  if (flag_convertendian)
                    for (int j1 = 0; j1 < nbytes; j1++) temp2[j1] = temp1[nbytes-1-j1];
                  else
                    for (int j1 = 0; j1 < nbytes; j1++) temp2[j1] = temp1[j1];

                  // Reinterpret as float or double
                  if (nbytes == 4)
                    {
                      buf = reinterpret_cast<float&>(temp2);
                      *btpos += 4UL;
                    }
                  else
                    {
                      buf = reinterpret_cast<double&>(temp2);
                      *btpos += 8UL;
                    }

                  // 1D index 
                  unsigned long int j = (unsigned long int)(ix) + 
                                        (unsigned long int)(nx)*(unsigned long int)(iy) +
                                        (unsigned long int)(nx)*(unsigned long int)(ny)*(unsigned long int)(iz) +
                                        (unsigned long int)(nx)*(unsigned long int)(ny)*(unsigned long int)(nz)*(unsigned long int)(ih);

                  // Fill in the data array
                  fld[j] = buf;
                }
              }
            }
          }
      }
    else  // idinv != 0
      {
        // Read and store in the reverse order
        for (int ix = 0; ix < nx; ix++)
          {
          for (int iy = 0; iy < ny; iy++)
            {
            for (int iz = 0; iz < nz; iz++)
              {
              for (int ih = 0; ih < nh; ih++)
                {
                  // Read data element from file
                  inputfile.read(reinterpret_cast<char*>(temp1), nbytes);

                  // Endian conversion
                  if (flag_convertendian)
                    for (int j1 = 0; j1 < nbytes; j1++) temp2[j1] = temp1[nbytes-1-j1];
                  else
                    for (int j1 = 0; j1 < nbytes; j1++) temp2[j1] = temp1[j1];

                  // Reinterpret as float or double
                  if (nbytes == 4)
                    {
                      buf = reinterpret_cast<float&>(temp2);
                      *btpos += 4UL;
                    }
                  else
                    {
                      buf = reinterpret_cast<double&>(temp2);
                      *btpos += 8UL;
                    }

                  // 1D index 
                  unsigned long int j = (unsigned long int)(ix) + 
                                        (unsigned long int)(nx)*(unsigned long int)(iy) +
                                        (unsigned long int)(nx)*(unsigned long int)(ny)*(unsigned long int)(iz) +
                                        (unsigned long int)(nx)*(unsigned long int)(ny)*(unsigned long int)(nz)*(unsigned long int)(ih);

                  // Fill in the data array
                  fld[j] = buf;
                }
              }
            }
          }
      }

    // If fortran sequential access file, read the record length again and discard the value.
    if (ifiletype == 0) // Fortran sequential with 4-byte record length
      {
        *btpos += 4UL;
        inputfile.read(reinterpret_cast<char*>(foo), 4); 
      }
    else if (ifiletype == 1) // Fortran sequential with 8-byte record length
      { 
        *btpos += 8UL;
        inputfile.read(reinterpret_cast<char*>(foo), 8); 
      }

    // Close input file
    inputfile.close();

    // Exit if cannot read
    if ( inputfile.fail() )
      {
        // Display error message
        cout << "Cannot read from " << filename << endl;
        throw std::exception();
      }
}


/* Write unsigned char type data set */
void write_field_gen_enc( const char *filename, unsigned char *fld, unsigned long int ntot_enc )
{
    ofstream outputfile;
    outputfile.open(filename, ios::binary|ios::out|ios::app);
    outputfile.write(reinterpret_cast<char*>(fld), ntot_enc);
    outputfile.close();
}


/* Read unsigned char type data set */
void read_field_gen_enc( ifstream &inputfile, unsigned char *fld, unsigned long int ntot_enc )
{
    inputfile.read(reinterpret_cast<char*>(fld), ntot_enc);
}


/* Write double type data set */
void write_field_gen_raw( const char *filename, int nbytes, double *fld, unsigned long int ntot )
{
    // I/O variable declarations
    ofstream outputfile;
    unsigned char *temp2;
    double buf;

    // Check nbytes value
    if ( (nbytes != 4) && (nbytes != 8) ) 
      {
        // Display error message
        cout << "Generic input nbytes must be equal to 4 or 8" << endl;
        throw std::exception();
      }

    // File open
    outputfile.open(filename, ios::binary|ios::out|ios::app);

    // 1d loop for all elements of the array
    for (int j = 0; j < ntot; j++)
      {
        // Load one element from the data array
        buf = fld[j];

        // Reinterpret as float or double
        if (nbytes == 4) 
          {
            float foo = float(buf);
            float *buf4 = &foo;
            temp2 = reinterpret_cast<unsigned char*>(buf4);
          }
        else
          {
            double *buf8 = &buf;
            temp2 = reinterpret_cast<unsigned char*>(buf8);
          }

        // Write data element in the file
        outputfile.write(reinterpret_cast<char*>(temp2), nbytes);
      }

    // File close
    outputfile.close();
}


/* Read double type data set */
void read_field_gen_raw( ifstream &inputfile, int nbytes, double *fld, unsigned long int ntot )
{
    // I/O variable declarations
    unsigned char temp1[nbytes];
    double buf;

    // Check nbytes value
    if ( (nbytes != 4) && (nbytes != 8) ) 
      {
        // Display error message
        cout << "Generic input nbytes must be equal to 4 or 8" << endl;
        throw std::exception();
      }

    // 1d loop for all elements of the array
    for (int j = 0; j < ntot; j++)
      {
        // Read data element from file
        inputfile.read(reinterpret_cast<char*>(temp1), nbytes);

        // Reinterpret as float or double
        if (nbytes == 4) 
          buf = reinterpret_cast<float&>(temp1);
        else 
          buf = reinterpret_cast<double&>(temp1);

        // Fill in the data array
        fld[j] = buf;
      }
}


/* Write a regular record in the encoding header file */
void write_header_gen_enc( const char *filename, int idset, int *nbytes, unsigned char *recl, int *nx, int *ny, int *nz, int *nh, int *idinv, int *icomp, double *tol_base, double *tolabs, double *midval, double *halfspanval, unsigned char *wlev, unsigned char *nlay, unsigned long int *ntot_enc, double *deps_vec, double *minval_vec, unsigned long int *len_enc_vec )
{
   // Open file
   ofstream fs;
   fs.open(filename, ofstream::out|ofstream::app);

   // Append with a new dataset
   fs << " -----" << endl;
   fs << idset << endl;
   // Write name reminder string
   fs << " nbytes; recl; nx; ny; nz; nh; idinv; icomp;";
   if (*icomp) fs << " tol_base; tolabs; midval; halfspanval; wlev; nlay; ntot_enc;";
   if (*ntot_enc > 0) fs << " deps_vec(1:nlay); minval_vec(1:nlay); len_enc_vec(1:nlay)" << endl; else fs << endl;
   // Write values
   fs << *nbytes << endl;
   for (int j = 0; j < 8; j++) fs << std::hex << static_cast<unsigned>(recl[j]) << " ";
   fs << std::dec << endl;
   fs << *nx << endl;
   fs << *ny << endl;
   fs << *nz << endl;
   fs << *nh << endl;
   fs << *idinv << endl;
   fs << *icomp << endl;
   // Only write the following if this field will be compressed
   if (*icomp > 0)
     {
     fs << setprecision(numeric_limits<long double>::digits10 + 1) << *tol_base << endl;
     fs << setprecision(numeric_limits<long double>::digits10 + 1) << *tolabs << endl;
     fs << setprecision(numeric_limits<long double>::digits10 + 1) << *midval << endl;
     fs << setprecision(numeric_limits<long double>::digits10 + 1) << *halfspanval << endl;
     fs << static_cast<unsigned>(*wlev) << endl;
     fs << static_cast<unsigned>(*nlay) << endl;
     fs << *ntot_enc << endl;
     // Only write the following for non-trivial datasets
     if (*ntot_enc > 0)
       {
         for (int j=0; j<*nlay; j++) 
           fs << setprecision(numeric_limits<long double>::digits10 + 1) << deps_vec[j] << " ";
         fs << endl;
         for (int j=0; j<*nlay; j++) 
           fs << setprecision(numeric_limits<long double>::digits10 + 1) << minval_vec[j] << " ";
         fs << endl;
         for (int j=0; j<*nlay; j++) 
           fs << len_enc_vec[j] << " ";
         fs << endl;
       }
     }

   // Close file
   fs.close();
}

/* Read a regular record from the encoding header file */
void read_header_gen_enc( ifstream &fs, int idset, int *nbytes, unsigned char *recl, int *nx, int *ny, int *nz, int *nh, int *idinv, int *icomp, double *tol_base, double *tolabs, double *midval, double *halfspanval, unsigned char *wlev, unsigned char *nlay, unsigned long int *ntot_enc, double *deps_vec, double *minval_vec, unsigned long int *len_enc_vec )
{
   string str; 
   // Skip 1 line
   getline(fs, str);
   // Read and check dataset id
   int idset1;
   fs >> idset1;
   if (idset1 != idset) 
     {
       // Display error message
       cout << "Encoding header file read error" << endl;
       cout << "Reading field " << idset << ", found field " << idset1 << endl;
       throw std::exception();
     }
   // Skip remaining characters in the line
   getline(fs, str);
   // Skip 1 line
   getline(fs, str);
   // Read values
   fs >> *nbytes;
   for (int j = 0; j < 8; j++) 
     {
       int buf;
       fs >> std::hex >> buf;
       recl[j] = buf;
     }
   fs >> std::dec;
   getline(fs, str);
   fs >> *nx;
   fs >> *ny;
   fs >> *nz;
   fs >> *nh;
   fs >> *idinv;
   fs >> *icomp;
   // Only read the following if this field is be compressed
   if (*icomp > 0)
     {
       fs >> *tol_base;
       fs >> *tolabs;
       fs >> *midval;
       fs >> *halfspanval;
       int buf;
       fs >> buf;
       *wlev = buf;
       fs >> buf;
       *nlay = buf;
       fs >> *ntot_enc;
       // Skip the line end
       getline(fs, str);
       // Only read the following for non-trivial datasets
       if (*ntot_enc > 0)
         {
           for (int j=0; j<*nlay; j++) 
             fs >> deps_vec[j];
           getline(fs, str);
           for (int j=0; j<*nlay; j++) 
             fs >> minval_vec[j];
           getline(fs, str);
           for (int j=0; j<*nlay; j++) 
             fs >> len_enc_vec[j];
           getline(fs, str);
         }
     }
   else getline(fs, str);
 
   // Print out values from the header file
   cout << "  tolabs; midval; halfspanval; wlev; nlay; ntot_enc;";
   if (*ntot_enc > 0) cout << " deps_vec(1:nlay); minval_vec(1:nlay); len_enc_vec(1:nlay)" << endl; else cout << endl;
   cout << "  " << *tolabs << " " << *midval << " " << *halfspanval << " " << static_cast<unsigned>(*wlev) << " " << static_cast<unsigned>(*nlay) << " " << *ntot_enc << endl;
   if (*ntot_enc > 0)
     {
       cout << "  ";
       for (int j=0; j<*nlay; j++) 
         cout << deps_vec[j] << " ";
       cout << endl;
       cout << "  ";
       for (int j=0; j<*nlay; j++) 
         cout << minval_vec[j] << " ";
       cout << endl;
       cout << "  ";
       for (int j=0; j<*nlay; j++) 
         cout << len_enc_vec[j] << " ";
       cout << endl;
     }
}

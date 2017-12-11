/*
    ctrl_aux.cpp : This file is part of WaveRange CFD data compression utility

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

#include <stdlib.h>
#include <string>
#include <cstring>
#include <iostream>
#include <fstream>
#include <exception>
#include <iomanip>
#include <limits>

#include "../core/defs.h"
#include "ctrl_aux.h"

using namespace std;


/* Read control file for restart files */
void read_control_file( const char *filename, int &nx, int &ny, int &nz, int &nprocx, int &nprocy, char dsettab[NDSMAX][256], int &ndset )
{
    const int szpar = 50;
    int nproc, ipar, npar = 0, expect_val = 0, pos = 0, state = 0;
    char c, partab[szpar][256], valtab[szpar][256], str[256];
    ifstream is;
    
    // Open the control file
    is.open(filename);
    if (!is) 
      {
        // Display error message
        cout << "Unable to open namelist file" << endl;
        throw std::exception();
      }
    
    // Fill in the parameter table
    while (is.get(c))
      {
        switch (c) {
        // These symbols are separators
        case '\n':
        case '&':
        case ' ':
        case '\'':
        case ',':
          {
            if (pos == 0) break;
            str[pos] = '\0';

            switch (state) {
            // Reading parameter name
            case 1:
              {
                // List of relevant parameter names
                if ( (strcmp(str,"nx") == 0) ||
                     (strcmp(str,"ny") == 0) ||
                     (strcmp(str,"nr") == 0) ||
                     (strcmp(str,"nproc") == 0) || 
                     (strcmp(str,"dim_size") == 0) || 
                     (strcmp(str,"var") == 0) || 
                     (strcmp(str,"rec") == 0) )
                  {
                    strcpy(partab[npar],str);
                    expect_val = 1;
                    //cout << "name = " << partab[npar] << endl;
                  }
                break;
              }
            // Reading parameter value
            case 2:
              {
                if (expect_val == 1)
                  {
                    strcpy(valtab[npar],str);
                    expect_val = 0;
                    //cout << "value = " << valtab[npar] << endl;
                    npar++;
                  }
              }
            }

            state = 0;
            pos = 0;
            break;
          }
        // Switch to parameter value reading state
        case '=':
          {
            state = 2;
            pos = 0;
            break;
          }
        // Reading parameter name or value in this case
        default:
          {
            if (state != 2) state = 1;
            str[pos++] = c;
          }
        }
      }

    // Close the control file
    is.close();

    // Assign values of the parameter variables
    // nx
    for (ipar = 0; ipar < npar; ipar++)
      if (strcmp(partab[ipar],"nx") == 0) nx = atoi(valtab[ipar]);
    // ny
    for (ipar = 0; ipar < npar; ipar++)
      if (strcmp(partab[ipar],"ny") == 0) ny = atoi(valtab[ipar]); 
    // nr (=nz)
    for (ipar = 0; ipar < npar; ipar++)
      if (strcmp(partab[ipar],"nr") == 0) nz = atoi(valtab[ipar]); 
    // nproc
    for (ipar = 0; ipar < npar; ipar++)
      if (strcmp(partab[ipar],"nproc") == 0) nproc = atoi(valtab[ipar]);    
    // dim_size (=nprocx)
    for (ipar = 0; ipar < npar; ipar++)
      if (strcmp(partab[ipar],"dim_size") == 0) nprocx = atoi(valtab[ipar]);    
    // nprocy 
    nprocy = nproc / nprocx;

    // Look for the first "var" record
    ipar = 0;
    while (strcmp(partab[ipar],"var") != 0) ipar++;
    // Fill in the dataset name table in the order of records in the file
    ndset = 0;
    while (ipar < npar)
      {
        strcpy(dsettab[atoi(valtab[ipar+1])-1],valtab[ipar]);
        //cout << "fld = " << dsettab[atoi(valtab[ipar+1])-1] << " rec = " << atoi(valtab[ipar+1])-1 << endl;
        ipar += 2;
        ndset++;
      }
}


/* Read control file for regular output in GRADS format */
void read_control_file_grads( const char *filename, int &nx, int &ny, int &nz, int &nt, char *dsetname )
{
    const int szpar = 50;
    int ipar, npar = 0, pos = 0, state = 1;
    char c, partab[szpar][256], valtab[szpar][256], str[256];
    ifstream is;
    
    // Open the control file
    is.open(filename);
    if (!is) 
      {
        // Display error message
        cout << "Unable to open namelist file" << endl;
        throw std::exception();
      }
    
    // Fill in the parameter table
    while (is.get(c))
      {
        // 3 possible states: 0 - skip; 1 - read name; 2 - read value
        switch (c) {
        // Separators
        case '\n':
        case '^':
        case ' ':
          {
            if (pos > 0)
              {
                str[pos] = '\0';

                switch (state) {
                // Reading parameter name
                case 1:
                  {
                    // List of relevant parameter names
                    if ( (strcmp(str,"DSET") == 0) ||
                         (strcmp(str,"XDEF") == 0) ||
                         (strcmp(str,"YDEF") == 0) ||
                         (strcmp(str,"ZDEF") == 0) || 
                         (strcmp(str,"TDEF") == 0) )
                      {
                        strcpy(partab[npar],str);
                        state = 2;
                        //cout << "name = " << partab[npar] << endl;
                      }
                    break;
                  }
                // Reading parameter value
                case 2:
                  {
                    if ( pos > 0 )
                      {
                        strcpy(valtab[npar],str);
                        state = 0;
                        //cout << "value = " << valtab[npar] << endl;
                        npar++;
                      }
                    break;
                  }
                }
                pos = 0;
              }
            // New parameter name is expected at new line
            if (c == '\n') state = 1;
            break;
          }
        // Reading parameter name or value in this case
        default:
          {
            str[pos++] = c;
          }
        }
      }

    // Close the control file
    is.close();

    // Assign values of the parameter variables
    // DSET (=dsetname)
    for (ipar = 0; ipar < npar; ipar++)
      if (strcmp(partab[ipar],"DSET") == 0) strcpy(dsetname,valtab[ipar]);    
    // XDEF (=nx)
    for (ipar = 0; ipar < npar; ipar++)
      if (strcmp(partab[ipar],"XDEF") == 0) nx = atoi(valtab[ipar]);
    // YDEF (=ny)
    for (ipar = 0; ipar < npar; ipar++)
      if (strcmp(partab[ipar],"YDEF") == 0) ny = atoi(valtab[ipar]); 
    // ZDEF (=nz)
    for (ipar = 0; ipar < npar; ipar++)
      if (strcmp(partab[ipar],"ZDEF") == 0) nz = atoi(valtab[ipar]); 
    // TDEF (=nt)
    for (ipar = 0; ipar < npar; ipar++)
      if (strcmp(partab[ipar],"TDEF") == 0) nt = atoi(valtab[ipar]);    
}


/* Write a field to an unformatted fortran binary file */
void write_field_mssg( const char *filename, int flag_convertendian, int nbytes, int idset, int nx, int ny, int nz, int nxloc, int nyloc, int ixst, int iyst, double *fld )
{
    // I/O variable declarations
    ofstream outputfile;
    unsigned char temp1[nbytes];
    double buf;

    // Check nbytes value
    if ( (nbytes != 4) && (nbytes != 8) ) 
      {
        // Display error message
        cout << "MSSG input nbytes must be equal to 4 or 8" << endl;
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

    // Write in fortran order
    for (int iz = 0; iz < nz; iz++)
      {
      for (int iy = iyst; iy < iyst+nyloc; iy++)
        {
        for (int ix = ixst; ix < ixst+nxloc; ix++)
          {
            // 1D index 
            unsigned long int j = (unsigned long int)(ix) + 
                                  (unsigned long int)(nx)*(unsigned long int)(iy) +
                                  (unsigned long int)(nx)*(unsigned long int)(ny)*(unsigned long int)(iz);

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
              for (int j = 0; j < nbytes; j++) temp1[j] = temp2[nbytes-1-j];

            // Write data element in the file
            outputfile.write(reinterpret_cast<char*>(temp1), nbytes);
          }
        }
      }

    // Close output file
    outputfile.close();
}


/* Read a field from an unformatted fortran binary file */
void read_field_mssg( const char *filename, int flag_convertendian, int nbytes, int idset, int nx, int ny, int nz, int nxloc, int nyloc, int ixst, int iyst, double *fld )
{
    // I/O variable declarations
    ifstream inputfile;
    unsigned char temp1[nbytes], temp2[nbytes];
    double buf;

    // Check nbytes value
    if ( (nbytes != 4) && (nbytes != 8) ) 
      {
        // Display error message
        cout << "MSSG input nbytes must be equal to 4 or 8" << endl;
        throw std::exception();
      }

    // Print out file name and dataset id 
    //cout << "Input data file name: " << filename << endl; 
    //cout << "Dataset id: " << idset << endl;  

    // Open input file
    inputfile.open(filename, ios::in|ios::binary);

    // Skip preceding datasets 
    for (int j = 0; j < idset; j++ )
      for (int iz = 0; iz < nz; iz++)
        for (int iy = iyst; iy < iyst+nyloc; iy++)
          for (int ix = ixst; ix < ixst+nxloc; ix++)
            inputfile.read(reinterpret_cast<char*>(temp1), nbytes);

    // Read in fortran order
    for (int iz = 0; iz < nz; iz++)
      {
      for (int iy = iyst; iy < iyst+nyloc; iy++)
        {
        for (int ix = ixst; ix < ixst+nxloc; ix++)
          {
            // Read data element from file
            inputfile.read(reinterpret_cast<char*>(temp1), nbytes);

            // Endian conversion
            if (flag_convertendian)
              for (int j = 0; j < nbytes; j++) temp2[j] = temp1[nbytes-1-j];

            // Reinterpret as float or double
            if (nbytes == 4)
              buf = reinterpret_cast<float&>(temp2);
            else
              buf = reinterpret_cast<double&>(temp2);

            // 1D index 
            unsigned long int j = (unsigned long int)(ix) + 
                                  (unsigned long int)(nx)*(unsigned long int)(iy) +
                                  (unsigned long int)(nx)*(unsigned long int)(ny)*(unsigned long int)(iz);

            // Fill in the data array
            fld[j] = buf;
          }
        }
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
void write_field_mssg_enc( const char *filename, unsigned char *fld, unsigned long int ntot_enc )
{
    ofstream outputfile;
    outputfile.open(filename, ios::binary|ios::out|ios::app);
    outputfile.write(reinterpret_cast<char*>(fld), ntot_enc);
    outputfile.close();
}


/* Read unsigned char type data set */
void read_field_mssg_enc( ifstream &inputfile, unsigned char *fld, unsigned long int ntot_enc )
{
    inputfile.read(reinterpret_cast<char*>(fld), ntot_enc);
}


/* Write a regular record in the encoding header file */
void write_header_mssg_enc( const char *filename, int idset, const char *dsetname, double *tolabs, double *midval, double *halfspanval, unsigned char *wlev, unsigned char *nlay, unsigned long int *ntot_enc, double *deps_vec, double *minval_vec, unsigned long int *len_enc_vec )
{
   // Open file
   ofstream fs;
   fs.open(filename, ofstream::out|ofstream::app);

   // Append with a new dataset
   fs << " -----" << endl;
   fs << idset+1 << endl;
   fs << " Data set name = " << dsetname << endl;
   // Write name reminder string
   fs << " tolabs; midval; halfspanval; wlev; nlay; ntot_enc;";
   if (*ntot_enc > 0) fs << " deps_vec(1:nlay); minval_vec(1:nlay); len_enc_vec(1:nlay)" << endl; else fs << endl;
   // Write values
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

   // Close file
   fs.close();
}

/* Read a regular record from the encoding header file */
void read_header_mssg_enc( ifstream &fs, int idset, char *dsetname, double *tolabs, double *midval, double *halfspanval, unsigned char *wlev, unsigned char *nlay, unsigned long int *ntot_enc, double *deps_vec, double *minval_vec, unsigned long int *len_enc_vec )
{
   string str; 
   // Skip 1 line
   getline(fs, str);
   // Read and check dataset id
   int idset1;
   fs >> idset1;
   if (idset1 != idset+1) 
     {
       // Display error message
       cout << "Encoding header file does not match with the control file" << endl;
       cout << "idset+1 = " << idset+1 << " idset1 = " << idset1 << endl;
       throw std::exception();
     }
   // Skip remaining characters in the line
   getline(fs, str);
   // Read dataset name
   getline(fs, str);
   str.erase(0,17);
   strcpy (dsetname, str.c_str());
   // Skip 1 line
   getline(fs, str);
   // Read values
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

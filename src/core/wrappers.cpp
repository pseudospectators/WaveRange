/*
    wrappers.cpp : This file is part of WaveRange CFD data compression utility

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
#include <float.h>

#include <iostream>
#include <exception>
#include <memory>

extern "C" {
#include "../rangecod/port.h"
#include "../rangecod/rangecod.h"
#include "../waveletcdf97_3d/waveletcdf97_3d.h"
}

#include "../core/defs.h"
#include "wrappers.h"

using namespace std;


/* Calculate local precision */
double lcl_prec(int nx, int ny, int nz, int jx, int jy, int jz, int mx, int my, int mz, double *cutoffvec)
{
    // Cartesian coordinates of the local precision cutoff block
    int kx = int(double(jx)/double(nx)*double(mx));
    int ky = int(double(jy)/double(ny)*double(my));
    int kz = int(double(jz)/double(nz)*double(mz));

    // Evaluate the local precision
    return cutoffvec[kx+mx*ky+mx*my*kz];
}


/* Encoding subroutine with wavelet transform and range coding */ 
void encoding_wrap(int nx, int ny, int nz, double *fld_1d, int wtflag, int mx, int my, int mz, double *cutoffvec, double& tolabs, double& midval, double& halfspanval, unsigned char& wlev, unsigned char& nlay, unsigned long int& ntot_enc, double *deps_vec, double *minval_vec, unsigned long int *len_enc_vec, unsigned char *data_enc)
{
    /* Wavelet decomposition */
    // Print wavelet decomposition status
    cout << "Wavelet decomposition..." << endl;

    // Total number of elements in the input array
    unsigned long int ntot = (unsigned long int)(nx)*(unsigned long int)(ny)*(unsigned long int)(nz);

    // Number of elements in the local cutoff array
    unsigned int mtot = mx*my*mz;

    // Wavelet transform depth, hardcoded, see header file; zero if no transform required
    if (wtflag) wlev = WAV_LVL; else wlev = 0;

    // Find the minimum and maximum values
    double minval = fld_1d[0];
    double maxval = fld_1d[0];
    for (unsigned long int j = 0; j < ntot; j++)
      {
        double fld1_tmp = fld_1d[j];
        if (fld1_tmp < minval) minval = fld1_tmp;
        if (fld1_tmp > maxval) maxval = fld1_tmp;
      }

    // Find the middle value and the half-span of the data values
    halfspanval = (maxval-minval)/2;
    midval = minval+halfspanval;

    // If the half-span is close to zero, encoding is impossible and unnecessary
    if (halfspanval <= 2*DBL_MIN)
      {
        // Encoded data array is empty and not used
        ntot_enc = 0;
        nlay = 0;
        tolabs = 0;

        // Exit from the subroutine
        return;
      }

    // Apply wavelet transform
    waveletcdf97_3d(nx,ny,nz,int(wlev),fld_1d);

    /* Range encoding */
    // Print encoding statement
    cout << "Range encoding..." << endl;

    // Alphabet size
    int q = 256;

    // Allocate the quantized input and output vectors
    unsigned char *fld_q = new unsigned char[ntot];
    unsigned char *enc_q = new unsigned char[(SAFETY_BUFFER_FACTOR+1UL)*ntot]; // Encoded array may be longer than the original

    // Output quantized data array length
    unsigned long int len_out_q = 0;

    // Output vector counter
    unsigned long int jtot = 0;

    // Byte layer counter
    unsigned char ilay = 0;

    // Minimum cutoff
    double tolrel = cutoffvec[0];
    for (unsigned int k=1; k<mtot; k++) if (cutoffvec[k] < tolrel) tolrel = cutoffvec[k];

    // Absolute tolerance
    tolabs = tolrel * max(fabs(minval),fabs(maxval));

    // Apply a correction for round-off errors in wavelet transform
    tolabs /= WAV_ACC_COEF;

    // Iteration break flag set to false by default
    unsigned char brflag = 0;

    // Quantize and encode all byte layers
    while (1)
    {
        // Calculate min and max
        maxval = fld_1d[0];
        minval = fld_1d[0];
        for(unsigned long int j = 1; j < ntot; j++)
          {     
                double fld1_tmp = fld_1d[j];
                if (fld1_tmp < minval) minval = fld1_tmp;
                if (fld1_tmp > maxval) maxval = fld1_tmp;
          }

        // Store the minimum value
        minval_vec[ilay] = minval;

        // Print min and max
        cout << "min=" << minval << " max=" << maxval << endl;

        // Quantize for a q-letter alphabet
        double deps = (maxval-minval)/(double)(q-1);

        // Impose the desired accuracy of the least significant bit plane
        if (deps < tolabs) 
        {
          deps = tolabs;
          brflag = 1;
        }

        // Save the quantization interval size
        deps_vec[ilay] = deps;

        // Two branches depending on the local precision mask activated or not
        if (mtot > 1)
          {
          // Loop for all jp in physical space
          for(unsigned long int jp = 0; jp < ntot; jp++)
            {   
              // Calculate 3D index in wavelet space
              int l, jwx, jwy, jwz;
              ind_p2w_3d( wlev, nx, ny, nz, jp%nx, (jp/nx)%ny, jp/nx/ny, &l, &jwx, &jwy, &jwz);

              // Uniform cutoff by default
              double precmask = tolabs;

              // Load precision cutoff function in physical space, only applied on smallest detail coefficients
              if (l <= LOC_CUTOFF_LVL) 
                precmask = tolabs/tolrel * lcl_prec(nx, ny, nz, jp%nx, (jp/nx)%ny, jp/nx/ny, mx, my, mz, cutoffvec);

              // 1D index in wavelet space
              unsigned long int jw = (unsigned long int)(jwx) + 
                                     (unsigned long int)(nx)*(unsigned long int)(jwy) +
                                     (unsigned long int)(nx)*(unsigned long int)(ny)*(unsigned long int)(jwz);
  
              // For every element, quantize depending on the local precision mask 
              if (maxval-minval < precmask) 
                { 
                  // Set to zero if the full range is less than the local precision
                  // Note that fld_q is in wavelet space
                  fld_q[jw] = 0;
                  fld_1d[jw] = minval; 
                } 
              else 
                { 
                  // Set to a value between 0 and 256 if the full range is greater than the local precision
                  fld_q[jw] = round( ( fld_1d[jw] - minval ) / deps );
                }
            }
          }
        // If the local precision mask is not activated, use faster loop
        else
          {
          // Loop for all jp in physical space
          for(unsigned long int jp = 0; jp < ntot; jp++)
            {   
              // Set to a value between 0 and 256 if the full range is greater than the local precision
              fld_q[jp] = round( ( fld_1d[jp] - minval ) / deps );
            }
          }   	 
    
        // The following only executes if all quantized data are non-zero
        // Print the resolution for this layer
        cout << "ilay=" << static_cast<unsigned>(ilay) << " deps=" << deps << endl;

        // Residual field
        for(unsigned long int j = 0; j < ntot; j++)
          fld_1d[j] = fld_1d[j] - ( fld_q[j]*deps + minval );

        // Check quantized data bounds
        unsigned char iminval = fld_q[0];
        unsigned char imaxval = fld_q[0];
        for(unsigned long int j = 1; j < ntot; j++)
          {
             unsigned char fld1_q_tmp = fld_q[j];
             if (fld1_q_tmp<iminval) iminval = fld1_q_tmp;
             if (fld1_q_tmp>imaxval) imaxval = fld1_q_tmp;
          }
        cout << "imin=" << static_cast<unsigned>(iminval) << " imax=" << static_cast<unsigned>(imaxval) << " med=" << fld_q[ntot/2UL]*deps + minval << endl;

        // Encode
        range_encode(fld_q,ntot,enc_q,len_out_q);

        // Store the encoded data
        len_enc_vec[ilay] = len_out_q;
        for(unsigned long int j = 0; j < len_out_q; j++) 
          {
            // Store the encoded data element
            data_enc[jtot++] = enc_q[j];

            // Check for overflow
            if (jtot > SAFETY_BUFFER_FACTOR*NLAYMAX*ntot)
              {
                cout << "Error: encoded array is too large. Use larger SAFETY_BUFFER_FACTOR" << endl;
                throw std::exception();
              }
          }

        // Print the encoded and original data set size
        cout << "len_out_q=" << len_out_q << " ntot=" << ntot << endl;

        // Update layer index
        ilay ++;

        // Stop when the required tolerance is reached for all samples
        if ( brflag ) 
          {
             // Stop iterations
             break;
          }
    }

    // Number of byte layers
    nlay = ilay;

    // Total encoded data array length
    ntot_enc = jtot;

    // Deallocate memory
    delete [] fld_q;
    delete [] enc_q;
}


/* Decoding subroutine with range decoding and inverse wavelet transform*/ 
void decoding_wrap(int nx, int ny, int nz, double *fld_1d, double& tolabs, double& midval, double& halfspanval, unsigned char& wlev, unsigned char& nlay, unsigned long int& ntot_enc, double *deps_vec, double *minval_vec, unsigned long int *len_enc_vec, unsigned char *data_enc)
{
    // Total number of elements
    unsigned long int ntot = (unsigned long int)(nx)*(unsigned long int)(ny)*(unsigned long int)(nz);

    // Case of trivial data
    if (ntot_enc == 0)
      {
        // Reconstruct
        for(unsigned long int j = 0; j < ntot; j++) fld_1d[j] = midval;

        // Exit 
        return;
      }

    /* Range decoding */
    // Print decoding status
    cout << "Range decoding..." << endl;

    // Allocate the quantized input and output vectors
    unsigned char *dec_q = new unsigned char[ntot];
    unsigned char *enc_q = new unsigned char[(SAFETY_BUFFER_FACTOR+1UL)*ntot]; // Encoded array may be longer than the original

    // Cumulative field
    for(unsigned long int j = 0; j < ntot; j++) fld_1d[j] = 0;

    // Input vector counter
    unsigned long int jtot = 0;

    // Output the decoded output vector   
    for (unsigned char ilay = 0; ilay < nlay; ilay++)
    {
        // Print layer index
        cout << "ilay=" << static_cast<unsigned>(ilay) << endl;

        // Parameters for reconstruction
        double deps = deps_vec[ilay];
        double minval = minval_vec[ilay];

        // Get the encoded data
        unsigned long int len_out_q = len_enc_vec[ilay];
        for(unsigned long int j = 0; j < len_out_q; j++) enc_q[j] = data_enc[jtot++];

        // Decode
        range_decode(enc_q,len_out_q,dec_q,ntot);

        // Check quantized data bounds
        unsigned char iminval = dec_q[0];
        unsigned char imaxval = dec_q[0];
        for(unsigned long int j = 1; j < ntot; j++)
          {
             if (dec_q[j]<iminval) iminval = dec_q[j];
             if (dec_q[j]>imaxval) imaxval = dec_q[j];
          }
        cout << "imin=" << static_cast<unsigned>(iminval) << " imax=" << static_cast<unsigned>(imaxval) <<  " med=" << dec_q[ntot/2UL]*deps + minval << endl;

        // Cumulative field
        for(unsigned long int j = 0; j < ntot; j++)
          fld_1d[j] = fld_1d[j] + ( dec_q[j]*deps + minval );
    }

    /* Wavelet reconstruction */
    // Print wavelet reconstruction status
    cout << "Wavelet reconstruction..." << endl;

    // Inverse wavelet transform if the data is non-trivial
    waveletcdf97_3d(nx,ny,nz,-int(wlev),fld_1d);

    // Deallocate memory
    delete [] enc_q;
    delete [] dec_q;
}


/* Use the range encoder to code an array fld_q */ 
void range_encode(unsigned char *fld_q, unsigned long int ntot, unsigned char *enc_q, unsigned long int& len_out_q)
{   freq counts[257], blocksize, i;
    int buffer[BLOCKSIZE];
    unsigned long int j;
    unsigned long int pos_in = 0;
    unsigned long int pos_out = 0;

    // Allocate a range coder object
    rangecoder *rc = (rangecoder*)malloc(sizeof(rangecoder));

    // Initialize data buffer
    init_databuf(rc,2*BLOCKSIZE+1000);

    // Start up the range coder, first byte 0, no header
    start_encoding(rc,0,0);

    // Coding: loop for all blocks of the input vector
    while (1)
    {
        // Put data in a buffer
        for (blocksize = 0; blocksize < BLOCKSIZE; blocksize++) {
          buffer[blocksize] = fld_q[pos_in];
          pos_in++; 
          if (pos_in>ntot) break;
        }

        // Block start marker
        encode_freq(rc,1,1,2);

        // Get the statistics 
        countblock(buffer,blocksize,counts);

        // Write the statistics.
        // Cant use putchar or other since we are after start of the rangecoder 
        // as you can see the rangecoder doesn't care where probabilities come 
        // from, it uses a flat distribution of 0..0xffff in encode_short. 
        for(i=0; i<256; i++)
            encode_short(rc,counts[i]);

        // Store in counters[i] the number of all bytes < i, so sum up 
        counts[256] = blocksize;
        for (i=256; i; i--)
            counts[i-1] = counts[i]-counts[i-1];

        // Output the encoded symbols 
        for(i=0; i<blocksize; i++) {
            int ch = buffer[i];
            encode_freq(rc,counts[ch+1]-counts[ch],counts[ch],counts[256]);
        }

        // Copy data from a buffer to the output vector
        for (j = 0; j < rc->datapos; j++) {
          enc_q[pos_out++] = rc->databuf[j];
        }

        // Reset the buffer    
        rc->datapos = 0;

        // Terminate if no more data
        if (blocksize<BLOCKSIZE) break;
    }

    // Flag absence of next block by a bit 
    encode_freq(rc,1,0,2);

    // Finalize the encoder 
    done_encoding(rc);

    // Copy data from a buffer to the output vector
    for (j = 0; j < rc->datapos; j++) {
      enc_q[pos_out++] = rc->databuf[j];
    }

    // True length of the encoded array
    len_out_q = pos_out;

    // Deallocate data buffer
    free_databuf(rc);

    // Deallocate range coder object
    free(rc);
}


/* Decode the range-encoded data enc_q */
void range_decode(unsigned char *enc_q, unsigned long int len_out_q, unsigned char *dec_q, unsigned long int ntot)
{   freq counts[257], blocksize, i, cf, symbol, middle, first, last;

    // Allocate a range coder object
    rangecoder *rc = (rangecoder*)malloc(sizeof(rangecoder));

    // Initialize data buffer
    rc->help = 0;
    unsigned long int pos_out = 0;
    init_databuf(rc,len_out_q);
    for(unsigned long int j = 0; j < len_out_q; j++) rc->databuf[j] = enc_q[j];
    rc->datalen = len_out_q;
    rc->datapos = 0;

    // Start decoding
    if (start_decoding(rc) != 0)
    {   fprintf(stderr,"could not successfully open input data\n");
        exit(1);
    }

    // Decoding: loop for all blocks of the input vector
    while (cf = decode_culfreq(rc,2))
    {   // Read the beginning of the block
        decode_update(rc,1,1,2);

        // Read frequencies
        readcounts(rc,counts);

        // Figure out blocksize by summing counts; also use counts as in encoder 
        blocksize = 0;
        for (i=0; i<256; i++)
        {   freq tmp = counts[i];
            counts[i] = blocksize;
            blocksize += tmp;
        }
        counts[256] = blocksize;

        for (i=0; i<blocksize; i++)
        {   // Decode frequency
            cf = decode_culfreq(rc,blocksize);
            // Figure out nearst symbol using a binary search
            first = 0;
            last = 256;
            middle = (first+last)/2;
            while (first <= last)
            {
              if(counts[middle] < cf)
              {
                first = middle + 1;
              }
              else if(counts[middle] == cf)
              {
                break;
              }
              else
              {
                last = middle - 1;
              }
              middle = (first + last)/2;
            }
            if(first > last)
            {
              middle = last;
            }
            // If some symbols have zero frequency, skip them
            for (symbol=middle; counts[symbol+1]<=cf; symbol++);
            // Update the decoder
            decode_update(rc, counts[symbol+1]-counts[symbol],counts[symbol],blocksize);
            // Store the decoded element
            dec_q[pos_out++] = symbol;
        }
    }

    // Finalize decoding 
    done_decoding(rc);

    // Deallocate data buffer
    free_databuf(rc);

    // Deallocate range coder object
    free(rc);
}


/*
    waveletcdf97_3d.c : This file is part of WaveRange CFD data compression utility

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

    Modified from the Matlab code by Pascal Getreuer (http://www.getreuer.info/home/waveletcdf97/)
*/

#include "../core/defs.h"
#include "waveletcdf97_3d.h"

/* Three-dimensional wavelet transform using CDF9/7 wavelets */
void waveletcdf97_3d(int N1in, int N2in, int N3in, int lvlin, double *X)
{
  // Lifting filter coefficients
  static const double lfc[4] = {-1.5861343420693648, -0.0529801185718856, 0.8829110755411875, 0.4435068520511142};

  // Scale factor
  static const double scl = 1.1496043988602418;

  // Indexes
  int k;
  unsigned long int i1, i2, i3, i, M1, M2, M3, M, N, Q;

  // Temporary vectors
  double *V, *V0, *V1;

  // Extrapolation coefficients for odd-sized arrays
  double ext[3];
  ext[0] = -2*lfc[0]*lfc[1]*lfc[2]/(1+2*lfc[1]*lfc[2]);
  ext[1] = -2*lfc[1]*lfc[2]/(1+2*lfc[1]*lfc[2]);
  ext[2] = -2*(lfc[0]+lfc[2]+3*lfc[0]*lfc[1]*lfc[2])/(1+2*lfc[1]*lfc[2]);

  // Initialize data size and level to input values
  unsigned long int N1 = (unsigned long int)(N1in);
  unsigned long int N2 = (unsigned long int)(N2in);
  unsigned long int N3 = (unsigned long int)(N3in);
  int lvl = lvlin;

  if (lvl >= 0)   
    // Forward transform
    {
      for (k = 1; k <= lvl; k++)
        {
          // Low-pass filtered vector length
          M1 = (N1/2UL) + ( (N1%2UL) > 0UL ? 1UL : 0UL );
          M2 = (N2/2UL) + ( (N2%2UL) > 0UL ? 1UL : 0UL );
          M3 = (N3/2UL) + ( (N3%2UL) > 0UL ? 1UL : 0UL );

          // Transform along the FIRST direction
          // At least two elements are required
          if (N1 > 1UL)
            {         
              // Array extents in the corresponding direction
              N = N1;
              M = M1;

              // Allocate temporary vectors
              V = malloc(N*sizeof(double));
              V0 = malloc(M*sizeof(double));
              V1 = malloc(M*sizeof(double));

              // Loop over the remaining two directions
              for (i3 = 0; i3 < N3; i3++)
                for (i2 = 0; i2 < N2; i2++)
                  {
                    // Place data elements in a contiguous vector
                    for (i1 = 0; i1 < N1; i1++) V[i1] = X[i1+N1in*i2+N1in*N2in*i3];

                    // Initialize low-pass and high-pass filtered vectors
                    for (i = 0; i < M; i++)
                      {
                        V0[i] = V[2UL*i];
                        if (2UL*i+1UL<N) V1[i] = V[2UL*i+1UL];
                      }

                    // Extrapolate if the vector length is odd
                    if (N%2UL) V1[M-1UL] = V0[M-2UL]*ext[0] + V1[M-2UL]*ext[1] + V0[M-1UL]*ext[2];

                    // Apply lifting stage 1
                    for (i = 0; i < M-1UL; i++) V1[i] += lfc[0]*(V0[i+1UL]+V0[i]);
                    V1[M-1UL] += lfc[0]*2*V0[M-1UL];

                    // Apply lifting stage 2
                    V0[0] += lfc[1UL]*2*V1[0UL];
                    for (i = 1; i < M; i++) V0[i] += lfc[1]*(V1[i]+V1[i-1UL]);

                    // Apply lifting stage 3
                    for (i = 0; i < M-1UL; i++) V1[i] += lfc[2]*(V0[i+1UL]+V0[i]);
                    V1[M-1UL] += lfc[2]*2*V0[M-1UL];

                    // Apply lifting stage 4
                    V0[0] += lfc[3]*2*V1[0];
                    for (i = 1; i < M; i++) V0[i] += lfc[3]*(V1[i]+V1[i-1UL]);

                    // Concatenate low-pass and high-pass vectors
                    for (i = 0; i < M; i++)
                      {
                        V[i] = V0[i]*scl;
                        if (2UL*i+1UL<N) V[i+M] = V1[i]/scl;
                      }

                    // Substitute the result in the 3D array
                    for (i1 = 0; i1 < N1; i1++) X[i1+N1in*i2+N1in*N2in*i3] = V[i1];
                  }

              // Deallocate arrays
              free(V1);
              free(V0);
              free(V); 
            }

          // Transform along the SECOND direction
          // At least two elements are required
          if (N2 > 1UL)
            {         
              // Array extents in the corresponding direction
              N = N2;
              M = M2;

              // Allocate temporary vectors
              V = malloc(N*sizeof(double));
              V0 = malloc(M*sizeof(double));
              V1 = malloc(M*sizeof(double));

              // Loop over the remaining two directions
              for (i3 = 0; i3 < N3; i3++)
                for (i1 = 0; i1 < N1; i1++)
                  {
                    // Place data elements in a contiguous vector
                    for (i2 = 0; i2 < N2; i2++) V[i2] = X[i1+N1in*i2+N1in*N2in*i3];

                    // Initialize low-pass and high-pass filtered vectors
                    for (i = 0; i < M; i++)
                      {
                        V0[i] = V[2UL*i];
                        if (2UL*i+1UL<N) V1[i] = V[2UL*i+1UL];
                      }

                    // Extrapolate if the vector length is odd
                    if (N%2UL) V1[M-1UL] = V0[M-2UL]*ext[0] + V1[M-2UL]*ext[1] + V0[M-1UL]*ext[2];

                    // Apply lifting stage 1
                    for (i = 0; i < M-1UL; i++) V1[i] += lfc[0]*(V0[i+1UL]+V0[i]);
                    V1[M-1] += lfc[0]*2*V0[M-1UL];

                    // Apply lifting stage 2
                    V0[0] += lfc[1]*2*V1[0];
                    for (i = 1; i < M; i++) V0[i] += lfc[1]*(V1[i]+V1[i-1UL]);

                    // Apply lifting stage 3
                    for (i = 0; i < M-1UL; i++) V1[i] += lfc[2]*(V0[i+1UL]+V0[i]);
                    V1[M-1] += lfc[2]*2*V0[M-1UL];

                    // Apply lifting stage 4
                    V0[0] += lfc[3]*2*V1[0];
                    for (i = 1; i < M; i++) V0[i] += lfc[3]*(V1[i]+V1[i-1UL]);

                    // Concatenate low-pass and high-pass vectors
                    for (i = 0; i < M; i++)
                      {
                        V[i] = V0[i]*scl;
                        if (2UL*i+1UL<N) V[i+M] = V1[i]/scl;
                      }

                    // Substitute the result in the 3D array
                    for (i2 = 0; i2 < N2; i2++) X[i1+N1in*i2+N1in*N2in*i3] = V[i2];
                  }

              // Deallocate arrays
              free(V1);
              free(V0);
              free(V); 
            }

          // Transform along the THIRD direction
          // At least two elements are required
          if (N3 > 1UL)
            {         
              // Array extents in the corresponding direction
              N = N3;
              M = M3;

              // Allocate temporary vectors
              V = malloc(N*sizeof(double));
              V0 = malloc(M*sizeof(double));
              V1 = malloc(M*sizeof(double));

              // Loop over the remaining two directions
              for (i2 = 0; i2 < N2; i2++)
                for (i1 = 0; i1 < N1; i1++)
                  {
                    // Place data elements in a contiguous vector
                    for (i3 = 0; i3 < N3; i3++) V[i3] = X[i1+N1in*i2+N1in*N2in*i3];

                    // Initialize low-pass and high-pass filtered vectors
                    for (i = 0; i < M; i++)
                      {
                        V0[i] = V[2UL*i];
                        if (2UL*i+1UL<N) V1[i] = V[2UL*i+1UL];
                      }

                    // Extrapolate if the vector length is odd
                    if (N%2UL) V1[M-1UL] = V0[M-2UL]*ext[0] + V1[M-2UL]*ext[1] + V0[M-1UL]*ext[2];

                    // Apply lifting stage 1
                    for (i = 0; i < M-1UL; i++) V1[i] += lfc[0]*(V0[i+1UL]+V0[i]);
                    V1[M-1UL] += lfc[0]*2*V0[M-1UL];

                    // Apply lifting stage 2
                    V0[0] += lfc[1]*2*V1[0];
                    for (i = 1; i < M; i++) V0[i] += lfc[1]*(V1[i]+V1[i-1UL]);

                    // Apply lifting stage 3
                    for (i = 0; i < M-1UL; i++) V1[i] += lfc[2]*(V0[i+1UL]+V0[i]);
                    V1[M-1UL] += lfc[2]*2*V0[M-1UL];

                    // Apply lifting stage 4
                    V0[0] += lfc[3]*2*V1[0];
                    for (i = 1; i < M; i++) V0[i] += lfc[3]*(V1[i]+V1[i-1UL]);

                    // Concatenate low-pass and high-pass vectors
                    for (i = 0; i < M; i++)
                      {
                        V[i] = V0[i]*scl;
                        if (2UL*i+1UL<N) V[i+M] = V1[i]/scl;
                      }

                    // Substitute the result in the 3D array
                    for (i3 = 0; i3 < N3; i3++) X[i1+N1in*i2+N1in*N2in*i3] = V[i3];
                  }

              // Deallocate arrays
              free(V1);
              free(V0);
              free(V); 
            }

          // Assign the subset array extents for the next iteration
          N1 = M1;
          N2 = M2;
          N3 = M3;
        }
    }
  else           
    // Inverse transform
    {
      for (k = 1+lvl; k <= 0; k++)
        {
          unsigned long int pow2k = 1;
          int kk;
          for (kk = 1; kk <= -k; kk++) pow2k *= 2UL;
          M1 = (N1/pow2k) + ( (N1%pow2k) > 0UL ? 1UL : 0UL );
          M2 = (N2/pow2k) + ( (N2%pow2k) > 0UL ? 1UL : 0UL );
          M3 = (N3/pow2k) + ( (N3%pow2k) > 0UL ? 1UL : 0UL );

          // Inverse transform along the THIRD direction
          // At least two elements are required
          if (M3 > 1UL)
            {         
              // Array extents in the corresponding direction
              M = M3;
              Q = (M/2UL) + ( (M%2UL) > 0UL ? 1UL : 0UL );

              // Allocate temporary vectors
              V = malloc(M*sizeof(double));
              V0 = malloc(Q*sizeof(double));
              V1 = malloc(Q*sizeof(double));

              // Loop over the remaining two directions
              for (i2 = 0; i2 < M2; i2++)
                for (i1 = 0; i1 < M1; i1++)
                  {
                    // Place data elements in a contiguous vector
                    for (i3 = 0; i3 < M3; i3++) V[i3] = X[i1+N1in*i2+N1in*N2in*i3];

                    // Initialize low-pass and high-pass filtered vectors
                    for (i = 0; i < Q; i++) V0[i] = V[i]/scl;
                    for (i = 0; i < M-Q; i++) V1[i] = V[i+Q]*scl;
                    if (M%2UL) V1[Q-1UL] = 0; 

                    // Apply lifting stage 1
                    V0[0] -= lfc[3]*2*V1[0];
                    for (i = 1; i < Q; i++) V0[i] -= lfc[3]*(V1[i]+V1[i-1UL]);

                    // Apply lifting stage 2
                    for (i = 0; i < Q-1UL; i++) V1[i] -= lfc[2]*(V0[i+1UL]+V0[i]);
                    V1[Q-1UL] -= lfc[2]*2*V0[Q-1UL];

                    // Apply lifting stage 3
                    V0[0] -= lfc[1]*2*V1[0];
                    for (i = 1; i < Q; i++) V0[i] -= lfc[1]*(V1[i]+V1[i-1UL]);

                    // Apply lifting stage 4
                    for (i = 0; i < Q-1UL; i++) V1[i] -= lfc[0]*(V0[i+1UL]+V0[i]);
                    V1[Q-1UL] -= lfc[0]*2*V0[Q-1UL];

                    // Concatenate low-pass and high-pass vectors
                    for (i = 0; i < Q; i++)
                      {
                        V[2UL*i] = V0[i];
                        if (2UL*i+1UL<M) V[2UL*i+1UL] = V1[i];
                      }

                    // Substitute the result in the 3D array
                    for (i3 = 0; i3 < M3; i3++) X[i1+N1in*i2+N1in*N2in*i3] = V[i3];
                  }

              // Deallocate arrays
              free(V1);
              free(V0);
              free(V); 
            }

          // Inverse transform along the SECOND direction
          // At least two elements are required
          if (M2 > 1UL)
            {         
              // Array extents in the corresponding direction
              M = M2;
              Q = (M/2UL) + ( (M%2UL) > 0UL ? 1UL : 0UL );

              // Allocate temporary vectors
              V = malloc(M*sizeof(double));
              V0 = malloc(Q*sizeof(double));
              V1 = malloc(Q*sizeof(double));

              // Loop over the remaining two directions
              for (i3 = 0; i3 < M3; i3++)
                for (i1 = 0; i1 < M1; i1++)
                  {
                    // Place data elements in a contiguous vector
                    for (i2 = 0; i2 < M2; i2++) V[i2] = X[i1+N1in*i2+N1in*N2in*i3];

                    // Initialize low-pass and high-pass filtered vectors
                    for (i = 0; i < Q; i++) V0[i] = V[i]/scl;
                    for (i = 0; i < M-Q; i++) V1[i] = V[i+Q]*scl;
                    if (M%2UL) V1[Q-1UL] = 0; 

                    // Apply lifting stage 1
                    V0[0] -= lfc[3]*2*V1[0];
                    for (i = 1; i < Q; i++) V0[i] -= lfc[3]*(V1[i]+V1[i-1UL]);

                    // Apply lifting stage 2
                    for (i = 0; i < Q-1UL; i++) V1[i] -= lfc[2]*(V0[i+1UL]+V0[i]);
                    V1[Q-1] -= lfc[2]*2*V0[Q-1];

                    // Apply lifting stage 3
                    V0[0] -= lfc[1]*2*V1[0];
                    for (i = 1; i < Q; i++) V0[i] -= lfc[1]*(V1[i]+V1[i-1UL]);

                    // Apply lifting stage 4
                    for (i = 0; i < Q-1UL; i++) V1[i] -= lfc[0]*(V0[i+1UL]+V0[i]);
                    V1[Q-1UL] -= lfc[0]*2*V0[Q-1UL];

                    // Concatenate low-pass and high-pass vectors
                    for (i = 0; i < Q; i++)
                      {
                        V[2UL*i] = V0[i];
                        if (2UL*i+1UL<M) V[2UL*i+1UL] = V1[i];
                      }

                    // Substitute the result in the 3D array
                    for (i2 = 0; i2 < M2; i2++) X[i1+N1in*i2+N1in*N2in*i3] = V[i2];
                  }

              // Deallocate arrays
              free(V1);
              free(V0);
              free(V); 
            }

          // Inverse transform along the FIRST direction
          // At least two elements are required
          if (M1 > 1UL)
            {         
              // Array extents in the corresponding direction
              M = M1;
              Q = (M/2UL) + ( (M%2UL) > 0UL ? 1UL : 0UL );

              // Allocate temporary vectors
              V = malloc(M*sizeof(double));
              V0 = malloc(Q*sizeof(double));
              V1 = malloc(Q*sizeof(double));

              // Loop over the remaining two directions
              for (i3 = 0; i3 < M3; i3++)
                for (i2 = 0; i2 < M2; i2++)
                  {
                    // Place data elements in a contiguous vector
                    for (i1 = 0; i1 < M1; i1++) V[i1] = X[i1+N1in*i2+N1in*N2in*i3];

                    // Initialize low-pass and high-pass filtered vectors
                    for (i = 0; i < Q; i++) V0[i] = V[i]/scl;
                    for (i = 0; i < M-Q; i++) V1[i] = V[i+Q]*scl;
                    if (M%2UL) V1[Q-1UL] = 0; 

                    // Apply lifting stage 1
                    V0[0] -= lfc[3]*2*V1[0];
                    for (i = 1; i < Q; i++) V0[i] -= lfc[3]*(V1[i]+V1[i-1UL]);

                    // Apply lifting stage 2
                    for (i = 0; i < Q-1UL; i++) V1[i] -= lfc[2]*(V0[i+1UL]+V0[i]);
                    V1[Q-1UL] -= lfc[2]*2*V0[Q-1UL];

                    // Apply lifting stage 3
                    V0[0] -= lfc[1]*2*V1[0];
                    for (i = 1; i < Q; i++) V0[i] -= lfc[1]*(V1[i]+V1[i-1UL]);

                    // Apply lifting stage 4
                    for (i = 0; i < Q-1UL; i++) V1[i] -= lfc[0]*(V0[i+1UL]+V0[i]);
                    V1[Q-1UL] -= lfc[0]*2*V0[Q-1UL];

                    // Concatenate low-pass and high-pass vectors
                    for (i = 0; i < Q; i++)
                      {
                        V[2UL*i] = V0[i];
                        if (2UL*i+1UL<M) V[2UL*i+1UL] = V1[i];
                      }

                    // Substitute the result in the 3D array
                    for (i1 = 0; i1 < M1; i1++) X[i1+N1in*i2+N1in*N2in*i3] = V[i1];
                  }

              // Deallocate arrays
              free(V1);
              free(V0);
              free(V); 
            }
        }
    }
}



/* Convert 3D index from physical space to wavelet space */
void ind_p2w_3d( int lvlin, int N1in, int N2in, int N3in, int i1in, int i2in, int i3in, int *lvl, int *i1, int *i2, int *i3 )
{
  // Indexes
  int k, M1, M2, M3;

  // Initialize data size and level to input values
  int N1 = N1in;
  int N2 = N2in;
  int N3 = N3in;
  *lvl = 0;
  *i1 = i1in;
  *i2 = i2in;
  *i3 = i3in;

  int chlvl = 0;

  if (lvlin >= 0)   
    // Forward transform
    {
      for (k = 1; k <= lvlin; k++)
        {
          // Low-pass filtered vector length
          M1 = (N1/2) + ( (N1%2) > 0 ? 1 : 0 );
          M2 = (N2/2) + ( (N2%2) > 0 ? 1 : 0 );
          M3 = (N3/2) + ( (N3%2) > 0 ? 1 : 0 );

          // Transform along the FIRST direction
          // At least two elements are required
          if (N1 > 1)
            {         
              // If the actual point index is inside the low-pass quadrant
              if ( (*i3 < N3) && (*i2 < N2) && (*i1 < N1) )
                  {
                    // Update the index
                    if (*i1%2) *i1 = *i1/2+M1; else *i1 /= 2;

                    // Update level flag
                    chlvl = 1;
                  }
            }

          // Transform along the SECOND direction
          // At least two elements are required
          if (N2 > 1)
            {         
              // If the actual point index is inside the low-pass quadrant
              if ( (*i3 < N3) && (*i2 < N2) && (*i1 < N1) )
                  {
                    // Update the index
                    if (*i2%2) *i2 = *i2/2+M2; else *i2 /= 2;

                    // Update level flag
                    chlvl = 1;
                  }
            }

          // Transform along the THIRD direction
          // At least two elements are required
          if (N3 > 1)
            {         
              // If the actual point index is inside the low-pass quadrant
              if ( (*i3 < N3) && (*i2 < N2) && (*i1 < N1) )
                  {
                    // Update the index
                    if (*i3%2) *i3 = *i3/2+M3; else *i3 /= 2;

                    // Update level flag
                    chlvl = 1;
                  }
            }

          // Assign the subset array extents for the next iteration
          N1 = M1;
          N2 = M2;
          N3 = M3;

          // Update the actual level
          if (chlvl) *lvl += 1;
        }
    }
}


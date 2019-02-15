/*
    waveletcdf97_3d.h : This file is part of WaveRange CFD data compression utility

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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Three-dimensional wavelet transform using CDF9/7 wavelets */
void waveletcdf97_3d(int N1in, int N2in, int N3in, int lvlin, double *X);

/* Convert 3D index from physical space to wavelet space */
void ind_p2w_3d(int lvlin, int N1in, int N2in, int N3in, int i1in, int i2in, int i3in, int *lvl, int *i1, int *i2, int *i3);

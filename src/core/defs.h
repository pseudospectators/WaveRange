/*
    defs.h : This file is part of WaveRange CFD data compression utility

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

/***** Constant parameters *****/
/* Version of the coder. Format: XYYZZ, where X is MSSG coder ver., YY is FluSI coder ver., ZZ is other */
#define CODER_VERSION 31400
/* Range coder block size. Must be less than 1<<16 */
#define BLOCKSIZE 60000
/* Maximum number of bit planes is 8 for 64-bit real type */
#define NLAYMAX 8UL
/* Uniform cutoff flag (0: false; 1: true) */
#define UNIFORM_CUTOFF 1
/* Number of top levels of wavelet decomposition subject to local (degraded precision) cutoff */
#define LOC_CUTOFF_LVL 1
/* Downscaling block size for non-uniform cutoff */
#define DS_BLOCK 16
/* Encoded array must be no larger than the input data size multiplied by this factor */
#define SAFETY_BUFFER_FACTOR 1UL
/* Wavelet reconstruction roundoff error correction coefficient = Linf_error/Linf_tolerance */
#define WAV_ACC_COEF 1.75
/* Maximum depth of wavelet transform */
#define WAV_LVL 4
/* Maximum number of datasets in a restart file */
#define NDSMAX 50
/* Number of digits in MSSG output file name extension */
#define MSSG_FILE_DIG 4
/* Number of elements in MSSG time record */
#define MSSG_TIME_REC_LEN 15
/* Relative tolerance of mask (indicator) function compression */
#define MSSG_MASK_TOLREL 0.126
/* Mask (indicator) function threshold value */
#define MSSG_MASK_THRESHOLD -1e-20

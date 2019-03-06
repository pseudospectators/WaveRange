!
!    example_fort.f90 : This file is part of WaveRange CFD data compression utility
!
!    Copyright (C) 2017-2019  Dmitry Kolomenskiy
!    Copyright (C) 2017-2019  Ryo Onishi
!    Copyright (C) 2017-2019  JAMSTEC
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.
!  
!    Reference:
!    doc/cfdproc2017.pdf
!    Dmitry Kolomenskiy, Ryo Onishi and Hitoshi Uehara "Wavelet-Based Compression of CFD Big Data"
!    Proceedings of the 31st Computational Fluid Dynamics Symposium, Kyoto, December 12-14, 2017
!    Paper No. C08-1
!
!    This work is supported by the FLAGSHIP2020, MEXT within the priority study4 
!    (Advancement of meteorological and global environmental predictions utilizing 
!    observational “Big Data”).
!

program example_fort

   ! This program is an example of calling WaveRange encoder and decoder 
   ! functions from a Fortran program. The original floating-point data
   ! are contained in a 3D array fld_ini(1:nx,1:ny,1:nz).
   ! After compression using encoding_wrap_f, all necessary information for 
   ! reconstruction is stored in the following variables: nx, ny, nz,   
   ! midval, halfspanval, wlev, nlay, ntot_enc, deps_vec(1:nlaymax), 
   ! minval_vec(1:nlaymax), len_enc_vec(1:nlaymax), data_enc(1:ntot_enc_max),
   ! where the suitable array bounds nlaymax and ntot_enc_max are
   ! determined by calling setup_wr_f.
   ! The reconstructed 3D array fld_rec(1:nx,1:ny,1:nz) is calculated
   ! by calling decoding_wrap_f. It differs from fld_ini by approximately
   ! tolrel*maxval(dabs(fld_ini))

   implicit none
   byte :: wlev, nlay
   byte, allocatable :: data_enc(:) 
   integer :: ix, iy, iz, nx, ny, nz, nlaymax, wtflag
   integer*8 :: ntot_enc_max,ntot_enc
   integer*8, allocatable :: len_enc_vec(:)
   double precision :: tolrel, midval, halfspanval, tolabs, &
           err_linf_abs, refval, err_linf_rel
   double precision, allocatable :: deps_vec(:), minval_vec(:), &
           fld_ini(:,:,:), fld_rec(:,:,:), fld_tmp(:,:,:)

   ! Print a message on startup
   print *, "FORTRAN EXAMPLE: SETTING UP THE PARAMETERS"
   
   ! Set up the dataset size
   nx = 64
   ny = 64
   nz = 64

   ! Set up the relative tolerance of lossy compression
   tolrel = 1e-6

   ! Wavelet transform (0:deactivate; 1:activate)
   wtflag = 1

   ! Determine the maximum number of bit planes nlaymax
   ! and the maximum length of the encoded data array ntot_enc_max
   call setup_wr_f(nx,ny,nz,nlaymax,ntot_enc_max)

   ! Allocate memory for all allocatable arrays
   allocate(fld_ini(nx,ny,nz),fld_rec(nx,ny,nz),fld_tmp(nx,ny,nz), &
           deps_vec(nlaymax),minval_vec(nlaymax), &
           len_enc_vec(nlaymax),data_enc(ntot_enc_max))

   ! Fill the 3D field with some values
   do iz = 1,nz
     do iy = 1,ny
       do ix = 1,nx
         fld_ini(ix,iy,iz) = 10.0d0* &
            dsin(dble(ix-1)/dble(nx)) * &
            dsin(dble(iy-1)/dble(ny))**2 * &
            dcos(dble(iz-1)/dble(nz))
       end do
     end do
   end do

   ! Print a message on compression
   print *, "FORTRAN EXAMPLE: COMPRESSING"

   ! The encoder overwrites the elements of its input array, therefore,
   ! the initial data elements are copied in a temporary array
   fld_tmp(:,:,:) = fld_ini(:,:,:)

   ! Apply data compression, using the following input data:
   ! nx,ny,nz,fld,wtflag,tolrel
   call encoding_wrap_f(nx,ny,nz,fld_tmp,wtflag,tolrel,tolabs, &
           midval,halfspanval,wlev,nlay,ntot_enc,deps_vec, &
           minval_vec,len_enc_vec,data_enc)

   ! Here, the output variables can be written in a file:
   ! tolabs,midval,halfspanval,wlev,nlay,ntot_enc,deps_vec,minval_vec,
   ! len_enc_vec,data_enc
   ! It may also be useful to include nx,ny,nz,nlaymax,ntot_enc_max

   ! ... add I/O statements here ...

   ! Print a message on reconstruction
   print *, "FORTRAN EXAMPLE: RECONSTRUCTING"

   ! Reconstruct the 3D field from the compressed format. Here, the input is:
   ! nx,ny,nz,midval,halfspanval,wlev,nlay,ntot_enc,deps_vec,minval_vec,
   ! len_enc_vec,data_enc
   ! The output is fld_rec
   call decoding_wrap_f(nx,ny,nz,fld_rec,midval,halfspanval,wlev, &
           nlay,ntot_enc,deps_vec,minval_vec,len_enc_vec,data_enc)

   ! Compare the original and the reconstructed field
   err_linf_abs = maxval(dabs(fld_rec-fld_ini))
   refval = maxval(dabs(fld_ini))
   err_linf_rel = err_linf_abs / refval
   print *, "FORTRAN EXAMPLE: Absolute error = ", err_linf_abs, &
           " Relative error = ", err_linf_rel, &
           " as normalized by ", refval

   ! Deallocate all allocatable arrays
   deallocate(fld_ini,fld_rec,fld_tmp,deps_vec,minval_vec, &
           len_enc_vec,data_enc)

   ! Print a message on exit
   print *, "FORTRAN EXAMPLE: EXIT"

end program example_fort


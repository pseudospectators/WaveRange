!
!    create_in_field.f90 : This file is part of WaveRange CFD data compression utility
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

program create_in_field

   ! This program produces an example input file for compression using 
   ! the 'generic' interface of WaveRange. It fills two three-dimensional 
   ! arrays with nontrivial real values, writes them in a single file 
   ! using unformatted output, then writes one more real number in the end 
   ! of the file before closing it. 

   implicit none
   integer :: ix, iy, iz, ifld, nx1, ny1, nz1, nx2, ny2, nz2, fln
   double precision, allocatable :: fld1(:,:,:), fld2(:,:,:)
   real :: fld3

   ! Print a message on startup
   print *, "DATASET GENERATOR: SETTING UP THE PARAMETERS"
   
   ! Set up the data size
   nx1 = 32
   ny1 = 32
   nz1 = 32
   nx2 = 64
   ny2 = 64
   nz2 = 64
   
   ! Allocate memory for all allocatable arrays
   allocate(fld1(nx1,ny1,nz1),fld2(nx2,ny2,nz2))

   ! Fill the 3D fields with some values
   do iz = 1,nz1
     do iy = 1,ny1
       do ix = 1,nx1
         fld1(ix,iy,iz) = 10.0d0* &
            dsin(dble(ix-1)/dble(nx1)) * &
            dsin(dble(iy-1)/dble(ny1))**2 * &
            dcos(dble(iz-1)/dble(nz1))
       end do
     end do
   end do
   do iz = 1,nz2
     do iy = 1,ny2
       do ix = 1,nx2
         fld2(ix,iy,iz) = 10.0d0* &
            dcos(dble(ix-1)/dble(nx2)) * &
            dsin(dble(iy-1)/dble(ny2))**2 * &
            dsin(dble(iz-1)/dble(nz2))
       end do
     end do
   end do
   fld3 = 3.14
   
   ! Binary file output
   print *, "DATASET GENERATOR: WRITING IN AN UNFORMATTED FILE"
   fln = 90
   open(fln,file='data.bin',form='unformatted')
   write(fln) (((fld1(ix,iy,iz),ix=1,nx1),iy=1,ny1),iz=1,nz1)
   write(fln) (((fld2(ix,iy,iz),ix=1,nx2),iy=1,ny2),iz=1,nz2)
   write(fln) fld3
   close(fln)   
   
   ! Deallocate all allocatable arrays
   deallocate(fld1,fld2)

   ! Print a message on exit
   print *, "DATASET GENERATOR: EXIT"

end program create_in_field

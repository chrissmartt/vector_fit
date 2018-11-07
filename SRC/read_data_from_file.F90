! This file forms part of the VECTOR_FIT project
!
! Simple Implementation of the Vector fitting process
!
! B. Gustavsen and A. Semlyen, "Rational Approximation 
! of Frequency Domain Responses by Vector Fitting", IEEE 
! Trans. Power Delivery, vol. 14, No. 3, July 1999, pp. 1052-1061.
!
! Copyright (C) 2018 University of Nottingham
!
! VECTOR_FIT is free software: you can redistribute it and/or modify it under the 
! terms of the GNU General Public License as published by the Free Software 
! Foundation, either version 3 of the License, or (at your option) any later 
! version.
! 
! VECTOR_FIT is distributed in the hope that it will be useful, but 
! WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
! or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License 
! for more details.
! 
! A copy of the GNU General Public License version 3 can be found in the 
! file COPYING.txt in the root or at <http://www.gnu.org/licenses/>.
! 
! VECTOR_FIT uses the EISPACK library. EISPACK is subject to 
! the GNU Lesser General Public License. A copy of the GNU Lesser General Public 
! License version can be found in the file COPPYING.LESSER.txt 
! or at <http://www.gnu.org/licenses/>.
! 
! The University of Nottingham can be contacted at: ggiemr@nottingham.ac.uk
!
! Author C Smartt
!
SUBROUTINE read_data_from_file(filename,read_loop,nf,fmin,fmax,s,f)

IMPLICIT NONE

character*256 :: filename
integer :: nf
complex*16 :: f(nf)
complex*16 :: s(nf)
real*8 :: fmin,fmax

! local variables

integer :: read_loop
real*8 ::f_in,real_value_in,imag_value_in

real*8,parameter     :: pi=3.1415926535897931d0
complex*16,parameter :: j=(0d0,1d0)

! START

! if read_loop=1 then count the number of samples
! if read_loop=2 then read the data into the frequency and function value arrays

! READ FILES IN TWO STAGES, FIRST GET THE NUMBER OF SAMPLES THEN ALLOCATE ARRAYS AND READ THE DATA
     
  OPEN(unit=10,file=filename,status='OLD',err=9000)
  nf=0
  fmin=1D30
  fmax=0d0
   
10 CONTINUE

    read(10,*,end=1000)f_in,real_value_in,imag_value_in
   
    nf=nf+1

    if (read_loop.eq.2) then

      s(nf)=j*2d0*pi*f_in
      f(nf)=dcmplx(real_value_in,imag_value_in)
      
      fmin=min(fmin,f_in)
      fmax=max(fmax,f_in)
     
    end if ! read_loop.EQ.2
   
    GOTO 10  ! read next line of data
   
1000 CONTINUE

  CLOSE(unit=10)   

  RETURN
  
9000 write(*,*)'Error opening file:',trim(filename)
  STOP
  
END SUBROUTINE read_data_from_file

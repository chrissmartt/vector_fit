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
SUBROUTINE calculate_error(nf,f,f_fit,s,error)

IMPLICIT NONE

integer :: nf
complex*16 :: f(nf)
complex*16 :: f_fit(nf)
complex*16 :: s(nf)

real*8     :: error

! local variables

integer :: n

! START


! Calculate the error in this approximation

  error=0d0

  do n=1,nf
      
    error=error+abs((f(n)-f_fit(n))**2)
  
  end do ! next frequency
  
  error=sqrt(error/nf)

END SUBROUTINE calculate_error

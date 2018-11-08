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
SUBROUTINE write_coefficients(order,a,c,c_hat,d,h,fmin,fmax,nf)

IMPLICIT NONE

integer :: order
complex*16 :: a(order)
complex*16 :: c(order)
complex*16 :: c_hat(order)
real*8 :: d
real*8 :: h
real*8 :: fmin,fmax
integer :: nf

! local variables

integer :: i
real*8 ::wnorm
real*8,parameter     :: pi=3.1415926535897931d0

! START

  write(*,*)
  write(*,*)'d=',real(d)
  write(*,*)'h=',real(h)
  write(*,*)'               pole (a)                          residue (c)'
  do i=1,order
    write(*,*)cmplx(a(i)),cmplx(c(i))
  end do 
  write(*,*)
  
  open(unit=50,file='Vfit.filter')
  
  wnorm=2d0*pi*fmax

  write(50,'(A)')'Vfit filter output'
  write(50,*)order,' # order'
  write(50,*)wnorm,' # wnorm'
  write(50,*)d,' # d'
  write(50,*)h*wnorm,' # h'
  write(50,*)'                       pole (a)                                         residue (c)'
  do i=1,order
    write(50,*)dble(a(i)/wnorm),imag(a(i)/wnorm),dble(c(i)/wnorm),imag(c(i)/wnorm)
  end do 
  
  write(50,*)'            wnorm_min             wnorm_max                  nw'
  write(50,*)fmin/fmax, 1d0 ,nf
  
  close(unit=50)
  

END SUBROUTINE write_coefficients

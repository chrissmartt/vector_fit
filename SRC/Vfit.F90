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

PROGRAM Vfit


IMPLICIT NONE

! parameters to set up a particular test case

integer :: test_case
integer :: order
logical :: include_const_term
logical :: include_s_term
logical :: log_freq_samples
integer :: nf
logical :: log_freq_starting_poles
logical :: complex_starting_poles
integer :: max_iterations

! General solution parameters

real*8  :: min_error=1D-6
!logical :: verbose=.TRUE.
logical :: verbose=.FALSE.
logical :: use_SVD_invert=.FALSE.

! Numerical constants

real*8,parameter     :: pi=3.1415926535897931d0
complex*16,parameter :: j=(0d0,1d0)

! Solution variables

real*8 :: freq
real*8 :: fmin,fmax,fstep_poles
real*8 :: log_fmin,log_fmax,log_fstep_poles

complex*16,allocatable :: f(:)
complex*16,allocatable :: s(:)

complex*16,allocatable :: f_fit(:)

complex*16,allocatable :: new_a(:)
complex*16,allocatable :: a(:)
complex*16,allocatable :: c(:)
complex*16,allocatable :: c_hat(:)
real*8 :: d
real*8 :: h

real*8     :: error

integer :: n
integer :: iteration
integer :: npairs,pole
real*8 :: alpha,beta

integer :: i,read_loop

character*256 :: filename

! START

write(*,*)'Enter the filename for the complex frequency domain data or 1-4 to run one of the internal test cases'
read(*,'(A)')filename
if (filename(1:1).Eq.'1') then
  test_case=1
else if (filename(1:1).Eq.'2') then
  test_case=2
else if (filename(1:1).Eq.'3') then
  test_case=3
else if (filename(1:1).Eq.'4') then
  test_case=4
else
  test_case=0
end if

if(verbose) write(*,*)'CALLED Vfit'

if    (test_case.EQ.1) then
  order=3
  log_freq_samples=.TRUE.
  nf=101
  include_const_term=.TRUE.
  include_s_term=.TRUE.
  log_freq_starting_poles=.TRUE.
  complex_starting_poles=.FALSE.
  max_iterations=1
  ALLOCATE( f(nf) )
  ALLOCATE( s(nf) )
  ALLOCATE( f_fit(nf) )
  ALLOCATE( new_a(order) )
  ALLOCATE( a(order) )
  ALLOCATE( c(order) )
  ALLOCATE( c_hat(order) )  
  CALL set_test_set_1(nf,fmin,fmax,s,f,log_freq_samples)
  
else if (test_case.EQ.2) then
  order=18
  log_freq_samples=.FALSE.
  nf=100
  include_const_term=.TRUE.
  include_s_term=.TRUE.
  log_freq_starting_poles=.FALSE.
  complex_starting_poles=.TRUE.
  max_iterations=1
  ALLOCATE( f(nf) )
  ALLOCATE( s(nf) )
  ALLOCATE( f_fit(nf) )
  ALLOCATE( new_a(order) )
  ALLOCATE( a(order) )
  ALLOCATE( c(order) )
  ALLOCATE( c_hat(order) )  
  CALL   set_test_set_2(nf,fmin,fmax,s,f,log_freq_samples)
  
else if (test_case.EQ.3) then
  order=5
  log_freq_samples=.FALSE.
  nf=100
  include_const_term=.TRUE.
  include_s_term=.TRUE.
  log_freq_starting_poles=.FALSE.
  complex_starting_poles=.TRUE.
  max_iterations=10
  ALLOCATE( f(nf) )
  ALLOCATE( s(nf) )
  ALLOCATE( f_fit(nf) )
  ALLOCATE( new_a(order) )
  ALLOCATE( a(order) )
  ALLOCATE( c(order) )
  ALLOCATE( c_hat(order) )  
  CALL   set_test_set_3(nf,fmin,fmax,s,f,log_freq_samples)
  
else if (test_case.EQ.4) then
  order=3
  log_freq_samples=.FALSE.
  nf=500
  include_const_term=.TRUE.
  include_s_term=.TRUE.
  log_freq_starting_poles=.FALSE.
  complex_starting_poles=.TRUE.
  max_iterations=10
  ALLOCATE( f(nf) )
  ALLOCATE( s(nf) )
  ALLOCATE( f_fit(nf) )
  ALLOCATE( new_a(order) )
  ALLOCATE( a(order) )
  ALLOCATE( c(order) )
  ALLOCATE( c_hat(order) )  
  CALL   set_test_set_4(nf,fmin,fmax,s,f,log_freq_samples)
  
else

! Read complex frequency domain data from a file
  read_loop=1
  CALL read_data_from_file(filename,read_loop,nf,fmin,fmax,s,f)
  ALLOCATE( f(nf) )
  ALLOCATE( s(nf) )
  ALLOCATE( f_fit(nf) )
  
  read_loop=2
  CALL read_data_from_file(filename,read_loop,nf,fmin,fmax,s,f)
  
  write(*,*)'Number of frequency domain samples read=',nf
  write(*,*)'fmin=',fmin,' fmax=',fmax
  
  write(*,*)'Enter the order of the vector fit solution'
  read(*,*)order
  ALLOCATE( new_a(order) )
  ALLOCATE( a(order) )
  ALLOCATE( c(order) )
  ALLOCATE( c_hat(order) )  
  
  write(*,*)'Enter the maximum number of iterations of the vector fit solution'
  read(*,*)max_iterations
  
  include_const_term=.TRUE.
  include_s_term=.TRUE.
  log_freq_starting_poles=.FALSE.
  complex_starting_poles=.TRUE.
  
end if

! Generate the staring poles

log_fmin=log10(fmin)
log_fmax=log10(fmax)

if (.NOT.complex_starting_poles) then

! Real starting poles
  
  if(order.GT.1) then
    fstep_poles=(fmax-fmin)/(order-1)
    log_fstep_poles=(log_fmax-log_fmin)/(order-1)
  else
    fstep_poles=0d0
    log_fstep_poles=0d0
  end if
   
  do i=1,order
  
    if (log_freq_starting_poles) then
      freq=10d0**( (log_fmin+(i-1)*log_fstep_poles) )
    else
      freq=fmin+(i-1)*fstep_poles
    end if
  
    a(i)=-2d0*pi*freq
  
  end do
  
else

! Complex starting poles
  npairs=order/2
  
  if (npairs.GT.1) then
    fstep_poles=(fmax-fmin)/(npairs-1)
    log_fstep_poles=(log_fmax-log_fmin)/(npairs-1)
  else
    fstep_poles=0d0
    log_fstep_poles=0d0
  end if
  
  pole=0
  do i=1,npairs
  
    if (npairs.GT.1) then
  
      if (log_freq_starting_poles) then
        freq=10d0**( (log_fmin+(i-1)*log_fstep_poles) )
      else
        freq=fmin+(i-1)*fstep_poles
      end if
      
    else
  
      if (log_freq_starting_poles) then
        freq=10d0**( (log_fmin+log_fmax)/2d0 )
      else
        freq=(fmax+fmin)/2d0
      end if
      
    end if
    
    beta=2d0*pi*freq
    alpha=-beta/100d0
    pole=pole+1
    a(pole)=alpha+j*beta
    pole=pole+1
    a(pole)=alpha-j*beta
  
  end do

  if (pole.NE.order) then
! we have an odd number of poles so we must include one real pole
    alpha=0.75*2d0*pi*fmin   
    pole=pole+1 
    a(pole)=-alpha
  end if

end if

  write(*,*)
  write(*,*)'        initial poles (w)             intitial poles (f)           '
  do i=1,order
    write(*,*)cmplx(a(i)),cmplx(a(i)/(2d0*pi))
  end do 
  write(*,*)

! Iteration Loop

do iteration=1,max_iterations

  if(verbose) write(*,*)'Iteration',iteration,' of ',max_iterations
! Pole calculation process

  if(verbose) write(*,*)'CALLING update_poles'
  CALL update_poles(order,a,new_a,c,c_hat,d,h,nf,f,s,include_const_term,include_s_term,use_SVD_invert,verbose)

  a(:)=new_a(:)
  
  if(verbose) write(*,*)'CALLING calculate_residues'
  CALL calculate_residues(order,a,new_a,c,c_hat,d,h,nf,f,s,include_const_term,include_s_term,use_SVD_invert,verbose)
         
! Calculate the new fit function

  if(verbose) write(*,*)'CALLING calculate_fit_function'
  CALL calculate_fit_function(order,a,c,d,h,nf,f_fit,s)
      
! Calculate the mean square error between the original funcition, f, and the function fit

  if(verbose) write(*,*)'CALLING calculate_error'
  CALL calculate_error(nf,f,f_fit,s,error)
  
  write(*,*)'Iteration:',iteration,' error=',error
  
!  CALL write_coefficients(order,a,c,c_hat,d,h,fmin,fmax)

  if (error.LT.min_error) then
    write(*,*)'Exiting after ',iteration,' iterations'
    EXIT
  end if

! Next iteration
end do

! Write the poles, zeros and error to screen

if(verbose) write(*,*)'CALLING write_coefficients'
CALL write_coefficients(order,a,c,c_hat,d,h,fmin,fmax)

! write filter functions for plotting

if(verbose) write(*,*)'CALLING write_data_to_plot'
CALL write_data_to_plot(nf,f,f_fit,s)

! Finish up

DEALLOCATE( f )
DEALLOCATE( s )
DEALLOCATE( f_fit )
DEALLOCATE( a )
DEALLOCATE( new_a )
DEALLOCATE( c )
DEALLOCATE( c_hat )

END PROGRAM Vfit

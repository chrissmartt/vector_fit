SUBROUTINE set_test_set_1(nf,fmin,fmax,s,f,log_freq_samples)

IMPLICIT NONE

integer :: nf
complex*16 :: f(nf)
complex*16 :: s(nf)
real*8 :: fmin,fmax
logical :: log_freq_samples

! local variables

real*8 :: fstep,freq
real*8 :: log_fmin,log_fmax,log_fstep

integer :: n

real*8,parameter     :: pi=3.1415926535897931d0
complex*16,parameter :: j=(0d0,1d0)

! START

fmin=1D0
fmax=1D4

if (log_freq_samples) then
  log_fmin=log10(fmin)
  log_fmax=log10(fmax)
  log_fstep=(log_fmax-log_fmin)/(nf-1)
else
  fstep=(fmax-fmin)/(nf-1)
end if

do n=1,nf

  if (log_freq_samples) then
    freq=10d0**( (log_fmin+(n-1)*log_fstep) )
  else
    freq=fmin+(n-1)*fstep
  end if

! s=jw
  s(n)=2d0*pi*j*freq

! evaluate function of s
  f(n)=  2d0/(s(n)+5d0) + (30d0+j*40d0)/(s(n)-(-100d0+j*500d0))+ (30d0-j*40d0)/(s(n)-(-100d0-j*500d0)) + 0.5d0
  
end do


END SUBROUTINE set_test_set_1
!
! ___________________________________________________
!
!
SUBROUTINE set_test_set_2(nf,fmin,fmax,s,f,log_freq_samples)

! Example from section 4.1 of vector fit paper

IMPLICIT NONE

integer :: nf
complex*16 :: f(nf)
complex*16 :: s(nf)
real*8 :: fmin,fmax
logical :: log_freq_samples

! local variables

real*8 :: fstep,freq
real*8 :: log_fmin,log_fmax,log_fstep

integer :: n

real*8,parameter     :: pi=3.1415926535897931d0
complex*16,parameter :: j=(0d0,1d0)

! START


fmin=1D0
fmax=100D3

if (log_freq_samples) then
  log_fmin=log10(fmin)
  log_fmax=log10(fmax)
  log_fstep=(log_fmax-log_fmin)/(nf-1)
else
  fstep=(fmax-fmin)/(nf-1)
end if

do n=1,nf

  if (log_freq_samples) then
    freq=10d0**( (log_fmin+(n-1)*log_fstep) )
  else
    freq=fmin+(n-1)*fstep
  end if

! s=jw
  s(n)=2d0*pi*j*freq
  
! evaluate function of s
  f(n)=      0.2d0+2D-5*s(n)        &
        -2d0*pi*3000d0         / ( s(n) - (2d0*pi*( -4500d0)            ))  &
        -2d0*pi*83000d0        / ( s(n) - (2d0*pi*( -41000d0)           ))  &
        +2d0*pi*(-5d0+j*7d3)   / ( s(n) - (2d0*pi*(-100d0+j*5d3)  ))  &
        +2d0*pi*(-5d0-j*7d3)   / ( s(n) - (2d0*pi*(-100d0-j*5d3)  ))  &
        +2d0*pi*(-20d0+j*18d3) / ( s(n) - (2d0*pi*(-120d0+j*15d3) ))  &
        +2d0*pi*(-20d0-j*18d3) / ( s(n) - (2d0*pi*(-120d0-j*15d3) ))  &
        +2d0*pi*(6d3+j*45d3)   / ( s(n) - (2d0*pi*(-3d3+j*35d3) ))  &
        +2d0*pi*(6d3-j*45d3)   / ( s(n) - (2d0*pi*(-3d3-j*35d3) ))  &
        +2d0*pi*(40d0 +j*60d3) / ( s(n) - (2d0*pi*(-200d0+j*45d3) ))  &
        +2d0*pi*(40d0-j*60d3)  / ( s(n) - (2d0*pi*(-200d0-j*45d3) ))  &
        +2d0*pi*(90d0 +j*10d3) / ( s(n) - (2d0*pi*(-1500d0+j*45d3)))  &
        +2d0*pi*(90d0-j*10d3)  / ( s(n) - (2d0*pi*(-1500d0-j*45d3)))  &
        +2d0*pi*(5d4+j*80d3)   / ( s(n) - (2d0*pi*(-5d2+j*70d3) ))  &
        +2d0*pi*(5d4-j*80d3)   / ( s(n) - (2d0*pi*(-5d2-j*70d3) ))  &
        +2d0*pi*(1d3+j*45d3)   / ( s(n) - (2d0*pi*(-1d3+j*73d3) ))  &
        +2d0*pi*(1d3-j*45d3)   / ( s(n) - (2d0*pi*(-1d3-j*73d3) ))  &
        +2d0*pi*(-5d3+j*92d3)  / ( s(n) - (2d0*pi*(-2d3+j*90d3) ))  &
        +2d0*pi*(-5d3-j*92d3)  / ( s(n) - (2d0*pi*(-2d3-j*90d3) )) 
        
end do


END SUBROUTINE set_test_set_2
!
! ___________________________________________________
!
!
SUBROUTINE set_test_set_3(nf,fmin,fmax,s,f,log_freq_samples)

IMPLICIT NONE

integer :: nf
complex*16 :: f(nf)
complex*16 :: s(nf)
real*8 :: fmin,fmax
logical :: log_freq_samples

! local variables

real*8 :: fstep,freq
real*8 :: log_fmin,log_fmax,log_fstep

integer :: n

real*8,parameter     :: pi=3.1415926535897931d0
complex*16,parameter :: j=(0d0,1d0)

! START


fmin=1D0
fmax=3D2

if (log_freq_samples) then
  log_fmin=log10(fmin)
  log_fmax=log10(fmax)
  log_fstep=(log_fmax-log_fmin)/(nf-1)
else
  fstep=(fmax-fmin)/(nf-1)
end if

do n=1,nf

  if (log_freq_samples) then
    freq=10d0**( (log_fmin+(n-1)*log_fstep) )
  else
    freq=fmin+(n-1)*fstep
  end if

! s=jw
  s(n)=2d0*pi*j*freq
  
! evaluate function of s
  f(n)=  2d0/(s(n)+5d0)                   &
  + (30d0+j*40d0)/(s(n)-(-100d0+j*500d0))+ (30d0-j*40d0)/(s(n)-(-100d0-j*500d0))   &
  + (20d0+j*90d0)/(s(n)-(-150d0+j*800d0))+ (20d0-j*90d0)/(s(n)-(-150d0-j*800d0))   &
  + 0.5d0+2d-5*s(n) 
  
end do


END SUBROUTINE set_test_set_3
!
! ___________________________________________________
!
!
SUBROUTINE set_test_set_4(nf,fmin,fmax,s,f,log_freq_samples)

IMPLICIT NONE

integer :: nf
complex*16 :: f(nf)
complex*16 :: s(nf)
real*8 :: fmin,fmax
logical :: log_freq_samples

! local variables

real*8 :: fstep,freq
real*8 :: log_fmin,log_fmax,log_fstep

real*8 :: L,C,R

integer :: n

real*8,parameter     :: pi=3.1415926535897931d0
complex*16,parameter :: j=(0d0,1d0)

! START
L=1d-6
C=1d-6
R=10d0

fmin=1D3
fmax=1d7

if (log_freq_samples) then
  log_fmin=log10(fmin)
  log_fmax=log10(fmax)
  log_fstep=(log_fmax-log_fmin)/(nf-1)
else
  fstep=(fmax-fmin)/(nf-1)
end if

do n=1,nf

  if (log_freq_samples) then
    freq=10d0**( (log_fmin+(n-1)*log_fstep) )
  else
    freq=fmin+(n-1)*fstep
  end if

! s=jw
  s(n)=2d0*pi*j*freq
  
! evaluate function of s
  f(n)=    R+s(n)*L+1d0/(s(n)*C)   &
       + (2d7+j*9d5)/(s(n)-(-150000d0+j*2d7))+ (2d7-j*9d5)/(s(n)-(-150000d0-j*2d7))

end do


END SUBROUTINE set_test_set_4
!
! ___________________________________________________
!
!

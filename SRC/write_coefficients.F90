SUBROUTINE write_coefficients(order,a,c,c_hat,d,h,fmin,fmax)

IMPLICIT NONE

integer :: order
complex*16 :: a(order)
complex*16 :: c(order)
complex*16 :: c_hat(order)
real*8 :: d
real*8 :: h
real*8 :: fmin,fmax

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
  write(50,*)'               pole (a)                          residue (c)'
  do i=1,order
    write(50,*)dble(a(i)/wnorm),imag(a(i)/wnorm),dble(c(i)/wnorm),imag(c(i)/wnorm)
  end do 
  
  write(50,*)'wnorm_min   wnorm_max      nw'
  write(50,*)fmin/fmax, 1d0 ,200
  
  close(unit=50)
  

END SUBROUTINE write_coefficients

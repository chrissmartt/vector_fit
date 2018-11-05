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

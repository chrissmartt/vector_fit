SUBROUTINE calculate_fit_function(order,a,c,d,h,nf,f_fit,s)

IMPLICIT NONE

integer :: order
complex*16 :: a(order)
complex*16 :: c(order)
real*8 :: d
real*8 :: h

integer :: nf
complex*16 :: f_fit(nf)
complex*16 :: s(nf)

! local variables

integer :: n,i

! START

  do n=1,nf
  
! evaluate fitted function of s
    f_fit(n)= d+s(n)*h
  
    do i=1,order
      f_fit(n)= f_fit(n)+ c(i)/(s(n)-a(i))
    end do 
  
  end do ! next frequency


END SUBROUTINE calculate_fit_function

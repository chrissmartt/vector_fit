SUBROUTINE write_data_to_plot(nf,f,f_fit,s)

IMPLICIT NONE

integer :: nf
complex*16 :: f(nf)
complex*16 :: f_fit(nf)
complex*16 :: s(nf)

real*8     :: error

! local variables

real*8,parameter     :: pi=3.1415926535897931d0
complex*16,parameter :: j=(0d0,1d0)

integer :: n
real*8  :: frequency

! START


! open a file and write the frequency domain function and function fit data

  open(unit=10,file='function_data.fout')

  do n=1,nf
      
    frequency=dble(s(n)/(2d0*pi*j))
    
    write(10,8000)frequency,real(f(n)),imag(f(n)),abs(f(n)),real(f_fit(n)),imag(f_fit(n)),abs(f_fit(n))
    
8000 format(7ES16.6)
  
  end do ! next frequency
  
  close(unit=10)

END SUBROUTINE write_data_to_plot

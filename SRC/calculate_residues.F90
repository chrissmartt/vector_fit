SUBROUTINE calculate_residues(order,a,new_a,c,c_hat,d,h,nf,f,s,include_const_term,include_s_term,use_SVD_invert,verbose)

IMPLICIT NONE

! Variables passed to subroutine

integer :: order                ! Number of poles
complex*16 :: a(order)          ! complex poles
complex*16 :: new_a(order)      ! complex poles
complex*16 :: c(order)          ! complex residues of f(s)
complex*16 :: c_hat(order)      ! complex residues of sigma(s) approximation
real*8 :: d                     ! constatnt term
real*8 :: h                     ! s term

integer :: nf
complex*16 :: f(nf)
complex*16 :: s(nf)

logical :: include_const_term
logical :: include_s_term
logical :: use_SVD_invert
logical :: verbose

! local variables

complex*16,parameter :: j=(0d0,1d0)
real*8,parameter     :: small=1d-8

integer :: nr,nr2,nc   ! number of rows and columns in the matrix equation
integer :: row,col
integer :: i

integer :: pole_type(order)
real*8  :: test1,test2

! Complex solution matrices and vectors
complex*16,allocatable :: Amat(:,:)
complex*16,allocatable :: x(:)
complex*16,allocatable :: b(:)

! Real solution matrices and vectors
real*8,allocatable :: Armat(:,:)
real*8,allocatable :: xr(:)
real*8,allocatable :: br(:)

! matrix for zero calculation
real*8,allocatable     :: Hrmat(:,:)
real*8,allocatable     :: Azr(:,:)
real*8,allocatable     :: bzr(:)
real*8,allocatable     :: czr_hat(:)
complex*16,allocatable :: zeros(:)

! START

if(verbose) then
  write(*,*)'CALLED calculate_residues'
  write(*,*)'order=',order
  write(*,*)'nf=',nf
end if

nr=nf
nr2=2*nr
nc=order   ! number of unknows: d,h,c1..c_order
if (include_const_term) nc=nc+1
if (include_s_term) nc=nc+1

! allocate memory for the matrix equation
if(verbose) write(*,*)'Allocating memory: nr=',nr,' nc=',nc,' order=',order

ALLOCATE( Amat(nr,nc) )
ALLOCATE( x(nc) )
ALLOCATE( b(nr) )

ALLOCATE( Armat(nr2,nc) )
ALLOCATE( xr(nc) )
ALLOCATE( br(nr2) )

ALLOCATE( Hrmat(order,order) )
ALLOCATE( Azr(order,order) )
ALLOCATE( bzr(order) )
ALLOCATE( czr_hat(order) )
ALLOCATE( zeros(order) )

! sort the poles into either real poles or complex conjugate pairs
i=1

do while(i.LE.order)

! Check the imaginary part of the pole a(i). if its magnitude is 
! greater than some small number assume we have a complex pole pair

  if (abs(imag(a(i))).GT.small) then
! assume we have the first pole of a complex pole pair
    pole_type(i)=1
    
! check that there is another pole
    if (i+1.GT.order) then
      write(*,*)'Error in update_poles: no second pole of conjugate pair'    
      write(*,*)'pole i  :',a(i)
      STOP
    end if
    
! Check that the poles i and i+1 are complex conjugates
! test 1: Real(pole(i)-pole(i+1))=0

    test1=dble(a(i)-a(i+1))
    
! test 2: Imag(pole(i)+pole(i+1))=0
    test2=imag(a(i)+a(i+1))
    
    if ( (abs(test1).GT.small).OR.(abs(test2).GT.small) ) then
! one of the tests has failed
      write(*,*)'Error in update_poles: poles do not form a conjugate pair'
      write(*,*)'pole i  :',a(i)
      write(*,*)'pole i+1:',a(i+1)
      write(*,*)'test1 value=',test1,' should be <',small
      write(*,*)'test2 value=',test2,' should be <',small
      STOP
    end if

    i=i+1
    pole_type(i)=2
  else
! This is assumed to be a real pole
    pole_type(i)=0
  end if
  
  i=i+1

end do

if(verbose) then
  write(*,*)'Pole types:'
  do i=1,order
    write(*,*)'pole',i,' type ',pole_type(i)
  end do
end if

! set up the linear equations

do row=1,nf   ! one row for each frequency of interest

  col=0
  
! (sigma f)_fit terms. note the form used for complex pole pairs in Appendix A
  do i=1,order
    if (pole_type(i).EQ.0) then
      col=col+1
      Amat(row,col)=(1d0,0d0)/(s(row)-a(i))
    else if (pole_type(i).EQ.1) then
      col=col+1
      Amat(row,col)=(1d0,0d0)/(s(row)-a(i))+(1d0,0d0)/(s(row)-a(i+1))
    else if (pole_type(i).EQ.2) then
      col=col+1
      Amat(row,col)=j*((1d0,0d0)/(s(row)-a(i-1))-(1d0,0d0)/(s(row)-a(i)))
    end if
  end do
  
! d term
  if (include_const_term) then
    col=col+1
    Amat(row,col)=(1d0,0d0)
  end if
  
! h term
  if (include_s_term) then
    col=col+1
    Amat(row,col)=s(row)
  end if
     
! right hand side

  b(row)=f(row)

end do ! next row (frequency)

! work out the real equations to solve

do row=1,nf   ! one row for each frequency of interest

  do col=1,nc
    Armat(row   ,col)=dble(Amat(row,col))
    Armat(row+nr,col)=imag(Amat(row,col))
    br(row   )=dble(b(row))
    br(row+nr)=imag(b(row)) 
  end do
  
end do ! next row of complex equation set

if (verbose) then
  write(*,*)'A matrix, nrows=nf, ncols=order+2'
  do row=1,nr2
    write(*,'(I6,100ES12.2)')row,Armat(row,:)
  end do
end if

! solve the linear equations 

if (use_SVD_invert) then
  if (verbose) write(*,*)'CALLING svd_solve'
  CALL svd_solve(Armat,nr2,nc,xr,br,verbose)  
else
  if (verbose) write(*,*)'CALLING morse_penrose_solve'
  CALL morse_penrose_solve(Armat,nr2,nc,xr,br,verbose)  
end if

if(verbose) then
  write(*,*)'Solution vector xr:'
  do i=1,nc
    write(*,*)i,xr(i)
  end do
end if

! calculate the complex solution vector, x and hence the
! solution for the zeros of sigma, c_hat

row=0
  
! (sigma f)_fit terms. note the form used for complex pole pairs in Appendix A
do i=1,order
  if (pole_type(i).EQ.0) then
    row=row+1
    x(row)=xr(row)
  else if (pole_type(i).EQ.1) then
    row=row+1
    x(row)  =xr(row)+j*xr(row+1)
  else if (pole_type(i).EQ.2) then
    row=row+1
    x(row)  =xr(row-1)-j*xr(row)
  end if
end do
  
  ! d term
if (include_const_term) then
  row=row+1
  x(row)=xr(row)
end if
  
! h term
if (include_s_term) then
  row=row+1
  x(row)=xr(row)
end if
  
if(verbose) then
  write(*,*)'Solution vector x:'
  do i=1,nc
    write(*,*)i,x(i)
  end do
end if

! we now have x, so get the unknowns c,d,h in order from the solution vector
row=0
  
do i=1,order
  row=row+1
  c(i)=x(row)
end do

if (include_const_term) then
  row=row+1
  d=x(row)
else
  d=0d0
end if

if (include_s_term) then
  row=row+1
  h=x(row)
else
  h=0d0
end if

if (verbose) then
  write(*,*)'c vector'
  do row=1,order
    write(*,*)row,c(row)
  end do
  
  write(*,*)'d=',d
  
  write(*,*)'h=',h 
  
end if


! deallocate memory 

DEALLOCATE( Amat )
DEALLOCATE( x )
DEALLOCATE( b )
DEALLOCATE( Armat )
DEALLOCATE( xr )
DEALLOCATE( br )
DEALLOCATE( Hrmat )
DEALLOCATE( Azr )
DEALLOCATE( bzr )
DEALLOCATE( czr_hat )
DEALLOCATE( zeros )

RETURN

END SUBROUTINE calculate_residues

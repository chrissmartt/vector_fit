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

  SUBROUTINE dinvert_Gauss_Jordan(A,n,AI,matdim) 

IMPLICIT NONE

! invert matrix using Gauss Jordan method with pivoting
       
  integer	:: n,matdim
  real*8	:: A(matdim,matdim)
  real*8	:: AI(matdim,matdim)

! local variables
  		       
  integer	:: row,col,reduce_col,i
  
  real*8	:: max_element
  real*8	:: pivot_element
  integer	:: max_row
  
  integer	:: pivot_row
  
  integer	:: pivot_row_save(matdim)
  
  real*8	:: row_multiplier
  real*8	:: swap

! START
  
! copy A to AI 
  AI(1:n,1:n)= A(1:n,1:n)

  pivot_row_save(1:matdim)=0
  
! loop over columns of the matrix and reduce each column in turn to identity matrix column     
  do reduce_col=1,n
  
! find the largest element in this column and use as the pivot element
    max_element=0d0
    max_row=0
    do row=reduce_col,n
      if (abs(AI(row,reduce_col)).GT.max_element) then
        max_element=abs(AI(row,reduce_col))
	max_row=row
      end if
    end do  
    
    if (max_row.eq.0) then
! all elements are zero so singular matrix
      write(*,*)'Singular matrix found in dinvert_Gauss_Jordan'
      write(*,*)'Reduction col=',reduce_col
      STOP
    end if
    
    pivot_row=max_row
    pivot_row_save(reduce_col)=pivot_row
    
! swap pivot row with the row reduce_col

    if (pivot_row.ne.reduce_col) then
      do col=1,n
        swap=AI(reduce_col,col)
	AI(reduce_col,col)=AI(pivot_row,col)
	AI(pivot_row,col)=swap
      end do 
    end if
    
    pivot_row=reduce_col   
    pivot_element=AI(reduce_col,reduce_col)
    
! operate on pivot row    
    do col=1,n
      if (col.ne.reduce_col) then
        AI(pivot_row,col) = AI(pivot_row,col)/pivot_element
      else	
        AI(pivot_row,col) = 1d0/pivot_element
      end if
    end do

! operate on rows other than the pivot row   
    do row=1,n
    
      if (row.ne.pivot_row) then
      
        row_multiplier=AI(row,reduce_col)
	
        do col=1,n
          if (col.ne.reduce_col) then
            AI(row,col) = AI(row,col)- AI(pivot_row,col)*row_multiplier
          else	
            AI(row,reduce_col) =-AI(pivot_row,reduce_col)*row_multiplier
          end if
        end do
	
      end if ! not pivot row
    
    end do ! next row
    
  end do ! next column of the matrix to reduce
  
  do reduce_col=n,1,-1
  
    if (reduce_col.ne.pivot_row_save(reduce_col)) then
! rows were swapped so must swap the corresponding columns

      do row=1,n
        swap=AI(row,pivot_row_save(reduce_col))
        AI(row,pivot_row_save(reduce_col))=AI(row,reduce_col)
	AI(row,reduce_col)=swap
      end do
      
    end if
    
  end do

  RETURN
  
  END
!
!________________________________________________________________________
!
!
  SUBROUTINE morse_penrose_solve(A,nr,nc,x,b,verbose)

IMPLICIT NONE

! left generalised inverse matrix 
       
  integer	:: nr,nc
  real*8	:: A(nr,nc)
  real*8	:: x(nc)
  real*8	:: b(nr)
  
  logical :: verbose

! local variables

  real*8	:: AT(nc,nr)
  real*8	:: ATA(nc,nc)
  real*8	:: ATAI(nc,nc)
  real*8	:: AI(nc,nr)
  integer i,j
		       
! START

  if(verbose) write(*,*)'CALLED morse_penrose_solve'

  if (nr.LT.nc) then
    write(*,*)'Error in morse_penrose_solve, nr < nc'
    STOP
  end if
  
  AT=TRANSPOSE(A)
  ATA=MATMUL(AT,A)
    
  CALL dinvert_Gauss_Jordan(ATA,nc,ATAI,nc) 
  
  AI=MATMUL(ATAI,AT)
  
  x=MATMUL(AI,b)  
  
  RETURN
  
  END
!
!________________________________________________________________________
!
!
  SUBROUTINE svd_solve(Ain,nr,nc,x,b,verbose)

IMPLICIT NONE

! left generalised inverse matrix 
       
  integer	:: nr,nc
  real*8	:: Ain(nr,nc)
  real*8	:: x(nc)
  real*8	:: b(nr)
  
  logical :: verbose

! local variables

  real*8	:: A(nr,nr)
  real*8	:: AI(nr,nr)
  integer i,j
		       
! START

  if(verbose) write(*,*)'CALLED svd_solve'

  if (nr.LT.nc) then
    write(*,*)'Error in svd_solve, nr < nc'
    STOP
  end if
  
  A(:,:)=0d0
  
  do i=1,nr
    do j=1,nc
      A(i,j)=Ain(i,j)
    end do
  end do
   
  CALL dsvd_invert(A,nr,nc,AI,nr) 
     
  do i=1,nc
    x(i)=0d0
    do j=1,nr
      x(i)=x(i)+AI(i,j)*b(j)
    end do
  end do
   
  RETURN
  
  END
!
! __________________________________________________
!
!
       subroutine dsvd(A,ar,ac,matdim,U,GAMMA,VT) 

IMPLICIT NONE

! invert matrix using SVD implemented in EISPACK
       
       integer 		:: ar,ac,matdim
       real*8 		:: A(matdim,matdim)
       real*8 		:: GAMMA(matdim)
       real*8 		:: U(matdim,matdim),VT(matdim,matdim)
       
!      EISPACK ARGUMENTS

       integer*4	:: m,n
       integer*4	:: ierr
       logical 		:: matu
       logical 		:: matv
       REAL*8,allocatable :: A_local(:,:)
       REAL*8,allocatable :: W_local(:)
       REAL*8,allocatable :: U_local(:,:)
       REAL*8,allocatable :: V_local(:,:)
       
       integer		:: row,col

! START
       
! perform singular value decomposition on A
      
      M=ar
      N=ac
      
      matu=.TRUE.
      matv=.TRUE.
      
      ALLOCATE( A_local(M,N) )
      ALLOCATE( W_local(N) )
      ALLOCATE( U_local(M,N) )
      ALLOCATE( V_local(N,N) )
      
      do row=1,M
        do col=1,N
	
	  A_local(row,col)=A(row,col)
	  
	end do
      end do
      
      CALL svd ( m, n, A_local, W_local , matu, U_local, matv, V_local, ierr )
      
      GAMMA(:)=0d0
      U(:,:)=0d0
      VT(:,:)=0d0
      
      do row=1,N	
	GAMMA(row)=W_local(row)	  
      end do
      
      do row=1,M
        do col=1,N	
	  U(row,col)=U_local(row,col)	  
	end do
      end do
      
      do row=1,N
        do col=1,N	
	  VT(row,col)=V_local(col,row)	  
	end do
      end do
      
      DEALLOCATE( A_local )
      DEALLOCATE( W_local )
      DEALLOCATE( U_local )
      DEALLOCATE( V_local )
		      
      return
      end
!
! __________________________________________________
!
!
       subroutine dsvd_invert(A,ar,ac,AI,matdim) 

IMPLICIT NONE

! invert matrix using SVD implemented in EISPACK
       
       integer ar,ac,matdim
       real*8 A(matdim,matdim)
       real*8 AI(matdim,matdim)
                     
       DOUBLE PRECISION  GAMMA(matdim)
       DOUBLE PRECISION  GAMMA_I(matdim)
       real*8   V(matdim,matdim),VT(matdim,matdim)
       real*8   U(matdim,matdim),UT(matdim,matdim)
       
       real*8 tm2(matdim,matdim)
       
       integer row,col

! START
     
! zero matrices 
     
       do row=1,matdim
         do col=1,matdim
	   U(row,col)=0.0D0
	   UT(row,col)=0.0D0
	   V(row,col)=0.0D0
	   VT(row,col)=0.0D0
	   tm2(row,col)=0.0D0
	 end do
	 GAMMA(row)=0.0D0
	 GAMMA_I(row)=0.0D0
       end do

! perform singular value decomposition on A
       
       call dsvd(A,ar,ac,matdim,U,GAMMA,VT)
              
! calculate gamma inverse
       
       do row=1,ac
       
         if (abs(gamma(row)).ne.(0.0D0)) then
           gamma_i(row)=1.0D0/gamma(row)
         else 
           gamma_i(row)=0.0D0
	 end if
       end do
              
! calculate transposes of U and VT
!       UT=TRANSPOSE(U)
!       V=TRANSPOSE(VT)
       call dtranspose(U,ar,ac,UT,matdim)
       call dtranspose(VT,ac,ac,V,matdim)       
                
! pre multiply U transpose by gamma inverse

       do row=1,ac
         do col=1,ar
	   tm2(row,col)=UT(row,col)*gamma_i(row)
	 end do
       end do
                     
! pre multiply result by V to give ai, the inverse of a
                               	
!       AI=MATMUL(V,TM2)
       call dmatmul(V,ac,ac,tm2,ac,ar,AI,matdim)	 

       RETURN
      
       end
!
! __________________________________________________
!
!
       subroutine deig(A,ar,GAMMA,matdim) 

IMPLICIT NONE

! calculate the eigenvalues of A using EISPACK

       integer ar,matdim
       real*8 A(matdim,matdim)
       complex*16 GAMMA(matdim)
       
!      EISPACK ARGUMENTS

       integer*4	:: n
       integer*4	:: ierr
       integer 		:: matz
       REAL*8,allocatable :: A_local(:,:)
       REAL*8,allocatable :: WR_local(:)
       REAL*8,allocatable :: WI_local(:)
       REAL*8,allocatable :: Z_local(:,:)
       
       integer		:: row,col

! START
         
  N=ar
  
  matz=0   ! set to zero integer to calculate eigenvalues
  
  ALLOCATE( A_local(N,N) )
  ALLOCATE( WR_local(N) )
  ALLOCATE( WI_local(N) )
  ALLOCATE( Z_local(N,N) )

! Copy matrix  
  do row=1,N
    do col=1,N

      A_local(row,col)=A(row,col)
      
    end do
  end do
  
  CALL rg ( n, A_local, WR_local, WI_local , matz, Z_local, ierr )
  
  if (ierr.ne.0) then
    write(*,*)'Error in deig'
    STOP
  end if
          
  do row=1,n
    gamma(row)=dcmplx(Wr_local(row),Wi_local(row)) 
  end do		
      
  DEALLOCATE( A_local )
  DEALLOCATE( WR_local )
  DEALLOCATE( WI_local )
  DEALLOCATE( Z_local )
    		  
  return
  end
!_____________________________________________________________________
!
! NAME
!    dmatmul
!
! DESCRIPTION
!    multiply real*8 matrices
!
!
! SEE ALSO
!   
!
! HISTORY
!
!     started 18/05/05 CJS
!
! COMMENTS
!   
!
  subroutine dmatmul(a,ar,ac,b,br,bc,ans,matdim)

! Modules used

IMPLICIT NONE

!Argument list variables

  integer matdim
  real*8 a(matdim,matdim),b(matdim,matdim),       &
  	     ans(matdim,matdim)
  integer ar,ac,br,bc

!Local variables

  real*8 sum
  integer row,col,i
 
! START

  if (ac.ne.br) then
    print*,'matrix dimension error in dmatmul'
    stop
  end if
  do row=1,ar
    do col=1,bc 
      sum=0d0
      do i=1,ac
  	sum=sum+a(row,i)*b(i,col)
      end do
      ans(row,col)=sum
    end do
  end do

! END
 
  return
  end
!
!_____________________________________________________________________
!
! NAME
!    dmatdmul
!
! DESCRIPTION
!    post multiply a real*8 matrix by a diagonal matrix
!
!
! SEE ALSO
!   
!
! HISTORY
!
!     started 18/05/05 CJS
!
! COMMENTS
!   
!
 subroutine dmatdmul(a,ar,ac,d,dr,ans,matdim)

! Modules used

IMPLICIT NONE

!Argument list variables

 integer matdim
 real*8 a(matdim,matdim),d(matdim),ans(matdim,matdim)
 integer ar,ac,dr

!Local variables

 integer r,c
 
! START

 if (ac.ne.dr) then
   print*,'matrix dimension error in dmatdmul'
   stop
 end if
 do r=1,ar
   do c=1,dr
     ans(r,c)=a(r,c)*d(c)
   end do
 end do

! END
 
 return
 end 
!
!_____________________________________________________________________
!
! NAME
!    dmatvmul
!
! DESCRIPTION
!    multiply a vector by a real*8 matrix
!
!
! SEE ALSO
!   
!
! HISTORY
!
!     started 18/05/05 CJS
!
! COMMENTS
!   
!
 subroutine dmatvmul(a,ar,ac,v,vr,ans,matdim)

! Modules used

IMPLICIT NONE

!Argument list variables

 integer matdim
 real*8 a(matdim,matdim),v(matdim),ans(matdim)
 integer ar,ac,vr

!Local variables

 real*8 sum
 integer r,c
 
! START

 if (ac.ne.vr) then
   print*,'matrix dimension error'
   stop
 end if
 do r=1,ar
   sum=0d0
   do c=1,vr
     sum=sum+a(r,c)*v(c)
   end do
   ans(r)=sum
 end do

! END
 
 return
 end
!
!_____________________________________________________________________
!
! NAME
!    dtranspose
!
! DESCRIPTION
!    
!    calculate  transpose of a matrix
!
! SEE ALSO
!   
!
! HISTORY
!
!     started 18/05/05 CJS
!
! COMMENTS
!
!
 subroutine dtranspose(a,nr,nc,at,matdim)

! Modules used

IMPLICIT NONE

!Argument list variables

 integer matdim,nr,nc
 real*8 a(matdim,matdim),at(matdim,matdim)

!Local variables

 integer r,c
 
! START

 do r=1,nr
   do c=1,nc
    at(c,r)=a(r,c)
   end do
 end do
 
! END
 
 return
 end

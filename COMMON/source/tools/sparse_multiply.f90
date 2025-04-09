! ********************************************************************* !
! Subroutine to multiply two sparse matrices
! Copyright (C) 2017  Erica Galetti
! 
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! 
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.

! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
! ********************************************************************* !


SUBROUTINE sparse_multiply(A_size,A_indices,A_values,B_size,B_indices,B_values,C_size,C_indices,C_values,max_number_of_nonzeros)

IMPLICIT NONE
INTEGER, DIMENSION(2), INTENT(IN)				:: A_size,B_size
INTEGER, DIMENSION(:,:), INTENT(IN)				:: A_indices,B_indices
DOUBLE PRECISION, DIMENSION(:), INTENT(IN)			:: A_values,B_values
INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT)			:: C_size
INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT)		:: C_indices
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, INTENT(OUT)	:: C_values
INTEGER, INTENT(IN)						:: max_number_of_nonzeros
INTEGER								:: nn,rr,cc,ii,jj
INTEGER								:: number_of_nonzero_columns,number_of_nonzero_rows
INTEGER								:: column_number,row_number
INTEGER, DIMENSION(:), ALLOCATABLE				:: A_column_indices,B_row_indices
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE			:: A_column_values,B_row_values
DOUBLE PRECISION						:: value
INTEGER, DIMENSION(:,:), ALLOCATABLE				:: work_indices
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE			:: work_values
REAL								:: start,finish

  IF (A_size(2).NE.B_size(1)) THEN
    WRITE(*,*)
    WRITE(*,*)'*****************************************'
    WRITE(*,*)'ERROR in subroutine sparse_multiply'
    WRITE(*,*)'-----------------------------------------'
    WRITE(*,*)'The number of columns of A does not match'
    WRITE(*,*)'the number of rows of B!!'
    WRITE(*,*)'TERMINATING PROGRAM!!'
    WRITE(*,*)'*****************************************'
    WRITE(*,*)
    STOP
  END IF
  
  ALLOCATE(C_size(2))
  C_size=[A_size(1),B_size(2)]
  IF (max_number_of_nonzeros.EQ.0) THEN
    ALLOCATE(C_indices(C_size(1)*C_size(2),2))
    ALLOCATE(C_values(C_size(1)*C_size(2)))
  ELSE
    ALLOCATE(C_indices(max_number_of_nonzeros,2))
    ALLOCATE(C_values(max_number_of_nonzeros))
  END IF
  
  call cpu_time(start)
  
  nn=0
  ! Loop over all rows of A
  DO rr=1,A_size(1)
    
    ! Check if row "rr" of A is all zeros
    number_of_nonzero_columns=SIZE(PACK(A_indices(:,2),A_indices(:,1).EQ.rr))
    IF (number_of_nonzero_columns.EQ.0) CYCLE
    
    ! If row "rr" of A contains nonzero elements, get their column indices and values
    ALLOCATE(A_column_indices(number_of_nonzero_columns))
    ALLOCATE(A_column_values(number_of_nonzero_columns))
    A_column_indices=PACK(A_indices(:,2),A_indices(:,1).EQ.rr)
    A_column_values=PACK(A_values,A_indices(:,1).EQ.rr)
    
    ! Loop over all columns of B
    DO cc=1,B_size(2)
      
      ! Check if column "cc" of B is all zeros
      number_of_nonzero_rows=SIZE(PACK(B_indices(:,1),B_indices(:,2).EQ.cc))
      IF (number_of_nonzero_rows.EQ.0) CYCLE
      
      ! If column "cc" of B contains nonzero elements, get their row indices and values
      ALLOCATE(B_row_indices(number_of_nonzero_rows))
      ALLOCATE(B_row_values(number_of_nonzero_rows))
      B_row_indices=PACK(B_indices(:,1),B_indices(:,2).EQ.cc)
      B_row_values=PACK(B_values,B_indices(:,2).EQ.cc)
      
      ! Loop over nonzero elements, multiply matching columns & rows, and add
      ! (equivalent to dot product between row "r" and column "cc")
      value=0
      DO jj=1,number_of_nonzero_columns
        column_number=A_column_indices(jj)
        DO ii=1,number_of_nonzero_rows
          row_number=B_row_indices(ii)
            IF (column_number.EQ.row_number) THEN
              value=value+A_column_values(jj)*B_row_values(ii)
            END IF
        END DO
      END DO
      
      ! Save value (if not equal to zero)
      IF (value.NE.0) THEN
        nn=nn+1
        C_indices(nn,:)=[rr,cc]
        C_values(nn)=value
      END IF
      
      DEALLOCATE(B_row_indices,B_row_values)
      
    END DO
    
    DEALLOCATE(A_column_indices,A_column_values)
    
  END DO
  
  ! Keep only nonzero elements
  ALLOCATE(work_indices(nn,2),work_values(nn))
  work_indices=C_indices(1:nn,:)
  work_values=C_values(1:nn)
  DEALLOCATE(C_indices,C_values)
  ALLOCATE(C_indices(nn,2),C_values(nn))
  C_indices=work_indices
  C_values=work_values
  DEALLOCATE(work_indices,work_values)
  
  call cpu_time(finish)
  print '("Time = ",f6.3," seconds.")',finish-start
  
END SUBROUTINE

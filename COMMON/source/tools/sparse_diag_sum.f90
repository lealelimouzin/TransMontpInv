! ********************************************************************* !
! Subroutine to sum a sparse matrix and a diagonal matrix
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


SUBROUTINE sparse_diag_sum(A_size,A_indices,A_values,D_size,D_indices,D_values,C_size,C_indices,C_values)

USE m_mrgrnk
IMPLICIT NONE
INTEGER, DIMENSION(:), INTENT(IN)                               :: A_size
INTEGER, INTENT(IN)                                             :: D_size
INTEGER, DIMENSION(:,:), INTENT(IN)                             :: A_indices
INTEGER, DIMENSION(:), INTENT(INOUT)                            :: D_indices
DOUBLE PRECISION, DIMENSION(:), INTENT(IN)                      :: A_values
DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT)                   :: D_values
INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT)                 :: C_size
INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT)               :: C_indices
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, INTENT(OUT)        :: C_values
INTEGER                                                         :: nn
INTEGER, DIMENSION(:), ALLOCATABLE                              :: diag_order
REAL                                                            :: start,finish
LOGICAL                                                         :: verbose=.FALSE.
INTERFACE 
        SUBROUTINE MRGRNK(XDONT,IRNGT)
                INTEGER, DIMENSION(:), INTENT(IN)       :: XDONT
                INTEGER, DIMENSION(:), INTENT(OUT)      :: IRNGT
        END SUBROUTINE MRGRNK
END INTERFACE

  IF (A_size(1).NE.D_size.OR.A_size(2).NE.D_size) THEN
    WRITE(*,*)
    WRITE(*,*)'***********************************'
    WRITE(*,*)'ERROR in subroutine sparse_diag_sum'
    WRITE(*,*)'-----------------------------------'
    WRITE(*,*)'The size of A does not match'
    WRITE(*,*)'the size of D!!'
    WRITE(*,*)'TERMINATING PROGRAM!!'
    WRITE(*,*)'***********************************'
    WRITE(*,*)
    STOP
  END IF
  
  IF (ALLOCATED(C_size)) DEALLOCATE(C_size)
  IF (ALLOCATED(C_indices)) DEALLOCATE(C_indices)
  IF (ALLOCATED(C_values)) DEALLOCATE(C_values)
  
  ! Sort the entries of D along the diagonal
  ALLOCATE(diag_order(SIZE(D_values)))
  CALL MRGRNK(D_indices,diag_order)
  D_indices=D_indices(diag_order)
  D_values=D_values(diag_order)
  DEALLOCATE(diag_order)
  
  ! Allocate space
  ALLOCATE(C_size(2))
  C_size=A_size
  ALLOCATE(C_indices(SIZE(A_values),2))
  ALLOCATE(C_values(SIZE(A_values)))
  C_indices=A_indices
  C_values=A_values
  
  IF (verbose) CALL cpu_time(start)
   
  ! Loop over the nonzero elements of A    
  DO nn=1,SIZE(A_values)
    
    IF (A_indices(nn,1).EQ.A_indices(nn,2)) THEN ! if we are on a diagonal element
    
      C_values(nn)=C_values(nn)+D_values(A_indices(nn,1))
      
    END IF
    
  END DO
    
  IF (verbose) CALL cpu_time(finish)
  IF (verbose) PRINT '("Time for A+DIAG = ",f12.9," seconds.")',finish-start
  
END SUBROUTINE

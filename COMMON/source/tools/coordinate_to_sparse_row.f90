! ********************************************************************* !
! Convert a sparse matrix from coordinate to compressed row format
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


SUBROUTINE coordinate_to_sparse_row(M_size,M_indices,M_values,M_rpointer)

! This subroutine rearranges the entries of a sparse matrix M (given in coordinate
! format) so that they are given into order of rows and then order of columns.
! It also creates a row pointer array which allows operations with other
! subroutines requiring sparse matrices in compressed row format.
!
! Erica Galetti, November 2016
! erica.galetti@ed.ac.uk

USE m_mrgrnk
IMPLICIT NONE
INTEGER, DIMENSION(:), INTENT(IN)                       :: M_size
INTEGER, DIMENSION(:,:), INTENT(INOUT)                  :: M_indices
DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT)           :: M_values
INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT)         :: M_rpointer
INTEGER, DIMENSION(:), ALLOCATABLE                      :: column_order,row_order,worki
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE             :: workd
INTEGER                                                 :: ii,jj,ini,fin,number_of_uniques

! The size of the row pointer array is the number of rows + 1
IF (ALLOCATED(M_rpointer)) DEALLOCATE(M_rpointer)
ALLOCATE(M_rpointer(M_size(1)+1))

! The last entry of the row pointer array is the number of nonzero
! elements in M + 1
M_rpointer(SIZE(M_rpointer))=SIZE(M_values)+1

! Find the sorted order of the row indices of M
ALLOCATE(row_order(SIZE(M_values)))
CALL I_MRGRNK(M_indices(:,1),row_order)

! Rearrange M into a sorted order
M_indices(:,1)=M_indices(row_order,1)
M_indices(:,2)=M_indices(row_order,2)
M_values=M_values(row_order)

! Create row pointer array
DO jj=1,M_size(1)
        M_rpointer(jj)=MINLOC(ABS(M_indices(:,1)-jj),1)
END DO

! Sort entries of each row by column
DO jj=1,M_size(1)
        ini=M_rpointer(jj)
        fin=M_rpointer(jj+1)-1
        ALLOCATE(column_order(SIZE( (/ (ii,ii=ini,fin) /) )))
        ALLOCATE(worki(SIZE( (/ (ii,ii=ini,fin) /) )))
        ALLOCATE(workd(SIZE( (/ (ii,ii=ini,fin) /) )))
        worki=M_indices(ini:fin,2)
        workd=M_values(ini:fin)
        CALL I_MRGRNK(worki,column_order)
        M_indices(ini:fin,2)=worki(column_order)
        M_values(ini:fin)=workd(column_order)
        DEALLOCATE(column_order,worki,workd)
END DO

! Free up memory
IF (ALLOCATED(column_order)) DEALLOCATE(column_order)
IF (ALLOCATED(row_order)) DEALLOCATE(row_order)
IF (ALLOCATED(worki)) DEALLOCATE(worki)
IF (ALLOCATED(workd)) DEALLOCATE(workd)

END SUBROUTINE
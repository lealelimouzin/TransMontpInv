! ********************************************************************* !
! This file is part of module "modelling_2_5D" in modelling_2_5D.f90.
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


SUBROUTINE calculate_A()

USE modelling_2_5D, ONLY: G_size,G_indices,G_cpointer,G_values,x_number_of_cells,z_number_of_cells,nax,naz,&
                          DS_size,DS_indices,D_rpointer,DS_values,A_size,A_indices,A_cpointer,A_values
IMPLICIT NONE
INTEGER                         :: ii,jj,nn,row,column,ini,fin
INTEGER                         :: number_of_nonzeros,entry_number,nx,nz,cpoint
DOUBLE PRECISION                :: value_0,value_1u,value_1l,value_nxu,value_nxl
INTEGER, ALLOCATABLE            :: DS_column_indices(:),G_row_indices(:)
DOUBLE PRECISION, ALLOCATABLE   :: DS_column_values(:),G_row_values(:)
INTERFACE
        SUBROUTINE get_G_values(G_indices,G_values,G_cpointer,column,G_row_indices,G_row_values)
                INTEGER, DIMENSION(:,:), INTENT(IN)                             :: G_indices
                DOUBLE PRECISION, DIMENSION(:), INTENT(IN)                      :: G_values
                INTEGER, DIMENSION(:), INTENT(IN)                               :: G_cpointer
                INTEGER, INTENT(IN)                                             :: column
                INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT)                 :: G_row_indices
                DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, INTENT(OUT)        :: G_row_values
        END SUBROUTINE
END INTERFACE
INTERFACE
        SUBROUTINE get_DS_values(DS_indices,DS_values,DS_rpointer,row,DS_column_indices,DS_column_values)
                INTEGER, DIMENSION(:,:), INTENT(IN)                             :: DS_indices
                DOUBLE PRECISION, DIMENSION(:), INTENT(IN)                      :: DS_values
                INTEGER, DIMENSION(:), INTENT(IN)                               :: DS_rpointer
                INTEGER, INTENT(IN)                                             :: row
                INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT)                 :: DS_column_indices
                DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, INTENT(OUT)        :: DS_column_values
        END SUBROUTINE
END INTERFACE
INTERFACE
        SUBROUTINE multiply(DS_column_indices,DS_column_values,G_row_indices,G_row_values,value)
                INTEGER, DIMENSION(:), INTENT(IN)                               :: DS_column_indices,G_row_indices
                DOUBLE PRECISION, DIMENSION(:), INTENT(IN)                      :: DS_column_values,G_row_values
                DOUBLE PRECISION, INTENT(OUT)                                   :: value
        END SUBROUTINE
END INTERFACE

! To calculate A=-D*S*G we can exploit the fact that we know its structure: it is a sparse matrix with nonzero
! entries along diagonal 0 (the main diagonal), along diagonals +1 and -1 (the first upper and lower diagonals),
! and along diagonals +nx and -nx. Since we also know the structure of D*S and G, we can easily work out the value
! of each nonzero entry of A.

nx=x_number_of_cells
nz=z_number_of_cells

! Loop over the the diagonal entries
 nn=0
 entry_number=0
 cpoint=0
DO jj=1,z_number_of_cells
        DO ii=1,x_number_of_cells
                
                nn=nn+1 ! column number of output matrix
                
                CALL get_G_values(G_indices,G_values,G_cpointer,nn,G_row_indices,G_row_values)
                
                value_0=0
                value_1u=0
                value_1l=0
                value_nxu=0
                value_nxl=0
                
                IF (jj.EQ.1) THEN
                        IF (ii.EQ.1) THEN
                                
                                CALL get_DS_values(DS_indices,DS_values,D_rpointer,nn,DS_column_indices,DS_column_values)
                                value_0=SUM(DS_column_values*G_row_values)
                                
                                !value_1u=value_1u
                                
                                CALL get_DS_values(DS_indices,DS_values,D_rpointer,nn+1,DS_column_indices,DS_column_values)
                                CALL multiply(DS_column_indices,DS_column_values,G_row_indices,G_row_values,value_1l)
                                
                                !value_nxu=value_nxu
                                
                                CALL get_DS_values(DS_indices,DS_values,D_rpointer,nn+nx,DS_column_indices,DS_column_values)
                                CALL multiply(DS_column_indices,DS_column_values,G_row_indices,G_row_values,value_nxl)
                                
                        ELSE IF (ii.GT.1.AND.ii.LT.x_number_of_cells) THEN
                                
                                CALL get_DS_values(DS_indices,DS_values,D_rpointer,nn,DS_column_indices,DS_column_values)
                                value_0=SUM(DS_column_values*G_row_values)
                                
                                CALL get_DS_values(DS_indices,DS_values,D_rpointer,nn-1,DS_column_indices,DS_column_values)
                                CALL multiply(DS_column_indices,DS_column_values,G_row_indices,G_row_values,value_1u)
                                
                                CALL get_DS_values(DS_indices,DS_values,D_rpointer,nn+1,DS_column_indices,DS_column_values)
                                CALL multiply(DS_column_indices,DS_column_values,G_row_indices,G_row_values,value_1l)
                                
                                !value_nxu=value_nxu
                                
                                CALL get_DS_values(DS_indices,DS_values,D_rpointer,nn+nx,DS_column_indices,DS_column_values)
                                CALL multiply(DS_column_indices,DS_column_values,G_row_indices,G_row_values,value_nxl)
                                
                        ELSE IF (ii.EQ.x_number_of_cells) THEN
                                
                                CALL get_DS_values(DS_indices,DS_values,D_rpointer,nn,DS_column_indices,DS_column_values)
                                value_0=SUM(DS_column_values*G_row_values)
                                
                                CALL get_DS_values(DS_indices,DS_values,D_rpointer,nn-1,DS_column_indices,DS_column_values)
                                CALL multiply(DS_column_indices,DS_column_values,G_row_indices,G_row_values,value_1u)
                                
                                !value_1l=value_1l
                                
                                !value_nxu=value_nxu
                                
                                CALL get_DS_values(DS_indices,DS_values,D_rpointer,nn+nx,DS_column_indices,DS_column_values)
                                CALL multiply(DS_column_indices,DS_column_values,G_row_indices,G_row_values,value_nxl)
                                
                        END IF
                        
                ELSE IF (jj.GT.1.AND.jj.LT.z_number_of_cells) THEN
                        
                        IF (ii.EQ.1) THEN
                                
                                CALL get_DS_values(DS_indices,DS_values,D_rpointer,nn,DS_column_indices,DS_column_values)
                                value_0=SUM(DS_column_values*G_row_values)
                                
                                !value_1u=value_1u
                                
                                CALL get_DS_values(DS_indices,DS_values,D_rpointer,nn+1,DS_column_indices,DS_column_values)
                                CALL multiply(DS_column_indices,DS_column_values,G_row_indices,G_row_values,value_1l)
                                
                                CALL get_DS_values(DS_indices,DS_values,D_rpointer,nn-nx,DS_column_indices,DS_column_values)
                                CALL multiply(DS_column_indices,DS_column_values,G_row_indices,G_row_values,value_nxu)
                                
                                CALL get_DS_values(DS_indices,DS_values,D_rpointer,nn+nx,DS_column_indices,DS_column_values)
                                CALL multiply(DS_column_indices,DS_column_values,G_row_indices,G_row_values,value_nxl)
                                
                        ELSE IF (ii.GT.1.AND.ii.LT.x_number_of_cells) THEN
                                
                                CALL get_DS_values(DS_indices,DS_values,D_rpointer,nn,DS_column_indices,DS_column_values)
                                value_0=SUM(DS_column_values*G_row_values)
                                
                                CALL get_DS_values(DS_indices,DS_values,D_rpointer,nn-1,DS_column_indices,DS_column_values)
                                CALL multiply(DS_column_indices,DS_column_values,G_row_indices,G_row_values,value_1u)
                                
                                CALL get_DS_values(DS_indices,DS_values,D_rpointer,nn+1,DS_column_indices,DS_column_values)
                                CALL multiply(DS_column_indices,DS_column_values,G_row_indices,G_row_values,value_1l)
                                
                                CALL get_DS_values(DS_indices,DS_values,D_rpointer,nn-nx,DS_column_indices,DS_column_values)
                                CALL multiply(DS_column_indices,DS_column_values,G_row_indices,G_row_values,value_nxu)
                                
                                CALL get_DS_values(DS_indices,DS_values,D_rpointer,nn+nx,DS_column_indices,DS_column_values)
                                CALL multiply(DS_column_indices,DS_column_values,G_row_indices,G_row_values,value_nxl)
                                
                        ELSE IF (ii.EQ.x_number_of_cells) THEN
                                
                                CALL get_DS_values(DS_indices,DS_values,D_rpointer,nn,DS_column_indices,DS_column_values)
                                value_0=SUM(DS_column_values*G_row_values)
                                
                                CALL get_DS_values(DS_indices,DS_values,D_rpointer,nn-1,DS_column_indices,DS_column_values)
                                CALL multiply(DS_column_indices,DS_column_values,G_row_indices,G_row_values,value_1u)
                                
                                !value_1l=value_1l
                                
                                CALL get_DS_values(DS_indices,DS_values,D_rpointer,nn-nx,DS_column_indices,DS_column_values)
                                CALL multiply(DS_column_indices,DS_column_values,G_row_indices,G_row_values,value_nxu)
                                
                                CALL get_DS_values(DS_indices,DS_values,D_rpointer,nn+nx,DS_column_indices,DS_column_values)
                                CALL multiply(DS_column_indices,DS_column_values,G_row_indices,G_row_values,value_nxl)
                                
                        END IF
                        
                ELSE IF (jj.EQ.z_number_of_cells) THEN
                        
                        IF (ii.EQ.1) THEN
                                
                                CALL get_DS_values(DS_indices,DS_values,D_rpointer,nn,DS_column_indices,DS_column_values)
                                value_0=SUM(DS_column_values*G_row_values)
                                
                                !value_1u=value_1u
                                
                                CALL get_DS_values(DS_indices,DS_values,D_rpointer,nn+1,DS_column_indices,DS_column_values)
                                CALL multiply(DS_column_indices,DS_column_values,G_row_indices,G_row_values,value_1l)
                                
                                CALL get_DS_values(DS_indices,DS_values,D_rpointer,nn-nx,DS_column_indices,DS_column_values)
                                CALL multiply(DS_column_indices,DS_column_values,G_row_indices,G_row_values,value_nxu)
                                
                                !value_nxl=value_nxl
                                
                        ELSE IF (ii.GT.1.AND.ii.LT.x_number_of_cells) THEN
                                
                                CALL get_DS_values(DS_indices,DS_values,D_rpointer,nn,DS_column_indices,DS_column_values)
                                value_0=SUM(DS_column_values*G_row_values)
                                
                                CALL get_DS_values(DS_indices,DS_values,D_rpointer,nn-1,DS_column_indices,DS_column_values)
                                CALL multiply(DS_column_indices,DS_column_values,G_row_indices,G_row_values,value_1u)
                                
                                CALL get_DS_values(DS_indices,DS_values,D_rpointer,nn+1,DS_column_indices,DS_column_values)
                                CALL multiply(DS_column_indices,DS_column_values,G_row_indices,G_row_values,value_1l)
                                
                                CALL get_DS_values(DS_indices,DS_values,D_rpointer,nn-nx,DS_column_indices,DS_column_values)
                                CALL multiply(DS_column_indices,DS_column_values,G_row_indices,G_row_values,value_nxu)
                                
                                !value_nxl=value_nxl
                                
                        ELSE IF (ii.EQ.x_number_of_cells) THEN
                                
                                CALL get_DS_values(DS_indices,DS_values,D_rpointer,nn,DS_column_indices,DS_column_values)
                                value_0=SUM(DS_column_values*G_row_values)
                                
                                CALL get_DS_values(DS_indices,DS_values,D_rpointer,nn-1,DS_column_indices,DS_column_values)
                                CALL multiply(DS_column_indices,DS_column_values,G_row_indices,G_row_values,value_1u)
                                
                                !value_1l=value_1l
                                
                                CALL get_DS_values(DS_indices,DS_values,D_rpointer,nn-nx,DS_column_indices,DS_column_values)
                                CALL multiply(DS_column_indices,DS_column_values,G_row_indices,G_row_values,value_nxu)
                                
                                !value_nxl=value_nxl
                                
                        END IF
                END IF
                
                A_cpointer(nn)=cpoint+1
                
                IF (value_nxu.NE.0) THEN
                        entry_number=entry_number+1
                        cpoint=cpoint+1
                        A_indices(entry_number,:)=(/ nn-nx,nn /)
                        A_values(entry_number)=-value_nxu
                END IF
                
                IF (value_1u.NE.0) THEN
                        entry_number=entry_number+1
                        cpoint=cpoint+1
                        A_indices(entry_number,:)=(/ nn-1,nn /)
                        A_values(entry_number)=-value_1u
                END IF
                
                IF (value_0.NE.0) THEN
                        entry_number=entry_number+1
                        cpoint=cpoint+1
                        A_indices(entry_number,:)=(/ nn,nn /)
                        A_values(entry_number)=-value_0
                END IF
                
                IF (value_1l.NE.0) THEN
                        entry_number=entry_number+1
                        cpoint=cpoint+1
                        A_indices(entry_number,:)=(/ nn+1,nn /)
                        A_values(entry_number)=-value_1l
                END IF
                
                IF (value_nxl.NE.0) THEN
                        entry_number=entry_number+1
                        cpoint=cpoint+1
                        A_indices(entry_number,:)=(/ nn+nx,nn /)
                        A_values(entry_number)=-value_nxl
                END IF
                
                IF (ALLOCATED(G_row_indices)) DEALLOCATE(G_row_indices)
                IF (ALLOCATED(G_row_values)) DEALLOCATE(G_row_values)
                IF (ALLOCATED(DS_column_indices)) DEALLOCATE(DS_column_indices)
                IF (ALLOCATED(DS_column_values)) DEALLOCATE(DS_column_values)
                
        END DO
END DO
A_cpointer(nn+1)=cpoint+1

END SUBROUTINE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE get_G_values(G_indices,G_values,G_cpointer,column,G_row_indices,G_row_values)

IMPLICIT NONE
INTEGER, DIMENSION(:,:), INTENT(IN)                             :: G_indices
DOUBLE PRECISION, DIMENSION(:), INTENT(IN)                      :: G_values
INTEGER, DIMENSION(:), INTENT(IN)                               :: G_cpointer
INTEGER, INTENT(IN)                                             :: column
INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT)                 :: G_row_indices
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, INTENT(OUT)        :: G_row_values
INTEGER                                                         :: ini,fin,ii

IF (ALLOCATED(G_row_indices)) DEALLOCATE(G_row_indices)
IF (ALLOCATED(G_row_values)) DEALLOCATE(G_row_values)

ini=G_cpointer(column) ! index of first element in column "column" of G
fin=G_cpointer(column+1)-1 ! index of last element in column "column" of G

ALLOCATE(G_row_indices(SIZE( (/ (ii,ii=ini,fin) /) ))) ! get indices of nonzero rows of G in column "column"
G_row_indices=G_indices(ini:fin,1)
ALLOCATE(G_row_values(SIZE( (/ (ii,ii=ini,fin) /) ))) ! get values of nonzero rows of G in column "column"
G_row_values=G_values(ini:fin)

END SUBROUTINE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE get_DS_values(DS_indices,DS_values,DS_rpointer,row,DS_column_indices,DS_column_values)

IMPLICIT NONE
INTEGER, DIMENSION(:,:), INTENT(IN)                             :: DS_indices
DOUBLE PRECISION, DIMENSION(:), INTENT(IN)                      :: DS_values
INTEGER, DIMENSION(:), INTENT(IN)                               :: DS_rpointer
INTEGER, INTENT(IN)                                             :: row
INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT)                 :: DS_column_indices
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, INTENT(OUT)        :: DS_column_values
INTEGER                                                         :: ini,fin,ii

IF (ALLOCATED(DS_column_indices)) DEALLOCATE(DS_column_indices)
IF (ALLOCATED(DS_column_values)) DEALLOCATE(DS_column_values)

ini=DS_rpointer(row) ! index of first element in row "row" of DS
fin=DS_rpointer(row+1)-1 ! index of last element in row "row" of DS

ALLOCATE(DS_column_indices(SIZE( (/ (ii,ii=ini,fin) /) ))) ! get indices of nonzero columns of DS in row "row"
DS_column_indices=DS_indices(ini:fin,2)
ALLOCATE(DS_column_values(SIZE( (/ (ii,ii=ini,fin) /) ))) ! get values of nonzero columns of DS in row "row"
DS_column_values=DS_values(ini:fin)

END SUBROUTINE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE multiply(DS_column_indices,DS_column_values,G_row_indices,G_row_values,value)

IMPLICIT NONE
INTEGER, DIMENSION(:), INTENT(IN)                               :: DS_column_indices,G_row_indices
DOUBLE PRECISION, DIMENSION(:), INTENT(IN)                      :: DS_column_values,G_row_values
DOUBLE PRECISION, INTENT(OUT)                                   :: value
INTEGER                                                         :: ii,jj

value=0
DO ii=1,SIZE(DS_column_indices)
        DO jj=1,SIZE(G_row_indices)
                IF (DS_column_indices(ii).EQ.G_row_indices(jj)) THEN
                        value=value+DS_column_values(ii)*G_row_values(jj)
                END IF
        END DO
END DO

END SUBROUTINE

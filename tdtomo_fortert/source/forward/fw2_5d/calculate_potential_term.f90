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


SUBROUTINE calculate_potential_term()

USE modelling_2_5D, ONLY: potential_locations,x_number_of_cells,z_number_of_cells,&
                          potential_term__size,potential_term__indices,&
                          potential_term__values,potential_term__cpointer
USE acquisition, ONLY: number_of_measurements
IMPLICIT NONE
INTEGER, ALLOCATABLE                            :: P1_size(:),P2_size(:)
INTEGER, ALLOCATABLE                            :: P1_indices(:,:),P2_indices(:,:)
DOUBLE PRECISION, ALLOCATABLE                   :: P1_values(:),P2_values(:)
INTEGER, ALLOCATABLE                            :: P1_cpointer(:),P2_cpointer(:)
INTEGER                                         :: ini,fin,nn,ii,jj
INTEGER, DIMENSION(:,:), ALLOCATABLE            :: worki
INTEGER, ALLOCATABLE                            :: workc(:)
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE     :: workd
INTEGER, DIMENSION(:), ALLOCATABLE              :: rows
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE     :: values
INTEGER, ALLOCATABLE                            :: rows_unique(:)
DOUBLE PRECISION                                :: value
INTERFACE
        SUBROUTINE interpolate_electrodes(electrode_list,Q_size,Q_indices,Q_values,Q_cpointer)
                DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN)                    :: electrode_list
                INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT)                 :: Q_size
                INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT)               :: Q_indices
                DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, INTENT(OUT)        :: Q_values
                INTEGER, ALLOCATABLE, INTENT(OUT)                               :: Q_cpointer(:)
        END SUBROUTINE
END INTERFACE
INTERFACE
        SUBROUTINE coordinate_to_sparse_column(M_size,M_indices,M_values,M_cpointer)
                INTEGER, DIMENSION(:), INTENT(IN)               :: M_size
                INTEGER, DIMENSION(:,:), INTENT(INOUT)          :: M_indices
                DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT)   :: M_values
                INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: M_cpointer
        END SUBROUTINE
END INTERFACE
INTERFACE
        SUBROUTINE trim_integer_array(integers_in,integers_out)
                INTEGER, DIMENSION(:), INTENT(IN)               :: integers_in
                INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: integers_out
        END SUBROUTINE
END INTERFACE

IF (ALLOCATED(potential_term__size)) DEALLOCATE(potential_term__size)
IF (ALLOCATED(potential_term__indices)) DEALLOCATE(potential_term__indices)
IF (ALLOCATED(potential_term__values)) DEALLOCATE(potential_term__values)
IF (ALLOCATED(potential_term__cpointer)) DEALLOCATE(potential_term__cpointer)

! Calculate term from positive potential electrodes
CALL interpolate_electrodes(potential_locations(:,1:2),P1_size,P1_indices,P1_values,P1_cpointer)

! Calculate term from negative potential electrodes
CALL interpolate_electrodes(potential_locations(:,3:4),P2_size,P2_indices,P2_values,P2_cpointer)

! Merge the P1 and P2 terms
ALLOCATE(worki(SIZE(P1_values)+SIZE(P2_values),2))
ALLOCATE(workd(SIZE(P1_values)+SIZE(P2_values)))
ini=1
fin=SIZE(P1_values)
worki(ini:fin,:)=P1_indices
workd(ini:fin)=P1_values
ini=fin+1
fin=ini-1+SIZE(P2_values)
worki(ini:fin,:)=P2_indices
workd(ini:fin)=-P2_values

! Sort the elements by column and row and convert to compressed column format
CALL coordinate_to_sparse_column((/x_number_of_cells*z_number_of_cells,number_of_measurements/),&
                                 worki,workd,workc)

! Deallocate the P1 and P2 arrays - no longer needed
DEALLOCATE(P1_size,P2_size,P1_indices,P2_indices,P1_values,P2_values,P1_cpointer,P2_cpointer)

! Sum elements with the same indices in case they are present
ALLOCATE(potential_term__size(2))
potential_term__size=(/x_number_of_cells*z_number_of_cells,number_of_measurements/)
ALLOCATE(potential_term__indices(SIZE(P1_values)+SIZE(P2_values),2))
ALLOCATE(potential_term__values(SIZE(P1_values)+SIZE(P2_values)))
nn=0
DO ii=1,number_of_measurements
        ini=workc(ii)
        fin=workc(ii+1)-1
        IF (ini.EQ.fin) THEN
                nn=nn+1
                potential_term__indices(nn,:)=worki(ini,:)
                potential_term__values(nn)=workd(ini)
        ELSE IF (ini.NE.fin) THEN
                ALLOCATE(rows(SIZE( (/ (jj,jj=ini,fin) /) )))
                ALLOCATE(values(SIZE( (/ (jj,jj=ini,fin) /) )))
                rows=worki(ini:fin,1)
                values=workd(ini:fin)
                CALL trim_integer_array(rows,rows_unique)
                IF (SIZE(rows).EQ.SIZE(rows_unique)) THEN
                        DO jj=1,SIZE(rows)
                                IF (values(jj).NE.0) THEN
                                        nn=nn+1
                                        potential_term__indices(nn,:)=(/rows(jj),ii/)
                                        potential_term__values(nn)=values(jj)
                                END IF
                        END DO
                ELSE
                        DO jj=1,SIZE(rows_unique)
                                value=SUM(PACK(values,rows.EQ.rows_unique(jj)))
                                IF (value.NE.0) THEN
                                        nn=nn+1
                                        potential_term__indices(nn,:)=(/rows_unique(jj),ii/)
                                        potential_term__values(nn)=value
                                END IF
                        END DO
                END IF
                DEALLOCATE(rows,values)
        END IF
END DO
DEALLOCATE(rows_unique)
DEALLOCATE(worki,workd,workc)
ALLOCATE(worki(nn,2),workd(nn))
worki=potential_term__indices(1:nn,:)
workd=potential_term__values(1:nn)
DEALLOCATE(potential_term__indices,potential_term__values)
ALLOCATE(potential_term__indices(nn,2),potential_term__values(nn))
potential_term__indices=worki
potential_term__values=workd
DEALLOCATE(worki,workd)

! Sort the elements by column and row and convert to compressed column format
CALL coordinate_to_sparse_column(potential_term__size,potential_term__indices,potential_term__values,potential_term__cpointer)

END SUBROUTINE


SUBROUTINE trim_integer_array(integers_in,integers_out)

USE m_unirnk
IMPLICIT NONE
INTEGER, DIMENSION(:), INTENT(IN)               :: integers_in
INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: integers_out
INTEGER                                         :: number_of_uniques
INTEGER, DIMENSION(:), ALLOCATABLE              :: worki

ALLOCATE(worki(SIZE(integers_in)))
worki=integers_in

CALL I_UNIRNK(integers_in,worki,number_of_uniques)

IF (ALLOCATED(integers_out)) DEALLOCATE(integers_out)
ALLOCATE(integers_out(number_of_uniques))
integers_out=integers_in(worki(1:number_of_uniques))

DEALLOCATE(worki)

END SUBROUTINE

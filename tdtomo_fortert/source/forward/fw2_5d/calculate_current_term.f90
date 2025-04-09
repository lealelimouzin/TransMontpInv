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


SUBROUTINE calculate_current_term()

USE modelling_2_5D, ONLY: current_locations,current_term,number_of_sources,cell_areas,x_number_of_cells,z_number_of_cells
IMPLICIT NONE
INTEGER, ALLOCATABLE            :: Q_size(:)
INTEGER, ALLOCATABLE            :: Q_indices(:,:)
DOUBLE PRECISION, ALLOCATABLE   :: Q_values(:)
INTEGER, ALLOCATABLE            :: Q_cpointer(:)
INTEGER                         :: ii,jj,nn,row,column
INTERFACE
        SUBROUTINE interpolate_electrodes(electrode_list,Q_size,Q_indices,Q_values,Q_cpointer)
                DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN)                    :: electrode_list
                INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT)                 :: Q_size
                INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT)               :: Q_indices
                DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, INTENT(OUT)        :: Q_values
                INTEGER, ALLOCATABLE, INTENT(OUT)                               :: Q_cpointer(:)
        END SUBROUTINE
END INTERFACE

IF (ALLOCATED(current_term)) DEALLOCATE(current_term)

ALLOCATE(current_term(x_number_of_cells*z_number_of_cells,number_of_sources))
 current_term=0

! Calculate term from positive current electrodes
CALL interpolate_electrodes(current_locations(:,1:2),Q_size,Q_indices,Q_values,Q_cpointer)
DO nn=1,SIZE(Q_values)
        row=Q_indices(nn,1)
        column=Q_indices(nn,2)
        current_term(row,column)=current_term(row,column)+Q_values(nn)
END DO

DEALLOCATE(Q_size,Q_indices,Q_values,Q_cpointer)

! Calculate term from negative current electrodes
CALL interpolate_electrodes(current_locations(:,3:4),Q_size,Q_indices,Q_values,Q_cpointer)
DO nn=1,SIZE(Q_values)
        row=Q_indices(nn,1)
        column=Q_indices(nn,2)
        current_term(row,column)=current_term(row,column)-Q_values(nn)
END DO

! Divide by cell area
DO nn=1,number_of_sources
        current_term(:,nn)=current_term(:,nn)/RESHAPE(cell_areas,(/x_number_of_cells*z_number_of_cells/))
END DO

END SUBROUTINE
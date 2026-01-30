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


SUBROUTINE calculate_data(resistance_data,info)

USE modelling_2_5D, ONLY: number_of_coefficients,k_parameters,g_parameters,&
                          current_term,number_of_sources,current_electrode_numbers,&
                          potential_term__size,potential_term__indices,&
                          potential_term__values,potential_term__cpointer,&
                          x_number_of_cells,z_number_of_cells,resistivity,&
                          A_size,A_indices,A_cpointer,A_values,&
                          potential_field,potential_temp,&
                          C_size,C_indices,C_values,L_size,L_indices,L_values,&
                          factors,dgssv_initialized,dgssv_do_clear
USE acquisition, ONLY: number_of_measurements
IMPLICIT NONE
DOUBLE PRECISION, DIMENSION(:), INTENT(OUT)     :: resistance_data
INTEGER, INTENT(OUT)                            :: info
! INTEGER*8                                       :: factors
INTEGER                                         :: ii,jj,nn,ini,fin
INTEGER, DIMENSION(:), ALLOCATABLE              :: rows
INTEGER                                         :: column
REAL                                            :: start,finish
LOGICAL                                         :: verbose=.FALSE.

IF (.NOT.ALLOCATED(potential_field)) ALLOCATE(potential_field(x_number_of_cells*z_number_of_cells,number_of_sources))
potential_field=0

IF (.NOT.ALLOCATED(potential_temp)) ALLOCATE(potential_temp(x_number_of_cells*z_number_of_cells,number_of_sources))

 C_values=RESHAPE(1/resistivity,(/x_number_of_cells*z_number_of_cells/))

DO nn=1,number_of_coefficients
        
        ! Sum A and (k_parameters(nn)**2)*C (diagonal)
        L_size=A_size
        L_indices=A_indices
        L_values=A_values
        DO ii=1,SIZE(A_values)
                IF (A_indices(ii,1).EQ.A_indices(ii,2)) THEN ! if we are on a diagonal element
                        L_values(ii)=L_values(ii)+(k_parameters(nn)**2)*C_values(A_indices(ii,1))
                END IF
        END DO
        
        potential_temp=0.5*current_term
        
        ! Initialise the linear solver from SuperLU
        IF (.NOT.dgssv_initialized) THEN
                IF (verbose) CALL cpu_time(start)
                CALL c_fortran_dgssv(1,L_size(1),SIZE(L_values),number_of_sources,&
                                L_values,L_indices(:,1),A_cpointer,potential_temp,L_size(2),factors,info)
                IF (verbose) WRITE(*,*)'info=',info
                IF (verbose) CALL cpu_time(finish)
                IF (verbose) PRINT '("Time for dgssv_1 = ",f12.9," seconds.")',finish-start
                dgssv_initialized=.TRUE.
        END IF
        
        ! Solve the linear system using SuperLU
        IF (verbose) CALL cpu_time(start)
        CALL c_fortran_dgssv(2,L_size(1),SIZE(L_values),number_of_sources,&
                               L_values,L_indices(:,1),A_cpointer,potential_temp,L_size(2),factors,info)
        IF (verbose) WRITE(*,*)'info=',info
        IF (verbose) CALL cpu_time(finish)
        IF (verbose) PRINT '("Time for dgssv_2 = ",f12.9," seconds.")',finish-start
        potential_field=potential_field+g_parameters(nn)*potential_temp
        
        ! Free the memory
        IF (dgssv_do_clear) THEN
                IF (verbose) CALL cpu_time(start)
                CALL c_fortran_dgssv(3,L_size(1),SIZE(L_values),number_of_sources,&
                                L_values,L_indices(:,1),A_cpointer,potential_temp,L_size(2),factors,info)
                IF (verbose) WRITE(*,*)'info=',info
                IF (verbose) CALL cpu_time(finish)
                IF (verbose) PRINT '("Time for dgssv_3 = ",f12.9," seconds.")',finish-start
        END IF
        
END DO

DO nn=1,number_of_measurements
        
        ini=potential_term__cpointer(nn)
        fin=potential_term__cpointer(nn+1)-1
        
        column=current_electrode_numbers(nn)
        
        resistance_data(nn)=SUM(potential_field(potential_term__indices(ini:fin,1),column)*potential_term__values(ini:fin))
        
END DO

END SUBROUTINE
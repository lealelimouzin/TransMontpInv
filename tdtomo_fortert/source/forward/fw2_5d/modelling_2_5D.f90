! ********************************************************************* !
! Forward modelling module for electrical resistivity in 2.5D.
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


MODULE modelling_2_5D

IMPLICIT NONE
! These parameters are needed for modelling
INTEGER, DIMENSION(:), ALLOCATABLE, SAVE                :: current_electrode_numbers
INTEGER, SAVE                                           :: number_of_sources
DOUBLE PRECISION, ALLOCATABLE, SAVE                     :: current_term(:,:)
INTEGER, ALLOCATABLE, SAVE                              :: potential_term__size(:)
INTEGER, ALLOCATABLE, SAVE                              :: potential_term__indices(:,:)
DOUBLE PRECISION, ALLOCATABLE, SAVE                     :: potential_term__values(:)
INTEGER, ALLOCATABLE, SAVE                              :: potential_term__cpointer(:)
INTEGER, DIMENSION(:), ALLOCATABLE, SAVE                :: D_size,G_size,DS_size
INTEGER, DIMENSION(:,:), ALLOCATABLE, SAVE              :: D_indices,G_indices,DS_indices
INTEGER, ALLOCATABLE, SAVE                              :: D_rpointer(:),G_cpointer(:)
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, SAVE       :: D_values,G_values,DS_values
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, SAVE     :: cell_areas
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, SAVE       :: x_spacings,z_spacings
INTEGER, SAVE                                           :: x_number_of_cells,z_number_of_cells
INTEGER, SAVE                                           :: number_of_coefficients
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, SAVE       :: k_parameters,g_parameters
INTEGER, SAVE                                           :: Nx,Nz,nax,naz
INTEGER, DIMENSION(:,:), ALLOCATABLE, SAVE              :: GRDax,GRDaz
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, SAVE       :: x_centers,z_centers
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, SAVE     :: resistivity
INTEGER, SAVE                                           :: S_size,C_size
INTEGER, DIMENSION(:), ALLOCATABLE, SAVE                :: S_indices,C_indices
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, SAVE       :: S_values,C_values
INTEGER, SAVE                                           :: number_of_nonzeros
INTEGER, DIMENSION(:), ALLOCATABLE, SAVE                :: A_size,A_cpointer
INTEGER, DIMENSION(:,:), ALLOCATABLE, SAVE              :: A_indices
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, SAVE       :: A_values
INTEGER, DIMENSION(:), ALLOCATABLE, SAVE                :: L_size
INTEGER, DIMENSION(:,:), ALLOCATABLE, SAVE              :: L_indices
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, SAVE       :: L_values
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, SAVE     :: potential_field,potential_temp
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, SAVE     :: resistivity_on_x_face,resistivity_on_z_face
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, SAVE     :: x_face_size,z_face_size
INTEGER, DIMENSION(:), ALLOCATABLE, SAVE                :: x_rows,x_columns,z_rows,z_columns
! These parameters are for set-up
CHARACTER (LEN=30)                                      :: fw2_5d_file
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE           :: current_locations,potential_locations
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE             :: x_edges,z_edges
INTEGER                                                 :: x_number_of_edges,z_number_of_edges
INTEGER                                                 :: np
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE           :: Dx,Dz
INTEGER, DIMENSION(:,:), ALLOCATABLE                    :: GRDp
LOGICAL                                                 :: verbose=.FALSE.
! These parameters are used for the call to subroutine c_fortran_dgssv
INTEGER*8, SAVE                                         :: factors
LOGICAL, SAVE                                           :: dgssv_initialized=.FALSE., dgssv_do_clear=.FALSE.


CONTAINS


SUBROUTINE initialize_fw2_5d()

IMPLICIT NONE
INTERFACE
        SUBROUTINE get_data()
        END SUBROUTINE
END INTERFACE
INTERFACE
        SUBROUTINE read_grid2D()
        END SUBROUTINE
END INTERFACE
INTERFACE
        SUBROUTINE calculate_divergence2D()
        END SUBROUTINE
END INTERFACE
INTERFACE
        SUBROUTINE calculate_gradient2D()
        END SUBROUTINE
END INTERFACE
INTERFACE
        SUBROUTINE get_fourier_coefficients()
        END SUBROUTINE
END INTERFACE
INTERFACE
        SUBROUTINE calculate_current_term()
        END SUBROUTINE
END INTERFACE
INTERFACE
        SUBROUTINE calculate_potential_term()
        END SUBROUTINE
END INTERFACE

IF (verbose) WRITE(*,*)'Formatting data...'
CALL get_data()

IF (verbose) WRITE(*,*)'Reading grid...'
CALL read_grid2D()

IF (verbose) WRITE(*,*)'Calculating divergence operator...'
CALL calculate_divergence2D()

IF (verbose) WRITE(*,*)'Calculating gradient operator...'
CALL calculate_gradient2D()

IF (verbose) WRITE(*,*)'Determining Fourier coefficients...'
CALL get_fourier_coefficients()
IF (verbose) WRITE(*,*)k_parameters
IF (verbose) WRITE(*,*)g_parameters

IF (verbose) WRITE(*,*)'Calculating current term...'
CALL calculate_current_term()

IF (verbose) WRITE(*,*)'Calculating potential term...'
CALL calculate_potential_term()

DEALLOCATE(current_locations,potential_locations)
DEALLOCATE(x_edges,z_edges)
DEALLOCATE(Dx,Dz,GRDp)

END SUBROUTINE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE run_fw2_5d(resistance_data,info)

IMPLICIT NONE
DOUBLE PRECISION, DIMENSION(:), INTENT(OUT)     :: resistance_data
INTEGER                                         :: info
INTEGER                                         :: ii,jj
INTERFACE
        SUBROUTINE average_conductivity_on_face()
        END SUBROUTINE
END INTERFACE
INTERFACE
        SUBROUTINE calculate_DS()
        END SUBROUTINE
END INTERFACE
INTERFACE
        SUBROUTINE sparse_diag_multiply(A_size,A_indices,A_values,D_size,D_indices,D_values,C_size,C_indices,C_values,A_cpointer)
                INTEGER, DIMENSION(:), INTENT(IN)                               :: A_size
                INTEGER, INTENT(IN)                                             :: D_size
                INTEGER, DIMENSION(:,:), INTENT(IN)                             :: A_indices
                INTEGER, DIMENSION(:), INTENT(INOUT)                            :: D_indices
                DOUBLE PRECISION, DIMENSION(:), INTENT(IN)                      :: A_values
                DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT)                   :: D_values
                INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT)                 :: C_size
                INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT)               :: C_indices
                DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, INTENT(OUT)        :: C_values
                INTEGER, DIMENSION(:), OPTIONAL, INTENT(IN)                     :: A_cpointer
        END SUBROUTINE
END INTERFACE
INTERFACE
        SUBROUTINE calculate_A()
        END SUBROUTINE
END INTERFACE
INTERFACE
        SUBROUTINE calculate_data(resistance_data,info)
                DOUBLE PRECISION, DIMENSION(:), INTENT(OUT)     :: resistance_data
                INTEGER, INTENT(OUT)                            :: info
        END SUBROUTINE
END INTERFACE

IF (verbose) WRITE(*,*)'Calculating average conductivity on each face...'
CALL average_conductivity_on_face()

IF (verbose) WRITE(*,*)'Calculating D*S...'
CALL calculate_DS()

IF (verbose) WRITE(*,*)'Calculating A...'
CALL calculate_A()

IF (verbose) WRITE(*,*)'Calculating data...'
CALL calculate_data(resistance_data,info)

END SUBROUTINE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE clear_memory_fw2_5d()

IMPLICIT NONE
INTEGER         :: info=0
REAL            :: start,finish

IF (verbose) CALL cpu_time(start)

CALL c_fortran_dgssv(3,0,0,0,0,0,0,0,0,factors,info)

IF (verbose) WRITE(*,*)'info=',info
IF (verbose) CALL cpu_time(finish)
IF (verbose) PRINT '("Time for dgssv_3 = ",f12.9," seconds.")',finish-start

END SUBROUTINE

END MODULE

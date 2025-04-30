! ********************************************************************* !
! Forward modelling module (part of TERT package)
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


MODULE forward

USE input_files, ONLY: modelling_method,forward_file,path_fwdFrac
USE chains, ONLY: rho_mat,opt_parameters,x_min,x_max,z_min,z_max
USE acquisition, ONLY: number_of_measurements
USE string_conversion

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE forward_modelling(number_of_cells,voronoi_model,modelled_data,info,chain_number,iteration)

IMPLICIT NONE
INTEGER, INTENT(IN)                                                             :: number_of_cells,chain_number,iteration
DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN)                                    :: voronoi_model
DOUBLE PRECISION, DIMENSION(number_of_measurements), INTENT(OUT)                :: modelled_data
INTEGER, INTENT(OUT)                                                            :: info
CHARACTER(100)                                                                  :: nametmp
CHARACTER(200)                                                                  :: namedir


IF (modelling_method.EQ.1) THEN
                
        CALL forward_with_fw2_5d(number_of_cells,voronoi_model(1:number_of_cells,[1,3,4]),modelled_data,info)

ELSEIF (modelling_method.EQ.2) THEN

        ! Go to the chain directory
        nametmp = "results_chains/"//"res_chain_"//TRIM(int2str(chain_number))
        namedir = ADJUSTL(ADJUSTR(path_fwdFrac)//nametmp)

        IF (iteration.GT.0) THEN
                CALL write2DFN(number_of_cells,voronoi_model,namedir)
        END IF

        CALL forward_with_ert_ddp_2_5d(namedir,chain_number,modelled_data,info)

ELSE      
        WRITE(*,*)
        WRITE(*,*)'********************************************'
        WRITE(*,*)'ERROR in subroutine forward_modelling'
        WRITE(*,*)'--------------------------------------------'
        WRITE(*,*)'Modelling method ',TRIM(int2str(modelling_method)),' does not exist!'
        WRITE(*,*)'Check input file mksamples.in'
        WRITE(*,*)'TERMINATING PROGRAM!!'
        WRITE(*,*)'********************************************'
        WRITE(*,*)
        STOP
END IF

END SUBROUTINE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Modelling method 1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE forward_with_fw2_5d(number_of_cells,voronoi_model,resistance_data,info)

USE modelling_2_5D
USE voronoi_operations
IMPLICIT NONE
INTEGER, INTENT(IN)                                     :: number_of_cells
DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN)            :: voronoi_model
DOUBLE PRECISION, DIMENSION(:), INTENT(OUT)             :: resistance_data
INTEGER, INTENT(OUT)                                    :: info
INTERFACE
        SUBROUTINE create_resistivity_grid2D(ncell,voronoi_model)
                INTEGER, INTENT(IN)                                     :: ncell
                DOUBLE PRECISION, DIMENSION(ncell,3), INTENT(IN)        :: voronoi_model
        END SUBROUTINE
END INTERFACE

IF (.NOT.ALLOCATED(cell_areas)) THEN
        
        ! Define input file name for forward modelling
        fw2_5d_file=forward_file
        
        ! Initialise the modelling grid and various modelling parameters
        CALL initialize_fw2_5d()
        
END IF

! Create 2D resistivity model
CALL create_resistivity_grid2D(number_of_cells,voronoi_model)

! Run forward modelling
CALL run_fw2_5d(resistance_data,info)

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Modelling method 2 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE write2DFN(number_of_fractures,fracture_model,namedir)
        IMPLICIT NONE 
        INTEGER,INTENT(IN)                                      :: number_of_fractures
        DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN)            :: fracture_model
        DOUBLE PRECISION                                        :: length,x1,z,x2,sigma,b,z1,z2,angle,aperture
        CHARACTER(200)                                          :: namedir
        CHARACTER(400)                                          :: filenametmp
        DOUBLE PRECISION, DIMENSION(6)                          :: fracture
        DOUBLE PRECISION, DIMENSION(5) :: bounds
        INTEGER                                                 :: J 
        
        filenametmp = ADJUSTL(ADJUSTR(namedir)//TRIM("/input/testDFN/DFN.txt"))

        IF (opt_parameters.EQ.0) THEN ! No additional parameter to be inverted
                open (unit=40, file=filenametmp, status='replace', action='write')
                        write(40,'(A14)') TRIM('DETERMINISTIC1')
                        write(40,*) int2str(number_of_fractures)
                        DO J=1,number_of_fractures
                                length = fracture_model(J,4)
                                x1 = fracture_model(J,1) - length/2
                                z = fracture_model(J,3)
                                x2 = fracture_model(J,1) + length/2
                                b = 0.03
                                sigma = 2
                                fracture = [x1, z, x2, z, b, sigma]
                                write(40,'(F6.3,1X,F6.3,1X,F6.3,1X,F6.3,1X,F4.2,1X,F3.1)') fracture
                        END DO
                close(40)
        ELSE ! Add the inversion of either the angle or the aperture
                open (unit=40, file=filenametmp, status='replace', action='write')
                        write(40,'(A14)') TRIM('DETERMINISTIC1')
                        write(40,*) int2str(number_of_fractures)
                        DO J=1,number_of_fractures
                                IF (opt_parameters.EQ.1) THEN !angle
                                        ! Read the parameters determined in chains.f90
                                        angle = fracture_model(J,5)
                                        length = fracture_model(J,4)
                                        x1 = fracture_model(J,1) - COS(angle)*length/2
                                        z1 = fracture_model(J,3) - SIN(angle)*length/2
                                        x2 = fracture_model(J,1) + COS(angle)*length/2
                                        z2 = fracture_model(J,3) + SIN(angle)*length/2
                                        b = 0.03
                                        ! Check if the fracture is inside the domain
                                        bounds = check_edges(x1,x2,z1,z2,angle)
                                        x1 = bounds(1) ; x2 = bounds(2)
                                        z1 = bounds(3) ; z2 = bounds(4)
                                        angle = bounds(5)
                                ELSE ! aperture
                                        length = fracture_model(J,4)
                                        x1 = fracture_model(J,1) - length/2
                                        z1 = fracture_model(J,3)
                                        x2 = fracture_model(J,1) + length/2
                                        z2 = fracture_model(J,3)
                                        b = fracture_model(J,5)
                                END IF
                                sigma = 2
                                fracture = [x1, z1, x2, z2, b, sigma]
                                ! Store the DFN in a text file
                                write(40,'(F6.3,1X,F6.3,1X,F6.3,1X,F6.3,1X,F4.2,1X,F3.1)') fracture
                        END DO
                close(40)
        END IF

END SUBROUTINE

SUBROUTINE forward_with_ert_ddp_2_5d(namedir,chain_number,modelled_data,info)
        IMPLICIT NONE
        INTEGER, INTENT(IN)                                     :: chain_number
        INTEGER                                                 :: CSTAT,I
        CHARACTER(200)                                          :: namedir,cdtmp
        CHARACTER(300)                                          :: textfileDFN, textfileEPM,cwd
        CHARACTER(300)                                          :: textfile
        DOUBLE PRECISION, DIMENSION(:), INTENT(OUT)             :: modelled_data
        DOUBLE PRECISION, DIMENSION(number_of_measurements)     :: modelled_data_epm, modelled_data_dfn
        INTEGER, INTENT(OUT)                                    :: info

        ! Go to the result directory
        CALL CHDIR (namedir)
        
        ! Run forward model
        CALL EXECUTE_COMMAND_LINE ("python3 ../../main.py",CMDSTAT=CSTAT)

        IF ((CSTAT > 0).OR.(CSTAT < 0)) THEN
                PRINT *, "Command execution failed when running Python in forward.f90"
        END IF

        ! Store current apparent resistivity model (ratio of simulated data)
        textfileDFN = ADJUSTL(ADJUSTR(namedir)//TRIM("/results_fwd/rho_a_DFN.txt"))
        open (unit=88, file=textfileDFN, status='old', action='read')
                do I=1,number_of_measurements
                        read(88,*) modelled_data_dfn(I) 
                end do
        close(88)

        textfileEPM = ADJUSTL(ADJUSTR(namedir)//TRIM("/results_fwd/rho_a_EPM.txt"))
        open (unit=89, file=textfileEPM, status='old', action='read')
                do I=1,number_of_measurements
                        read(89,*) modelled_data_epm(I) 
                end do
        close(89)
        
        modelled_data = modelled_data_dfn-(modelled_data_epm-rho_mat) ! ratio of simulated data

        ! Go back to the bin directory
        cdtmp = "../../../../../bin/"
        CALL CHDIR (TRIM(cdtmp))

        info = 0
        
END SUBROUTINE

RECURSIVE FUNCTION check_edges(x1,x2,z1,z2,angle) result(bounds)
        DOUBLE PRECISION, INTENT(IN):: x1,x2,z1,z2,angle
        DOUBLE PRECISION :: x1_new,x2_new,z1_new,z2_new,angle_new
        DOUBLE PRECISION,DIMENSION(5):: bounds
        
        x1_new = x1 ; x2_new = x2 ;z1_new = z1 ; z2_new = z2 ; angle_new = angle 

        IF (x1.LT.x_min) THEN
                x1_new = x_min 
                angle_new = ATAN( (z2_new-z1_new) / (x2_new-x1_new))
                bounds(:) = check_edges(x1_new,x2_new,z1_new,z2_new,angle_new)
        ELSEIF (x2.LT.x_min) THEN
                x2_new = x_min 
                angle_new = ATAN( (z1_new-z2_new) / (x1_new-x2_new))
                bounds(:) = check_edges(x1_new,x2_new,z1_new,z2_new,angle_new)
        ELSEIF (x1.GT.x_max) THEN
                x1_new = x_max
                angle_new = ATAN( (z1_new-z2_new) / (x1_new-x2_new))
                bounds(:) = check_edges(x1_new,x2_new,z1_new,z2_new,angle_new)
        ELSEIF (x2.GT.x_max) THEN
                x2_new = x_max 
                angle_new = ATAN( (z2_new-z1_new) / (x2_new-x1_new))
                bounds(:) = check_edges(x1_new,x2_new,z1_new,z2_new,angle_new)      
        ELSEIF (z1.LT.z_min) THEN
                z1_new = z_min
                angle_new = ATAN( (z2_new-z1_new) / (x2_new-x1_new))
                bounds(:) = check_edges(x1_new,x2_new,z1_new,z2_new,angle_new)
        ELSEIF (z2.LT.z_min) THEN
                z2_new = z_min
                angle_new = ATAN( (z2_new-z1_new) / (x2_new-x1_new))
                bounds(:) = check_edges(x1_new,x2_new,z1_new,z2_new,angle_new)
        ELSEIF (z1.GT.z_max) THEN
                z1_new = z_max
                angle_new = ATAN( (z2_new-z1_new) / (x2_new-x1_new))
                bounds(:) = check_edges(x1_new,x2_new,z1_new,z2_new,angle_new)
        ELSEIF (z2.GT.z_max) THEN
                z2_new = z_max
                angle_new = ATAN( (z2_new-z1_new) / (x2_new-x1_new))
                bounds(:) = check_edges(x1_new,x2_new,z1_new,z2_new,angle_new)
        ELSE
                bounds(:) = [x1_new,x2_new,z1_new,z2_new,angle_new]
        END IF
END FUNCTION check_edges


END MODULE

! ********************************************************************* !
! Subroutine to initialise Voronoi models in mksamples.f90
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

! ********************************************************************* !
! Scripts adapted to inversion of preferential pathways in TransMontpInv
! package.
! Copyright (C) 2026 Léa Lelimouzin
! ********************************************************************* !

SUBROUTINE initialize_models()

USE chains
USE random_generator
USE debug
USE input_files, ONLY: modelling_method,path_fwdFrac
USE string_conversion
IMPLICIT NONE
INTEGER                 :: nn,CSTAT,I
INTEGER                 :: random_so_far
DOUBLE PRECISION        :: random_deviate
REAL,DIMENSION(5,6)     :: tempo
CHARACTER(100)          :: namefile,nametmp
CHARACTER(300)          :: namedir
CHARACTER(400)          :: filename
LOGICAL                 :: Existent=.FALSE.

IF (debug_mode) WRITE(debug_unit,*)'--- Chain ',TRIM(chain_name),': initializing models ---'

IF (run_number.EQ.1) THEN

        voronoi_current=0

        !----- If Voronoi tesselation (fw2_5d) -----
        IF (modelling_method.EQ.1) THEN
                ! Initial random number of cells 
                ncell_current=INT(ncell_min+unirand()*(ncell_max-ncell_min))
        
                ! Place the cell nuclei at random locations and assign them random values
                DO nn=1,ncell_current
                        voronoi_current(nn,1)=x_min+unirand()*(x_max-x_min)             
                        voronoi_current(nn,2)=y_min+unirand()*(y_max-y_min)             
                        voronoi_current(nn,3)=z_min+unirand()*(z_max-z_min)              
                        voronoi_current(nn,4)=value_min+unirand()*(value_max-value_min) 
                END DO
        
        !----- If fracture representation (ert_ddp_2_5d) -----
        ELSE IF (modelling_method.EQ.2) THEN
                ! Create a new directory
                nametmp = "results_chains/"//"res_chain_"//TRIM(int2str(chain_number))
                namedir = ADJUSTL(ADJUSTR(path_fwdFrac)//nametmp)

                CALL EXECUTE_COMMAND_LINE ("mkdir -p "//namedir)

                CALL CHDIR (namedir)

                CALL EXECUTE_COMMAND_LINE ("python3 ../../initialize_fracmodels.py",CMDSTAT=CSTAT)

                IF ((CSTAT > 0).OR.(CSTAT < 0)) THEN
                        PRINT *, "Command execution failed when running Python in initialize_models.f90"
                END IF

                ! Fill up current arrays
                filename = ADJUSTL(ADJUSTR(namedir)//TRIM("/input/testDFN/DFN.txt"))
                INQUIRE(FILE=filename, EXIST=Existent)
                IF (EXistent) THEN
                        open (unit=15, file=filename, status='old', action='read')
                        read(15,*)
                        read(15,*) ncell_current
                                do I=1,ncell_current
                                        read(15,*) tempo(I,:) ! fractures from DFN.txt
                                end do
                        close(15)
                        DO nn=1,ncell_current
                                voronoi_current(nn,1) = tempo(nn,1) + ABS(tempo(nn,3)-tempo(nn,1))/2 ! x centre
                                voronoi_current(nn,3) = tempo(nn,2) + ABS(tempo(nn,4)-tempo(nn,2))/2 ! z centre
                                !voronoi_current(nn,4) = ABS(tempo(nn,3)-tempo(nn,1)) ! length
                                voronoi_current(nn,4) = SQRT( (tempo(nn,3)-tempo(nn,1))**2 + (tempo(nn,4)-tempo(nn,2))**2 ) ! length
                                IF (opt_parameters.EQ.1) THEN 
                                        !voronoi_current(nn,5) = 0 ! angle
                                        voronoi_current(nn,5) = ACOS( ABS(voronoi_current(nn,1)-tempo(nn,1)) &
                                                                * 2/voronoi_current(nn,4) ) ! angle
                                ELSEIF (opt_parameters.EQ.2) THEN
                                        voronoi_current(nn,5) = 0.03 ! aperture
                                END IF
                        end do
                END IF
        END IF

        ! Initialise noise parameters
        IF (noise_type.EQ.0) THEN
                DO nn=1,number_of_datasets
                        sigma_current(nn)=sigma_min+unirand()*(sigma_max-sigma_min)
                END DO
        ELSE IF (noise_type.EQ.1.OR.noise_type.EQ.2.OR.noise_type.EQ.5) THEN
                DO nn=1,number_of_datasets
                        a_current(nn)=a_min+unirand()*(a_max-a_min)
                        b_current(nn)=b_min+unirand()*(b_max-b_min)
                END DO
        ELSE IF (noise_type.EQ.3.OR.noise_type.EQ.6) THEN
                DO nn=1,number_of_datasets
                        !lambda_current(nn)=lambda_min+unirand()*(lambda_max-lambda_min)
                        lambda_current(nn)=1 ! 18/07/25
                END DO
        END IF

        ! Initialise acceptance ratios
        acceptance_total=0
        acceptance_birth=0
        acceptance_death=0
        acceptance_move=0
        acceptance_move_main=0
        acceptance_move_delayed=0
        acceptance_value=0
        acceptance_value_main=0
        acceptance_value_delayed=0
        acceptance_noise=0
        IF (noise_type.EQ.1.OR.noise_type.EQ.2.OR.noise_type.EQ.5) THEN
                acceptance_a=0
                acceptance_b=0
        END IF
        acceptance_tempering=0
        IF (opt_parameters.GT.0) THEN
                acceptance_opt=0
                acceptance_opt_main=0
                acceptance_opt_delayed=0
        END IF
        
        ! Initialise the counter for accepted samples
        accepted_in_run=0
        accepted_so_far=0
        
        ! Initialise the overall sample counter
        samples_so_far=0
        samples_at_start=0
        
        ! Initialise temperatures
        IF (parallel_tempering.EQ.1) THEN
                temperature_values=temperature_values0
                temperature_indices=temperature_indices0
        END IF
        
ELSE IF (run_number.GT.1) THEN

        ! Open file containing last model of previous run
        OPEN(UNIT=80,FILE=TRIM("results/")//TRIM(last_tmp),STATUS='OLD',FORM='UNFORMATTED')
        
        !-------- MODELLING METHOD 1----------
        IF (modelling_method.EQ.1) THEN
                ! Read Voronoi tessellation
                READ(80)ncell_current
                ncell_proposed=ncell_current
                DO nn=1,ncell_current
                        READ(80)voronoi_current(nn,1),voronoi_current(nn,2),voronoi_current(nn,3),voronoi_current(nn,4)
                        voronoi_proposed(nn,:)=voronoi_current(nn,:)
                        voronoi_delayed(nn,:)=voronoi_current(nn,:)
                END DO
        !-------------------------------------
        ELSE IF (modelling_method.EQ.2) THEN 
                ! Read fractures or preferential pathways
                READ(80)ncell_current
                ncell_proposed=ncell_current
                IF (opt_parameters.EQ.0) THEN 
                        DO nn=1,ncell_current
                                READ(80)voronoi_current(nn,1),voronoi_current(nn,2),voronoi_current(nn,3),voronoi_current(nn,4)
                                voronoi_proposed(nn,:)=voronoi_current(nn,:)
                                voronoi_delayed(nn,:)=voronoi_current(nn,:)
                        END DO
                ELSE
                        DO nn=1,ncell_current
                                READ(80)voronoi_current(nn,1),voronoi_current(nn,2),voronoi_current(nn,3),voronoi_current(nn,4),&
                                voronoi_current(nn,5)
                                voronoi_proposed(nn,:)=voronoi_current(nn,:)
                                voronoi_delayed(nn,:)=voronoi_current(nn,:)
                        END DO
                END IF
        END IF
        
        ! Read noise parameters
        IF (noise_type.EQ.0) THEN
                DO nn=1,number_of_datasets
                        READ(80)sigma_current(nn)
                END DO
        ELSE IF (noise_type.EQ.1.OR.noise_type.EQ.2.OR.noise_type.EQ.5) THEN
                DO nn=1,number_of_datasets
                        READ(80)a_current(nn),b_current(nn)
                END DO
        ELSE IF (noise_type.EQ.3.OR.noise_type.EQ.6) THEN
                DO nn=1,number_of_datasets
                        READ(80)lambda_current(nn)
                END DO
        END IF
        
        ! Read acceptance ratios
        READ(80)accepted_so_far
        accepted_so_far_other=accepted_so_far
        accepted_in_run=0
        READ(80)acceptance_total(1),acceptance_total(2),acceptance_total(3)
        READ(80)acceptance_birth(1),acceptance_birth(2),acceptance_birth(3)
        READ(80)acceptance_death(1),acceptance_death(2),acceptance_death(3)
        READ(80)acceptance_move(1),acceptance_move(2),acceptance_move(3)
        READ(80)acceptance_value(1),acceptance_value(2),acceptance_value(3)
        READ(80)acceptance_move_main(1),acceptance_move_main(2),acceptance_move_main(3)
        READ(80)acceptance_move_delayed(1),acceptance_move_delayed(2),acceptance_move_delayed(3)
        READ(80)acceptance_value_main(1),acceptance_value_main(2),acceptance_value_main(3)
        READ(80)acceptance_value_delayed(1),acceptance_value_delayed(2),acceptance_value_delayed(3)
        READ(80)acceptance_noise(1),acceptance_noise(2),acceptance_noise(3)
        IF (noise_type.EQ.1.OR.noise_type.EQ.2.OR.noise_type.EQ.5) THEN
                READ(80)acceptance_a(1),acceptance_a(2),acceptance_a(3)
                READ(80)acceptance_b(1),acceptance_b(2),acceptance_b(3)
        END IF
        READ(80)acceptance_tempering(1),acceptance_tempering(2),acceptance_tempering(3)
        IF (opt_parameters.GT.0) THEN
                READ(80)acceptance_opt(1),acceptance_opt(2),acceptance_opt(3) 
                READ(80)acceptance_opt_main(1),acceptance_opt_main(2),acceptance_opt_main(3)
                READ(80)acceptance_opt_delayed(1),acceptance_opt_delayed(2),acceptance_opt_delayed(3)
        END IF
        
        ! Read number of random numbers used so far and update counter
        READ(80)random_so_far
        DO nn=1,random_so_far
                random_deviate=unirand()
        END DO
        
        ! Read number of samples completed so far
        READ(80)samples_so_far
        samples_at_start=samples_so_far
        
        ! Read modelled data from last accepted model
        DO nn=1,number_of_measurements
                READ(80)data_current(nn)
        END DO
        
        ! Read temperatures
        IF (parallel_tempering.EQ.1) THEN
                DO nn=1,number_of_temperatures
                        READ(80)temperature_indices(nn),temperature_values(nn)
                END DO
        END IF
        
        CLOSE(80)
        
END IF

IF (debug_mode) WRITE(debug_unit,*)'--- Chain ',TRIM(chain_name),': initialized models ---'

END SUBROUTINE


SUBROUTINE save_first_model()

USE chains
USE debug
IMPLICIT NONE
INTEGER         :: nn
 
 IF (debug_mode) WRITE(debug_unit,*)'--- Chain ',TRIM(chain_name),': saving first model ---'
 
 voronoi_0=voronoi_current
 chain_ncell(1)=ncell_current
 DO nn=1,number_of_datasets
        IF (noise_type.EQ.0) THEN
                chain_noise(1,nn)=sigma_current(nn)
        ELSE IF (noise_type.EQ.1.OR.noise_type.EQ.2.OR.noise_type.EQ.5) THEN
                chain_noise(1,nn)=a_current(nn)
                chain_noise(1,nn+number_of_datasets)=b_current(nn)
        ELSE IF (noise_type.EQ.3.OR.noise_type.EQ.6) THEN
                chain_noise(1,nn)=lambda_current(nn)
        END IF
 END DO
 chain_misfits(1)=misfit_current
 
 IF (debug_mode) WRITE(debug_unit,*)'--- Chain ',TRIM(chain_name),': saved first model ---'
 
END SUBROUTINE
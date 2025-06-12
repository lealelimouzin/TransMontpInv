! ********************************************************************* !
! Module containing basic Markov chain parameters and subroutines
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
! The module has been modified to support the inversion of thin
! preferential pathways such as fractures.
! Copyrigth (C) 2025 Lea Lelimouzin
!
! The same licence as presented above applied to this adaptation of the
! code.
! ********************************************************************* ! 


MODULE chains

USE acquisition, ONLY: number_of_measurements,number_of_datasets
USE voronoi_operations
USE random_generator
USE debug
USE string_conversion

!USE stdlib_stats_distribution_uniform, ONLY: unidist => rvs_uniform
IMPLICIT NONE
! Prior
LOGICAL                                                 :: inside_prior,inside_prior_delayed
INTEGER, SAVE                                           :: ncell_min,ncell_max
DOUBLE PRECISION, SAVE                                  :: value_min,value_max
DOUBLE PRECISION, SAVE                                  :: x_min,x_max,y_min,y_max,z_min,z_max
INTEGER, SAVE                                           :: noise_type
DOUBLE PRECISION, SAVE                                  :: a_min,a_max,b_min,b_max
DOUBLE PRECISION, SAVE                                  :: lambda_min,lambda_max
DOUBLE PRECISION, SAVE                                  :: sigma_min,sigma_max
DOUBLE PRECISION, SAVE                                  :: sigma_mat,rho_mat
INTEGER,SAVE                                            :: opt_parameters
DOUBLE PRECISION,SAVE                                   :: opt_min,opt_max 
! Proposals
LOGICAL, SAVE                                           :: use_delayed_value,use_delayed_position
LOGICAL, SAVE                                           :: use_delayed_opt
DOUBLE PRECISION, SAVE                                  :: value_change,value_change_delayed
DOUBLE PRECISION, SAVE                                  :: x_change,x_change_delayed
DOUBLE PRECISION, SAVE                                  :: y_change,y_change_delayed
DOUBLE PRECISION, SAVE                                  :: z_change,z_change_delayed
DOUBLE PRECISION, SAVE                                  :: birth_change
DOUBLE PRECISION, SAVE                                  :: a_change,b_change
DOUBLE PRECISION, SAVE                                  :: lambda_change
DOUBLE PRECISION, SAVE                                  :: sigma_change
INTEGER, SAVE                                           :: likelihood_type
DOUBLE PRECISION                                        :: l_change
DOUBLE PRECISION, SAVE                                  :: opt_change,opt_change_delayed
! Models
INTEGER, SAVE                                           :: number_of_parameters
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE           :: voronoi_current,voronoi_proposed,voronoi_delayed,voronoi_0
INTEGER                                                 :: ncell_current,ncell_proposed,cell_number
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE             :: a_current,a_proposed,b_current,b_proposed
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE             :: lambda_current,lambda_proposed
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE             :: sigma_current,sigma_proposed
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE             :: noise_current,noise_proposed,noise_delayed
INTEGER                                                 :: dataset_number
LOGICAL                                                 :: delayed_value,delayed_position,delayed_opt
LOGICAL                                                 :: change_a,change_b,change_sigma,change_lambda
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE             :: data_current,data_proposed,data_delayed
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE             :: noise_x_current,noise_x_proposed,noise_x_delayed
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE           :: residuals
! Acceptance parameters
DOUBLE PRECISION                                        :: alpha_component,alpha,alpha_delayed
DOUBLE PRECISION                                        :: likelihood_current,likelihood_proposed,likelihood_delayed
DOUBLE PRECISION                                        :: likelihood_normalization
DOUBLE PRECISION                                        :: misfit_current,misfit_proposed,misfit_delayed
! Markov chain parameters
INTEGER                                                 :: step
INTEGER                                                 :: chain_number,number_of_chains
CHARACTER (LEN=30)                                      :: chain_name
INTEGER, SAVE                                           :: run_number
CHARACTER (LEN=30), SAVE                                :: run_name
DOUBLE PRECISION, SAVE                                  :: max_runtime
INTEGER, SAVE                                           :: number_of_samples,max_number_of_samples
INTEGER                                                 :: samples_in_run,samples_so_far,samples_at_start
INTEGER, SAVE                                           :: display_interval
INTEGER, SAVE                                           :: no_data
! Parallel tempering
INTEGER, SAVE                                           :: parallel_tempering,jump_type,tempering_start
INTEGER, SAVE                                           :: number_of_temperatures,number_of_1s,number_of_non1s
DOUBLE PRECISION                                        :: temperature
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, SAVE       :: temperature_values,temperature_values0
INTEGER, DIMENSION(:), ALLOCATABLE, SAVE                :: temperature_indices,temperature_indices0
DOUBLE PRECISION                                        :: likelihood_T1,likelihood_T2
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, SAVE       :: noise_T1,noise_T2
LOGICAL                                                 :: swap_temperatures
! Chain evolution
LOGICAL                                                 :: birth,death,move,value,noise,option
INTEGER, DIMENSION(:,:), ALLOCATABLE                    :: chain_steps
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE           :: chain_values
INTEGER, DIMENSION(:), ALLOCATABLE                      :: chain_indices
INTEGER, DIMENSION(:), ALLOCATABLE                      :: chain_ncell
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE           :: chain_noise
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE             :: chain_misfits
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE           :: chain_acceptance
INTEGER, DIMENSION(:,:), ALLOCATABLE                    :: chain_temperatures
CHARACTER (LEN=30)                                      :: last_model_file,last_tmp
INTEGER                                                 :: accepted_so_far,accepted_in_run,accepted_so_far_other
INTEGER, DIMENSION(:), ALLOCATABLE                      :: accepted_so_far_all
DOUBLE PRECISION, DIMENSION(3)                          :: acceptance_total
DOUBLE PRECISION, DIMENSION(3)                          :: acceptance_birth,acceptance_death,acceptance_move,acceptance_value
DOUBLE PRECISION, DIMENSION(3)                          :: acceptance_opt
DOUBLE PRECISION, DIMENSION(3)                          :: acceptance_move_main,acceptance_move_delayed
DOUBLE PRECISION, DIMENSION(3)                          :: acceptance_value_main,acceptance_value_delayed
DOUBLE PRECISION, DIMENSION(3)                          :: acceptance_opt_main,acceptance_opt_delayed
DOUBLE PRECISION, DIMENSION(3)                          :: acceptance_noise,acceptance_a,acceptance_b
DOUBLE PRECISION, DIMENSION(3)                          :: acceptance_tempering
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE             :: acceptance_all
! Pi
DOUBLE PRECISION, PARAMETER                             :: pi=DBLE(4)*ATAN(DBLE(1))
! Output files
INTEGER, SAVE                                           :: samples_unit=50
INTEGER, SAVE                                           :: ncell_unit=51
INTEGER, SAVE                                           :: noise_unit=52
INTEGER, SAVE                                           :: misfits_unit=53
INTEGER, SAVE                                           :: acceptance_unit=54
INTEGER, SAVE                                           :: temperatures_unit=55
INTEGER, SAVE                                           :: residuals_unit=56

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initializes the Markov chains by allocating all necessary arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE initialize_chains()
        IMPLICIT NONE

        IF (debug_mode) WRITE(debug_unit,*)'--- Chain ',TRIM(chain_name),': initializing chains ---'

        IF (opt_parameters.EQ.0) THEN 
                ALLOCATE(acceptance_all(13))
        ELSE
                ALLOCATE(acceptance_all(16))
        END IF

        ALLOCATE(voronoi_current(ncell_max,number_of_parameters),voronoi_proposed(ncell_max,number_of_parameters),&
                voronoi_delayed(ncell_max,number_of_parameters),voronoi_0(ncell_max,number_of_parameters))
        ALLOCATE(a_current(number_of_datasets),a_proposed(number_of_datasets))
        ALLOCATE(b_current(number_of_datasets),b_proposed(number_of_datasets))
        ALLOCATE(lambda_current(number_of_datasets),lambda_proposed(number_of_datasets))
        ALLOCATE(sigma_current(number_of_datasets),sigma_proposed(number_of_datasets))
        ALLOCATE(noise_current(number_of_measurements),noise_proposed(number_of_measurements),noise_delayed(number_of_measurements))

        ALLOCATE(data_current(number_of_measurements),data_proposed(number_of_measurements),data_delayed(number_of_measurements))
        IF (noise_type.EQ.1.OR.noise_type.EQ.2.OR.noise_type.EQ.5) THEN
                ALLOCATE(noise_x_current(number_of_measurements),&
                        noise_x_proposed(number_of_measurements),&
                        noise_x_delayed(number_of_measurements))
        END IF

        IF (run_number.EQ.1) THEN
                samples_in_run=number_of_samples+1
        ELSE IF (run_number.GT.1) THEN
                samples_in_run=number_of_samples
        END IF

        ALLOCATE(residuals(number_of_measurements,samples_in_run))

        ALLOCATE(chain_steps(samples_in_run,2))
        chain_steps=0

        ALLOCATE(chain_values(samples_in_run,number_of_parameters))
        chain_values=0

        ALLOCATE(chain_indices(samples_in_run))
        chain_indices=0

        ALLOCATE(chain_ncell(samples_in_run))
        chain_ncell=0

        IF (noise_type.EQ.1.OR.noise_type.EQ.2.OR.noise_type.EQ.5) THEN
                ALLOCATE(chain_noise(samples_in_run,2*number_of_datasets))
                chain_noise=0
        ELSE IF (noise_type.EQ.0.OR.noise_type.EQ.3.OR.noise_type.EQ.6) THEN
                ALLOCATE(chain_noise(samples_in_run,number_of_datasets))
                chain_noise=0
        END IF
        
        ALLOCATE(chain_misfits(samples_in_run))
        chain_misfits=0
        IF (opt_parameters.EQ.0) THEN 
                IF (parallel_tempering.EQ.0) THEN
                        ALLOCATE(chain_acceptance(samples_in_run,6))
                ELSE IF (parallel_tempering.EQ.1) THEN
                        ALLOCATE(chain_acceptance(samples_in_run,7))
                END IF
        ELSE
                IF (parallel_tempering.EQ.0) THEN
                        ALLOCATE(chain_acceptance(samples_in_run,7))
                ELSE IF (parallel_tempering.EQ.1) THEN
                        ALLOCATE(chain_acceptance(samples_in_run,8))
                END IF
        END IF
        chain_acceptance=0

        ALLOCATE(accepted_so_far_all(number_of_chains))

        ! Name of file containing last model of run
        last_model_file='last.tmp'
        last_tmp=TRIM(last_model_file) // TRIM(chain_name)

        IF (debug_mode) WRITE(debug_unit,*)'--- Chain ',TRIM(chain_name),': initialized chains ---'
END SUBROUTINE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Define Markov chain step type (1=cell birth, 2=cell death, 3=cell move, 4=value change, 5=noise change,
! 6=angle or aperture change)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE step_type()
        IMPLICIT NONE
        birth=.FALSE.
        death=.FALSE.
        move=.FALSE.
        value=.FALSE.
        noise=.FALSE.
        option=.FALSE.
        IF (step.EQ.1) birth=.TRUE.
        IF (step.EQ.2) death=.TRUE.
        IF (step.EQ.3) move=.TRUE.
        IF (step.EQ.4) value=.TRUE.
        IF (step.EQ.5) noise=.TRUE.
        IF (step.EQ.6) option=.TRUE.
END SUBROUTINE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Adds a Voronoi cell to the current model
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE add_cell(modelling_method)
        IMPLICIT NONE
        INTEGER          :: modelling_method
        REAL             :: x1,x2

        IF (ncell_current.EQ.ncell_max) THEN
                inside_prior=.FALSE.
                RETURN
        END IF

        ncell_proposed=ncell_current+1 ! number of cells increases by 1

        ! Create random voronoi cell or a fracture
        voronoi_proposed(ncell_proposed,1)=x_min+(x_max-x_min)*unirand() ! x centre of the fracture
        voronoi_proposed(ncell_proposed,2)=y_min+(y_max-y_min)*unirand() ! y centre of the fracture (if 3D) 
        voronoi_proposed(ncell_proposed,3)=z_min+(z_max-z_min)*unirand() ! z centre (depth) of the fracture
        
        ! Set the other parameters
        IF (modelling_method.EQ.1) THEN ! If fw_2_5d
                CALL find_closest_nucleus(voronoi_proposed(ncell_proposed,1:3),&
                                        voronoi_current(1:ncell_current,1:3),&
                                        ncell_current,&
                                        cell_number)
                voronoi_proposed(ncell_proposed,4)=voronoi_current(cell_number,4)&
                                                +birth_change*gausrand() 
        ELSEIF (modelling_method.EQ.2) THEN ! If ert_ddp_2_5d
                ! Length                            
                voronoi_proposed(ncell_proposed,4)=powerrand(value_min,value_max) 
                ! Option
                IF (opt_parameters.GT.0) THEN
                        voronoi_proposed(ncell_proposed,5) = opt_min + unirand()*(opt_max-opt_min)  
                END IF
                ! Check if the fracture is still in the domain
                CALL check_frac_inside_domain(voronoi_proposed(ncell_proposed,:))
        END IF

        CALL check_inside_prior()

        IF (inside_prior) THEN
                voronoi_proposed(ncell_proposed+1:ncell_max,:)=0
                IF (modelling_method.EQ.1) THEN ! If fw_2_5d
                        alpha_component=LOG(birth_change*SQRT(2*pi))&
                                +((voronoi_proposed(ncell_proposed,4)-voronoi_current(cell_number,4))**2/(2*birth_change**2))
                END IF
        END IF

        write(*,*)'ADD at chain ', TRIM(int2str(chain_number))
        write(*,*)'ADD voronoi_proposed = ', voronoi_proposed(ncell_proposed,:) 

END SUBROUTINE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Delete a Voronoi cell / fracture from the current model
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE delete_cell(modelling_method)
        IMPLICIT NONE
        INTEGER         :: closest,modelling_method

        IF (ncell_current.EQ.ncell_min) THEN
                inside_prior=.FALSE.
                RETURN
        END IF

        ncell_proposed=ncell_current-1

        ! Randomly choose a cell to delete
        cell_number=CEILING(unirand()*ncell_current) 
        IF (cell_number.EQ.0) cell_number=1

        ! Delete the cell
        voronoi_proposed(cell_number,:)=voronoi_current(ncell_current,:) ! Replace the deleted cell with the last cell
        voronoi_proposed(ncell_current:ncell_max,:)=0                    ! set everything to zero below the proposed number of cells

        ! Additional component to alpha calculation if fw_2_5d
        IF (modelling_method.EQ.1) THEN
                CALL find_closest_nucleus(voronoi_current(cell_number,1:3),&
                                        voronoi_proposed(1:ncell_proposed,1:3),&
                                        ncell_proposed,&
                                        closest)
                alpha_component=LOG(1/(birth_change*SQRT(2*pi)))&
                        +(-(voronoi_proposed(closest,4)-voronoi_current(cell_number,4))**2/(2*birth_change**2))     
        END IF

        write(*,*)'DELETE at chain ', TRIM(int2str(chain_number))

END SUBROUTINE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Moves a Voronoi cell / fracture from the current model to a new location
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE move_cell(modelling_method)
        IMPLICIT NONE
        INTEGER :: modelling_method
        REAL    :: x1,x2

        ! Randomly choose a cell to move
        cell_number=CEILING(unirand()*ncell_current) 
        IF (cell_number.EQ.0) cell_number=1

        ! Perturb the location of the chosen cell nucleus
        voronoi_proposed(cell_number,1)=voronoi_current(cell_number,1)&
                                        +gausrand()*(x_max-x_min)*x_change/100
        voronoi_proposed(cell_number,2)=voronoi_current(cell_number,2)&
                                        +gausrand()*(y_max-y_min)*y_change/100
        voronoi_proposed(cell_number,3)=voronoi_current(cell_number,3)&
                                        +gausrand()*(z_max-z_min)*z_change/100

        ! If ert_ddp_2_5d: check if the fracture is still in the domain
        IF (modelling_method.EQ.2) THEN
                CALL check_frac_inside_domain(voronoi_proposed(cell_number,:))  
        END IF

        write(*,*) 'MOVE at chain ', TRIM(int2str(chain_number))
        write(*,*) 'MOVE voronoi_proposed = ', voronoi_proposed(cell_number,:) 

        CALL check_inside_prior()
END SUBROUTINE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Change the resistivity value (voronoi cell) / length value (fracture)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE change_value(modelling_method)
        IMPLICIT NONE
        INTEGER :: modelling_method
        REAL    :: x1,x2

        cell_number=CEILING(unirand()*ncell_current) ! randomly choose a cell to perturb
        IF (cell_number.EQ.0) cell_number=1

        IF (modelling_method.EQ.1) THEN ! If fw_2_5d
                voronoi_proposed(cell_number,4)=voronoi_current(cell_number,4)&
                                                +gausrand()*value_change
        ELSE IF (modelling_method.EQ.2) THEN ! If ert_ddp_2_5d
                l_change=(gausrand()*value_change)*powerrand_std(value_min,value_max)
                voronoi_proposed(cell_number,4)=voronoi_current(cell_number,4)+l_change
                ! Check if the fracture is still in the domain
                CALL check_frac_inside_domain(voronoi_proposed(cell_number,:))
        END IF

        write(*,*) 'VALUE at chain ', TRIM(int2str(chain_number))
        write(*,*) 'VALUE voronoi_proposed = ', voronoi_proposed(cell_number,:) 

        CALL check_inside_prior()
END SUBROUTINE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Change data noise parameters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE change_noise()
        IMPLICIT NONE
        DOUBLE PRECISION        :: random_deviate
        
        IF (number_of_datasets.GT.1) THEN
                dataset_number=CEILING(unirand()*number_of_datasets)
                IF (dataset_number.EQ.0) dataset_number=1
                IF (dataset_number.GT.number_of_datasets) dataset_number=number_of_datasets
        ELSE
                dataset_number=1
        END IF

        IF (noise_type.EQ.0) THEN
                change_sigma=.TRUE.
                sigma_proposed(dataset_number)=sigma_current(dataset_number)&
                                        +gausrand()*sigma_change
        ELSE IF (noise_type.EQ.1.OR.noise_type.EQ.2.OR.noise_type.EQ.5) THEN
                random_deviate=unirand()
                IF (random_deviate.LT.0.5) THEN
                        change_a=.TRUE.
                        a_proposed(dataset_number)=a_current(dataset_number)&
                                                +gausrand()*a_change
                ELSE
                        change_b=.TRUE.
                        b_proposed(dataset_number)=b_current(dataset_number)&
                                                +gausrand()*b_change
                END IF
        ELSE IF (noise_type.EQ.3.OR.noise_type.EQ.6) THEN
                change_lambda=.TRUE.
                lambda_proposed(dataset_number)=lambda_current(dataset_number)&
                                                +gausrand()*lambda_change
        END IF

        write(*,*) 'NOISE at chain ', TRIM(int2str(chain_number))

        CALL check_inside_prior()
END SUBROUTINE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Option: change the angle or aperture of a fracture, only if ert_ddp_2_5d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE change_option()
        IMPLICIT NONE
        
        ! Randomly choose a cell to be perturbed
        cell_number=CEILING(unirand()*ncell_current) 
        IF (cell_number.EQ.0) cell_number=1

        ! Change option
        voronoi_proposed(cell_number,5)=voronoi_current(cell_number,5)+gausrand()*(opt_max-opt_min)*opt_change/100.

        ! Check is the fracture is inside the domain
        CALL check_frac_inside_domain(voronoi_proposed(cell_number,:))

        write(*,*)'OPTION ', TRIM(int2str(opt_parameters)), ' at chain ', TRIM(int2str(chain_number))
        write(*,*)'OPTION voronoi_proposed = ', voronoi_proposed(cell_number,:)

        CALL check_inside_prior()
END SUBROUTINE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DELAYED: move a Voronoi cell / fracture in the current model to a new location
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE move_cell_delayed(modelling_method)
        IMPLICIT NONE
        INTEGER :: modelling_method

        delayed_position=.TRUE.

        ! Perturb the location of the chosen cell nucleus
        voronoi_delayed(cell_number,1)=voronoi_current(cell_number,1)&
                                +gausrand()*(x_max-x_min)*x_change_delayed/100
        voronoi_delayed(cell_number,2)=voronoi_current(cell_number,2)&
                                +gausrand()*(y_max-y_min)*y_change_delayed/100
        voronoi_delayed(cell_number,3)=voronoi_current(cell_number,3)&
                                +gausrand()*(z_max-z_min)*z_change_delayed/100

        ! If ert_ddp_2_5d: check if the fracture is still in the domain
        IF (modelling_method.EQ.2) THEN 
                CALL check_frac_inside_domain(voronoi_delayed(cell_number,:))
        END IF

        write(*,*)'MOVE DELAYED at chain ', TRIM(int2str(chain_number))
        write(*,*)'MOVE DELAYED voronoi_proposed is ', voronoi_delayed(cell_number,:) 

        CALL check_inside_prior()
END SUBROUTINE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DELAYED: electrical resistivity value (voronoi cell) or length (fracture)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE change_value_delayed(modelling_method)
        IMPLICIT NONE
        INTEGER :: modelling_method

        delayed_value=.TRUE.

        IF (modelling_method.EQ.1) THEN ! If fw_2_5d
                voronoi_delayed(cell_number,4)=voronoi_current(cell_number,4)&
                                                +gausrand()*value_change_delayed
        ELSE IF (modelling_method.EQ.2) THEN ! If ert_ddp_2_5d
                ! Change length
                l_change=(gausrand()*value_change_delayed)*powerrand_std(value_min,value_max)
                voronoi_delayed(cell_number,4)=voronoi_delayed(cell_number,4)+l_change
                ! Check if the fracture is still in the domain
                CALL check_frac_inside_domain(voronoi_delayed(cell_number,:))
        END IF

        write(*,*)'VALUE DELAYED at chain ', TRIM(int2str(chain_number))
        write(*,*)'VALUE DELAYED voronoi_proposed is ', voronoi_delayed(cell_number,:) 

        CALL check_inside_prior()
END SUBROUTINE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DELAYED: change the angle or aperture of a fracture, only if ert_ddp_2_5d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE change_option_delayed()
        IMPLICIT NONE

        delayed_opt=.TRUE.
        
        ! Randomly choose a cell to perturb
        cell_number=CEILING(unirand()*ncell_current) 
        IF (cell_number.EQ.0) cell_number=1
        
        ! Change option value
        voronoi_delayed(cell_number,5)=voronoi_current(cell_number,5)+gausrand()*(opt_max-opt_min)*opt_change_delayed/100.

        ! Check is the fracture is inside the domain
        CALL check_frac_inside_domain(voronoi_delayed(cell_number,:))

        write(*,*)'OPTION DELAYED ', TRIM(int2str(opt_parameters)), ' at chain ', TRIM(int2str(chain_number))
        write(*,*)'OPTION DELAYED voronoi_delayed = ', voronoi_delayed(cell_number,1:4), voronoi_delayed(cell_number,5)*180./pi 

        CALL check_inside_prior()       
END SUBROUTINE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Check if the model parameters respect the prior limits
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE check_inside_prior()
        IMPLICIT NONE

        IF (birth) THEN
                IF (voronoi_proposed(ncell_proposed,1).LT.x_min.OR.&
                voronoi_proposed(ncell_proposed,1).GT.x_max) inside_prior=.FALSE.
                IF (voronoi_proposed(ncell_proposed,2).LT.y_min.OR.&
                voronoi_proposed(ncell_proposed,2).GT.y_max) inside_prior=.FALSE.
                IF (voronoi_proposed(ncell_proposed,3).LT.z_min.OR.&
                voronoi_proposed(ncell_proposed,3).GT.z_max) inside_prior=.FALSE.
                IF (voronoi_proposed(ncell_proposed,4).LT.value_min.OR.&
                voronoi_proposed(ncell_proposed,4).GT.value_max) inside_prior=.FALSE.
        ELSE IF (move) THEN
                IF (.NOT.delayed_position) THEN
                        IF (voronoi_proposed(cell_number,1).LT.x_min.OR.&
                        voronoi_proposed(cell_number,1).GT.x_max) inside_prior=.FALSE.
                        IF (voronoi_proposed(cell_number,2).LT.y_min.OR.&
                        voronoi_proposed(cell_number,2).GT.y_max) inside_prior=.FALSE.
                        IF (voronoi_proposed(cell_number,3).LT.z_min.OR.&
                        voronoi_proposed(cell_number,3).GT.z_max) inside_prior=.FALSE.
                ELSE IF (delayed_position) THEN
                        IF (voronoi_delayed(cell_number,1).LT.x_min.OR.&
                        voronoi_delayed(cell_number,1).GT.x_max) inside_prior_delayed=.FALSE.
                        IF (voronoi_delayed(cell_number,2).LT.y_min.OR.&
                        voronoi_delayed(cell_number,2).GT.y_max) inside_prior_delayed=.FALSE.
                        IF (voronoi_delayed(cell_number,3).LT.z_min.OR.&
                        voronoi_delayed(cell_number,3).GT.z_max) inside_prior_delayed=.FALSE.
                END IF
        ELSE IF (value) THEN
                IF (.NOT.delayed_value) THEN
                        IF (voronoi_proposed(cell_number,4).LT.value_min.OR.&
                        voronoi_proposed(cell_number,4).GT.value_max) inside_prior=.FALSE.
                ELSE IF (delayed_value) THEN
                        IF (voronoi_delayed(cell_number,4).LT.value_min.OR.&
                        voronoi_delayed(cell_number,4).GT.value_max) inside_prior_delayed=.FALSE.
                END IF
        ELSE IF (noise) THEN
                IF (change_sigma) THEN
                        IF (sigma_proposed(dataset_number).LT.sigma_min.OR.&
                        sigma_proposed(dataset_number).GT.sigma_max) inside_prior=.FALSE.
                ELSE IF (change_a) THEN
                        IF (a_proposed(dataset_number).LT.a_min.OR.&
                        a_proposed(dataset_number).GT.a_max) inside_prior=.FALSE.
                ELSE IF (change_b) THEN
                        IF (b_proposed(dataset_number).LT.b_min.OR.&
                        b_proposed(dataset_number).GT.b_max) inside_prior=.FALSE.
                ELSE IF (change_lambda) THEN
                        IF (lambda_proposed(dataset_number).LT.lambda_min.OR.&
                        lambda_proposed(dataset_number).GT.lambda_max) inside_prior=.FALSE.
                END IF
        ELSE IF (option) THEN
                IF (.NOT.delayed_opt) THEN
                        IF (voronoi_proposed(cell_number,5).LT.opt_min.OR.&
                        voronoi_proposed(cell_number,5).GT.opt_max) inside_prior=.FALSE.
                ELSE IF (delayed_opt) THEN
                        IF (voronoi_delayed(cell_number,5).LT.opt_min.OR.&
                        voronoi_delayed(cell_number,5).GT.opt_max) inside_prior_delayed=.FALSE.
                END IF
        END IF
END SUBROUTINE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Re-initialize the models
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE reset_proposed_models()
        IMPLICIT NONE
        inside_prior=.TRUE.
        inside_prior_delayed=.TRUE.
        ncell_proposed=ncell_current
        voronoi_proposed=voronoi_current
        voronoi_delayed=voronoi_current
        a_proposed=a_current
        b_proposed=b_current
        lambda_proposed=lambda_current
        sigma_proposed=sigma_current
        noise_proposed=noise_current
        noise_delayed=noise_current
        delayed_position=.FALSE.
        delayed_value=.FALSE.
        delayed_opt=.FALSE.
        change_sigma=.FALSE.
        change_a=.FALSE.
        change_b=.FALSE.
        change_lambda=.FALSE.
        alpha_component=0
        cell_number=0
        data_proposed=data_current
        data_delayed=data_current
        likelihood_proposed=likelihood_current
        likelihood_delayed=likelihood_current
        misfit_proposed=misfit_current
        misfit_delayed=misfit_current
        IF (noise_type.EQ.1.OR.noise_type.EQ.2.OR.noise_type.EQ.5) THEN
                noise_x_proposed=noise_x_current
                noise_x_delayed=noise_x_current
        END IF
END SUBROUTINE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutine to check if the proposed fracture is inside the boundaries of the domain
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE check_frac_inside_domain(fracture)
        IMPLICIT NONE
        DOUBLE PRECISION, DIMENSION(:) :: fracture 
        DOUBLE PRECISION               :: length,x1,x2,z1,z2,angle
        INTEGER(KIND=4), PARAMETER     :: dim_num = 2
        INTEGER(KIND=4)                :: ival
        REAL(KIND=8)                   :: p1(dim_num)
        REAL(KIND=8)                   :: p2(dim_num)
        REAL(KIND=8)                   :: pa(dim_num),pa_old(dim_num)
        REAL(KIND=8)                   :: pb(dim_num),pb_old(dim_num)
        
        ! Define the domain
        p1 = (/ x_min,z_min /)
        p2 = (/ x_max,z_max /)

        ! Define fracture bounds
        IF (opt_parameters.EQ.1) THEN 
                angle = fracture(5)
        ELSE
                angle = 0
        END IF
        length = fracture(4)
        x1 = fracture(1) - COS(angle)*length/2
        z1 = fracture(3) - SIN(angle)*length/2
        x2 = fracture(1) + COS(angle)*length/2
        z2 = fracture(3) + SIN(angle)*length/2

        ! Check if the fracture is inside the domain
        pa = (/ x1,z1 /) ; pa_old = pa
        pb = (/ x2,z2 /) ; pb_old = pb

        CALL box_segment_clip_2d(p1,p2,pa,pb,ival)

        !WRITE ( *, '(2x,4f8.4)' ) pa_old(1:dim_num), pb_old(1:dim_num)
        IF ( ival == -1 ) THEN
                WRITE ( *, '(a)' ) '  Line is outside the box.'
        !ELSE IF ( ival == 0 ) THEN
        !        WRITE ( *, '(a)' ) '  Line is inside the box.'
        !ELSE IF ( ival == 1 ) THEN
        !        WRITE ( *, '(a)' ) '  Line is partly outside the box.'
                !WRITE ( *, '(2x,2f8.4)' ) pa(1:dim_num)
        !ELSE IF ( ival == 2 ) THEN
        !        WRITE ( *, '(a)' ) '  Line is partly outside the box.'
                !WRITE ( *, '(2x,16x,2f8.4)' ) pb(1:dim_num)
        !ELSE IF ( ival == 3 ) THEN
        !        WRITE ( *, '(a)' ) '  Line is partly outside the box.'
                !WRITE ( *, '(2x,4f8.4)' ) pa(1:dim_num), pb(1:dim_num)
        END IF

        IF (ival.EQ.1.OR.ival.EQ.2.OR.ival.EQ.3) THEN
                ! Update fracture
                IF (pa(1)<pb(1)) THEN 
                        fracture(1) = pa(1) + abs( (pb(1)-pa(1))/2 ) !x
                ELSE
                        fracture(1) = pa(1) - abs( (pb(1)-pa(1))/2 ) !x
                END IF
                IF (pa(2)<pb(2)) THEN
                        fracture(3) = pa(2) + abs( (pb(2)-pa(2))/2 )!z
                ELSE
                        fracture(3) = pa(2) - abs( (pb(2)-pa(2))/2 )!z
                END IF
                fracture(4) = (pb(1)-pa(1))/COS(angle) !l
                if (fracture(4).LT.0.) THEN
                        WRITE(*,*)"WARNING: negative length"
                END IF
        END IF
END SUBROUTINE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutine to clip a line according to a rectangular domain
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE box_segment_clip_2d ( p1, p2, pa, pb, ival )

!*****************************************************************************80
!! BOX_SEGMENT_CLIP_2D uses a box to clip a line segment in 2D.
!
!  Discussion:
!
!    A box in 2D is a rectangle with sides aligned on coordinate
!    axes.  It can be described by its low and high corners, P1 and P2
!    as the set of points P satisfying:
!
!      P1(1:2) <= P(1:2) <= P2(1:2).
!
!    A line segment is the finite portion of a line that lies between
!    two points.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    21 March 2011
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) P1(2), P2(2), the low and high corners of the box.
!
!    Input/output, real ( kind = 8 ) PA(2), PB(2); on input, the endpoints 
!    of a line segment.  On output, the endpoints of the portion of the
!    line segment that lies inside the box.  However, if no part of the
!    initial line segment lies inside the box, the output value is the
!    same as the input value.
!
!    Output, integer ( kind = 4 ) IVAL:
!    -1, no part of the line segment is within the box.
!     0, no clipping was necessary.
!     1, PA was clipped.
!     2, PB was clipped.
!     3, PA and PB were clipped.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 2

  logical ( kind = 4 ) clip_a
  logical ( kind = 4 ) clip_b
  integer ( kind = 4 ) ival
  real ( kind = 8 ) p1(dim_num)
  real ( kind = 8 ) p2(dim_num)
  real ( kind = 8 ) pa(dim_num)
  real ( kind = 8 ) pb(dim_num)
  real ( kind = 8 ) q(dim_num)

  clip_a = .false.
  clip_b = .false.
!
!  Require that XMIN <= X.
!
  if ( pa(1) < p1(1) .and. pb(1) < p1(1) ) then
    ival = -1
    return
  end if

  if ( pa(1) < p1(1) .and. p1(1) <= pb(1) ) then
    q(1) = p1(1)
    q(2) = pa(2) + ( pb(2) - pa(2) ) * ( q(1) - pa(1) ) / ( pb(1) - pa(1) )
    pa(1:2) = q(1:2)
    clip_a = .true.
  else if ( p1(1) <= pa(1) .and. pb(1) < p1(1) ) then
    q(1) = p1(1)
    q(2) = pa(2) + ( pb(2) - pa(2) ) * ( q(1) - pa(1) ) / ( pb(1) - pa(1) )
    pb(1:2) = q(1:2)
    clip_b = .true.
  end if
!
!  Require that X <= XMAX.
!
  if ( p2(1) < pa(1) .and. p2(1) < pb(1) ) then
    ival = -1
    return
  end if

  if ( p2(1) < pa(1) .and. pb(1) <= p2(1) ) then
    q(1) = p2(1)
    q(2) = pa(2) + ( pb(2) - pa(2) ) * ( q(1) - pa(1) ) / ( pb(1) - pa(1) )
    pa(1:2) = q(1:2)
    clip_a = .true.
  else if ( pa(1) <= p2(1) .and. p2(1) < pb(1) ) then
    q(1) = p2(1)
    q(2) = pa(2) + ( pb(2) - pa(2) ) * ( q(1) - pa(1) ) / ( pb(1) - pa(1) )
    pb(1:2) = q(1:2)
    clip_b = .true.
  end if
!
!  Require that YMIN <= Y.
!
  if ( pa(2) < p1(2) .and. pb(2) < p1(2) ) then
    ival = -1
    return
  end if

  if ( pa(2) < p1(2) .and. p1(2) <= pb(2) ) then
    q(2) = p1(2)
    q(1) = pa(1) + ( pb(1) - pa(1) ) * ( q(2) - pa(2) ) / ( pb(2) - pa(2) )
    pa(1:2) = q(1:2)
    clip_a = .true.
  else if ( p1(2) <= pa(2) .and. pb(2) < p1(2) ) then
    q(2) = p1(2)
    q(1) = pa(1) + ( pb(1) - pa(1) ) * ( q(2) - pa(2) ) / ( pb(2) - pa(2) )
    pb(1:2) = q(1:2)
    clip_b = .true.
  end if
!
!  Require that Y <= YMAX.
!
  if ( p2(2) < pa(2) .and. p2(2) < pb(2) ) then
    ival = -1
    return
  end if

  if ( p2(2) < pa(2) .and. pb(2) <= p2(2) ) then
    q(2) = p2(2)
    q(1) = pa(1) + ( pb(1) - pa(1) ) * ( q(2) - pa(2) ) / ( pb(2) - pa(2) )
    pa(1:2) = q(1:2)
    clip_a = .true.
  else if ( pa(2) <= p2(2) .and. p2(2) < pb(2) ) then
    q(2) = p2(2)
    q(1) = pa(1) + ( pb(1) - pa(1) ) * ( p2(2) - pa(2) ) / ( pb(2) - pa(2) )
    pb(1:2) = q(1:2)
    clip_b = .true.
  end if

  ival = 0

  if ( clip_a ) then
    ival = ival + 1
  end if

  if ( clip_b ) then
    ival = ival + 2
  end if

  return
END SUBROUTINE

END MODULE
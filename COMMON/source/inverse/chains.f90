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
DOUBLE PRECISION                                        :: l_change ! needed when changing length of a fracture
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
DOUBLE PRECISION, DIMENSION(3)                          :: acceptance_move_main,acceptance_move_delayed
DOUBLE PRECISION, DIMENSION(3)                          :: acceptance_value_main,acceptance_value_delayed
DOUBLE PRECISION, DIMENSION(3)                          :: acceptance_opt_main,acceptance_opt_delayed
DOUBLE PRECISION, DIMENSION(3)                          :: acceptance_noise,acceptance_a,acceptance_b
DOUBLE PRECISION, DIMENSION(3)                          :: acceptance_tempering
DOUBLE PRECISION, DIMENSION(3)                          :: acceptance_opt
!DOUBLE PRECISION, DIMENSION(13)                         :: acceptance_all
DOUBLE PRECISION, DIMENSION(16)                         :: acceptance_all
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
SUBROUTINE initialize_chains()
! Initializes the Markov chains by allocating all necessary arrays
IMPLICIT NONE

IF (debug_mode) WRITE(debug_unit,*)'--- Chain ',TRIM(chain_name),': initializing chains ---'

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
SUBROUTINE step_type()
! Define Markov chain step type (1=cell birth, 2=cell death, 3=cell move, 4=value change, 5=noise change)
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
SUBROUTINE add_cell(modelling_method)
! Adds a Voronoi cell to the current model
IMPLICIT NONE
INTEGER          :: modelling_method
REAL             :: x1,x2
DOUBLE PRECISION :: mean_currentlength

IF (ncell_current.EQ.ncell_max) THEN
        inside_prior=.FALSE.
        RETURN
END IF

ncell_proposed=ncell_current+1 ! number of cells increases by 1

! Create random voronoi_current(cell_location for the new cell nucleus / new fracture centre)
voronoi_proposed(ncell_proposed,1)=x_min+(x_max-x_min)*unirand() ! x centre of the fracture
voronoi_proposed(ncell_proposed,2)=y_min+(y_max-y_min)*unirand() ! y centre of the fracture (if 3D) 
voronoi_proposed(ncell_proposed,3)=z_min+(z_max-z_min)*unirand() ! z centre (depth) of the fracture


!---------- If fw_2.5d ----------
IF (modelling_method.EQ.1) THEN
        ! Create velocity for the new cell by finding the velocity at the location of the new
        ! nucleus in the previous tessellation and perturbing it
        CALL find_closest_nucleus(voronoi_proposed(ncell_proposed,1:3),&
                                voronoi_current(1:ncell_current,1:3),&
                                ncell_current,&
                                cell_number)

        voronoi_proposed(ncell_proposed,4)=voronoi_current(cell_number,4)&
                                        +birth_change*gausrand() 

!---------- If ert_ddp_2_5d ----------
ELSEIF (modelling_method.EQ.2) THEN
        !voronoi_proposed(ncell_proposed,4)=voronoi_current(cell_number,4)&
        !                                +birth_change*powerrand(value_min)
                                        
        voronoi_proposed(ncell_proposed,4)=powerrand(value_min,value_max)
        
        ! Check if the fracture edge values are in the domain boundaries
        write(*,*) 'Add frac at chain ', TRIM(int2str(chain_number)), ': l (bf check) = ', voronoi_proposed(ncell_proposed,4) 
        x1 = voronoi_proposed(ncell_proposed,1) - voronoi_proposed(ncell_proposed,4)/2
        x2 = voronoi_proposed(ncell_proposed,1) + voronoi_proposed(ncell_proposed,4)/2
        IF (x1.LT.x_min) THEN 
                voronoi_proposed(ncell_proposed,4) = abs(x2-x_min)                                 ! recalculate the length
                voronoi_proposed(ncell_proposed,1) = x_min + voronoi_proposed(ncell_proposed,4)/2  ! recalculate the centre
        ELSE IF (x2.GT.x_max) THEN
                voronoi_proposed(ncell_proposed,4) = abs(x_max-x1)                                 ! recalculate the length
                voronoi_proposed(ncell_proposed,1) = x_max - voronoi_proposed(ncell_proposed,4)/2  ! recalculate the centre
        END IF  
        write(*,*) 'Add frac at chain ', TRIM(int2str(chain_number)), ': l (after check) = ', voronoi_proposed(ncell_proposed,4) 

        IF (opt_parameters.GT.0) THEN 
                voronoi_proposed(ncell_proposed,5) = opt_min+(opt_max-opt_min)*unirand()
        END IF
        write(*,*) 'Add frac at chain ', TRIM(int2str(chain_number)), ': angle/aperture = ', voronoi_proposed(ncell_proposed,5) 

END IF

CALL check_inside_prior()

IF (inside_prior) THEN
        voronoi_proposed(ncell_proposed+1:ncell_max,:)=0 ! set everything to zero below the proposed number of cells
        ! Calculate term needed to evaluate the acceptance ratio later

        !---------- If fw_2.5d ----------
        IF (modelling_method.EQ.1) THEN
                alpha_component=LOG(birth_change*SQRT(2*pi))&
                        +((voronoi_proposed(ncell_proposed,4)-voronoi_current(cell_number,4))**2/(2*birth_change**2))
        
        !---------- If ert_ddp_2_5d ----------
        !ELSEIF (modelling_method.EQ.2) THEN

                !alpha_component=1/powerrand_k(value_min,value_max)
                !write(*,*) 'alpha component (add) at chain ', TRIM(int2str(chain_number)), '= ', alpha_component

        END IF

END IF

END SUBROUTINE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE delete_cell(modelling_method)
! Delete a Voronoi cell from the current model
IMPLICIT NONE
INTEGER         :: closest,modelling_method

IF (ncell_current.EQ.ncell_min) THEN
        inside_prior=.FALSE.
        RETURN
END IF
write(*,*) 'delete a cell at chain ', TRIM(int2str(chain_number))

ncell_proposed=ncell_current-1 ! number of cells decreases by 1

cell_number=CEILING(unirand()*ncell_current) ! randomly choose a cell to delete
IF (cell_number.EQ.0) cell_number=1

!---------- If fw_2.5d ----------
IF (modelling_method.EQ.1) THEN
        ! Replace the deleted cell with the last cell
        voronoi_proposed(cell_number,:)=voronoi_current(ncell_current,:)
        voronoi_proposed(ncell_current:ncell_max,:)=0 ! set everything to zero below the proposed number of cells

        ! Find the closest nucleus in the new tessellation to the one that was deleted
        CALL find_closest_nucleus(voronoi_current(cell_number,1:3),&
                                voronoi_proposed(1:ncell_proposed,1:3),&
                                ncell_proposed,&
                                closest)
        ! Calculate term needed to evaluate the acceptance ratio later
        alpha_component=LOG(1/(birth_change*SQRT(2*pi)))&
                +(-(voronoi_proposed(closest,4)-voronoi_current(cell_number,4))**2/(2*birth_change**2))

!---------- If ert_ddp_2_5d ----------
ELSEIF (modelling_method.EQ.2) THEN
        
        IF (cell_number.EQ.1) THEN
                voronoi_proposed(1:ncell_proposed,:)=voronoi_current(2:ncell_current,:)
        ELSE IF (cell_number.EQ.ncell_current) THEN
                voronoi_proposed(1:ncell_proposed,:)=voronoi_current(1:ncell_proposed,:)
        ELSE
                voronoi_proposed(1:cell_number-1,:)=voronoi_current(1:cell_number-1,:)
                voronoi_proposed(cell_number:ncell_proposed,:)=voronoi_current(cell_number+1:ncell_current,:)
        END IF
        
        ! Calculate term needed to evaluate the acceptance ratio later
        !alpha_component= powerrand_k(value_min,value_max)
        !write(*,*) 'alpha component (delete) at chain ', TRIM(int2str(chain_number)), '= ', alpha_component
END IF

END SUBROUTINE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE move_cell(modelling_method)
! Moves a Voronoi cell / fracture in the current model to a new location
IMPLICIT NONE
INTEGER :: modelling_method
REAL :: x1,x2
write(*,*) 'move a cell at chain ', TRIM(int2str(chain_number))

cell_number=CEILING(unirand()*ncell_current) ! randomly choose a cell to move

IF (cell_number.EQ.0) cell_number=1

! Perturb the location of the chosen cell nucleus
voronoi_proposed(cell_number,1)=voronoi_current(cell_number,1)&
                                +gausrand()*(x_max-x_min)*x_change/100
voronoi_proposed(cell_number,2)=voronoi_current(cell_number,2)&
                                +gausrand()*(y_max-y_min)*y_change/100
voronoi_proposed(cell_number,3)=voronoi_current(cell_number,3)&
                                +gausrand()*(z_max-z_min)*z_change/100

IF (modelling_method.EQ.2) THEN
        ! Check if the fracture edge values are in the domain boundaries
        write(*,*) 'Move frac at chain ', TRIM(int2str(chain_number)), ': l (bf check) = ', voronoi_proposed(ncell_proposed,4) 
        x1 = voronoi_proposed(cell_number,1) - voronoi_proposed(cell_number,4)/2
        x2 = voronoi_proposed(cell_number,1) + voronoi_proposed(cell_number,4)/2
        IF (x1.LT.x_min) THEN 
                voronoi_proposed(cell_number,4) = abs(x2-x_min)                              ! recalculate the length
                voronoi_proposed(cell_number,1) = x_min + voronoi_proposed(cell_number,4)/2  ! recalculate the centre
        ELSE IF (x2.GT.x_max) THEN
                voronoi_proposed(cell_number,4) = abs(x_max-x1)                              ! recalculate the length
                voronoi_proposed(cell_number,1) = x_max - voronoi_proposed(cell_number,4)/2  ! recalculate the centre
        END IF   
        write(*,*) 'Move frac at chain ', TRIM(int2str(chain_number)), ': l (after check) = ', voronoi_proposed(ncell_proposed,4) 
END IF

CALL check_inside_prior()

END SUBROUTINE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE change_value(modelling_method)

IMPLICIT NONE
INTEGER :: modelling_method
REAL    :: x1,x2
write(*,*) 'change length at chain ', TRIM(int2str(chain_number))

cell_number=CEILING(unirand()*ncell_current) ! randomly choose a cell to perturb
IF (cell_number.EQ.0) cell_number=1

!---------- If fw_2.5d ----------
IF (modelling_method.EQ.1) THEN
        ! Perturb the value of the chosen cell
        voronoi_proposed(cell_number,4)=voronoi_current(cell_number,4)&
                                        +gausrand()*value_change
!---------- If ert_ddp_2_5d ----------
ELSE IF (modelling_method.EQ.2) THEN
        !voronoi_proposed(cell_number,4)=voronoi_current(cell_number,4)&
        !                                +powerrand(value_min)*value_change
        !l_change=gausrand()*powerrand_std(value_min,value_max)
        ! proposotion CC 10/04/2025
        l_change=(gausrand()*value_change)*powerrand_std(value_min,value_max)
        ! ou
        ! l_change=(unirand()-0.5)*2.*powerrand_std(value_min,value_max)*value_change
        ! WARNING UNIT VALUE CHANGE: %
        voronoi_proposed(cell_number,4)=voronoi_current(cell_number,4)+l_change
        !IF (voronoi_current(cell_number,4).LE.l_change) THEN
        !        voronoi_proposed(cell_number,4)=voronoi_current(cell_number,4)+l_change
        !ELSE
        !        IF (unirand().LE.0.5) THEN 
        !                voronoi_proposed(cell_number,4)=voronoi_current(cell_number,4)-l_change 
        !        ELSE
        !                voronoi_proposed(cell_number,4)=voronoi_current(cell_number,4)+l_change
        !        END IF
        !END IF
        !voronoi_proposed(cell_number,4)=voronoi_current(cell_number,4)&
        !                                +unirand()*powerrand_std(value_min,value_max)

        ! Check if the fracture edge values are in the domain boundaries
        x1 = voronoi_proposed(cell_number,1) - voronoi_proposed(cell_number,4)/2
        x2 = voronoi_proposed(cell_number,1) + voronoi_proposed(cell_number,4)/2
        write(*,*) 'Change value at chain ', TRIM(int2str(chain_number)), ': l (bf check) = ', voronoi_proposed(ncell_proposed,4) 
        IF (x1.LT.x_min) THEN 
                voronoi_proposed(cell_number,4) = abs(x2-x_min)                              ! recalculate the length
                voronoi_proposed(cell_number,1) = x_min + voronoi_proposed(cell_number,4)/2  ! recalculate the centre
        ELSE IF (x2.GT.x_max) THEN
                voronoi_proposed(cell_number,4) = abs(x_max-x1)                              ! recalculate the length
                voronoi_proposed(cell_number,1) = x_max - voronoi_proposed(cell_number,4)/2  ! recalculate the centre
        END IF    
        write(*,*) 'Change value at chain ', TRIM(int2str(chain_number)), ': l (after check) = ', voronoi_proposed(ncell_proposed,4) 
END IF

CALL check_inside_prior()

END SUBROUTINE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE change_noise()
! Changes data noise parameters
IMPLICIT NONE
DOUBLE PRECISION        :: random_deviate
write(*,*) 'change noise at chain ', TRIM(int2str(chain_number))

IF (number_of_datasets.GT.1) THEN
        dataset_number=CEILING(unirand()*number_of_datasets) ! randomly choose a dataset to perturb
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

CALL check_inside_prior()

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Option not imlemented for fw2_5d
SUBROUTINE change_option()
        ! Changes data option parameters
        IMPLICIT NONE
        DOUBLE PRECISION        :: random_deviate

        write(*,*) 'change optional parameter ', TRIM(int2str(opt_parameters)), ' at chain ', TRIM(int2str(chain_number))
        
        cell_number=CEILING(unirand()*ncell_current) ! randomly choose a cell to perturb
        IF (cell_number.EQ.0) cell_number=1

        voronoi_proposed(cell_number,5)=voronoi_current(cell_number,5)+gausrand()*opt_change
        write(*,*)'angle change_option',voronoi_proposed(cell_number,5)

        CALL check_inside_prior()
        
END SUBROUTINE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE move_cell_delayed(modelling_method)
! Moves a Voronoi cell in the current model to a new location
IMPLICIT NONE
INTEGER :: modelling_method
REAL :: x1,x2
write(*,*) 'move cell delayed at chain ', TRIM(int2str(chain_number))

delayed_position=.TRUE.

! Perturb the location of the chosen cell nucleus
voronoi_delayed(cell_number,1)=voronoi_current(cell_number,1)&
                               +gausrand()*(x_max-x_min)*x_change_delayed/100
voronoi_delayed(cell_number,2)=voronoi_current(cell_number,2)&
                               +gausrand()*(y_max-y_min)*y_change_delayed/100
voronoi_delayed(cell_number,3)=voronoi_current(cell_number,3)&
                               +gausrand()*(z_max-z_min)*z_change_delayed/100

IF (modelling_method.EQ.2) THEN
        ! Check if the fracture edge values are in the domain boundaries
        x1 = voronoi_delayed(cell_number,1) - voronoi_delayed(cell_number,4)/2
        x2 = voronoi_delayed(cell_number,1) + voronoi_delayed(cell_number,4)/2
        IF (x1.LT.x_min) THEN 
                voronoi_delayed(cell_number,4) = abs(x2-x_min)                             ! recalculate the length
                voronoi_delayed(cell_number,1) = x_min + voronoi_delayed(cell_number,4)/2  ! recalculate the centre
        ELSE IF (x2.GT.x_max) THEN
                voronoi_delayed(cell_number,4) = abs(x_max-x1)                             ! recalculate the length
                voronoi_delayed(cell_number,1) = x_max - voronoi_delayed(cell_number,4)/2  ! recalculate the centre
        END IF   
END IF

CALL check_inside_prior()

END SUBROUTINE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE change_value_delayed(modelling_method)

IMPLICIT NONE
INTEGER :: modelling_method
REAL    :: x1,x2
write(*,*) 'change value delayed at chain ', TRIM(int2str(chain_number))

delayed_value=.TRUE.

!---------- If fw_2.5d ----------
IF (modelling_method.EQ.1) THEN
        ! Perturb the value of the chosen cell
        voronoi_delayed(cell_number,4)=voronoi_current(cell_number,4)&
                                        +gausrand()*value_change_delayed
!---------- If ert_ddp_2_5d ----------
ELSE IF (modelling_method.EQ.2) THEN

        !l_change=gausrand()*powerrand_std(value_min,value_max)
        l_change=(gausrand()*value_change_delayed)*powerrand_std(value_min,value_max)
        voronoi_delayed(cell_number,4)=voronoi_delayed(cell_number,4)+l_change
        
        !voronoi_delayed(cell_number,4)=voronoi_current(cell_number,4)&
        !                                +powerrand(value_min,value_max)*value_change_delayed
                                        
        ! Check if the fracture edge values are in the domain boundaries
        x1 = voronoi_delayed(cell_number,1) - voronoi_delayed(cell_number,4)/2
        x2 = voronoi_delayed(cell_number,1) + voronoi_delayed(cell_number,4)/2
        IF (x1.LT.x_min) THEN 
                voronoi_delayed(cell_number,4) = abs(x2-x_min)                             ! recalculate the length
                voronoi_delayed(cell_number,1) = x_min + voronoi_delayed(cell_number,4)/2  ! recalculate the centre
        ELSE IF (x2.GT.x_max) THEN
                voronoi_delayed(cell_number,4) = abs(x_max-x1)                             ! recalculate the length
                voronoi_delayed(cell_number,1) = x_max - voronoi_delayed(cell_number,4)/2  ! recalculate the centre
        END IF   
END IF

CALL check_inside_prior()

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Option not imlemented for fw2_5d
SUBROUTINE change_option_delayed()
        ! Changes data option parameters
        IMPLICIT NONE
        DOUBLE PRECISION        :: random_deviate

        delayed_opt=.TRUE.

        write(*,*) 'change optional parameter delayed ', TRIM(int2str(opt_parameters)), ' at chain ', TRIM(int2str(chain_number))
        
        cell_number=CEILING(unirand()*ncell_current) ! randomly choose a cell to perturb
        IF (cell_number.EQ.0) cell_number=1

        voronoi_proposed(cell_number,5)=voronoi_current(cell_number,5)+gausrand()*opt_change_delayed

        CALL check_inside_prior()
        
END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE check_inside_prior()
! Check that the proposed model falls inside the prior
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

END MODULE
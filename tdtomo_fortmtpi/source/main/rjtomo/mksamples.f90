! ********************************************************************* !
! Program to generate Voronoi samples using the rj-McMC algorithm
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
! Scripts adapted to inversion of preferential pathways:
! Copyright (C) 2026 Léa Lelimouzin
! ********************************************************************* !


PROGRAM mksamples

USE forward
USE acquisition
USE input_files
USE chains
USE random_generator
USE string_conversion
USE noise_calculation
USE alpha_calculation
USE chain_update
USE tempering
USE debug
USE MPI
IMPLICIT NONE
REAL                    :: timer_start,timer_end
DOUBLE PRECISION        :: seconds_elapsed,hours_elapsed,minutes_elapsed
INTEGER                 :: sleep_time,sleep_total
INTEGER                 :: ii,jj,kk,nn
INTEGER                 :: iteration
DOUBLE PRECISION        :: random_deviate
INTEGER                 :: info,last_accepted
LOGICAL                 :: length_ok
LOGICAL                 :: accepted,delayed,accepted_main,accepted_delayed
LOGICAL                 :: time_is_up,temper
CHARACTER(LEN=70)       :: message
! For MPI
INTEGER                 :: status(MPI_STATUS_SIZE)
INTEGER                 :: number_of_cores,core_number,ierror,rr
INTERFACE
        SUBROUTINE print_status(start_end,acceptances,iteration,number_of_chains)
                INTEGER, INTENT(IN)                                     :: start_end
                DOUBLE PRECISION, DIMENSION(:), INTENT(IN), OPTIONAL    :: acceptances
                INTEGER, INTENT(IN), OPTIONAL                           :: iteration,number_of_chains
        END SUBROUTINE
END INTERFACE
INTERFACE
        SUBROUTINE check_output_files(length_ok)
                LOGICAL, INTENT(OUT)    :: length_ok
        END SUBROUTINE
END INTERFACE
INTERFACE
        SUBROUTINE initialize_models()
        END SUBROUTINE
END INTERFACE
INTERFACE
        SUBROUTINE write_debug_info_mk(info,iteration,random_deviate,accepted,delayed)
                INTEGER, INTENT(IN)                     :: info
                INTEGER, INTENT(IN), OPTIONAL           :: iteration
                DOUBLE PRECISION, INTENT(IN), OPTIONAL  :: random_deviate
                LOGICAL, INTENT(IN), OPTIONAL           :: accepted
                LOGICAL, INTENT(IN), OPTIONAL           :: delayed
        END SUBROUTINE
END INTERFACE
INTERFACE
        SUBROUTINE save_first_model()
        END SUBROUTINE
END INTERFACE
INTERFACE
        SUBROUTINE write_last_model()
        END SUBROUTINE
END INTERFACE
INTERFACE
        SUBROUTINE write_last_model_tmp(iteration)
                INTEGER, INTENT(IN)           :: iteration
        END SUBROUTINE
END INTERFACE
INTERFACE
        SUBROUTINE write_chains_to_file_tmp(iteration)
                INTEGER, INTENT(IN)           :: iteration
        END SUBROUTINE
END INTERFACE
INTERFACE
        SUBROUTINE write_chains_to_file()
        END SUBROUTINE
END INTERFACE
INTERFACE
        SUBROUTINE choose_step(iteration,step)
                INTEGER, INTENT(IN)     :: iteration
                INTEGER, INTENT(OUT)    :: step
        END SUBROUTINE
END INTERFACE
INTERFACE
        SUBROUTINE write_inversion_parameters()
        END SUBROUTINE
END INTERFACE


! Start timer
CALL cpu_time(timer_start)
time_is_up=.FALSE.

! Start iteration counter
iteration=0

! Start parallelisation
CALL MPI_INIT(ierror)
CALL MPI_COMM_SIZE(MPI_COMM_WORLD,number_of_cores,ierror)
CALL MPI_COMM_RANK(MPI_COMM_WORLD,core_number,ierror)
!  core_number=0

! Get chain number
!  chain_number=1
 chain_number=core_number+1
 number_of_chains=number_of_cores
 WRITE(chain_name,*)chain_number
 chain_name=TRIM(ADJUSTL(chain_name))

IF (debug_mode) CALL initialize_debug(chain_name)

IF (core_number.EQ.0) CALL print_status(1)

WRITE(message,*)'Hello from rj-McMC tomography and chain ',TRIM(chain_name),' on core ',TRIM(int2str(core_number)),'!'
IF (core_number.EQ.0) THEN
        WRITE(*,*)message
        DO ii=1,number_of_cores-1
                CALL MPI_RECV(message,70,MPI_CHARACTER,ii,1,MPI_COMM_WORLD,status,ierror)
                WRITE(*,*)message
        END DO
ELSE
        CALL MPI_SSEND(message,70,MPI_CHARACTER,0,1,MPI_COMM_WORLD,ierror)
END IF

! Get run number from command line
CALL getarg(1,run_name) 
CALL check_run_number()

! Read mksamples.in
CALL read_mksamples_in()
IF (input_error) THEN
        CALL MPI_FINALIZE(ierror)
        STOP
END IF

! Set up sleep time between chain status display
sleep_total=0
IF (add_sleep.EQ.1) THEN
        sleep_time=1
ELSE IF (add_sleep.EQ.0) THEN
        sleep_time=0
END IF

! Read data files
IF (debug_mode) WRITE(debug_unit,*)'--- Chain ',TRIM(chain_name),': reading data ---'
CALL read_data(0)
IF (debug_mode) WRITE(debug_unit,*)'--- Chain ',TRIM(chain_name),': read data ---'

! Initialise random number generator
IF (debug_mode) WRITE(debug_unit,*)'--- Chain ',TRIM(chain_name),': initializing random numbers ---'
CALL initialize_random_generator(chain_number)
IF (debug_mode) WRITE(debug_unit,*)'--- Chain ',TRIM(chain_name),': initialized random numbers ---'

! Initialise Markov chains and allocate necessary arrays
CALL initialize_chains()

! Initialise models
IF (run_number.EQ.1) THEN
        info=1
        DO WHILE (info.NE.0)
                CALL initialize_models()
                IF (no_data.EQ.0) CALL forward_modelling(ncell_current,voronoi_current,data_current,info,chain_number,iteration)
                IF (no_data.EQ.1) info=0
        END DO
        !residuals(:,1)=observed_data-LOG(data_current) 
        residuals(:,1)=observed_data-data_current 
ELSE IF (run_number.GT.1) THEN
        info=0
        CALL initialize_models()
        IF (samples_so_far.GE.max_number_of_samples) THEN
                IF (core_number.EQ.0) THEN
                        WRITE(*,*)
                        WRITE(*,*)'********************************************'
                        WRITE(*,*)'Maximum number of samples (',TRIM(int2str(max_number_of_samples)),')'
                        WRITE(*,*)'has been reached.'
                        WRITE(*,*)'TERMINATING PROGRAM!!'
                        WRITE(*,*)'********************************************'
                        WRITE(*,*)
                END IF
                CALL MPI_FINALIZE(ierror)
                STOP
        END IF
        IF (core_number.EQ.0) THEN
                accepted_so_far_all(1)=accepted_so_far
                DO ii=1,number_of_cores-1
                        CALL MPI_RECV(accepted_so_far_other,1,MPI_INTEGER,ii,1,MPI_COMM_WORLD,status,ierror)
                        accepted_so_far_all(ii+1)=accepted_so_far_other
                END DO
        ELSE
                CALL MPI_SSEND(accepted_so_far_other,1,MPI_INTEGER,0,1,MPI_COMM_WORLD,ierror)
        END IF

        IF (core_number.EQ.0) THEN
                CALL CHDIR(TRIM("results/"))
                CALL check_output_files(length_ok)
                CALL CHDIR("../")
        END IF
        CALL MPI_BCAST(length_ok,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierror)
        IF (.NOT.length_ok) THEN
                CALL MPI_FINALIZE(ierror)
                STOP
        END IF
END IF

! Initialise tempering
IF (parallel_tempering.EQ.1) CALL initialize_tempering()

! Write temperature information if in debug mode
IF (debug_mode) THEN
        DO ii=0,number_of_cores-1
                CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)
                IF (ii.EQ.core_number) CALL write_debug_info_mk(1)
        END DO
        CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)
END IF

! Write message to screen in case inversion is run without data
IF (core_number.EQ.0.AND.no_data.EQ.1) THEN
        WRITE(*,*)
        WRITE(*,*)'*** RUNNING INVERSION WITHOUT DATA ***'
END IF

! Write initial point of chain to screen
IF (core_number.EQ.0) THEN
        WRITE(*,*)
        IF (run_number.EQ.1) WRITE(*,*)'STARTING CHAINS FROM SAMPLE ',TRIM(int2str(samples_at_start))
        IF (run_number.GT.1) WRITE(*,*)'RE-STARTING CHAINS FROM SAMPLE ',TRIM(int2str(samples_at_start+1))
        WRITE(*,*)
        WRITE(*,*)'Max run time is: ',TRIM(dble2str(max_runtime)),' hr (',&
                  TRIM(dble2str(max_runtime*DBLE(3600))),' sec)'
        WRITE(*,*)
END IF

! Write number of measurements to screen
IF (core_number.EQ.0) THEN
        WRITE(*,*)'Number of measurements: ',TRIM(int2str(number_of_measurements))
        WRITE(*,*)
END IF

! Write initial number of cells to screen
WRITE(message,*)'Initial number of cells on chain ',TRIM(chain_name),': ',TRIM(int2str(ncell_current))

IF (core_number.EQ.0) THEN
        WRITE(*,*)message
        DO ii=1,number_of_cores-1
                CALL MPI_RECV(message,70,MPI_CHARACTER,ii,1,MPI_COMM_WORLD,status,ierror)
                WRITE(*,*)message
        END DO
ELSE
        CALL MPI_SSEND(message,70,MPI_CHARACTER,0,1,MPI_COMM_WORLD,ierror)
END IF
IF (core_number.EQ.0) WRITE(*,*)

! Calculate data uncertainties
IF (no_data.EQ.0) CALL calculate_noise(0,data_current)

! Write initial model to screen if in debug mode
IF (debug_mode) THEN
        DO ii=0,number_of_cores-1
                CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)
                IF (ii.EQ.core_number) CALL write_debug_info_mk(2)
        END DO
CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)
END IF

! Calculate the initial likelihood and misfit
IF (no_data.EQ.0) CALL calculate_likelihood(number_of_measurements,observed_data,&
                                            data_current,noise_current,&
                                            likelihood_type,likelihood_current,&
                                            misfit_current)

write(*,*)'mksamples initialization likelihood_type: ',likelihood_type ! TO BE DELETED
write(*,*)'mksamples initialization likelihood_current: ',likelihood_current ! TO BE DELETED
write(*,*)'mksamples initialization misfit_current: ',misfit_current ! TO BE DELETED
write(*,*)'mksamples initialization data_current: ',data_current ! TO BE DELETED

! Write initial misfit to screen
WRITE(message,*)'Initial misfit on chain ',TRIM(chain_name),': ',TRIM(dble2str(misfit_current))
IF (core_number.EQ.0) THEN
        WRITE(*,*)message
        DO ii=1,number_of_cores-1
                CALL MPI_RECV(message,70,MPI_CHARACTER,ii,1,MPI_COMM_WORLD,status,ierror)
                WRITE(*,*)message
        END DO
ELSE
        CALL MPI_SSEND(message,70,MPI_CHARACTER,0,1,MPI_COMM_WORLD,ierror)
END IF
IF (core_number.EQ.0) WRITE(*,*)

! Write initial likelihood to screen if in debug mode
IF (debug_mode) THEN
        DO ii=0,number_of_cores-1
                CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)
                IF (ii.EQ.core_number) CALL write_debug_info_mk(3,iteration)
        END DO
CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)
END IF

! Save the very first model of the chain
IF (run_number.EQ.1) CALL save_first_model()

! Start the sampling
DO WHILE (iteration.LT.number_of_samples.AND.&
          .NOT.time_is_up.AND.&
          samples_at_start+iteration.LT.max_number_of_samples)
        
        iteration=iteration+1

        CALL reset_proposed_models()
        accepted=.FALSE.
        delayed=.FALSE.
        accepted_main=.FALSE.
        accepted_delayed=.FALSE.
        temper=.FALSE.
        IF (parallel_tempering.EQ.1.AND.&
            samples_at_start+iteration.GE.tempering_start) temper=.TRUE.
        
        IF (debug_mode) CALL write_debug_info_mk(-1,iteration)
        
        ! Choose the type of step
        CALL choose_step(iteration,step)
        CALL step_type()

        ! Apply the step
        IF (birth) THEN
                CALL add_cell(modelling_method)
        ELSE IF (death) THEN
                CALL delete_cell(modelling_method)
        ELSE IF (move) THEN
                CALL move_cell(modelling_method)
        ELSE IF (value) THEN
                CALL change_value(modelling_method)
        ELSE IF (noise) THEN
                CALL change_noise()
        ELSE IF (option) THEN ! Option only implemented for ert_ddp_2_5d forward model
                CALL change_option()
        END IF
        

        ! Only go through this if we are inside the prior (otherwise model is automatically rejected)
        IF (inside_prior) THEN
                
                IF (.NOT.noise.AND.no_data.EQ.0) THEN ! if we are at a step involving a change in the model...
                        CALL forward_modelling(ncell_proposed,voronoi_proposed,data_proposed,info,chain_number,iteration) ! ...then model the data
                        IF (info.NE.0) THEN
                                WRITE(*,*)'********************************************'
                                WRITE(*,*)'Problem with forward modelling'
                                WRITE(*,*)'at iteration ',TRIM(int2str(iteration)),' on chain ',chain_name
                                WRITE(*,*)'REPEATING SAMPLE!!'
                                WRITE(*,*)'********************************************'
                                iteration=iteration-1
                                CYCLE
                        END IF
                END IF             

                IF (no_data.EQ.0) THEN
                        CALL calculate_noise(1,data_proposed)
                        CALL calculate_likelihood(number_of_measurements,observed_data,&
                                                  data_proposed,noise_proposed,&
                                                  likelihood_type,likelihood_proposed,&
                                                  misfit_proposed)
                        CALL calculate_normalization(number_of_measurements,data_proposed,noise_current,&
                                                     noise_proposed,likelihood_normalization)
                END IF
                
                CALL calculate_alpha(likelihood_current,likelihood_proposed,&
                                     likelihood_normalization,temper,alpha)
                
                write(*,*) 'NOISE proposed at chain ', TRIM(chain_name), ' =', noise_proposed(1)
                write(*,*) 'NOISE current at chain ', TRIM(chain_name), ' =', noise_current(1)
                write(*,*) '  Likelihood proposed at chain ', TRIM(chain_name), ' =', likelihood_proposed
                write(*,*) '  Likelihood current at chain ', TRIM(chain_name), ' =', likelihood_current 
                write(*,*) '  Likelihood normalization at chain ', TRIM(chain_name), ' =',likelihood_normalization

                ! Metropolis-Hastings acceptance criterion
                random_deviate=unirand()
                write(*,*) 'MH criterion at chain ',TRIM(chain_name), ' =', LOG(random_deviate)
                IF (LOG(random_deviate).LE.alpha) accepted=.TRUE.
                
                IF (move.AND.accepted) accepted_main=.TRUE.
                IF (value.AND.accepted) accepted_main=.TRUE.
                IF (option.AND.accepted) accepted_main=.TRUE. ! Option only implemented for ert_ddp_2_5d forward model
                write(*,*) '     Is the proposed model accepted at chain ', TRIM(chain_name), ' ? ', accepted
        END IF
        
        ! Write step to screen if in debug mode
        IF (debug_mode) THEN
                DO ii=0,number_of_cores-1
                        CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)
                        IF (ii.EQ.core_number) CALL write_debug_info_mk(4,iteration,random_deviate,accepted)
                END DO
        END IF

        ! Delayed rejection for "move" and "value" steps
        IF ((move.OR.value).AND..NOT.accepted) THEN
                
                delayed=.TRUE.
                IF (move) CALL move_cell_delayed(modelling_method)
                IF (value) CALL change_value_delayed(modelling_method)
                IF (option) CALL change_option_delayed() ! Option only implemented for ert_ddp_2_5d forward model
                
                ! Only go through this if we are inside the prior (otherwise model is automatically rejected)
                IF (inside_prior_delayed) THEN
                        
                        IF (no_data.EQ.0) THEN
                                
                                CALL forward_modelling(ncell_proposed,voronoi_delayed,data_delayed,info,chain_number,iteration)
                                IF (info.NE.0) THEN
                                        WRITE(*,*)'********************************************'
                                        WRITE(*,*)'Problem with forward modelling'
                                        WRITE(*,*)'at iteration ',TRIM(int2str(iteration)),' on chain ',chain_name
                                        WRITE(*,*)'REPEATING SAMPLE!!'
                                        WRITE(*,*)'********************************************'
                                        iteration=iteration-1
                                        CYCLE
                                END IF

                                CALL calculate_noise(2,data_delayed)
                                CALL calculate_likelihood(number_of_measurements,observed_data,data_delayed,noise_delayed,&
                                                          likelihood_type,likelihood_delayed,misfit_delayed)
                                CALL calculate_normalization(number_of_measurements,data_delayed,noise_current,&
                                                             noise_delayed,likelihood_normalization)
                                
                        END IF
                        
                        CALL calculate_alpha_delayed(likelihood_current,likelihood_proposed,&
                                                     likelihood_delayed,likelihood_normalization,&
                                                     temper,alpha_delayed)
                        
                        write(*,*) 'NOISE proposed DELAYED at chain ', TRIM(chain_name), '=', noise_proposed(1)
                        write(*,*) 'NOISE current DELAYED at chain ', TRIM(chain_name), '=', noise_current(1)
                        write(*,*) '  Likelihood proposed DELAYED at chain ', TRIM(chain_name), '=', likelihood_proposed
                        write(*,*) '  Likelihood current DELAYED at chain ', TRIM(chain_name), '=', likelihood_current 
                        write(*,*) '  Likelihood normalization DELAYED at chain ', TRIM(chain_name), '=',likelihood_normalization

                        ! Metropolis-Hastings acceptance criterion
                        random_deviate=unirand()
                        IF (LOG(random_deviate).LE.alpha_delayed) accepted=.TRUE.
                        write(*,*) 'MH acceptance criterion DELAYED at chain ',TRIM(chain_name),&
                                         &'=', LOG(random_deviate)
                        IF (accepted) accepted_delayed=.TRUE.
                        write(*,*) '     Is the poposed DELAYED model accepted at chain ', TRIM(chain_name), '? ', accepted_delayed
                END IF
                
        END IF
        
        ! Write step to screen if in debug mode
        IF (debug_mode) THEN
                DO ii=0,number_of_cores-1
                        CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)
                        IF (ii.EQ.core_number) CALL write_debug_info_mk(5,iteration,random_deviate,accepted,delayed)
                END DO
        END IF
        
        ! Update the Markov chain
        IF (inside_prior.AND.accepted.AND..NOT.delayed) THEN ! the model was accepted at the first attempt
                CALL update_models(1)
                last_accepted=iteration
        ELSE IF (inside_prior_delayed.AND.accepted.AND.delayed) THEN ! the model was accepted with delayed rejection
                CALL update_models(2)
                last_accepted=-iteration
        END IF
        CALL update_acceptance(accepted,accepted_main,accepted_delayed)
        CALL update_chain_arrays(accepted,iteration)
        
        ! Apply parallel tempering
        IF (temper) THEN
                ! Choose temperatures
                IF (chain_number.EQ.1) THEN
                        CALL choose_T1()
                        CALL choose_T2()
                END IF
                CALL MPI_BCAST(T1_chain,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
                CALL MPI_BCAST(T2_chain,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
                
                ! Get the likelihood and noise to be used
                IF (chain_number.EQ.T1_chain) THEN
                        likelihood_T1=likelihood_current
                        noise_T1=noise_current
                END IF
                !CALL MPI_BCAST(likelihood_T1,1,MPI_DOUBLE_PRECISION,T1_chain-1,MPI_COMM_WORLD,ierror)
                !CALL MPI_BCAST(noise_T1,number_of_measurements,MPI_DOUBLE_PRECISION,T1_chain-1,MPI_COMM_WORLD,ierror)
                IF (chain_number.EQ.T2_chain) THEN
                        likelihood_T2=likelihood_current
                        noise_T2=noise_current
                END IF
                !CALL MPI_BCAST(likelihood_T2,1,MPI_DOUBLE_PRECISION,T2_chain-1,MPI_COMM_WORLD,ierror)
                !CALL MPI_BCAST(noise_T2,number_of_measurements,MPI_DOUBLE_PRECISION,T2_chain-1,MPI_COMM_WORLD,ierror)
                
                IF (T1_chain.NE.1) THEN
                        IF (chain_number.EQ.1) THEN
                                CALL MPI_RECV(likelihood_T1,1,MPI_DOUBLE_PRECISION,&
                                              T1_chain-1,1,MPI_COMM_WORLD,status,ierror)
                                CALL MPI_RECV(noise_T1,number_of_measurements,MPI_DOUBLE_PRECISION,&
                                              T1_chain-1,1,MPI_COMM_WORLD,status,ierror)
                        ELSE IF (chain_number.EQ.T1_chain) THEN
                                CALL MPI_SSEND(likelihood_T1,1,MPI_DOUBLE_PRECISION,&
                                               0,1,MPI_COMM_WORLD,ierror)
                                CALL MPI_SSEND(noise_T1,number_of_measurements,MPI_DOUBLE_PRECISION,&
                                               0,1,MPI_COMM_WORLD,ierror)
                        END IF
                END IF
                IF (T2_chain.NE.1) THEN
                        IF (chain_number.EQ.1) THEN
                                CALL MPI_RECV(likelihood_T2,1,MPI_DOUBLE_PRECISION,&
                                              T2_chain-1,1,MPI_COMM_WORLD,status,ierror)
                                CALL MPI_RECV(noise_T2,number_of_measurements,MPI_DOUBLE_PRECISION,&
                                              T2_chain-1,1,MPI_COMM_WORLD,status,ierror)
                        ELSE IF (chain_number.EQ.T2_chain) THEN
                                CALL MPI_SSEND(likelihood_T2,1,MPI_DOUBLE_PRECISION,&
                                               0,1,MPI_COMM_WORLD,ierror)
                                CALL MPI_SSEND(noise_T2,number_of_measurements,MPI_DOUBLE_PRECISION,&
                                               0,1,MPI_COMM_WORLD,ierror)
                        END IF
                END IF
                
                ! Check whether the temperatures are swapping
                IF (chain_number.EQ.1) THEN
                        IF (inside_prior_delayed) THEN 
                                CALL do_tempering(data_delayed)
                        ELSE
                                CALL do_tempering(data_proposed)
                        END IF
                        CALL update_temperatures(iteration)
                END IF
                !IF (chain_number.EQ.1) THEN
                !        CALL do_tempering()
                !        CALL update_temperatures(iteration)
                !END IF
                CALL MPI_BCAST(temperature_values,number_of_temperatures,MPI_DOUBLE_PRECISION,&
                               0,MPI_COMM_WORLD,ierror)
                CALL MPI_BCAST(temperature_indices,number_of_temperatures,MPI_INT,&
                               0,MPI_COMM_WORLD,ierror)
                temperature=temperature_values(chain_number)
                
                ! Write temperature information if in debug mode
                IF (debug_mode) THEN
                        DO ii=0,number_of_cores-1
                                CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)
                                IF (ii.EQ.core_number) CALL write_debug_info_mk(6,iteration)
                        END DO
                        CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)
                END IF
                
        END IF
        
        ! Write chain status to screen
        IF (MOD(iteration,display_interval).EQ.0) THEN
                DO ii=0,number_of_cores-1
                        CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)
                        IF (ii.EQ.core_number) CALL print_status(3,acceptance_all,iteration,number_of_cores)
                        ! Sleeping the CPU allows the status of each chain to be printed in order
                        IF (add_sleep.EQ.1) THEN
                                CALL sleep(sleep_time)
                                sleep_total=sleep_total+sleep_time
                        END IF
                END DO
                IF (core_number.EQ.0) CALL print_status(4,acceptance_all)
        END IF
        
        CALL cpu_time(timer_end)
        seconds_elapsed=DBLE(timer_end-timer_start+REAL(sleep_total))
        hours_elapsed=FLOOR(seconds_elapsed/DBLE(3600))
        minutes_elapsed=FLOOR((seconds_elapsed-hours_elapsed*DBLE(3600))/DBLE(60))
        
        IF (MOD(iteration,display_interval).EQ.0.AND.core_number.EQ.0) THEN
                WRITE(*,*)'TIME ELAPSED: ',TRIM(int2str(INT(hours_elapsed))),'h ',&
                                           TRIM(int2str(INT(minutes_elapsed))),'m ',&
                                           TRIM(int2str(INT(seconds_elapsed&
                                                          -(hours_elapsed*DBLE(3600)&
                                                          +minutes_elapsed*DBLE(60))))),'s'
                WRITE(*,*)'--------------------------------------'
        END IF
        IF (seconds_elapsed.GE.max_runtime*DBLE(3600)) time_is_up=.TRUE.
        CALL MPI_BCAST(time_is_up,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierror)
        
        IF (debug_mode) WRITE(debug_unit,*)'----- Time on chain ',TRIM(chain_name),': ',&
                             TRIM(dble2str(seconds_elapsed)),' seconds -----'
        
        IF (debug_mode) CALL write_debug_info_mk(-2,iteration)
        
        !CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)

        ! Write chains content to file every x iteration(s)
        IF (MOD(iteration,1).EQ.0) THEN
                CALL CHDIR( TRIM("results/") )
                CALL write_last_model_tmp(iteration)
                CALL write_chains_to_file_tmp(iteration)
                CALL CHDIR('../')
        END IF

        
END DO

IF (time_is_up.AND.core_number.EQ.0) THEN
        WRITE(*,*)
        WRITE(*,*)'TIME IS UP - ITERATION ',TRIM(int2str(iteration)),&
                  '/',TRIM(int2str(number_of_samples))
END IF

! Get final number of samples (in this run and total)
samples_in_run=iteration
samples_so_far=samples_at_start+samples_in_run

! Save last model, Markov chains and inversion parameters to file
IF (parallel_tempering.EQ.1) THEN
        CALL MPI_BCAST(acceptance_tempering,3,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierror)
END IF
CALL CHDIR( TRIM("results/") )
!CALL write_chains_to_file()
CALL write_last_model()
IF (core_number.EQ.0) CALL write_inversion_parameters()
CALL CHDIR('../')

! Save a copy of the results per run
CALL SYSTEM('scp -r results/ results_run_'//TRIM(int2str(run_number)))

! Stop timer
CALL cpu_time(timer_end)
seconds_elapsed=DBLE(timer_end-timer_start+REAL(sleep_total))
hours_elapsed=seconds_elapsed/DBLE(3600)

! Terminate parallelisation
CALL MPI_FINALIZE(ierror)

IF (core_number.EQ.0) WRITE(*,*)
IF (core_number.EQ.0) PRINT '("Elapsed time is ",f6.2," hours.")',hours_elapsed
IF (debug_mode) WRITE(debug_unit,*)'----- Random count ',TRIM(chain_name),': ',TRIM(int2str(random_counter)),' -----'
IF (core_number.EQ.0) CALL print_status(2)

IF (debug_mode) CALL close_debug(chain_name)

END PROGRAM
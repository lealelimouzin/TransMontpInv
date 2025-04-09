! ********************************************************************* !
! Program to process the outputs of mksamples.f90
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


PROGRAM procsamples

USE input_files
USE chains, ONLY: chain_number,chain_name,step,&
                  samples_unit,ncell_unit,noise_unit,misfits_unit,&
                  acceptance_unit,temperatures_unit,residuals_unit
USE processing
USE acquisition, ONLY: observed_data
USE string_conversion
USE forward
USE debug
IMPLICIT NONE
REAL                            :: timer_start,timer_end
INTEGER                         :: ii,jj,kk,nn,pp
INTEGER                         :: core_number=0
LOGICAL                         :: temper
INTEGER                         :: samples_in_ensemble
CHARACTER(LEN=30)               :: debug_ext
INTERFACE
        SUBROUTINE read_inversion_parameters()
        END SUBROUTINE
END INTERFACE
INTERFACE
        SUBROUTINE write_debug_info_pr(info,iteration)
                INTEGER, INTENT(IN)                     :: info
                INTEGER, INTENT(IN), OPTIONAL           :: iteration
        END SUBROUTINE
END INTERFACE
INTERFACE
        SUBROUTINE write_results_to_file(samples_in_ensemble)
                INTEGER, INTENT(IN)                     :: samples_in_ensemble
        END SUBROUTINE
END INTERFACE

! Start timer
CALL cpu_time(timer_start)

IF (debug_mode) debug_ext='_proc'
IF (debug_mode) CALL initialize_debug(debug_ext)

! Read procsamples.in
CALL read_procsamples_in()
IF (input_error) THEN
        STOP
END IF

! Read inversion parameters file
CALL read_inversion_parameters()

! Check Markov chain parameters for processing
CALL check_processing_parameters()
IF (input_error) THEN
        STOP
END IF

! Initialise arrays
CALL initialize_processing()

! Sort out tmperatures if parallel tempering was employed
IF (parallel_tempering.EQ.1) CALL arrange_temperatures()
IF (parallel_tempering.EQ.0) temper=.TRUE.
IF (debug_mode) CALL write_debug_info_pr(1)

! Start counter for number of samples in ensemble
samples_in_ensemble=0

! Loop over all chains (i.e. over all processors on which the chains were run)
WRITE(*,*)'Looping over chains...'
DO chain_number=1,number_of_chains
        
        ! Get chain number
        WRITE(chain_name,*)chain_number
        chain_name=TRIM(ADJUSTL(chain_name))
        
        ! Check whether the chain should be ignored
        IF (number_of_bad_chains.GT.0) THEN
                IF (ANY(bad_chains.EQ.chain_number)) THEN
                        WRITE(*,*)'   ...chain ',TRIM(chain_name),' ignored!'
                        CYCLE
                END IF
        END IF
        WRITE(*,*)'   ...processing chain ',TRIM(chain_name),' of ',TRIM(int2str(number_of_chains))
        
        CALL open_files()
        
        ! Read and discard burn-in samples
        iteration=0
        DO WHILE (iteration.LE.burn_in)
                IF (iteration.EQ.0) THEN ! If we are at the first sample of the chain...
                        ! 1) read and save the number of cells
                        READ(ncell_unit)pp,ncell_current
                        ncells_iter(iteration+1)=ncells_iter(iteration+1)+ncell_current
                        !write(*,*)'Sample',iteration,'has',ncell_current,'cells'
                        ncell_all_first(chain_number)=ncell_current
                        ! 2) read and save the Voronoi geometry
                        DO nn=1,ncell_current
                                READ(samples_unit)voronoi_current(nn,1),&
                                                  voronoi_current(nn,2),&
                                                  voronoi_current(nn,3),&
                                                  voronoi_current(nn,4)
                        END DO
                        voronoi_all_first(:,:,chain_number)=voronoi_current
                        ! 3) read the data noise parameters
                        IF (noise_type.NE.4) THEN
                                DO nn=1,number_of_datasets
                                        IF (noise_type.EQ.0) READ(noise_unit)pp,sigma_current(nn)
                                        IF (noise_type.EQ.1.OR.&
                                            noise_type.EQ.2.OR.&
                                            noise_type.EQ.5) READ(noise_unit)pp,&
                                                                             a_current(nn),&
                                                                             b_current(nn)
                                        IF (noise_type.EQ.3) READ(noise_unit)pp,lambda_current(nn)
                                END DO
                        ELSE IF (noise_type.EQ.4.AND.chain_number.EQ.1) THEN
                                DO nn=1,number_of_measurements
                                        READ(noise_unit)noise_current(nn)
                                END DO
                                CALL evaluate_posteriors(.FALSE.,.TRUE.)
                        END IF
                        ! 4) read the misfit
                        READ(misfits_unit)pp,misfit_current
                        misfit_iter(iteration+1)=misfit_iter(iteration+1)+misfit_current
                        ! 5) read the residuals
                        READ(residuals_unit)save_residuals
                        IF (save_residuals.EQ.1) THEN
                                DO nn=1,number_of_measurements
                                        READ(residuals_unit)residuals(nn,1)
                                        residuals(nn,2)=residuals(nn,1)/observed_data(nn)*100
                                END DO
                        END IF
                        IF (debug_mode) CALL write_debug_info_pr(2)
                ELSE ! ...otherwise...
                        ! 1) read the number of cells
                        READ(ncell_unit)pp,ncell_current
                        ncells_iter(iteration+1)=ncells_iter(iteration+1)+ncell_current
                        !write(*,*)'Sample',iteration,'has',ncell_current,'cells'
                        ! 2) read the Markov chain step
                        READ(samples_unit)pp,&
                                          step,&
                                          accept,&
                                          cell_number,&
                                          x_new,&
                                          y_new,&
                                          z_new,&
                                          value_new
                        !write(*,*)iteration,pp,step,accept,cell_number,x_new,y_new,z_new,value_new
                        CALL update_model()
                        ! 3) read the data noise parameters
                        IF (noise_type.NE.4) THEN
                                DO nn=1,number_of_datasets
                                        IF (noise_type.EQ.0) READ(noise_unit)pp,sigma_current(nn)
                                        IF (noise_type.EQ.1.OR.&
                                            noise_type.EQ.2.OR.&
                                            noise_type.EQ.5) READ(noise_unit)pp,&
                                                                             a_current(nn),&
                                                                             b_current(nn)
                                        IF (noise_type.EQ.3) READ(noise_unit)pp,lambda_current(nn)
                                END DO
                        END IF
                        ! 4) read the misfit
                        READ(misfits_unit)pp,misfit_current
                        misfit_iter(iteration+1)=misfit_iter(iteration+1)+misfit_current
                        ! 5) read the residuals
                        IF (save_residuals.EQ.1.AND.accept.EQ.1.AND.step.NE.5) THEN
                                DO nn=1,number_of_measurements
                                        READ(residuals_unit)residuals(nn,1)
                                        residuals(nn,2)=residuals(nn,1)/observed_data(nn)*100
                                END DO
                        END IF
                        IF (debug_mode) CALL write_debug_info_pr(4,iteration)
                END IF
                iteration=iteration+1
        END DO
        
        ! Start collecting samples after the burn-in period
        DO WHILE (iteration.LE.samples_so_far)
                ! 1) read the number of cells
                READ(ncell_unit)pp,ncell_current
                ncells_iter(iteration+1)=ncells_iter(iteration+1)+ncell_current
                !write(*,*)'Sample',iteration,'has',ncell_current,'cells'
                ! 2) read the Markov chain step
                READ(samples_unit)pp,&
                                  step,&
                                  accept,&
                                  cell_number,&
                                  x_new,&
                                  y_new,&
                                  z_new,&
                                  value_new
                !write(*,*)iteration,pp,step,accept,cell_number,x_new,y_new,z_new,value_new
                CALL update_model()
                ! 3) read the data noise parameters
                IF (noise_type.NE.4) THEN
                        DO nn=1,number_of_datasets
                                IF (noise_type.EQ.0) READ(noise_unit)pp,sigma_current(nn)
                                IF (noise_type.EQ.1.OR.&
                                    noise_type.EQ.2.OR.&
                                    noise_type.EQ.5) READ(noise_unit)pp,&
                                                                     a_current(nn),&
                                                                     b_current(nn)
                                IF (noise_type.EQ.3) READ(noise_unit)pp,lambda_current(nn)
                        END DO
                END IF
                ! 4) read the misfit
                READ(misfits_unit)pp,misfit_current
                misfit_iter(iteration+1)=misfit_iter(iteration+1)+misfit_current
                ! 5) read the residuals
                IF (save_residuals.EQ.1.AND.accept.EQ.1.AND.step.NE.5) THEN
                        DO nn=1,number_of_measurements
                                READ(residuals_unit)residuals(nn,1)
                                residuals(nn,2)=residuals(nn,1)/observed_data(nn)*100
                        END DO
                END IF
                IF (debug_mode) CALL write_debug_info_pr(4,iteration)
                
                ! If parallel tempering was used, only keep samples at target temperature T=1
                IF (parallel_tempering.EQ.1) THEN
                        IF (temperatures(iteration,chain_number).EQ.1) temper=.TRUE.
                        IF (temperatures(iteration,chain_number).NE.1) temper=.FALSE.
                END IF
                ! Add sample to posterior
                IF (MOD(iteration-burn_in,thinning).EQ.0.AND.temper) THEN
                        ! Number of samples in the final ensemble increases by 1
                        samples_in_ensemble=samples_in_ensemble+1
                        ! Turn Voronoi model into xyz model
                        CALL voronoi_to_xyz()
                        ! Add sample to posteriors
                        CALL evaluate_posteriors(.TRUE.,.TRUE.)
                        ! Add Voronoi model to the statistical moments
                        CALL evaluate_statistics()
                        ! Add Voronoi model to the node density array
                        CALL calculate_node_density()
                        ! Add residuals to residuals histogram
                        IF (save_residuals.EQ.1) CALL calculate_residuals()
                END IF
                iteration=iteration+1
        END DO
        ncell_all_last(chain_number)=ncell_current
        voronoi_all_last(:,:,chain_number)=voronoi_current
        
        CALL close_files()
        
END DO

! Evaluate various statistics
moment_1=moment_1/samples_in_ensemble
moment_2=moment_2/samples_in_ensemble
moment_3=moment_3/samples_in_ensemble
moment_4=moment_4/samples_in_ensemble
node_density=DBLE(node_count)/DBLE(samples_in_ensemble)
average=moment_1
stdev=SQRT(moment_2-moment_1**2)
WHERE (harmean.EQ.0)
        harmean=EPSILON(DBLE(0))
END WHERE
harmean=1/(harmean/DBLE(samples_in_ensemble))
root_mean_square=SQRT(moment_2)
skewness=(moment_3&
          -3*moment_2*average&
          +3*moment_1*average**2&
          -average**3)&
          /(SQRT(moment_2-average**2)**3)
kurtosis=(moment_4&
          -4*moment_3*average&
          +6*moment_2*average**2&
          -4*moment_1*average**3&
          +average**4)&
          /((moment_2-average**2)**2)&
          -3
CALL calculate_median_and_maximum(samples_in_ensemble)
CALL calculate_entropy(samples_in_ensemble)
ncells_iter=ncells_iter/(number_of_chains-number_of_bad_chains)
misfit_iter=misfit_iter/(number_of_chains-number_of_bad_chains)

! Save results to file
CALL write_results_to_file(samples_in_ensemble)

WRITE(*,*)
WRITE(*,*)'Total samples in ensemble: ',TRIM(int2str(samples_in_ensemble))

! Stop timer
CALL cpu_time(timer_end)
IF (core_number.EQ.0) WRITE(*,*)
IF (core_number.EQ.0) PRINT '("Elapsed time is ",f6.2," hours.")',(timer_end-timer_start)/DBLE(3600)

IF (debug_mode) CALL close_debug(debug_ext)

END PROGRAM
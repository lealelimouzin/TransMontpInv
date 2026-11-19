! ********************************************************************* !
! Module for reading input files (part of TERT package)
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

MODULE input_files

USE chains, ONLY: number_of_chains,number_of_measurements,number_of_parameters
CHARACTER (LEN=30), SAVE                        :: samples_file,ncell_file,noise_file,acceptance_file,&
                                                   residuals_file,misfits_file,temperatures_file
INTEGER, SAVE                                   :: save_residuals,add_sleep
INTEGER, SAVE                                   :: modelling_method
CHARACTER (LEN=30), SAVE                        :: forward_file
LOGICAL                                         :: input_error
CHARACTER (LEN=30), SAVE                        :: average_file,stdev_file,median_file,maximum_file,harmean_file,&
                                                   rms_file,entropy_file,skewness_file,kurtosis_file,ncellsi_file,misfiti_file,&
                                                   postncells_file,postnoise_file,postvalue_file,node_file,first_file,last_file,&
                                                   residuals_histo_file
INTEGER, SAVE                                   :: x_number_of_pixels,y_number_of_pixels,z_number_of_pixels,&
                                                   value_number_of_pixels,noise_number_of_pixels,residuals_number_of_pixels,&
                                                   vtx_or_xyz
DOUBLE PRECISION, SAVE                          :: x_pixel_size,y_pixel_size,z_pixel_size,max_abs_residual,max_percentage_residual
INTEGER, SAVE                                   :: burn_in,thinning,stop_sampling,number_of_bad_chains
INTEGER, DIMENSION(:), ALLOCATABLE, SAVE        :: bad_chains

CONTAINS


SUBROUTINE read_mksamples_in()
! Reads mksamples.in
USE acquisition, ONLY: electrode_file,data_file,datasets_file,number_of_datasets
USE chains, ONLY: number_of_datasets,ncell_min,ncell_max,value_min,value_max,x_min,x_max,y_min,y_max,z_min,z_max,&
                  noise_type,a_min,a_max,b_min,b_max,lambda_min,lambda_max,sigma_min,sigma_max,&
                  value_change,value_change_delayed,x_change,y_change,z_change,x_change_delayed,y_change_delayed,z_change_delayed,&
                  birth_change,a_change,b_change,lambda_change,sigma_change,likelihood_type,&
                  max_runtime,number_of_samples,max_number_of_samples,display_interval,no_data,&
                  parallel_tempering,jump_type,tempering_start,sigma_mat,rho_mat,number_of_parameters,opt_parameters,&
                  opt_min,opt_max,opt_change,opt_change_delayed,pi,powerlaw_a,aperture_fixed,conductivity_fixed
IMPLICIT NONE
INTEGER                 :: nn
CHARACTER(LEN=3)        :: section

input_error=.FALSE.
OPEN(UNIT=60,FILE='mksamples.in',STATUS='OLD')
READ(60,*)
READ(60,*)
READ(60,*)
READ(60,'(a30)')electrode_file
READ(60,'(a30)')data_file
READ(60,*)number_of_datasets
IF (number_of_datasets.LT.1) number_of_datasets=1
READ(60,'(a30)')datasets_file
READ(60,*)
READ(60,*)
READ(60,*)
READ(60,'(a30)')samples_file
READ(60,'(a30)')ncell_file
READ(60,'(a30)')noise_file
READ(60,'(a30)')misfits_file
READ(60,'(a30)')acceptance_file
READ(60,'(a30)')temperatures_file
READ(60,*)save_residuals
READ(60,'(a30)')residuals_file
READ(60,*)
READ(60,*)
READ(60,*)
READ(60,*)number_of_parameters
READ(60,*)
READ(60,*)
READ(60,*)
READ(60,*)ncell_min,ncell_max
READ(60,*)value_min,value_max
READ(60,*)powerlaw_a
READ(60,*)x_min,x_max
READ(60,*)y_min,y_max
READ(60,*)z_min,z_max
READ(60,*)noise_type
READ(60,*)a_min,a_max
READ(60,*)b_min,b_max
READ(60,*)lambda_min,lambda_max
READ(60,*)sigma_min,sigma_max
READ(60,*)sigma_mat
READ(60,*)rho_mat
! --- Option to inverse also on either the angle or the aperture of the fractures
READ(60,*)opt_parameters
IF (opt_parameters.EQ.1) THEN !angle
        READ(60,*)opt_min,opt_max
        opt_min=opt_min*pi/180 ! Convert degrees to radians
        opt_max=opt_max*pi/180 ! Convert degrees to radians
        READ(60,*)
        READ(60,*)aperture_fixed
ELSEIF (opt_parameters.EQ.2) THEN !aperture
        READ(60,*)
        READ(60,*)opt_min,opt_max
        READ(60,*)
ELSE
        READ(60,*)
        READ(60,*)
        READ(60,*)aperture_fixed
END IF
READ(60,*)conductivity_fixed
! ---
READ(60,*)
READ(60,*)
READ(60,*)
READ(60,*)value_change
READ(60,*)value_change_delayed
READ(60,*)x_change,y_change,z_change
READ(60,*)x_change_delayed,y_change_delayed,z_change_delayed
READ(60,*)birth_change
READ(60,*)a_change
READ(60,*)b_change
READ(60,*)lambda_change
READ(60,*)sigma_change
! --- Option to inverse also on either the angle or the aperture of the fractures
READ(60,*)opt_change
READ(60,*)opt_change_delayed
! ---
READ(60,*)
READ(60,*)
READ(60,*)
READ(60,*)likelihood_type
READ(60,*)
READ(60,*)
READ(60,*)
READ(60,*)modelling_method
READ(60,'(a30)')forward_file
READ(60,*)
READ(60,*)
READ(60,*)
READ(60,*)max_runtime
READ(60,*)number_of_samples
READ(60,*)max_number_of_samples
READ(60,*)display_interval
READ(60,*)add_sleep
READ(60,*)no_data
READ(60,*)
READ(60,*)
READ(60,*)
READ(60,*)parallel_tempering
IF (parallel_tempering.EQ.1) THEN
        CALL read_tempering_parameters() 
END IF
IF (input_error) section='(I)'
CLOSE(60)

IF (save_residuals.NE.0.AND.save_residuals.NE.1) THEN
        input_error=.TRUE.
        section='(B)'
END IF
IF (likelihood_type.NE.0.AND.likelihood_type.NE.1) THEN
        input_error=.TRUE.
        section='(F)'
END IF
IF (max_runtime.LE.0) THEN
        input_error=.TRUE.
        section='(H)'
END IF
IF (number_of_samples.LE.0.OR.max_number_of_samples.LT.0) THEN
        input_error=.TRUE.
        section='(H)'
END IF
IF (display_interval.LE.0) THEN
        input_error=.TRUE.
        section='(H)'
END IF
IF (add_sleep.NE.0.AND.add_sleep.NE.1) THEN
        input_error=.TRUE.
        section='(H)'
END IF
IF (no_data.NE.0.AND.no_data.NE.1) THEN
        input_error=.TRUE.
        section='(H)'
END IF
IF (parallel_tempering.NE.0.AND.parallel_tempering.NE.1) THEN
        input_error=.TRUE.
        section='(I)'
END IF
IF (jump_type.NE.0.AND.jump_type.NE.1.AND.jump_type.NE.2) THEN
        input_error=.TRUE.
        section='(I)'
END IF
IF (tempering_start.LT.0) THEN
        input_error=.TRUE.
        section='(I)'
END IF
IF (input_error) THEN
        WRITE(*,*)
        WRITE(*,*)'********************************************'
        WRITE(*,*)'ERROR in mksamples.in'
        WRITE(*,*)'--------------------------------------------'
        WRITE(*,*)'Check section ',section,' of mksamples.in.'
        WRITE(*,*)'TERMINATING PROGRAM!!'
        WRITE(*,*)'********************************************'
        WRITE(*,*)
END IF

IF (max_number_of_samples.EQ.0) max_number_of_samples=HUGE(INT(0))

END SUBROUTINE


SUBROUTINE read_tempering_parameters()
! Reads the parallel tempering parameters from the input file
USE chains, ONLY: parallel_tempering,jump_type,tempering_start,number_of_1s,number_of_non1s,number_of_temperatures,&
                  number_of_chains,temperature_values0,temperature_indices0,temperature_values,temperature_indices
USE m_mrgrnk
IMPLICIT NONE
INTEGER                                 :: nn
INTEGER, DIMENSION(:), ALLOCATABLE      :: temp_order
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: tempo
READ(60,*)jump_type
READ(60,*)tempering_start
READ(60,*)number_of_1s,number_of_non1s
IF (parallel_tempering.EQ.1) THEN
        number_of_temperatures=number_of_1s+number_of_non1s
        IF (number_of_temperatures.NE.number_of_chains) THEN
                WRITE(*,*)
                WRITE(*,*)'********************************************'
                WRITE(*,*)'ERROR in parallel tempering'
                WRITE(*,*)'--------------------------------------------'
                WRITE(*,*)'The sum of the number of temperatures =1 and'
                WRITE(*,*)'>1 must be equal to the number of parallel'
                WRITE(*,*)'chains!!'
                WRITE(*,*)'Check mksamples.in and mksamples command.'
                WRITE(*,*)'TERMINATING PROGRAM!!'
                WRITE(*,*)'********************************************'
                WRITE(*,*)
                input_error=.TRUE.
                RETURN
        END IF
        ALLOCATE(temperature_values(number_of_temperatures),temperature_values0(number_of_temperatures))
        ALLOCATE(temperature_indices(number_of_temperatures),temperature_indices0(number_of_temperatures))
        ALLOCATE(temp_order(number_of_non1s))
        temperature_values0(1:number_of_1s)=1

        IF (number_of_non1s.GT.2) THEN 
                ALLOCATE(tempo(number_of_non1s))
                READ(60,*)temperature_values0(number_of_1s+1)
                CALL linspace(from=2,to=temperature_values0(number_of_1s+1),&
                        &array=tempo)
                temperature_values0(number_of_1s+1:number_of_temperatures)=tempo
                DEALLOCATE(tempo)
                READ(60,*)
        ELSE
                DO nn=number_of_1s+1,number_of_temperatures
                        READ(60,*)temperature_values0(nn)
                END DO
                IF (number_of_non1s.EQ.1) THEN
                        READ(60,*)
                END IF
        END IF
        
        CALL D_MRGRNK(temperature_values0(number_of_1s+1:number_of_temperatures),temp_order)
        temperature_values0(number_of_1s+1:number_of_temperatures)=temperature_values0(number_of_1s+temp_order)
        DEALLOCATE(temp_order)
        temperature_indices0(1:number_of_1s)=(/(nn,nn=1,number_of_1s)/)
        temperature_indices0(number_of_1s+1:number_of_temperatures)=(/(nn,nn=number_of_1s+1,number_of_temperatures)/)
END IF

END SUBROUTINE


SUBROUTINE check_run_number()
! Checks the format of the run number stated as the first argument of ./mksamples
USE chains, ONLY: run_name,run_number
IMPLICIT NONE
IF (run_name.EQ.'') THEN
        WRITE(*,*)
        WRITE(*,*)'********************************************'
        WRITE(*,*)'ERROR reading run number'
        WRITE(*,*)'--------------------------------------------'
        WRITE(*,*)'The run number must be set from the command'
        WRITE(*,*)'line as the first argument of mksamples!!'
        WRITE(*,*)'Check argument list of ./mksamples'
        WRITE(*,*)'TERMINATING PROGRAM!!'
        WRITE(*,*)'********************************************'
        WRITE(*,*)
        STOP
END IF
READ(run_name,*)run_number
IF (run_number.LT.1) THEN
        WRITE(*,*)
        WRITE(*,*)'********************************************'
        WRITE(*,*)'ERROR reading run number'
        WRITE(*,*)'--------------------------------------------'
        WRITE(*,*)'The run number must be an integer >=1!!'
        WRITE(*,*)'Check argument list of ./mksamples'
        WRITE(*,*)'TERMINATING PROGRAM!!'
        WRITE(*,*)'********************************************'
        WRITE(*,*)
        STOP
END IF
END SUBROUTINE


SUBROUTINE read_procsamples_in()
! Reads procsamples.in
USE chains, ONLY: number_of_datasets,ncell_min,ncell_max,value_min,value_max,x_min,x_max,y_min,y_max,z_min,z_max,&
                  noise_type,a_min,a_max,b_min,b_max,lambda_min,lambda_max,sigma_min,sigma_max,&
                  value_change,value_change_delayed,x_change,y_change,z_change,x_change_delayed,y_change_delayed,z_change_delayed,&
                  birth_change,a_change,b_change,lambda_change,sigma_change,likelihood_type,&
                  max_runtime,number_of_samples,max_number_of_samples,display_interval,no_data,&
                  parallel_tempering,jump_type,tempering_start
IMPLICIT NONE
INTEGER                 :: nn
CHARACTER(LEN=3)        :: section

input_error=.FALSE.
OPEN(UNIT=60,FILE='procsamples.in',STATUS='OLD')
READ(60,*)
READ(60,*)
READ(60,*)
READ(60,'(a30)')samples_file
READ(60,'(a30)')ncell_file
READ(60,'(a30)')noise_file
READ(60,'(a30)')misfits_file
READ(60,'(a30)')acceptance_file
READ(60,'(a30)')temperatures_file
READ(60,'(a30)')residuals_file
READ(60,*)
READ(60,*)
READ(60,*)
READ(60,'(a30)')average_file
READ(60,'(a30)')stdev_file
READ(60,'(a30)')median_file
READ(60,'(a30)')maximum_file
READ(60,'(a30)')harmean_file
READ(60,'(a30)')rms_file
READ(60,'(a30)')entropy_file
READ(60,'(a30)')skewness_file
READ(60,'(a30)')kurtosis_file
READ(60,'(a30)')ncellsi_file
READ(60,'(a30)')misfiti_file
READ(60,'(a30)')postncells_file
READ(60,'(a30)')postnoise_file
READ(60,'(a30)')postvalue_file
READ(60,'(a30)')node_file
READ(60,'(a30)')first_file
READ(60,'(a30)')last_file
READ(60,'(a30)')residuals_histo_file
READ(60,*)
READ(60,*)
READ(60,*)
READ(60,*)vtx_or_xyz
READ(60,*)x_number_of_pixels,y_number_of_pixels,z_number_of_pixels
READ(60,*)x_pixel_size,y_pixel_size,z_pixel_size
READ(60,*)value_number_of_pixels
READ(60,*)noise_number_of_pixels
READ(60,*)residuals_number_of_pixels
READ(60,*)max_abs_residual,max_percentage_residual
max_abs_residual=ABS(max_abs_residual)
max_percentage_residual=ABS(max_percentage_residual)
READ(60,*)
READ(60,*)
READ(60,*)
READ(60,*)burn_in
READ(60,*)thinning
READ(60,*)stop_sampling
READ(60,*)number_of_bad_chains
IF (number_of_bad_chains.GT.0) THEN
        ALLOCATE(bad_chains(number_of_bad_chains))
        DO nn=1,number_of_bad_chains
                READ(60,*)bad_chains(nn)
        END DO
END IF
CLOSE(60)

IF (vtx_or_xyz.NE.0.AND.vtx_or_xyz.NE.1) THEN
        input_error=.TRUE.
        section='(C)'
END IF
IF (x_number_of_pixels.LE.0.OR.y_number_of_pixels.LE.0.OR.z_number_of_pixels.LE.0) THEN
        input_error=.TRUE.
        section='(C)'
END IF
IF (x_pixel_size.LE.0.OR.y_pixel_size.LE.0.OR.z_pixel_size.LE.0) THEN
        input_error=.TRUE.
        section='(C)'
END IF
IF (value_number_of_pixels.LE.0.OR.noise_number_of_pixels.LE.0) THEN
        input_error=.TRUE.
        section='(C)'
END IF
IF (residuals_number_of_pixels.LE.0) THEN
        input_error=.TRUE.
        section='(C)'
END IF
IF (burn_in.LT.0.OR.thinning.LE.0.OR.stop_sampling.EQ.0.OR.stop_sampling.LT.-1) THEN
        input_error=.TRUE.
        section='(D)'
END IF
IF (number_of_bad_chains.LT.0) THEN
        input_error=.TRUE.
        section='(D)'
END IF
IF (input_error) THEN
        WRITE(*,*)
        WRITE(*,*)'********************************************'
        WRITE(*,*)'ERROR in procsamples.in'
        WRITE(*,*)'--------------------------------------------'
        WRITE(*,*)'Check section ',section,' of procsamples.in.'
        WRITE(*,*)'TERMINATING PROGRAM!!'
        WRITE(*,*)'********************************************'
        WRITE(*,*)
END IF

END SUBROUTINE


SUBROUTINE linspace(from, to, array) ! Automatically associate a temperature to each Markov chain
        IMPLICIT NONE
        INTEGER, INTENT(IN)           :: from
        DOUBLE PRECISION, INTENT(IN)  :: to
        DOUBLE PRECISION, INTENT(OUT) :: array(:)
        DOUBLE PRECISION              :: range
        INTEGER                       :: n, i

        n=size(array)
        range=to-from
        IF (n==0) RETURN 
        IF (n==1) THEN
            array(1)=from
            RETURN
        END IF
    
        DO i=1,n
            array(i)=from+range*(i-1)/(n-1)
        END DO
END SUBROUTINE


END MODULE
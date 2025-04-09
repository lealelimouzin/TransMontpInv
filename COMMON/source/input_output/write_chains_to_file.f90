! ********************************************************************* !
! Subroutine to write Markov chains to file in mksamples.f90
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


SUBROUTINE write_chains_to_file()

USE chains
USE input_files, ONLY: samples_file,residuals_file,ncell_file,noise_file,&
                       acceptance_file,misfits_file,temperatures_file,&
                       save_residuals
IMPLICIT NONE
INTEGER                 :: pp,nn,ii
LOGICAL                 :: file_open
CHARACTER(LEN=3)        :: adv
INTEGER                 :: file_type_1=1 ! file of number of cells
INTEGER                 :: file_type_2=2 ! file of chain steps
INTEGER                 :: file_type_3=3 ! file of noise parameters
INTEGER                 :: file_type_4=4 ! file of misfits
INTEGER                 :: file_type_5=5 ! file of acceptance ratios
INTEGER                 :: file_type_6=6 ! file of temperatures
INTEGER                 :: file_type_7=7 ! file of residuals
        
IF (run_number.EQ.1) THEN
        
        OPEN(UNIT=samples_unit,&
             FILE=TRIM(samples_file) // TRIM(chain_name),&
             STATUS='REPLACE',&
             FORM='UNFORMATTED',&
             ACCESS='SEQUENTIAL')
        
        OPEN(UNIT=residuals_unit,&
             FILE=TRIM(residuals_file) // TRIM(chain_name),&
             STATUS='REPLACE',&
             FORM='UNFORMATTED',&
             ACCESS='SEQUENTIAL')
        
        OPEN(UNIT=ncell_unit,&
             FILE=TRIM(ncell_file) // TRIM(chain_name),&
             STATUS='REPLACE',&
             FORM='UNFORMATTED',&
             ACCESS='SEQUENTIAL')
        
        IF (noise_type.NE.4) OPEN(UNIT=noise_unit,&
                                  FILE=TRIM(noise_file) // TRIM(chain_name),&
                                  STATUS='REPLACE',&
                                  FORM='UNFORMATTED',&
                                  ACCESS='SEQUENTIAL')
        IF (noise_type.EQ.4.AND.chain_number.EQ.1) OPEN(UNIT=noise_unit,&
                                                        FILE=TRIM(noise_file),&
                                                        STATUS='REPLACE',&
                                                        FORM='UNFORMATTED',&
                                                        ACCESS='SEQUENTIAL')
        OPEN(UNIT=misfits_unit,&
             FILE=TRIM(misfits_file) // TRIM(chain_name),&
             STATUS='REPLACE',&
             FORM='UNFORMATTED',&
             ACCESS='SEQUENTIAL')
        
        OPEN(UNIT=acceptance_unit,&
             FILE=TRIM(acceptance_file) // TRIM(chain_name),&
             STATUS='REPLACE',&
             FORM='UNFORMATTED',&
             ACCESS='SEQUENTIAL')
        
        IF (parallel_tempering.EQ.1.AND.chain_number.EQ.1) OPEN(UNIT=temperatures_unit,&
                                                                FILE=TRIM(temperatures_file),&
                                                                STATUS='REPLACE',&
                                                                FORM='UNFORMATTED',&
                                                                ACCESS='SEQUENTIAL')
        
        DO pp=0,samples_in_run
                
                ! Number of cells
                IF (pp.EQ.0) WRITE(ncell_unit)file_type_1
                WRITE(ncell_unit)pp,chain_ncell(pp+1)
                
                ! Chain steps
                IF (pp.EQ.0) THEN
                        WRITE(samples_unit)file_type_2
                        WRITE(samples_unit)chain_ncell(pp+1)
                        DO nn=1,chain_ncell(pp+1)
                                WRITE(samples_unit)voronoi_0(nn,1),&
                                                   voronoi_0(nn,2),&
                                                   voronoi_0(nn,3),&
                                                   voronoi_0(nn,4)
                        END DO
                ELSE
                        WRITE(samples_unit)pp,chain_steps(pp+1,1),&
                                              chain_steps(pp+1,2),&
                                              chain_indices(pp+1),&
                                              chain_values(pp+1,1),&
                                              chain_values(pp+1,2),&
                                              chain_values(pp+1,3),&
                                              chain_values(pp+1,4)
                END IF
                
                ! Residuals
                IF (pp.EQ.0) THEN
                        WRITE(residuals_unit)file_type_7
                        WRITE(residuals_unit)save_residuals
                END IF
                IF (save_residuals.EQ.1) THEN
                        IF (pp.EQ.0) THEN
                                DO nn=1,number_of_measurements
                                        WRITE(residuals_unit)residuals(nn,pp+1)
                                END DO
                        ELSE
                                IF (chain_steps(pp+1,2).EQ.1.AND.chain_steps(pp+1,1).NE.5) THEN
                                        DO nn=1,number_of_measurements
                                                WRITE(residuals_unit)residuals(nn,pp+1)
                                        END DO
                                END IF
                        END IF
                END IF
                
                ! Noise parameters
                IF (pp.EQ.0) WRITE(noise_unit)file_type_3
                IF (pp.EQ.0) WRITE(noise_unit)noise_type
                IF (noise_type.NE.4) THEN
                        IF (pp.EQ.0) WRITE(noise_unit)number_of_datasets
                        DO nn=1,number_of_datasets
                                IF (noise_type.EQ.0.OR.&
                                    noise_type.EQ.3.OR.noise_type.EQ.6) WRITE(noise_unit)pp,&
                                                                      chain_noise(pp+1,nn)
                                IF (noise_type.EQ.1.OR.&
                                    noise_type.EQ.2.OR.&
                                    noise_type.EQ.5) WRITE(noise_unit)pp,&
                                                                      chain_noise(pp+1,nn),&
                                                                      chain_noise(pp+1,nn+number_of_datasets)
                        END DO
                ELSE IF (noise_type.EQ.4.AND.chain_number.EQ.1.AND.pp.EQ.0) THEN
                        WRITE(noise_unit)number_of_measurements
                        DO nn=1,number_of_measurements
                                WRITE(noise_unit)noise_current(nn)
                        END DO
                END IF
                
                ! Misfit
                IF (pp.EQ.0) WRITE(misfits_unit)file_type_4
                WRITE(misfits_unit)pp,chain_misfits(pp+1)
                
                ! Acceptance ratios
                IF (pp.EQ.0) WRITE(acceptance_unit)file_type_5
                WRITE(acceptance_unit)pp,chain_acceptance(pp+1,1),&
                                         chain_acceptance(pp+1,2),&
                                         chain_acceptance(pp+1,3),&
                                         chain_acceptance(pp+1,4),&
                                         chain_acceptance(pp+1,5),&
                                         chain_acceptance(pp+1,6)
                
                ! Temperatures
                IF (parallel_tempering.EQ.1.AND.chain_number.EQ.1) THEN
                        IF (pp.EQ.0) THEN
                                WRITE(temperatures_unit)file_type_6
                                WRITE(temperatures_unit)number_of_temperatures
                                DO nn=1,number_of_temperatures
                                        WRITE(temperatures_unit)temperature_values0(nn)
                                END DO
                        ELSE
                                IF (pp.GE.tempering_start) WRITE(temperatures_unit)pp,&
                                                                                   chain_temperatures(pp,1),&
                                                                                   chain_temperatures(pp,2),&
                                                                                   chain_temperatures(pp,3),&
                                                                                   chain_acceptance(pp+1,7)
                        END IF
                END IF
                
        END DO
        
ELSE IF (run_number.GT.1) THEN
        
        OPEN(UNIT=samples_unit,&
             FILE=TRIM(samples_file) // TRIM(chain_name),&
             STATUS='OLD',&
             FORM='UNFORMATTED',&
             ACCESS='SEQUENTIAL',&
             POSITION='APPEND')
        
        OPEN(UNIT=residuals_unit,&
             FILE=TRIM(residuals_file) // TRIM(chain_name),&
             STATUS='OLD',&
             FORM='UNFORMATTED',&
             ACCESS='SEQUENTIAL',&
             POSITION='APPEND')
        
        OPEN(UNIT=ncell_unit,&
             FILE=TRIM(ncell_file) // TRIM(chain_name),&
             STATUS='OLD',&
             FORM='UNFORMATTED',&
             ACCESS='SEQUENTIAL',&
             POSITION='APPEND')
        
        IF (noise_type.NE.4) OPEN(UNIT=noise_unit,&
                                  FILE=TRIM(noise_file) // TRIM(chain_name),&
                                  STATUS='OLD',&
                                  FORM='UNFORMATTED',&
                                  ACCESS='SEQUENTIAL',&
                                  POSITION='APPEND')
        
        OPEN(UNIT=misfits_unit,&
             FILE=TRIM(misfits_file) // TRIM(chain_name),&
             STATUS='OLD',&
             FORM='UNFORMATTED',&
             ACCESS='SEQUENTIAL',&
             POSITION='APPEND')
        
        OPEN(UNIT=acceptance_unit,&
             FILE=TRIM(acceptance_file) // TRIM(chain_name),&
             STATUS='OLD',&
             FORM='UNFORMATTED',&
             ACCESS='SEQUENTIAL',&
             POSITION='APPEND')
        
        IF (parallel_tempering.EQ.1.AND.chain_number.EQ.1) OPEN(UNIT=temperatures_unit,&
                                                                FILE=TRIM(temperatures_file),&
                                                                STATUS='OLD',&
                                                                FORM='UNFORMATTED',&
                                                                ACCESS='SEQUENTIAL',&
                                                                POSITION='APPEND')
        
        DO pp=1,samples_in_run
                
                ! Number of cells
                WRITE(ncell_unit)samples_at_start+pp,chain_ncell(pp)
                
                ! Chain steps
                WRITE(samples_unit)samples_at_start+pp,chain_steps(pp,1),&
                                                       chain_steps(pp,2),&
                                                       chain_indices(pp),&
                                                       chain_values(pp,1),&
                                                       chain_values(pp,2),&
                                                       chain_values(pp,3),&
                                                       chain_values(pp,4)
                
                ! Residuals
                IF (save_residuals.EQ.1) THEN
                        IF (chain_steps(pp,2).EQ.1.AND.chain_steps(pp,1).NE.5) THEN
                                DO nn=1,number_of_measurements
                                        WRITE(residuals_unit)residuals(nn,pp)
                                END DO
                        END IF
                END IF
                
                ! Noise parameters
                IF (noise_type.NE.4) THEN
                        DO nn=1,number_of_datasets
                                IF (noise_type.EQ.0.OR.&
                                    noise_type.EQ.3.OR.noise_type.EQ.6) WRITE(noise_unit)samples_at_start+pp,&
                                                                      chain_noise(pp,nn)
                                IF (noise_type.EQ.1.OR.&
                                    noise_type.EQ.2.OR.&
                                    noise_type.EQ.5) WRITE(noise_unit)samples_at_start+pp,&
                                                                      chain_noise(pp,nn),&
                                                                      chain_noise(pp,nn+number_of_datasets)
                        END DO
                END IF
                
                ! Misfit
                WRITE(misfits_unit)samples_at_start+pp,chain_misfits(pp)
                
                ! Acceptance ratios
                WRITE(acceptance_unit)samples_at_start+pp,chain_acceptance(pp,1),&
                                                          chain_acceptance(pp,2),&
                                                          chain_acceptance(pp,3),&
                                                          chain_acceptance(pp,4),&
                                                          chain_acceptance(pp,5),&
                                                          chain_acceptance(pp,6)
                
                ! Temperatures
                IF (parallel_tempering.EQ.1.AND.chain_number.EQ.1) THEN
                        IF (samples_at_start+pp.GE.tempering_start) WRITE(temperatures_unit)samples_at_start+pp,&
                                                                                            chain_temperatures(pp,1),&
                                                                                            chain_temperatures(pp,2),&
                                                                                            chain_temperatures(pp,3),&
                                                                                            chain_acceptance(pp,7)
                END IF
                
        END DO
        
END IF

CLOSE(samples_unit)
CLOSE(residuals_unit)
CLOSE(ncell_unit)
INQUIRE(noise_unit,OPENED=file_open)
IF (file_open) CLOSE(noise_unit)
CLOSE(misfits_unit)
CLOSE(acceptance_unit)
IF (parallel_tempering.EQ.1.AND.chain_number.EQ.1) CLOSE(temperatures_unit)

END SUBROUTINE

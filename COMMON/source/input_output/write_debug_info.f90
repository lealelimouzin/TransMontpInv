! ********************************************************************* !
! Subroutines to write mksamples/procsamples debug information
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


SUBROUTINE write_debug_info_mk(info,iteration,random_deviate,accepted,delayed)

USE chains
USE string_conversion
USE debug
IMPLICIT NONE
INTEGER, INTENT(IN)                     :: info
INTEGER, INTENT(IN), OPTIONAL           :: iteration
DOUBLE PRECISION, INTENT(IN), OPTIONAL  :: random_deviate
LOGICAL, INTENT(IN), OPTIONAL           :: accepted
LOGICAL, INTENT(IN), OPTIONAL           :: delayed
INTEGER                                 :: nn

IF (info.LT.0) THEN
        IF (info.EQ.-1) THEN
                WRITE(debug_unit,*)
                WRITE(debug_unit,*)'********* START ITERATION ',TRIM(int2str(iteration)),' *********'
        ELSE IF (info.EQ.-2) THEN
                WRITE(debug_unit,*)'********** END ITERATION ',TRIM(int2str(iteration)),' **********'
        END IF
        RETURN
END IF

IF (.NOT.PRESENT(delayed).OR.delayed) WRITE(debug_unit,*)'--------- START DEBUG INFO ---------'
IF (info.EQ.1) THEN ! write temperatures
        
        IF (parallel_tempering.EQ.1) THEN
                WRITE(debug_unit,*)'INITIAL TEMPERATURES FROM CHAIN ',TRIM(int2str(chain_number))
                DO nn=1,number_of_temperatures
                        WRITE(debug_unit,*)'Chain ',TRIM(int2str(nn)),': i=',TRIM(int2str(temperature_indices(nn))),&
                                                                      ', T=',TRIM(dble2str(temperature_values(nn)))
                END DO
        ELSE IF (parallel_tempering.EQ.0) THEN
                WRITE(debug_unit,*)'PARALLEL TEMPERING NOT IN USE'
        END IF
        
ELSE IF (info.EQ.2) THEN ! write initial Voronoi model
        
        WRITE(debug_unit,*)'INITIAL MODEL ON CHAIN ',TRIM(chain_name)
        WRITE(debug_unit,*)'Noise type: ',TRIM(int2str(noise_type))
        IF (noise_type.EQ.0) THEN
                DO nn=1,number_of_datasets
                        WRITE(debug_unit,*)'sigma ',TRIM(int2str(nn)),': ',TRIM(dble2str(sigma_current(nn)))
                END DO
        ELSE IF (noise_type.EQ.3.OR.noise_type.EQ.6) THEN
                DO nn=1,number_of_datasets
                        WRITE(debug_unit,*)'lambda ',TRIM(int2str(nn)),': ',TRIM(dble2str(lambda_current(nn)))
                END DO
        ELSE IF (noise_type.EQ.4) THEN
                DO nn=1,number_of_measurements
                        WRITE(debug_unit,*)'sigma ',TRIM(int2str(nn)),': ',TRIM(dble2str(noise_current(nn)))
                END DO
        ELSE IF (noise_type.EQ.1.OR.noise_type.EQ.2.OR.noise_type.EQ.5) THEN
                DO nn=1,number_of_datasets
                        WRITE(debug_unit,*)'a',TRIM(int2str(nn)),': ',TRIM(dble2str(a_current(nn))),&
                                  ', b',TRIM(int2str(nn)),': ',TRIM(dble2str(b_current(nn)))
                END DO
        END IF
        WRITE(debug_unit,*)'Number of cells: ',TRIM(int2str(ncell_current))
        WRITE(debug_unit,*)'Geometry:'
        DO nn=1,ncell_current
                WRITE(debug_unit,*)TRIM(int2str(nn)),' ',TRIM(dble2str(voronoi_current(nn,1))),&
                                            ' ',TRIM(dble2str(voronoi_current(nn,2))),&
                                            ' ',TRIM(dble2str(voronoi_current(nn,3))),&
                                            ' ',TRIM(dble2str(voronoi_current(nn,4)))
        END DO
        
ELSE IF (info.EQ.3) THEN ! write likelihoods
        
        WRITE(debug_unit,*)'LIKELIHOODS ON CHAIN ',TRIM(chain_name)
        WRITE(debug_unit,*)'Iteration: ',TRIM(int2str(iteration))
        WRITE(debug_unit,*)'Current: ',TRIM(dble2str(likelihood_current))
        WRITE(debug_unit,*)'Proposed: ',TRIM(dble2str(likelihood_proposed))
        WRITE(debug_unit,*)'Delayed: ',TRIM(dble2str(likelihood_delayed))
        
ELSE IF (info.EQ.4) THEN ! write steps - main
        
        WRITE(debug_unit,*)'STEP TYPE ON CHAIN ',TRIM(chain_name)
        WRITE(debug_unit,*)'Iteration: ',TRIM(int2str(iteration))
        IF (birth) THEN
                WRITE(debug_unit,*)TRIM(int2str(step)),': birth'
                WRITE(debug_unit,*)TRIM(int2str(ncell_proposed)),&
                      ' ',TRIM(dble2str(voronoi_proposed(ncell_proposed,1))),&
                      ' ',TRIM(dble2str(voronoi_proposed(ncell_proposed,2))),&
                      ' ',TRIM(dble2str(voronoi_proposed(ncell_proposed,3))),&
                      ' ',TRIM(dble2str(voronoi_proposed(ncell_proposed,4)))
        ELSE IF (death) THEN
                WRITE(debug_unit,*)TRIM(int2str(step)),': death'
                WRITE(debug_unit,*)TRIM(int2str(cell_number))
        ELSE IF (move.OR.value) THEN
                IF (move) WRITE(debug_unit,*)TRIM(int2str(step)),': move'
                IF (value) WRITE(debug_unit,*)TRIM(int2str(step)),': value'
                WRITE(debug_unit,*)TRIM(int2str(cell_number)),&
                      ' ',TRIM(dble2str(voronoi_proposed(cell_number,1))),&
                      ' ',TRIM(dble2str(voronoi_proposed(cell_number,2))),&
                      ' ',TRIM(dble2str(voronoi_proposed(cell_number,3))),&
                      ' ',TRIM(dble2str(voronoi_proposed(cell_number,4)))
        ELSE IF (noise) THEN
                WRITE(debug_unit,*)TRIM(int2str(step)),': noise'
                IF (noise_type.EQ.0) THEN
                        DO nn=1,number_of_datasets
                                WRITE(debug_unit,*)'sigma ',TRIM(int2str(nn)),': ',TRIM(dble2str(sigma_proposed(nn)))
                        END DO
                ELSE IF (noise_type.EQ.3.OR.noise_type.EQ.6) THEN
                        DO nn=1,number_of_datasets
                                WRITE(debug_unit,*)'lambda ',TRIM(int2str(nn)),': ',TRIM(dble2str(lambda_proposed(nn)))
                        END DO
                ELSE IF (noise_type.EQ.4) THEN
                        DO nn=1,number_of_measurements
                                WRITE(debug_unit,*)'sigma ',TRIM(int2str(nn)),': ',TRIM(dble2str(noise_proposed(nn)))
                        END DO
                ELSE IF (noise_type.EQ.1.OR.noise_type.EQ.2.OR.noise_type.EQ.5) THEN
                        DO nn=1,number_of_datasets
                                WRITE(debug_unit,*)'a',TRIM(int2str(nn)),': ',TRIM(dble2str(a_proposed(nn))),&
                                          ', b',TRIM(int2str(nn)),': ',TRIM(dble2str(b_proposed(nn)))
                        END DO
                END IF
        END IF
        IF (inside_prior) WRITE(debug_unit,*)'inside prior: YES'
        IF (.NOT.inside_prior) WRITE(debug_unit,*)'inside prior: NO'
        IF (inside_prior) THEN
                WRITE(debug_unit,*)'log(r)=',TRIM(dble2str(LOG(random_deviate)))
                WRITE(debug_unit,*)'alpha=',TRIM(dble2str(alpha))
                IF (LOG(random_deviate).LE.alpha) WRITE(debug_unit,*)'log(r)<alpha: YES'
                IF (.NOT.(LOG(random_deviate).LE.alpha)) WRITE(debug_unit,*)'log(r)<alpha: NO'
        END IF
        IF (accepted) WRITE(debug_unit,*)'accepted: YES'
        IF (.NOT.accepted) WRITE(debug_unit,*)'accepted: NO'
        
ELSE IF (info.EQ.5) THEN ! write steps - main
        
        IF (delayed) THEN
                
                WRITE(debug_unit,*)'STEP TYPE ON CHAIN ',TRIM(chain_name)
                WRITE(debug_unit,*)'Iteration: ',TRIM(int2str(iteration)),' delayed'
                IF (move) WRITE(debug_unit,*)TRIM(int2str(step)),': move'
                IF (value) WRITE(debug_unit,*)TRIM(int2str(step)),': value'
                WRITE(debug_unit,*)TRIM(int2str(cell_number)),&
                      ' ',TRIM(dble2str(voronoi_delayed(cell_number,1))),&
                      ' ',TRIM(dble2str(voronoi_delayed(cell_number,2))),&
                      ' ',TRIM(dble2str(voronoi_delayed(cell_number,3))),&
                      ' ',TRIM(dble2str(voronoi_delayed(cell_number,4)))
                IF (inside_prior_delayed) WRITE(debug_unit,*)'inside prior delayed: YES'
                IF (.NOT.inside_prior_delayed) WRITE(debug_unit,*)'inside prior delayed: NO'
                IF (inside_prior_delayed) THEN
                        WRITE(debug_unit,*)'log(r)=',TRIM(dble2str(LOG(random_deviate)))
                        WRITE(debug_unit,*)'alpha_delayed=',TRIM(dble2str(alpha_delayed))
                        IF (LOG(random_deviate).LE.alpha_delayed) WRITE(debug_unit,*)'log(r)<alpha_delayed: YES'
                        IF (.NOT.(LOG(random_deviate).LE.alpha_delayed)) WRITE(debug_unit,*)'log(r)<alpha_delayed: NO'
                END IF
                IF (accepted) WRITE(debug_unit,*)'accepted: YES'
                IF (.NOT.accepted) WRITE(debug_unit,*)'accepted: NO'
                
        ELSE
                
                RETURN
                
        END IF
        
ELSE IF (info.EQ.6) THEN ! temperature information
        
        WRITE(debug_unit,*)'TEMPERATURES FROM CHAIN ',TRIM(int2str(chain_number)),' (used in next iteration)'
        WRITE(debug_unit,*)'Iteration: ',TRIM(int2str(iteration))
        DO nn=1,number_of_temperatures
                WRITE(debug_unit,*)'Chain ',TRIM(int2str(nn)),': i=',TRIM(int2str(temperature_indices(nn))),&
                                                              ', T=',TRIM(dble2str(temperature_values(nn)))
        END DO
        
END IF
IF (.NOT.PRESENT(delayed).OR.delayed) WRITE(debug_unit,*)'---------- END DEBUG INFO ----------'

END SUBROUTINE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE write_debug_info_pr(info,iteration)

USE chains
USE debug
USE string_conversion
USE processing, ONLY: temperatures,accept,x_new,y_new,z_new,value_new
IMPLICIT NONE
INTEGER, INTENT(IN)                     :: info
INTEGER, INTENT(IN), OPTIONAL           :: iteration
INTEGER                                 :: nn,tt
LOGICAL                                 :: accepted

IF (info.EQ.1) THEN ! write temperatures
        
        WRITE(debug_unit,*)'------------------------------------'
        IF (parallel_tempering.EQ.1) THEN
                WRITE(debug_unit,*)'INVERSION TEMPERATURES:'
                DO nn=1,samples_so_far
                        WRITE(debug_unit,*)'********* Iteration ',TRIM(int2str(nn)),' *********'
                        DO tt=1,number_of_temperatures
                                WRITE(debug_unit,*)'   - Chain',TRIM(int2str(tt)),': T=',temperatures(nn,tt)
                        END DO
                END DO
                WRITE(debug_unit,*)'**************************************'
        ELSE IF (parallel_tempering.EQ.0) THEN
                WRITE(debug_unit,*)'PARALLEL TEMPERING NOT IN USE'
        END IF
        
ELSE IF (info.EQ.2) THEN ! write initial Voronoi model
        
        WRITE(debug_unit,*)'------------------------------------'
        WRITE(debug_unit,*)'INITIAL MODEL ON CHAIN ',TRIM(chain_name)
        WRITE(debug_unit,*)'Noise type: ',TRIM(int2str(noise_type))
        IF (noise_type.EQ.0) THEN
                DO nn=1,number_of_datasets
                        WRITE(debug_unit,*)'sigma ',TRIM(int2str(nn)),': ',TRIM(dble2str(sigma_current(nn)))
                END DO
        ELSE IF (noise_type.EQ.3.OR.noise_type.EQ.6) THEN
                DO nn=1,number_of_datasets
                        WRITE(debug_unit,*)'lambda ',TRIM(int2str(nn)),': ',TRIM(dble2str(lambda_current(nn)))
                END DO
        ELSE IF (noise_type.EQ.4) THEN
                DO nn=1,number_of_measurements
                        WRITE(debug_unit,*)'sigma ',TRIM(int2str(nn)),': ',TRIM(dble2str(noise_current(nn)))
                END DO
        ELSE IF (noise_type.EQ.1.OR.noise_type.EQ.2.OR.noise_type.EQ.5) THEN
                DO nn=1,number_of_datasets
                        WRITE(debug_unit,*)'a',TRIM(int2str(nn)),': ',TRIM(dble2str(a_current(nn))),&
                                  ', b',TRIM(int2str(nn)),': ',TRIM(dble2str(b_current(nn)))
                END DO
        END IF
        WRITE(debug_unit,*)'Number of cells: ',TRIM(int2str(ncell_current))
        WRITE(debug_unit,*)'Geometry:'
        DO nn=1,ncell_current
                WRITE(debug_unit,*)TRIM(int2str(nn)),' ',TRIM(dble2str(voronoi_current(nn,1))),&
                                            ' ',TRIM(dble2str(voronoi_current(nn,2))),&
                                            ' ',TRIM(dble2str(voronoi_current(nn,3))),&
                                            ' ',TRIM(dble2str(voronoi_current(nn,4)))
        END DO
        WRITE(debug_unit,*)'------------------------------------'
        
ELSE IF (info.EQ.4) THEN ! write steps
        
        IF (accept.EQ.1) accepted=.TRUE.
        IF (accept.EQ.0) accepted=.FALSE.
        CALL step_type()
        
        WRITE(debug_unit,*)'********* START ITERATION ',TRIM(int2str(iteration)),' *********'
        WRITE(debug_unit,*)'STEP TYPE ON CHAIN ',TRIM(chain_name)
        IF (birth) WRITE(debug_unit,*)TRIM(int2str(step)),': birth'
        IF (death) WRITE(debug_unit,*)TRIM(int2str(step)),': death'
        IF (move) WRITE(debug_unit,*)TRIM(int2str(step)),': move'
        IF (value) WRITE(debug_unit,*)TRIM(int2str(step)),': value'
        IF (noise) WRITE(debug_unit,*)TRIM(int2str(step)),': noise'
        IF (accepted) WRITE(debug_unit,*)'accepted: YES'
        IF (.NOT.accepted) WRITE(debug_unit,*)'accepted: NO'
        IF (accepted) THEN
                IF (birth) THEN
                        WRITE(debug_unit,*)TRIM(int2str(ncell_current)),&
                                       ' ',TRIM(dble2str(x_new)),&
                                       ' ',TRIM(dble2str(y_new)),&
                                       ' ',TRIM(dble2str(z_new)),&
                                       ' ',TRIM(dble2str(value_new))
                ELSE IF (death) THEN
                        WRITE(debug_unit,*)TRIM(int2str(cell_number))
                ELSE IF (move) THEN
                        WRITE(debug_unit,*)TRIM(int2str(cell_number)),&
                                       ' ',TRIM(dble2str(x_new)),&
                                       ' ',TRIM(dble2str(y_new)),&
                                       ' ',TRIM(dble2str(z_new)),&
                                       ' ',TRIM(dble2str(voronoi_current(cell_number,4)))
                ELSE IF (value) THEN
                        WRITE(debug_unit,*)TRIM(int2str(cell_number)),&
                                       ' ',TRIM(dble2str(voronoi_current(cell_number,1))),&
                                       ' ',TRIM(dble2str(voronoi_current(cell_number,2))),&
                                       ' ',TRIM(dble2str(voronoi_current(cell_number,3))),&
                                       ' ',TRIM(dble2str(value_new))
                ELSE IF (noise) THEN
                        IF (noise_type.EQ.0) THEN
                                DO nn=1,number_of_datasets
                                        WRITE(debug_unit,*)'sigma ',TRIM(int2str(nn)),': ',TRIM(dble2str(sigma_current(nn)))
                                END DO
                        ELSE IF (noise_type.EQ.3.OR.noise_type.EQ.6) THEN
                                DO nn=1,number_of_datasets
                                        WRITE(debug_unit,*)'lambda ',TRIM(int2str(nn)),': ',TRIM(dble2str(lambda_current(nn)))
                                END DO
                        ELSE IF (noise_type.EQ.4) THEN
                                DO nn=1,number_of_measurements
                                        WRITE(debug_unit,*)'sigma ',TRIM(int2str(nn)),': ',TRIM(dble2str(noise_current(nn)))
                                END DO
                        ELSE IF (noise_type.EQ.1.OR.noise_type.EQ.2.OR.noise_type.EQ.5) THEN
                                DO nn=1,number_of_datasets
                                        WRITE(debug_unit,*)'a',TRIM(int2str(nn)),': ',TRIM(dble2str(a_current(nn))),&
                                                         ', b',TRIM(int2str(nn)),': ',TRIM(dble2str(b_current(nn)))
                                END DO
                        END IF
                END IF
        END IF
        WRITE(debug_unit,*)'********** END ITERATION ',TRIM(int2str(iteration)),' **********'
        
END IF

END SUBROUTINE
! ********************************************************************* !
! Module for noise calculation in mksamples.f90
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


MODULE noise_calculation

USE acquisition, ONLY: observed_data,data_noise,dataset_associations,&
                       number_of_measurements,number_of_datasets,&
                       std_val
USE chains, ONLY: noise_type,a_current,a_proposed,b_current,b_proposed,&
                  lambda_current,lambda_proposed,sigma_current,sigma_proposed,&
                  noise_current,noise_proposed,noise_delayed,&
                  birth,death,move,value,noise,chain_name,sigma_mat
USE debug

CONTAINS

SUBROUTINE calculate_noise(current_proposed_delayed,modelled_data)

IMPLICIT NONE
INTEGER, INTENT(IN)                             :: current_proposed_delayed
DOUBLE PRECISION, DIMENSION(:), INTENT(IN)      :: modelled_data

IF (debug_mode) WRITE(debug_unit,*)'--- Chain ',TRIM(chain_name),': calculating noise ---'

IF (current_proposed_delayed.EQ.0) THEN
        
        IF (noise_type.EQ.0) THEN
                CALL noise_sigma(sigma_current,noise_current)
        ELSE IF (noise_type.EQ.1) THEN
                CALL noise_ab(a_current,b_current,modelled_data,noise_current)
        ELSE IF (noise_type.EQ.2) THEN
                CALL noise_ab(a_current,b_current,observed_data,noise_current)
        ELSE IF (noise_type.EQ.3) THEN
                CALL noise_lambda(lambda_current,data_noise,noise_current)
        ELSE IF (noise_type.EQ.4) THEN
                noise_current=data_noise
        ELSE IF (noise_type.EQ.6) THEN
                CALL noise_frac(lambda_current,data_noise,sigma_mat,noise_current)
        END IF
        
ELSE IF (current_proposed_delayed.EQ.1) THEN
        
        IF (noise) THEN
                IF (noise_type.EQ.0) THEN
                        CALL noise_sigma(sigma_proposed,noise_proposed)
                ELSE IF (noise_type.EQ.1) THEN
                        CALL noise_ab(a_proposed,b_proposed,modelled_data,noise_proposed)
                ELSE IF (noise_type.EQ.2) THEN
                        CALL noise_ab(a_proposed,b_proposed,observed_data,noise_proposed)
                ELSE IF (noise_type.EQ.3) THEN
                        CALL noise_lambda(lambda_proposed,data_noise,noise_proposed)
                ELSE IF (noise_type.EQ.6) THEN
                        CALL noise_frac(lambda_proposed,data_noise,sigma_mat,noise_proposed)
                END IF
        ELSE
                IF (noise_type.EQ.1) THEN
                        CALL noise_ab(a_current,b_current,modelled_data,noise_proposed)
                ELSE
                        noise_proposed=noise_current
                END IF
        END IF
        
ELSE IF (current_proposed_delayed.EQ.2) THEN ! delayed rejection only occurs for "move" and "value" steps
        
        IF (noise_type.EQ.1) THEN
                CALL noise_ab(a_current,b_current,modelled_data,noise_delayed)
        ELSE
                noise_delayed=noise_current
        END IF
        
END IF

IF (debug_mode) WRITE(debug_unit,*)'--- Chain ',TRIM(chain_name),': calculated noise ---'

END SUBROUTINE

! Noise type = 0
SUBROUTINE noise_sigma(sigma,estimated_noise) ! estimated_noise = noise_proposed
IMPLICIT NONE
DOUBLE PRECISION, DIMENSION(number_of_datasets), INTENT(IN)             :: sigma
DOUBLE PRECISION, DIMENSION(number_of_measurements), INTENT(OUT)        :: estimated_noise
INTEGER                                                                 :: ii,jj
DO ii=1,number_of_measurements
        jj=dataset_associations(ii)
        estimated_noise(ii)=sigma(jj)
END DO
END SUBROUTINE

! Noise type = 1,2
SUBROUTINE noise_ab(a,b,data,estimated_noise)
IMPLICIT NONE
DOUBLE PRECISION, DIMENSION(number_of_datasets), INTENT(IN)             :: a,b
DOUBLE PRECISION, DIMENSION(number_of_measurements), INTENT(IN)         :: data
DOUBLE PRECISION, DIMENSION(number_of_measurements), INTENT(OUT)        :: estimated_noise
INTEGER                                                                 :: ii,jj
DO ii=1,number_of_measurements
        jj=dataset_associations(ii)
        estimated_noise(ii)=SQRT(a(jj)**2*data(ii)**2+b(jj)**2)
END DO
END SUBROUTINE

! Noise type = 3
SUBROUTINE noise_lambda(lambda,data_uncertainty,estimated_noise) ! estimated_noise = noise_proposed
IMPLICIT NONE
DOUBLE PRECISION, DIMENSION(number_of_datasets), INTENT(IN)             :: lambda
DOUBLE PRECISION, DIMENSION(number_of_measurements), INTENT(IN)         :: data_uncertainty
DOUBLE PRECISION, DIMENSION(number_of_measurements), INTENT(OUT)        :: estimated_noise
INTEGER                                                                 :: ii,jj
DO ii=1,number_of_measurements
        jj=dataset_associations(ii)
        estimated_noise(ii)=lambda(jj)*data_uncertainty(ii)
END DO
END SUBROUTINE

! Noise type = 6
SUBROUTINE noise_frac(lambda,data_uncertainty,sigma_mat,estimated_noise) ! estimated_noise = noise_proposed
IMPLICIT NONE
DOUBLE PRECISION, DIMENSION(number_of_datasets), INTENT(IN)             :: lambda
DOUBLE PRECISION, DIMENSION(number_of_measurements), INTENT(IN)         :: data_uncertainty
DOUBLE PRECISION, INTENT(IN)                                            :: sigma_mat
DOUBLE PRECISION, DIMENSION(number_of_measurements), INTENT(OUT)        :: estimated_noise
INTEGER                                                                 :: ii,jj
DO ii=1,number_of_measurements
        jj=dataset_associations(ii)
        !estimated_noise(ii)=lambda(jj)*data_uncertainty(ii)
        !estimated_noise(ii)=lambda(jj)*(data_uncertainty(ii)+std_val(ii)) ! test 18/07/25
        estimated_noise(ii)=lambda(jj)*(data_uncertainty(ii)+sigma_mat)
END DO
END SUBROUTINE

END MODULE
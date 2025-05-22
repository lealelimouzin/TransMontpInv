! ********************************************************************* !
! Module to apply parallel tempering in mksamples.f90
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


MODULE tempering

USE chains, ONLY: parallel_tempering,jump_type,tempering_start,number_of_temperatures,&
                  number_of_1s,number_of_non1s,temperature_values,temperature_indices,&
                  temperature_values0,temperature_indices0,temperature,&
                  likelihood_T1,likelihood_T2,noise_T1,noise_T2,swap_temperatures,&
                  acceptance_tempering,acceptance_all,chain_temperatures,&
                  number_of_chains,number_of_measurements,number_of_samples,chain_number,&
                  chain_acceptance,run_number
USE random_generator
USE debug
USE string_conversion
IMPLICIT NONE
INTEGER                 :: T1_chain,T1_index
DOUBLE PRECISION        :: T1_value
INTEGER                 :: T2_chain,T2_index
DOUBLE PRECISION        :: T2_value


CONTAINS


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE initialize_tempering()

IMPLICIT NONE

IF (.NOT.(ALLOCATED(noise_T1))) ALLOCATE(noise_T1(number_of_measurements))
IF (.NOT.(ALLOCATED(noise_T2))) ALLOCATE(noise_T2(number_of_measurements))
ALLOCATE(chain_temperatures(number_of_samples,3))

temperature=temperature_values(chain_number)

END SUBROUTINE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE choose_T1()
! The first temperature is chosen randomly
IMPLICIT NONE
INTEGER         :: ii

T1_chain=0
DO WHILE (.NOT.(ANY((/(ii,ii=1,number_of_chains)/).EQ.T1_chain)))
        T1_chain=NINT(unirand()*(number_of_chains+1))
END DO
T1_value=temperature_values(T1_chain)
T1_index=temperature_indices(T1_chain)

write(*,*) 'Tempering : T1 is ', TRIM(int2str(T1_chain)), ' with a temperature of ', TRIM(dble2str(T1_value))

IF (debug_mode) WRITE(debug_unit,*)'----- T1: chain=',TRIM(int2str(T1_chain)),&
                          ', index=',TRIM(int2str(T1_index)),&
                          ', value=',TRIM(dble2str(T1_value)),' -----'

END SUBROUTINE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE choose_T2()
! Choose the second temperature depending on the choice made in mksamples.in
IMPLICIT NONE
INTEGER         :: ii

T2_chain=0
IF (jump_type.EQ.0) THEN ! choose a random temperature different from T1
        DO WHILE (.NOT.(ANY((/(ii,ii=1,number_of_chains)/).EQ.T2_chain)).OR.T1_value.EQ.T2_value)
                T2_chain=NINT(unirand()*(number_of_chains+1))
                T2_value=temperature_values(T2_chain)
                T2_index=temperature_indices(T2_chain)
        END DO
ELSE IF (jump_type.EQ.1) THEN ! choose a temperature randomly between the two neighbours of T1
        IF (T1_value.EQ.1) THEN ! if T1=1, choose the next temperature above it
                T2_value=temperature_values0(number_of_1s+1)
                T2_index=temperature_indices0(number_of_1s+1)
        ELSE IF (T1_value.EQ.temperature_values0(number_of_temperatures)) THEN ! if T1 is the maximum temperature, choose the first below it
                T2_value=temperature_values0(number_of_temperatures-1)
                T2_index=temperature_indices0(number_of_temperatures-1)
        ELSE ! in all other cases choose randomly between the two neighbours of T1
                IF (unirand().LT.0.5) THEN
                        ii=-1 ! one below
                ELSE
                        ii=1 ! one above
                END IF
                T2_value=temperature_values0(T1_index+ii)
                T2_index=temperature_indices0(T1_index+ii)
        END IF
        T2_chain=MINLOC(ABS(temperature_indices-T2_index),1)
ELSE IF (jump_type.EQ.2) THEN ! choose the nearest temperature to T1
        IF (T1_value.EQ.1) THEN ! if T1=1, choose the next temperature above it
                T2_value=temperature_values0(number_of_1s+1)
                T2_index=temperature_indices0(number_of_1s+1)
        ELSE IF (T1_value.EQ.temperature_values0(number_of_temperatures)) THEN ! if T1 is the maximum temperature, choose the first below it
                T2_value=temperature_values0(number_of_temperatures-1)
                T2_index=temperature_indices0(number_of_temperatures-1)
        ELSE ! in all other cases choose the nearest temperature between the two neighbours of T1
                ii=MINLOC(ABS(temperature_values0-T1_value),1,temperature_values0.NE.T1_value)
                T2_value=temperature_values0(ii)
                T2_index=temperature_indices0(ii)
        END IF
        T2_chain=MINLOC(ABS(temperature_indices-T2_index),1)
END IF

write(*,*) 'Tempering : T2 is ', TRIM(int2str(T2_chain)), ' with a temperature of ', TRIM(dble2str(T2_value))

IF (debug_mode) WRITE(debug_unit,*)'----- T2: chain=',TRIM(int2str(T2_chain)),&
                          ', index=',TRIM(int2str(T2_index)),&
                          ', value=',TRIM(dble2str(T2_value)),' -----'

END SUBROUTINE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE do_tempering(modelled_data)

USE alpha_calculation
IMPLICIT NONE
DOUBLE PRECISION, DIMENSION(:), INTENT(IN)  :: modelled_data
DOUBLE PRECISION        :: tempering_normalization
DOUBLE PRECISION        :: alpha_swap
DOUBLE PRECISION        :: random_deviate

swap_temperatures=.FALSE.

CALL calculate_normalization(number_of_measurements,modelled_data,noise_T1,noise_T2,tempering_normalization)

alpha_swap=MINVAL([DBLE(0),-likelihood_T1/T2_value&
                           +likelihood_T2/T2_value&
                           -likelihood_T2/T1_value&
                           +likelihood_T1/T1_value&
                           +(1/T1_value-1/T2_value)*tempering_normalization])

random_deviate=unirand()
IF (LOG(random_deviate).LE.alpha_swap) swap_temperatures=.TRUE.

IF (debug_mode) THEN
        IF (swap_temperatures) WRITE(debug_unit,*)'----- Swapping: YES ------'
        IF (.NOT.swap_temperatures) WRITE(debug_unit,*)'----- Swapping: NO ------'
END IF

END SUBROUTINE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE update_temperatures(iteration)

IMPLICIT NONE
INTEGER         :: iteration

IF (swap_temperatures) THEN
        temperature_values((/T1_chain,T2_chain/))=temperature_values((/T2_chain,T1_chain/))
        temperature_indices((/T1_chain,T2_chain/))=temperature_indices((/T2_chain,T1_chain/))
        chain_temperatures(iteration,:)=(/1,T1_chain,T2_chain/)
        acceptance_tempering(1)=acceptance_tempering(1)+1
ELSE
        chain_temperatures(iteration,:)=(/0,T1_chain,T2_chain/)
        acceptance_tempering(2)=acceptance_tempering(2)+1
END IF

IF (SUM(acceptance_tempering(1:2)).GT.0) acceptance_tempering(3)=100*acceptance_tempering(1)/SUM(acceptance_tempering(1:2))
acceptance_all(13)=acceptance_tempering(3)
IF (run_number.EQ.1) chain_acceptance(iteration+1,7)=acceptance_tempering(3)
IF (run_number.GT.1) chain_acceptance(iteration,7)=acceptance_tempering(3)

END SUBROUTINE

END MODULE
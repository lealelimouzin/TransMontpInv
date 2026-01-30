! ********************************************************************* !
! Subroutine to write last Markov chain model to file in mksamples.f90
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


SUBROUTINE write_last_model()

USE chains
USE random_generator
IMPLICIT NONE
INTEGER         :: nn

! Open file containing last model of current run
OPEN(UNIT=90,FILE=TRIM(last_tmp) // '.txt',STATUS='REPLACE') ! this is just for reference
OPEN(UNIT=80,FILE=last_tmp,STATUS='REPLACE',FORM='UNFORMATTED') ! this will be read in the following iteration

! Write Voronoi tessellation
WRITE(90,*)'Voronoi model'
WRITE(90,*)ncell_current
WRITE(80)ncell_current
DO nn=1,ncell_current
        IF (opt_parameters.EQ.0) THEN
                WRITE(90,*)nn,voronoi_current(nn,1),voronoi_current(nn,2),voronoi_current(nn,3),voronoi_current(nn,4)
                WRITE(80)voronoi_current(nn,1),voronoi_current(nn,2),voronoi_current(nn,3),voronoi_current(nn,4)
        ELSE
                WRITE(90,*)nn,voronoi_current(nn,1),voronoi_current(nn,2),voronoi_current(nn,3),voronoi_current(nn,4),&
                voronoi_current(nn,5)
                WRITE(80)voronoi_current(nn,1),voronoi_current(nn,2),voronoi_current(nn,3),voronoi_current(nn,4),&
                voronoi_current(nn,5)
        END IF
END DO

! Write noise parameters
WRITE(90,*)'Noise parameters'
IF (noise_type.EQ.0) THEN
        DO nn=1,number_of_datasets
                WRITE(90,*)sigma_current(nn)
                WRITE(80)sigma_current(nn)
        END DO
ELSE IF (noise_type.EQ.1.OR.noise_type.EQ.2.OR.noise_type.EQ.5) THEN
        DO nn=1,number_of_datasets
                WRITE(90,*)a_current(nn),b_current(nn)
                WRITE(80)a_current(nn),b_current(nn)
        END DO
ELSE IF (noise_type.EQ.3.OR.noise_type.EQ.6) THEN
        DO nn=1,number_of_datasets
                WRITE(90,*)lambda_current(nn)
                WRITE(80)lambda_current(nn)
        END DO
END IF

! Write acceptance ratios
WRITE(90,*)'Acceptance ratios'
WRITE(90,*)accepted_so_far+accepted_in_run
WRITE(90,*)acceptance_total(1),acceptance_total(2),acceptance_total(3)
WRITE(90,*)acceptance_birth(1),acceptance_birth(2),acceptance_birth(3)
WRITE(90,*)acceptance_death(1),acceptance_death(2),acceptance_death(3)
WRITE(90,*)acceptance_move(1),acceptance_move(2),acceptance_move(3)
WRITE(90,*)acceptance_value(1),acceptance_value(2),acceptance_value(3)
WRITE(90,*)acceptance_move_main(1),acceptance_move_main(2),acceptance_move_main(3)
WRITE(90,*)acceptance_move_delayed(1),acceptance_move_delayed(2),acceptance_move_delayed(3)
WRITE(90,*)acceptance_value_main(1),acceptance_value_main(2),acceptance_value_main(3)
WRITE(90,*)acceptance_value_delayed(1),acceptance_value_delayed(2),acceptance_value_delayed(3)
WRITE(90,*)acceptance_noise(1),acceptance_noise(2),acceptance_noise(3)
IF (noise_type.EQ.1.OR.noise_type.EQ.2.OR.noise_type.EQ.5) THEN
        WRITE(90,*)acceptance_a(1),acceptance_a(2),acceptance_a(3)
        WRITE(90,*)acceptance_b(1),acceptance_b(2),acceptance_b(3)
END IF
WRITE(90,*)acceptance_tempering(1),acceptance_tempering(2),acceptance_tempering(3)
IF (opt_parameters.GT.0) THEN 
        WRITE(90,*)acceptance_opt(1),acceptance_opt(2),acceptance_opt(3) ! 25/11/24
        WRITE(90,*)acceptance_opt_main(1),acceptance_opt_main(2),acceptance_opt_main(3)
        WRITE(90,*)acceptance_opt_delayed(1),acceptance_opt_delayed(2),acceptance_opt_delayed(3)
END IF
WRITE(80)accepted_so_far+accepted_in_run
WRITE(80)acceptance_total(1),acceptance_total(2),acceptance_total(3)
WRITE(80)acceptance_birth(1),acceptance_birth(2),acceptance_birth(3)
WRITE(80)acceptance_death(1),acceptance_death(2),acceptance_death(3)
WRITE(80)acceptance_move(1),acceptance_move(2),acceptance_move(3)
WRITE(80)acceptance_value(1),acceptance_value(2),acceptance_value(3)
WRITE(80)acceptance_move_main(1),acceptance_move_main(2),acceptance_move_main(3)
WRITE(80)acceptance_move_delayed(1),acceptance_move_delayed(2),acceptance_move_delayed(3)
WRITE(80)acceptance_value_main(1),acceptance_value_main(2),acceptance_value_main(3)
WRITE(80)acceptance_value_delayed(1),acceptance_value_delayed(2),acceptance_value_delayed(3)
WRITE(80)acceptance_noise(1),acceptance_noise(2),acceptance_noise(3)
IF (noise_type.EQ.1.OR.noise_type.EQ.2.OR.noise_type.EQ.5) THEN
        WRITE(80)acceptance_a(1),acceptance_a(2),acceptance_a(3)
        WRITE(80)acceptance_b(1),acceptance_b(2),acceptance_b(3)
END IF
WRITE(80)acceptance_tempering(1),acceptance_tempering(2),acceptance_tempering(3)
IF (opt_parameters.GT.0) THEN 
        WRITE(80)acceptance_opt(1),acceptance_opt(2),acceptance_opt(3)
        WRITE(80)acceptance_opt_main(1),acceptance_opt_main(2),acceptance_opt_main(3)
        WRITE(80)acceptance_opt_delayed(1),acceptance_opt_delayed(2),acceptance_opt_delayed(3)
END IF

! Write number of random numbers used so far
WRITE(90,*)'Random counter'
WRITE(90,*)random_counter
WRITE(80)random_counter

! Write number of samples completed so far
WRITE(90,*)'Total samples so far'
WRITE(90,*)samples_so_far
WRITE(80)samples_so_far

! Write modelled data from current model
WRITE(90,*)'Modelled data'
DO nn=1,number_of_measurements
        WRITE(90,*)data_current(nn)
        WRITE(80)data_current(nn)
END DO

! Write temperatures
IF (parallel_tempering.EQ.1) THEN
        WRITE(90,*)'Temperatures'
        DO nn=1,number_of_temperatures
                WRITE(90,*)temperature_indices(nn),temperature_values(nn)
                WRITE(80)temperature_indices(nn),temperature_values(nn)
        END DO
END IF

CLOSE(90)
CLOSE(80)

END SUBROUTINE
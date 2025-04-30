! ********************************************************************* !
! Subroutine to choose step type in mksamples.f90
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


SUBROUTINE choose_step(iteration,step)

USE chains, ONLY: noise_type,opt_parameters
USE random_generator
IMPLICIT NONE
INTEGER, INTENT(IN)     :: iteration
INTEGER, INTENT(OUT)    :: step
DOUBLE PRECISION        :: random_deviate

! Pick a random number
random_deviate=unirand()

IF (MOD(iteration,2).EQ.0) THEN
        IF (random_deviate.LT.DBLE(0.333)) THEN
                ! BIRTH - add a cell
                step=1
        ELSE IF (random_deviate.LT.DBLE(0.666)) THEN
                ! DEATH - remove a cell
                step=2
        ELSE
                ! MOVE - move a cell
                step=3
        END IF
ELSE
        IF (opt_parameters.GT.0) THEN 
                IF (noise_type.NE.4) THEN
                        IF (random_deviate.LT.DBLE(0.333)) THEN
                                ! VALUE - change the value of a cell
                                step=4
                        ELSE IF (random_deviate.LT.DBLE(0.666)) THEN
                                ! NOISE - change a data noise parameter
                                step=5
                        ELSE
                                ! OPTION - change either an angle or an aperture
                                step=6
                        END IF  
                ELSE IF (noise_type.EQ.4) THEN
                        IF (random_deviate.LT.DBLE(0.5)) THEN
                                ! VALUE - change the value of a cell
                                step=4
                        ELSE
                                ! OPTION - change either an angle or an aperture
                                step=6   
                        END IF
                END IF
        ELSE
                IF (noise_type.NE.4) THEN
                        IF (random_deviate.LT.DBLE(0.500)) THEN
                                ! VALUE - change the value of a cell
                                step=4
                        ELSE
                                ! NOISE - change a data noise parameter
                                step=5
                        END IF  
                ELSE IF (noise_type.EQ.4) THEN
                                ! VALUE - change the value of a cell
                                step=4
                END IF
        END IF
END IF
step=6
END SUBROUTINE

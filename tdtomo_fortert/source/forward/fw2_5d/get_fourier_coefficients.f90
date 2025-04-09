! ********************************************************************* !
! This file is part of module "modelling_2_5D" in modelling_2_5D.f90.
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


SUBROUTINE get_fourier_coefficients()

USE modelling_2_5D, ONLY: number_of_coefficients,k_parameters,g_parameters
IMPLICIT NONE
INTERFACE
        SUBROUTINE optimize_k_g()
        END SUBROUTINE
END INTERFACE

IF (ALLOCATED(k_parameters)) DEALLOCATE(k_parameters)
IF (ALLOCATED(g_parameters)) DEALLOCATE(g_parameters)

IF (number_of_coefficients.EQ.0) THEN
        ALLOCATE(k_parameters(4),g_parameters(4))
        k_parameters=(/0.0217102,0.2161121,1.0608400,5.0765870/)
        g_parameters=(/0.0463660,0.2365931,1.0382080,5.3648010/)
        number_of_coefficients=4
ELSE IF (number_of_coefficients.GT.0) THEN
        CALL optimize_k_g()
ELSE
        WRITE(*,*)
        WRITE(*,*)'********************************************'
        WRITE(*,*)'ERROR in subroutine get_fourier_coefficients'
        WRITE(*,*)'--------------------------------------------'
        WRITE(*,*)'The number of Fourier coefficients must'
        WRITE(*,*)'be an integer >=0!!'
        WRITE(*,*)'Check first line of input file forward.in'
        WRITE(*,*)'TERMINATING PROGRAM!!'
        WRITE(*,*)'********************************************'
        WRITE(*,*)
        STOP
END IF

END SUBROUTINE
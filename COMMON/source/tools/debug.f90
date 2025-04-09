! ********************************************************************* !
! Module for debugging in mksamples.f90 and procsamples.f90
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


MODULE debug

IMPLICIT NONE
LOGICAL, SAVE           :: debug_mode=.FALSE. ! .TRUE. runs the code in debug mode, .FALSE. in normal mode
INTEGER, SAVE           :: debug_unit=5000

 CONTAINS

SUBROUTINE initialize_debug(chain_name)

IMPLICIT NONE
CHARACTER (LEN=30)      :: chain_name

OPEN(UNIT=debug_unit,FILE='debug.txt' // TRIM(chain_name),STATUS='REPLACE',FORM='FORMATTED',ACCESS='SEQUENTIAL')

END SUBROUTINE


SUBROUTINE close_debug(chain_name)

IMPLICIT NONE
CHARACTER (LEN=30)      :: chain_name

CLOSE(debug_unit)

END SUBROUTINE

END MODULE
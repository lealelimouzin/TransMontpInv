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


SUBROUTINE calculate_DS()

USE modelling_2_5D, ONLY: D_size,D_indices,D_values,&
                          S_size,S_indices,S_values,&
                          DS_size,DS_indices,DS_values
IMPLICIT NONE
INTEGER                         :: nn,cc
DOUBLE PRECISION, DIMENSION(1)  :: value
REAL                            :: start,finish
LOGICAL                         :: verbose=.FALSE.

  IF (D_size(2).NE.S_size) THEN
    WRITE(*,*)
    WRITE(*,*)'*****************************************'
    WRITE(*,*)'ERROR in subroutine calculate_DS'
    WRITE(*,*)'-----------------------------------------'
    WRITE(*,*)'The number of columns of D does not match'
    WRITE(*,*)'the number of rows of S!!'
    WRITE(*,*)'TERMINATING PROGRAM!!'
    WRITE(*,*)'*****************************************'
    WRITE(*,*)
    STOP
  END IF
  
  ! Initial values
  DS_size=D_size
  DS_indices=D_indices
  DS_values=D_values
  
  IF (verbose) CALL cpu_time(start)
  
  ! Loop over the nonzero elements of the divergence array
  DO nn=1,SIZE(D_values)
    
    ! Get column number
    cc=D_indices(nn,2)
    
    ! Give warning message in case the entries of diagonal
    ! matrix S are not ordered, then get value of S in column cc
    IF (cc.NE.S_indices(cc)) THEN
      WRITE(*,*)
      WRITE(*,*)'**********************************'
      WRITE(*,*)'WARNING in subroutine calculate_DS'
      WRITE(*,*)'----------------------------------'
      WRITE(*,*)'The entries of diagonal matrix S'
      WRITE(*,*)'are not ordered!!'
      WRITE(*,*)'**********************************'
      WRITE(*,*)
      value=PACK(S_values,S_indices.EQ.cc)
    ELSE
      value(1)=S_values(cc)
    END IF
    
    ! Multiply
    DS_values(nn)=DS_values(nn)*value(1)
    
  END DO
    
  IF (verbose) CALL cpu_time(finish)
  IF (verbose) PRINT '("Time for A+B = ",f12.9," seconds.")',finish-start

END SUBROUTINE
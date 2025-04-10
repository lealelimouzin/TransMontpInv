! ********************************************************************* !
! Program to add random Gaussian noise to a resistance dataset
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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This program creates random Gaussian noise and adds it to !
! a resistance dataset.                                     !
! Noise can be crated independently for different datasets, !
! with the same level of noise for all readings or with a   !
! resistance-dependent level of noise.                      !
!                                                           !
! Erica Galetti, September 2017                             !
! erica.galetti@ed.ac.uk                                    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

PROGRAM add_noise

USE string_conversion
USE random_generator
IMPLICIT NONE
INTEGER :: noise_type,rand_seed,number_of_datasets
INTEGER :: nn,IOstatus,switch,switch2,dataset
INTEGER:: p1,p2,c1,c2,number_of_measurements
DOUBLE PRECISION :: resistance,uncertainty,error
CHARACTER(LEN=30) :: measurements_file,datasets_file
CHARACTER(LEN=30) :: measurements_with_noise_file,noise_file
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: noise

OPEN(UNIT=10,FILE='add_noise.in',STATUS='OLD')
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,1)measurements_file
READ(10,1)datasets_file
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,1)measurements_with_noise_file
READ(10,1)noise_file
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,*)noise_type
READ(10,*)rand_seed
READ(10,*)number_of_datasets
IF (noise_type.EQ.0.OR.noise_type.EQ.2) THEN
        ALLOCATE(noise(number_of_datasets,1))
ELSE IF (noise_type.EQ.1) THEN
        ALLOCATE(noise(number_of_datasets,2))
ELSE
        WRITE(*,*)'Invalid input: noise_type'
        WRITE(*,*)'TERMINATING PROGRAM!!'
        STOP
END IF
DO dataset=1,number_of_datasets
        IF (noise_type.EQ.0.OR.noise_type.EQ.2) THEN
                READ(10,*)noise(dataset,1)
        ELSE IF (noise_type.EQ.1) THEN
                READ(10,*)noise(dataset,1),noise(dataset,2)
        END IF
END DO
1   FORMAT(a30)
CLOSE(10)

CALL initialize_random_generator(rand_seed)

IF (number_of_datasets.EQ.1) THEN
        WRITE(*,*)'Creating Gaussian noise for 1 dataset'
ELSE IF (number_of_datasets.GT.1) THEN
        WRITE(*,*)'Creating Gaussian noise for ',&
                  TRIM(int2str(number_of_datasets)),' datasets'
ELSE
        WRITE(*,*)'Invalid input: number_of_datasets'
        WRITE(*,*)'TERMINATING PROGRAM!!'
        STOP 
END IF

! Open INPUT files
OPEN(UNIT=10,FILE=measurements_file,STATUS='OLD')
READ(10,*)number_of_measurements
IF (number_of_datasets.GT.1) OPEN(UNIT=20,FILE=datasets_file,STATUS='OLD')

! Open OUTPUT files
OPEN(UNIT=30,FILE=noise_file,STATUS='REPLACE')
OPEN(UNIT=50,FILE=measurements_with_noise_file,STATUS='REPLACE')
WRITE(50,*)number_of_measurements

! Create noisy data
DO nn=1,number_of_measurements
        
        ! Read from input files
        READ(10,*)switch,p1,p2,c1,c2,resistance
        IF (number_of_datasets.GT.1) THEN
                READ(20,*)dataset
        ELSE IF (number_of_datasets.EQ.1) THEN
                dataset=1
        END IF
        
        ! Calculate noise and write to file
        IF (switch.GT.0) THEN
                IF (noise_type.EQ.0) THEN
                        uncertainty=noise(dataset,1)
                ELSE IF (noise_type.EQ.1) THEN
                        uncertainty=SQRT(noise(dataset,1)**2*resistance**2+noise(dataset,2)**2)
                ELSE IF (noise_type.EQ.2) THEN
                        uncertainty=noise(dataset,1)/DBLE(100)*resistance
                END IF
                error=gausrand()*uncertainty
                WRITE(30,*)error,uncertainty
                WRITE(50,*)nn,p1,p2,c1,c2,resistance+error,uncertainty
        ELSE IF (switch.EQ.0) THEN
                WRITE(30,*)0,0
                WRITE(50,*)0,p1,p2,c1,c2,100.0,0.1
        ELSE
                WRITE(*,*)'Invalid switch in measurements'
                WRITE(*,*)'file: ',measurements_file
                WRITE(*,*)'TERMINATING PROGRAM!!' 
                STOP
        END IF
        
END DO
CLOSE(10)
IF (number_of_datasets.GT.1) CLOSE(20)
CLOSE(30)
CLOSE(50)

DEALLOCATE(noise)

END PROGRAM

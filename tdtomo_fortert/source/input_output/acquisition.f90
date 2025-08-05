! ********************************************************************* !
! Module for reading data files (part of TERT package)
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


MODULE acquisition

IMPLICIT NONE
CHARACTER (LEN=30), SAVE                                :: electrode_file,data_file,datasets_file
CHARACTER (LEN=30), SAVE                                :: std_file ! 07/07/25
INTEGER, SAVE                                           :: number_of_electrodes,number_of_measurements,number_of_datasets
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, SAVE     :: electrodes
INTEGER, DIMENSION(:,:), ALLOCATABLE, SAVE              :: quadrupoles
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, SAVE       :: observed_data,data_noise
INTEGER, DIMENSION(:), ALLOCATABLE, SAVE                :: dataset_associations
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, SAVE       :: std_val ! 07/07/25


CONTAINS


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE read_data(problem_type)

IMPLICIT NONE
INTEGER :: problem_type

CALL read_electrodes()
CALL read_measurements(problem_type)
CALL read_distances()

END SUBROUTINE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE read_electrodes()
! Reads the file containing electrode coordinates. This file is structured as follows:
!
!       line 1          number of electrodes in array
!       line 2-end      electrode coordinates given as:
!                               column 1        electrode number
!                               column 2        electrode x coordinate
!                               column 3        electrode z coordinate (if inversion is 2D) OR
!                                               electrode y coordinate (if inversion is 3D)
!                               column 4        electrode z coordinate (if inversion is 3D)

IMPLICIT NONE
INTEGER                 :: electrode_number
DOUBLE PRECISION        :: x_electrode,y_electrode,z_electrode
INTEGER                 :: nn,IOstatus
LOGICAL                 :: verbose=.FALSE.

IF (ALLOCATED(electrodes)) DEALLOCATE(electrodes)

OPEN(UNIT=10,FILE=electrode_file,STATUS='OLD')
READ(10,*)number_of_electrodes
ALLOCATE(electrodes(number_of_electrodes,3)) ! check whether a y coordinate is given (file has 4 columns)
DO nn=1,number_of_electrodes
        READ(10,*,IOSTAT=IOstatus)electrode_number,x_electrode,y_electrode,z_electrode
        electrodes(electrode_number,:)=[x_electrode,y_electrode,z_electrode]
END DO
IF (verbose) WRITE(*,*)'IOstatus=',IOstatus

IF (IOstatus.EQ.0) THEN ! if no error was encountered, a y coordinate is given (file has 4 columns)
        IF (verbose) WRITE(*,*)'Coordinates given for electrodes: x, y, z'
ELSEIF (IOstatus.NE.0) THEN ! in this case a y coordinate is not given (file has 3 columns)
        IF (verbose) WRITE(*,*)'Coordinates given for electrodes: x, z'
        REWIND(10)
        READ(10,*)
        DEALLOCATE(electrodes)
        ALLOCATE(electrodes(number_of_electrodes,2))
        DO nn=1,number_of_electrodes
                READ(10,*)electrode_number,x_electrode,z_electrode
                electrodes(electrode_number,:)=[x_electrode,z_electrode]
        END DO
END IF
IF (verbose) WRITE(*,*)'IOstatus=',IOstatus
CLOSE(10)

END SUBROUTINE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE read_measurements(problem_type)
! Reads the file containing resistance measurements. This file is structured as follows:
!
!       line 1          number of measurements
!       line 2-end      measurements given as:
!                               column 1        switch (>0 if measurement is to be included in the inversion;
!                                               =0 if measurement should not be used in the inversion)
!                               column 2-5      electrode numbers for P+ P- C+ C-
!                               column 6        measured resistance (in Ohm)
!                               column 7        uncertainty in measured resistance (in Ohm)

IMPLICIT NONE
INTEGER                                         :: problem_type
INTEGER                                         :: switch,p1,p2,c1,c2,dataset
DOUBLE PRECISION                                :: resistance,error
INTEGER                                         :: nn,valid_count
INTEGER, DIMENSION(:,:), ALLOCATABLE            :: worki
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE     :: workd

IF (ALLOCATED(worki)) DEALLOCATE(worki)
IF (ALLOCATED(workd)) DEALLOCATE(workd)
IF (ALLOCATED(quadrupoles)) DEALLOCATE(quadrupoles)
IF (ALLOCATED(observed_data)) DEALLOCATE(observed_data)
IF (ALLOCATED(data_noise)) DEALLOCATE(data_noise)
IF (ALLOCATED(dataset_associations)) DEALLOCATE(dataset_associations)

OPEN(UNIT=10,FILE=data_file,STATUS='OLD')
READ(10,*)number_of_measurements
ALLOCATE(quadrupoles(number_of_measurements,4))
ALLOCATE(observed_data(number_of_measurements))
ALLOCATE(data_noise(number_of_measurements))
ALLOCATE(dataset_associations(number_of_measurements))
IF (number_of_datasets.GT.1) THEN
        OPEN(UNIT=20,FILE=datasets_file,STATUS='OLD')
ELSE IF (number_of_datasets.EQ.1) THEN
        dataset_associations=1
END IF

valid_count=0
DO nn=1,number_of_measurements
        IF (problem_type.EQ.0) READ(10,*)switch,p1,p2,c1,c2,resistance,error ! inverse problem: we read the data
        IF (problem_type.EQ.1) READ(10,*)switch,p1,p2,c1,c2 ! forward problem: we only need the quadrupoles
        IF (number_of_datasets.GT.1.AND.problem_type.EQ.0) READ(20,*)dataset
        IF (switch.GT.0) THEN
                valid_count=valid_count+1
                quadrupoles(valid_count,:)=[p1,p2,c1,c2]
                IF (problem_type.EQ.0) THEN
                        observed_data(valid_count)=resistance
                        data_noise(valid_count)=error
                        IF (number_of_datasets.GT.1) dataset_associations(valid_count)=dataset
                END IF
        END IF
END DO
CLOSE(10)
IF (number_of_datasets.GT.1.AND.problem_type.EQ.0) CLOSE(20)
IF (valid_count.LT.number_of_measurements) THEN
        number_of_measurements=valid_count ! number of valid readings
        ! Resize quadrupoles to the number of valid readings
        ALLOCATE(worki(valid_count,4))
        worki=quadrupoles(1:valid_count,:)
        DEALLOCATE(quadrupoles)
        ALLOCATE(quadrupoles(valid_count,4))
        quadrupoles=worki
        DEALLOCATE(worki)
        IF (problem_type.EQ.0) THEN
                ! Resize observed_data to the number of valid readings
                ALLOCATE(workd(valid_count))
                workd=observed_data(1:valid_count)
                DEALLOCATE(observed_data)
                ALLOCATE(observed_data(valid_count))
                observed_data=workd
                ! Resize data_noise to the number of valid readings
                workd=data_noise(1:valid_count)
                DEALLOCATE(data_noise)
                ALLOCATE(data_noise(valid_count))
                data_noise=workd
                DEALLOCATE(workd)
                ! Resize dataset_associations to the number of valid readings
                ALLOCATE(worki(valid_count,1))
                worki(:,1)=dataset_associations(1:valid_count)
                DEALLOCATE(dataset_associations)
                ALLOCATE(dataset_associations(valid_count))
                dataset_associations=worki(:,1)
                DEALLOCATE(worki)
        END IF
END IF

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE read_distances()
! Reads the file containing electrode coordinatesdistances between electrodes of a same quadrupole. This file is structured as follows:
!
!       line 1          number of electrodes in array
!       line 2-end      distances 

IMPLICIT NONE
INTEGER                 :: electrode_number
DOUBLE PRECISION        :: val
INTEGER                 :: nn,IOstatus
LOGICAL                 :: verbose=.FALSE.

IF (ALLOCATED(std_val)) DEALLOCATE(std_val)


OPEN(UNIT=10,FILE=std_file,STATUS='OLD')
READ(10,*)number_of_electrodes

ALLOCATE(std_val(number_of_electrodes))

DO nn=1,number_of_electrodes
        READ(10,*,IOSTAT=IOstatus)val
        std_val(nn)=val
END DO
IF (verbose) WRITE(*,*)'IOstatus=',IOstatus

!IF (IOstatus.EQ.0) THEN ! if no error was encountered, a y coordinate is given (file has 4 columns)
!        IF (verbose) WRITE(*,*)'Coordinates given for electrodes: x, y, z'
!ELSEIF (IOstatus.NE.0) THEN ! in this case a y coordinate is not given (file has 3 columns)
!        IF (verbose) WRITE(*,*)'Coordinates given for electrodes: x, z'
!        REWIND(10)
!        READ(10,*)
!        DEALLOCATE(electrodes)
!        ALLOCATE(electrodes(number_of_electrodes,2))
!        DO nn=1,number_of_electrodes
!                READ(10,*)electrode_number,x_electrode,z_electrode
!                electrodes(electrode_number,:)=[x_electrode,z_electrode]
!        END DO
!END IF
IF (verbose) WRITE(*,*)'IOstatus=',IOstatus
CLOSE(10)

END SUBROUTINE

END MODULE
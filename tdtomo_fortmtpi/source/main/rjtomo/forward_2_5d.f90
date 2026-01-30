! ********************************************************************* !
! Program to model electrical resistance data in 2.5D
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


PROGRAM forward_2_5d

USE acquisition
USE modelling_2_5d
USE string_conversion
IMPLICIT NONE
INTEGER :: nn,info,i,j,ii,jj
INTEGER :: model_type,log10_rho,number_of_cells
CHARACTER(LEN=30) :: model_file,forward_file
CHARACTER(LEN=30) :: resistance_file,apparent_resistivity_file
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: voronoi
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: resistance,apparent_resistivity
LOGICAL :: verb=.FALSE.
INTERFACE
        SUBROUTINE create_resistivity_grid2D(ncell,voronoi_model)
                INTEGER, INTENT(IN)                                     :: ncell
                DOUBLE PRECISION, DIMENSION(ncell,3), INTENT(IN)        :: voronoi_model
        END SUBROUTINE
END INTERFACE
INTERFACE
        SUBROUTINE calculate_apparent_resistivity(electrodes,quadrupoles,resistance,apparent_resistivity)
                DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN)    :: electrodes
                INTEGER, DIMENSION(:,:), INTENT(IN)             :: quadrupoles
                DOUBLE PRECISION, DIMENSION(:), INTENT(IN)      :: resistance
                DOUBLE PRECISION, DIMENSION(:), INTENT(OUT)     :: apparent_resistivity
        END SUBROUTINE
END INTERFACE


! Read input file
OPEN(UNIT=10,FILE='forward_2_5d.in',STATUS='OLD')
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,1)electrode_file
READ(10,1)data_file
READ(10,*)model_type
READ(10,1)model_file
READ(10,*)log10_rho
READ(10,1)forward_file
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,1)resistance_file
READ(10,1)apparent_resistivity_file
1   FORMAT(a30)
CLOSE(10)

CALL read_data(1)
ALLOCATE(resistance(number_of_measurements))
ALLOCATE(apparent_resistivity(number_of_measurements))

IF (verb) THEN
        
        WRITE(*,*)'Number of electrodes: ',TRIM(int2str(number_of_electrodes))
        DO nn=1,number_of_electrodes
                WRITE(*,*)electrodes(nn,:)
        END DO
        WRITE(*,*)
        
        WRITE(*,*)'Number of measurements: ',TRIM(int2str(number_of_measurements))
        DO nn=1,number_of_measurements
                WRITE(*,*)nn,quadrupoles(nn,:)
        END DO
        WRITE(*,*)
        
END IF

! Define input file name for forward modelling
fw2_5d_file=forward_file

! Initialise the modelling grid and various modelling parameters
CALL initialize_fw2_5d()

! Read model
OPEN(UNIT=20,FILE=model_file,STATUS='OLD')
IF (model_type.EQ.1) THEN
        
        ! Read Voronoi cells
        READ(20,*)number_of_cells
        ALLOCATE(voronoi(number_of_cells,4))
        DO nn=1,number_of_cells
                READ(20,*)i,voronoi(nn,1),voronoi(nn,2),voronoi(nn,3),voronoi(nn,4)
        END DO

        IF (log10_rho.EQ.0) THEN
                voronoi(:,4)=LOG10(voronoi(:,4))
        ELSE IF (log10_rho.EQ.1) THEN
                voronoi(:,4)=voronoi(:,4)
        ELSE
                WRITE(*,*)'Invalid input: log10_rho'
                WRITE(*,*)'TERMINATING PROGRAM!!'
                STOP
        END IF
        
        ! Create 2D resistivity model
        CALL create_resistivity_grid2D(number_of_cells,voronoi(1:number_of_cells,[1,3,4]))
        
ELSE IF (model_type.EQ.0) THEN
        
        ! Read model from top-left to bottom-right corner
        DO ii=1,x_number_of_cells
                DO jj=1,z_number_of_cells
                        READ(20,*)resistivity(ii,jj)
                END DO
        END DO

        IF (log10_rho.EQ.0) THEN
                resistivity=resistivity
        ELSE IF (log10_rho.EQ.1) THEN
                resistivity=10**resistivity
        ELSE
                WRITE(*,*)'Invalid input: log10_rho'
                WRITE(*,*)'TERMINATING PROGRAM!!'
                STOP
        END IF
        
ELSE
        
        WRITE(*,*)'Invalid input: model_type'
        WRITE(*,*)'TERMINATING PROGRAM!!'
        STOP
        
END IF
CLOSE(20)

! Run forward modelling
CALL run_fw2_5d(resistance,info)
IF (info.NE.0) THEN
        WRITE(*,*)'Modelling error.'
        WRITE(*,*)'TERMINATING PROGRAM!!'
        STOP
END IF

CALL clear_memory_fw2_5d()

! Write resistances to file
OPEN(UNIT=50,FILE=resistance_file,STATUS='REPLACE')
WRITE(50,*)number_of_measurements
DO nn=1,number_of_measurements
        WRITE(50,'(5I6,2F18.8)')nn,quadrupoles(nn,1),quadrupoles(nn,2),&
                                   quadrupoles(nn,3),quadrupoles(nn,4),&
                                   resistance(nn),LOG10(resistance(nn))
END DO
CLOSE(50)

CALL calculate_apparent_resistivity(electrodes,quadrupoles,resistance,apparent_resistivity)

! Write apparent resistivities to file
OPEN(UNIT=50,FILE=apparent_resistivity_file,STATUS='REPLACE')
WRITE(50,*)number_of_measurements
DO nn=1,number_of_measurements
        WRITE(50,'(5I6,2F18.8)')nn,quadrupoles(nn,1),quadrupoles(nn,2),&
                                   quadrupoles(nn,3),quadrupoles(nn,4),&
                                   apparent_resistivity(nn),LOG10(apparent_resistivity(nn))
END DO
CLOSE(50)

END PROGRAM


SUBROUTINE calculate_apparent_resistivity(electrodes,quadrupoles,resistance,apparent_resistivity)

IMPLICIT NONE
DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN)    :: electrodes
INTEGER, DIMENSION(:,:), INTENT(IN)             :: quadrupoles
DOUBLE PRECISION, DIMENSION(:), INTENT(IN)      :: resistance
DOUBLE PRECISION, DIMENSION(:), INTENT(OUT)     :: apparent_resistivity
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE     :: geometric_factor
DOUBLE PRECISION, PARAMETER                     :: pi=DBLE(4)*ATAN(DBLE(1))

ALLOCATE(geometric_factor(SIZE(resistance)))

geometric_factor=1/ABS(electrodes(quadrupoles(:,1),1)-electrodes(quadrupoles(:,3),1)) ! P+C+
geometric_factor=geometric_factor-1/ABS(electrodes(quadrupoles(:,1),1)-electrodes(quadrupoles(:,4),1)) ! P+C-
geometric_factor=geometric_factor+1/ABS(electrodes(quadrupoles(:,2),1)-electrodes(quadrupoles(:,4),1)) ! P-C-
geometric_factor=geometric_factor-1/ABS(electrodes(quadrupoles(:,2),1)-electrodes(quadrupoles(:,3),1)) ! P-C+

geometric_factor=2*pi/geometric_factor

apparent_resistivity=resistance*geometric_factor

DEALLOCATE(geometric_factor)

END SUBROUTINE

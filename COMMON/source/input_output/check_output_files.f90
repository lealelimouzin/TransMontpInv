! ********************************************************************* !
! Subroutine to check consistency of output files from mksamples.f90
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


SUBROUTINE check_output_files(length_ok)

USE chains
USE string_conversion
USE acquisition, ONLY: number_of_measurements
USE input_files, ONLY: samples_file,residuals_file,ncell_file,noise_file,&
                       acceptance_file,misfits_file,temperatures_file,&
                       save_residuals
IMPLICIT NONE
LOGICAL, INTENT(OUT)                    :: length_ok
INTEGER                                 :: nn,IOstatus
INTEGER, DIMENSION(:,:), ALLOCATABLE    :: file_lengths
INTEGER, DIMENSION(:), ALLOCATABLE      :: ncell_initial
INTEGER                                 :: pp,accept,noise_count
CHARACTER(100) :: samples_file_test

ALLOCATE(file_lengths(number_of_chains,8))
file_lengths=0
ALLOCATE(ncell_initial(number_of_chains))
ncell_initial=0

length_ok=.TRUE.

DO nn=1,number_of_chains

        ! File of chain steps
        OPEN(UNIT=samples_unit,FILE=TRIM(samples_file) // TRIM(int2str(nn)),&
             STATUS='OLD',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
        READ(samples_unit) ! ignore line of file type
        READ(samples_unit)ncell_initial(nn)
        noise_count=0
        IOstatus=0
        DO WHILE (IOstatus.EQ.0)
                IF (file_lengths(nn,3).LT.ncell_initial(nn)) THEN
                        READ(samples_unit,IOSTAT=IOstatus)
                ELSE IF (file_lengths(nn,3).GE.ncell_initial(nn)) THEN
                        READ(samples_unit,IOSTAT=IOstatus)pp,step,accept
                        IF (IOstatus.EQ.0.AND.step.EQ.5.AND.accept.EQ.1) noise_count=noise_count+1
                END IF
                IF (IOstatus.EQ.0) THEN
                        file_lengths(nn,3)=file_lengths(nn,3)+1 ! length of file of chain steps (actual)
                ELSE
                        EXIT
                END IF
        END DO
        CLOSE(samples_unit)
        
        ! File of number of cells
        OPEN(UNIT=ncell_unit,FILE=TRIM(ncell_file) // TRIM(int2str(nn)),&
             STATUS='OLD',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
        READ(ncell_unit) ! ignore line of file type
        IOstatus=0
        DO WHILE (IOstatus.EQ.0)
                READ(ncell_unit,IOSTAT=IOstatus)
                IF (IOstatus.EQ.0) THEN
                        file_lengths(nn,2)=file_lengths(nn,2)+1 ! length of file of number of cells
                ELSE
                        EXIT
                END IF
        END DO
        file_lengths(nn,1)=ncell_initial(nn)+file_lengths(nn,2)-1 ! length of file of chain steps (expected)
        CLOSE(ncell_unit)
        
        ! File of noise
        IF (noise_type.NE.4)&
        OPEN(UNIT=noise_unit,FILE=TRIM(noise_file) // TRIM(int2str(nn)),&
             STATUS='OLD',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
        IF (noise_type.EQ.4)&
        OPEN(UNIT=noise_unit,FILE=TRIM(noise_file),&
             STATUS='OLD',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
        READ(noise_unit) ! ignore line of file type
        READ(noise_unit) ! ignore line of noise type
        READ(noise_unit) ! ignore line with number of datasets/readings
        IOstatus=0
        DO WHILE (IOstatus.EQ.0)
                READ(noise_unit,IOSTAT=IOstatus)
                IF (IOstatus.EQ.0) THEN
                        file_lengths(nn,4)=file_lengths(nn,4)+1 ! length of file of noise parameters
                ELSE
                        EXIT
                END IF
        END DO
        CLOSE(noise_unit)
        
        ! File of misfits
        OPEN(UNIT=misfits_unit,FILE=TRIM(misfits_file) // TRIM(int2str(nn)),&
             STATUS='OLD',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
        READ(misfits_unit) ! ignore line of file type
        IOstatus=0
        DO WHILE (IOstatus.EQ.0)
                READ(misfits_unit,IOSTAT=IOstatus)
                IF (IOstatus.EQ.0) THEN
                        file_lengths(nn,5)=file_lengths(nn,5)+1 ! length of file of number of misfits
                ELSE
                        EXIT
                END IF
        END DO
        CLOSE(misfits_unit)
        
        ! File of acceptance ratios
        OPEN(UNIT=acceptance_unit,FILE=TRIM(acceptance_file) // TRIM(int2str(nn)),&
             STATUS='OLD',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
        READ(acceptance_unit) ! ignore line of file type
        IOstatus=0
        DO WHILE (IOstatus.EQ.0)
                READ(acceptance_unit,IOSTAT=IOstatus)
                IF (IOstatus.EQ.0) THEN
                        file_lengths(nn,6)=file_lengths(nn,6)+1 ! length of file of number of acceptance ratios
                ELSE
                        EXIT
                END IF
        END DO
        CLOSE(acceptance_unit)
        
        ! File of temperatures
        IF (parallel_tempering.EQ.1) THEN
                OPEN(UNIT=temperatures_unit,FILE=TRIM(temperatures_file),&
                     STATUS='OLD',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
                READ(temperatures_unit) ! ignore line of file type
                READ(temperatures_unit) ! ignore line with number of tempratures
                IOstatus=0
                DO WHILE (IOstatus.EQ.0)
                        READ(temperatures_unit,IOSTAT=IOstatus)
                        IF (IOstatus.EQ.0) THEN
                                file_lengths(nn,7)=file_lengths(nn,7)+1 ! length of file of number of temperatures
                        ELSE
                                EXIT
                        END IF
                END DO
                CLOSE(temperatures_unit)
        END IF
        
        ! File of residuals
        IF (save_residuals.EQ.1) THEN
                OPEN(UNIT=residuals_unit,FILE=TRIM(residuals_file) // TRIM(int2str(nn)),&
                     STATUS='OLD',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
                READ(residuals_unit) ! ignore line of file type
                READ(residuals_unit) ! ignore line with save_residuals
                IOstatus=0
                DO WHILE (IOstatus.EQ.0)
                        READ(residuals_unit,IOSTAT=IOstatus)
                        IF (IOstatus.EQ.0) THEN
                                file_lengths(nn,8)=file_lengths(nn,8)+1 ! length of file of residuals
                        ELSE
                                EXIT
                        END IF
                END DO
                CLOSE(residuals_unit)
        END IF
        
        ! The length of samples_file should be given by
        ! the length of ncell_file + the initial number of cells - 1
        IF (file_lengths(nn,1).NE.file_lengths(nn,3)) THEN
                length_ok=.FALSE.
                EXIT
        END IF
        
        ! The length of ncell_file, misfits_file and acceptance_file
        ! should be the same
        IF (file_lengths(nn,2).NE.file_lengths(nn,5).OR.&
            file_lengths(nn,2).NE.file_lengths(nn,6).OR.&
            file_lengths(nn,5).NE.file_lengths(nn,6)) THEN
                length_ok=.FALSE.
                EXIT
        END IF
        
        IF (noise_type.NE.4) THEN
                ! The length of noise_file should be given by the length
                ! of ncell_file * number_of_datasets
                IF (file_lengths(nn,2)*number_of_datasets.NE.file_lengths(nn,4)) THEN
                        length_ok=.FALSE.
                        EXIT
                END IF
        ELSE IF (noise_type.EQ.4) THEN
                ! The length of noise_file should be equal to the
                ! number of readings
                IF (file_lengths(nn,4).NE.number_of_measurements) THEN
                        length_ok=.FALSE.
                        EXIT
                END IF
        END IF
        
        ! The length of ncell_file should be equal to the number
        ! of samples so far + 1
        IF (file_lengths(nn,2).NE.samples_so_far+1) THEN
                length_ok=.FALSE.
                !write(*,*)'l210 whith chain ', chain_number, ': ',length_ok ! 25/11/24
                !write(*,*)'l213 check output file_lengths: ',file_lengths(nn,2)
                !write(*,*)'l213 check output fsamples_so_far+1: ',samples_so_far+1
                EXIT
        END IF
        
        ! Compare the length of temperatures_file with the
        ! number of samples so far
        IF (parallel_tempering.EQ.1) THEN
                IF (tempering_start.LE.1) THEN
                        IF (file_lengths(nn,7).NE.samples_so_far+number_of_temperatures) THEN
                                length_ok=.FALSE.
                                EXIT
                        END IF
                ELSE IF (tempering_start.GT.1) THEN
                        IF (file_lengths(nn,7).NE.samples_so_far+number_of_temperatures-tempering_start+1) THEN
                                length_ok=.FALSE.
                                EXIT
                        END IF
                END IF
        END IF
        
        ! The length of residuals_file should be given by the
        ! number of readings * ( number of accepted samples - number of noise steps + 1)
        IF (save_residuals.EQ.1) THEN
                IF (file_lengths(nn,8).NE.(accepted_so_far_all(nn)-noise_count+1)*number_of_measurements) THEN
                        length_ok=.FALSE.
                        EXIT
                END IF
        END IF
        
END DO

IF (length_ok) THEN
        IF (ANY(file_lengths(:,2).NE.file_lengths(1,2))) length_ok=.FALSE.
END IF

IF (.NOT.length_ok) THEN
        WRITE(*,*)
        WRITE(*,*)'********************************************'
        WRITE(*,*)'Lengths of output files are inconsistent :-('
        WRITE(*,*)'TERMINATING PROGRAM!!'
        WRITE(*,*)'********************************************'
        WRITE(*,*)'Checked up to chain ',TRIM(int2str(nn)),', which showed'
        WRITE(*,*)'inconsistencies. See below for details.'
        WRITE(*,*)'--------------------------------------------'
        DO nn=1,number_of_chains
                WRITE(*,*)'CHAIN ',TRIM(int2str(nn))
                WRITE(*,*)TRIM(samples_file) // TRIM(int2str(nn)),' actual: ',TRIM(int2str(file_lengths(nn,3)))
                WRITE(*,*)TRIM(samples_file) // TRIM(int2str(nn)),' expected: ',TRIM(int2str(file_lengths(nn,1)))
                WRITE(*,*)TRIM(ncell_file) // TRIM(int2str(nn)),': ',TRIM(int2str(file_lengths(nn,2)))
                WRITE(*,*)TRIM(noise_file) // TRIM(int2str(nn)),' actual: ',TRIM(int2str(file_lengths(nn,4)))
                IF (noise_type.NE.4) WRITE(*,*)TRIM(noise_file) // TRIM(int2str(nn)),&
                                 ' expected: ',TRIM(int2str(file_lengths(nn,2)*number_of_datasets))
                IF (noise_type.EQ.4) WRITE(*,*)TRIM(noise_file) // TRIM(int2str(nn)),&
                                 ' expected: ',TRIM(int2str(number_of_measurements))
                WRITE(*,*)TRIM(misfits_file) // TRIM(int2str(nn)),': ',TRIM(int2str(file_lengths(nn,5)))
                WRITE(*,*)TRIM(acceptance_file) // TRIM(int2str(nn)),': ',TRIM(int2str(file_lengths(nn,6)))
                IF (parallel_tempering.EQ.1) WRITE(*,*)TRIM(temperatures_file),': ',TRIM(int2str(file_lengths(nn,7)))
                IF (save_residuals.EQ.1) WRITE(*,*)TRIM(residuals_file),' actual: ',TRIM(int2str(file_lengths(nn,8)))
                IF (save_residuals.EQ.1) WRITE(*,*)TRIM(residuals_file),&
                                     ' expected: ',TRIM(int2str((accepted_so_far_all(nn)-noise_count+1)*number_of_measurements))
                IF (nn.NE.number_of_chains) WRITE(*,*)'--------------------------------------------'
        END DO
        WRITE(*,*)'********************************************'
        WRITE(*,*)
ELSE
        WRITE(*,*)
        WRITE(*,*)'Lengths of output files are OK :-)'
END IF

END SUBROUTINE

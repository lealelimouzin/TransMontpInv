! ********************************************************************* !
! Program to convert a binary Markov chain file to text
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


PROGRAM chain2txt

IMPLICIT NONE
CHARACTER(LEN=30)                       :: file_name
INTEGER                                 :: file_type
INTEGER                                 :: noise_type,number_of_datasets
INTEGER                                 :: number_of_measurements,number_of_temperatures
INTEGER                                 :: save_residuals
INTEGER                                 :: nn,iter,ii,IOstatus
INTEGER                                 :: ncell
DOUBLE PRECISION                        :: dble1,dble2,dble3,dble4,dble5,dble6,dble7
DOUBLE PRECISION, DIMENSION(4)          :: dbles
INTEGER                                 :: int1,int2,int3
CHARACTER(LEN=*), PARAMETER             :: FMT51="(1I10,1I8)"
CHARACTER(LEN=*), PARAMETER             :: FMT50v="(4F22.10)"
CHARACTER(LEN=*), PARAMETER             :: FMT50s="(1I10,2I2,1I8)"
CHARACTER(LEN=*), PARAMETER             :: FMT50d="(1F22.10)"
CHARACTER(LEN=*), PARAMETER             :: FMT50i="(1I4)"
CHARACTER(LEN=*), PARAMETER             :: FMT54="(1I10,6F14.6)"
CHARACTER(LEN=*), PARAMETER             :: FMT55="(4I10,1F14.6)"
CHARACTER(LEN=3)                        :: adv

! Get file name from command line
CALL getarg(1,file_name)

! Open binary file and read file type
OPEN(UNIT=10,FILE=TRIM(file_name),&
     STATUS='OLD',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
READ(10)file_type
! Open text file
OPEN(UNIT=20,FILE=TRIM(file_name) // '.txt',&
     STATUS='REPLACE',FORM='FORMATTED',ACCESS='SEQUENTIAL')

! Read different entries depending on the file type
IOstatus=0
IF (file_type.EQ.1) THEN ! file of number of cells
        DO WHILE (IOstatus.EQ.0)
                READ(10,IOSTAT=IOstatus)iter,ncell
                IF (IOstatus.NE.0) EXIT
                WRITE(20,FMT51)iter,ncell
        END DO
ELSE IF (file_type.EQ.2) THEN ! file of chain steps
        READ(10)ncell
        DO nn=1,ncell
                READ(10)dble1,dble2,dble3,dble4
                WRITE(20,FMT50v)dble1,dble2,dble3,dble4
        END DO
        DO WHILE (IOstatus.EQ.0)
                READ(10,IOSTAT=IOstatus)iter,int1,int2,int3,dbles(1),dbles(2),dbles(3),dbles(4)
                IF (IOstatus.NE.0) EXIT
                WRITE(20,FMT50s,ADVANCE='NO')iter,int1,int2,int3
                DO nn=1,4
                        IF (nn.LT.4) adv="NO"
                        IF (nn.EQ.4) adv="YES"
                        IF (dbles(nn).NE.0) THEN
                                WRITE(20,FMT50d,ADVANCE=adv)dbles(nn)
                        ELSE
                                WRITE(20,FMT50i,ADVANCE=adv)0
                        END IF
                END DO
        END DO
ELSE IF (file_type.EQ.3) THEN ! file of noise parameters
        READ(10)noise_type
        IF (noise_type.EQ.4) THEN
                READ(10)number_of_measurements
                WRITE(20,*)number_of_measurements
                DO nn=1,number_of_measurements
                        READ(10)dble1
                        WRITE(20,*)dble1
                END DO
        ELSE
                READ(10)number_of_datasets
                DO WHILE (IOstatus.EQ.0)
                        DO nn=1,number_of_datasets
                                IF (noise_type.EQ.0.OR.noise_type.EQ.3.OR.noise_type.EQ.6) THEN
                                        READ(10,IOSTAT=IOstatus)iter,dble1
                                        IF (IOstatus.NE.0) EXIT
                                        WRITE(20,*)iter,dble1
                                ELSE IF (noise_type.EQ.1.OR.noise_type.EQ.2.OR.noise_type.EQ.5) THEN
                                        READ(10,IOSTAT=IOstatus)iter,dble1,dble2
                                        IF (IOstatus.NE.0) EXIT
                                        WRITE(20,*)iter,dble1,dble2
                                END IF
                        END DO
                        IF (IOstatus.NE.0) EXIT
                END DO
        END IF
ELSE IF (file_type.EQ.4) THEN ! file of misfits
        DO WHILE (IOstatus.EQ.0)
                READ(10,IOSTAT=IOstatus)iter,dble1
                IF (IOstatus.NE.0) EXIT
                WRITE(20,*)iter,dble1
        END DO
ELSE IF (file_type.EQ.5) THEN ! file of acceptance ratios
        DO WHILE (IOstatus.EQ.0)
                READ(10,IOSTAT=IOstatus)iter,dble1,dble2,dble3,dble4,dble5,dble6
                IF (IOstatus.NE.0) EXIT
                WRITE(20,FMT54)iter,dble1,dble2,dble3,dble4,dble5,dble6
        END DO
ELSE IF (file_type.EQ.6) THEN ! file of temperatures
        READ(10)number_of_temperatures
        WRITE(20,*)number_of_temperatures
        DO nn=1,number_of_temperatures
                READ(10)dble1
                WRITE(20,*)dble1
        END DO
        DO WHILE (IOstatus.EQ.0)
                READ(10,IOSTAT=IOstatus)iter,int1,int2,int3,dble7
                IF (IOstatus.NE.0) EXIT
                WRITE(20,FMT55)iter,int1,int2,int3,dble7
        END DO
ELSE IF (file_type.EQ.7) THEN ! file of residuals
        READ(10)save_residuals
        IF (save_residuals.EQ.1) THEN
                DO WHILE (IOstatus.EQ.0)
                        READ(10,IOSTAT=IOstatus)dble1
                        IF (IOstatus.NE.0) EXIT
                        WRITE(20,*)dble1
                END DO
        END IF
END IF

CLOSE(10)
CLOSE(20)

END PROGRAM
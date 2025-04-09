! ********************************************************************* !
! Program to calculate value posteriors from the procsamples.f90 outputs
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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This code reads the file of value posteriors produced by     !
! procsamples and writes text files of posterior distributions !
! at discrete [x y z] locations.                               !
!                                                              !
! Erica Galetti, March 2017                                    !
! erica.galetti@ed.ac.uk                                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

PROGRAM value_post
USE input_files, ONLY: postvalue_file,x_number_of_pixels,y_number_of_pixels,z_number_of_pixels,value_number_of_pixels
USE processing, ONLY: x_pixels,y_pixels,z_pixels,value_pixels
USE string_conversion
IMPLICIT NONE
INTEGER :: ii,jj,kk,ll,nn
DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE :: posterior_value
INTEGER :: number_of_points
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: points
DOUBLE PRECISION :: x_point,y_point,z_point
CHARACTER (LEN=30) :: output_file


OPEN(UNIT=10,FILE='value_3D_post.in',STATUS='OLD')
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,1)postvalue_file
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,1)output_file
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,*)number_of_points
ALLOCATE(points(number_of_points,3))
DO ii=1,number_of_points
        READ(10,*)points(ii,1),points(ii,2),points(ii,3)
END DO
1   FORMAT(a30)
CLOSE(10)

OPEN(UNIT=20,FILE=postvalue_file,STATUS='OLD')
READ(20,*)x_number_of_pixels
READ(20,*)y_number_of_pixels
READ(20,*)z_number_of_pixels
READ(20,*)value_number_of_pixels
ALLOCATE(posterior_value(z_number_of_pixels,x_number_of_pixels,y_number_of_pixels,value_number_of_pixels))
ALLOCATE(x_pixels(x_number_of_pixels))
ALLOCATE(y_pixels(y_number_of_pixels))
ALLOCATE(z_pixels(z_number_of_pixels))
ALLOCATE(value_pixels(value_number_of_pixels))
DO ii=1,x_number_of_pixels
        READ(20,*)x_pixels(ii)
END DO
DO jj=1,y_number_of_pixels
        READ(20,*)y_pixels(jj)
END DO
DO kk=1,z_number_of_pixels
        READ(20,*)z_pixels(kk)
END DO
DO ll=1,value_number_of_pixels
        READ(20,*)value_pixels(ll)
END DO
DO jj=1,y_number_of_pixels
        DO ii=1,x_number_of_pixels
                DO kk=1,z_number_of_pixels
                        DO ll=1,value_number_of_pixels
                                READ(20,*)posterior_value(kk,ii,jj,ll)
                        END DO
                END DO
        END DO
END DO
CLOSE(20)

DO nn=1,number_of_points

        OPEN(UNIT=54,FILE=TRIM(output_file) // '_point_' // TRIM(int2str(nn)) // '.out',STATUS='REPLACE')
        
        x_point=points(nn,1)
        y_point=points(nn,2)
        z_point=points(nn,3)
        
        ii=MINLOC(ABS(x_pixels-x_point),DIM=1)
        jj=MINLOC(ABS(y_pixels-y_point),DIM=1)
        kk=MINLOC(ABS(z_pixels-z_point),DIM=1)
        
        DO ll=1,value_number_of_pixels
                WRITE(54,*)value_pixels(ll),posterior_value(kk,ii,jj,ll)
        END DO
        
        CLOSE(54)
                
END DO

DEALLOCATE(x_pixels,y_pixels,z_pixels,value_pixels,posterior_value,points)

END PROGRAM value_post

! ********************************************************************* !
! Subroutine to write a grid file in vtx format in procsamples.f90
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


SUBROUTINE write_vtx(file_unit,map)

USE processing, ONLY: x_number_of_pixels,y_number_of_pixels,z_number_of_pixels,&
                      x_pixel_width,y_pixel_width,z_pixel_width
USE chains, ONLY: x_min,x_max,y_min,y_max,z_min,z_max
IMPLICIT NONE
INTEGER, INTENT(IN)                             :: file_unit
DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(IN)  :: map
INTEGER                                         :: dimensions
INTEGER                                         :: ii,jj,kk
INTEGER                                         :: i,j,k
INTEGER                                         :: v_number_of_pixels,h_number_of_pixels
DOUBLE PRECISION                                :: v_max,h_min,v_pixel_width,h_pixel_width

IF (x_number_of_pixels.EQ.1.OR.&
    y_number_of_pixels.EQ.1.OR.&
    z_number_of_pixels.EQ.1) THEN
        dimensions=2
ELSE
        dimensions=3
END IF

IF (dimensions.EQ.2) THEN
        
        IF (x_number_of_pixels.EQ.1) THEN ! yz plane
                
                v_number_of_pixels=z_number_of_pixels
                v_max=z_max
                v_pixel_width=z_pixel_width
                h_number_of_pixels=y_number_of_pixels
                h_min=y_min
                h_pixel_width=y_pixel_width
                
        ELSE IF (y_number_of_pixels.EQ.1) THEN ! xz plane
                
                v_number_of_pixels=z_number_of_pixels
                v_max=z_max
                v_pixel_width=z_pixel_width
                h_number_of_pixels=x_number_of_pixels
                h_min=x_min
                h_pixel_width=x_pixel_width
                
        ELSE IF (z_number_of_pixels.EQ.1) THEN ! xy plane
                
                v_number_of_pixels=y_number_of_pixels
                v_max=y_max
                v_pixel_width=y_pixel_width
                h_number_of_pixels=x_number_of_pixels
                h_min=x_min
                h_pixel_width=x_pixel_width
                
        END IF
        
        WRITE(file_unit,*)v_number_of_pixels,h_number_of_pixels
        WRITE(file_unit,'(3f14.8)')v_max,h_min
        WRITE(file_unit,'(3f14.8)')v_pixel_width,h_pixel_width
        
        WRITE(file_unit,'(1X)')
        
        DO i=0,h_number_of_pixels+1
                
                ii=i
                IF (i.EQ.0) ii=1
                IF (i.EQ.h_number_of_pixels+1) ii=h_number_of_pixels
                
                DO j=v_number_of_pixels+1,0,-1
                        
                        jj=j
                        IF (j.EQ.0) jj=1
                        IF (j.EQ.v_number_of_pixels+1) jj=v_number_of_pixels
                        
                        IF (x_number_of_pixels.EQ.1) WRITE(file_unit,*)map(jj,1,ii) ! yz plane
                        IF (y_number_of_pixels.EQ.1) WRITE(file_unit,*)map(jj,ii,1) ! xz plane
                        IF (z_number_of_pixels.EQ.1) WRITE(file_unit,*)map(1,ii,jj) ! xy plane
                        
                END DO
                
                WRITE(file_unit,'(1X)')
                
        END DO
        
        WRITE(file_unit,'(1X)')
        
ELSE IF (dimensions.EQ.3) THEN
        
        WRITE(file_unit,*)y_number_of_pixels,x_number_of_pixels,z_number_of_pixels
        WRITE(file_unit,'(3f14.8)')y_max,x_min,z_min
        WRITE(file_unit,'(3f14.8)')y_pixel_width,x_pixel_width,z_pixel_width
        WRITE(file_unit,'(1X)')
        
        DO k=0,z_number_of_pixels+1
                
                IF (dimensions.EQ.2) THEN
                        IF (k.EQ.0.OR.k.EQ.z_number_of_pixels+1) CYCLE
                END IF
                kk=k
                IF (k.EQ.0) kk=1
                IF (k.EQ.z_number_of_pixels+1) kk=z_number_of_pixels
                
                DO i=0,x_number_of_pixels+1
                        
                        ii=i
                        IF (i.EQ.0) ii=1
                        IF (i.EQ.x_number_of_pixels+1) ii=x_number_of_pixels
                        
                        DO j=y_number_of_pixels+1,0,-1
                                
                                jj=j
                                IF (j.EQ.0) jj=1
                                IF (j.EQ.y_number_of_pixels+1) jj=y_number_of_pixels
                                
                                WRITE(file_unit,*)map(kk,ii,jj)
                                
                        END DO
                        
                        WRITE(file_unit,'(1X)')
                        
                END DO
                
                WRITE(file_unit,'(1X)')
                
        END DO
        
END IF

END SUBROUTINE
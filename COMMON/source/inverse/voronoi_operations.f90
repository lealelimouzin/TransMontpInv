! ********************************************************************* !
! Module for Voronoi-tessellation related operations
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


MODULE voronoi_operations

CONTAINS

SUBROUTINE find_closest_nucleus(point,voronoi_model,number_of_cells,nucleus)
! Finds the closest nucleus of a Voronoi tessellation to a point.
!
! INPUT:        - point                 1D array of [x,y,z] coordinates of point
!               - voronoi model         2D array where each line contains the
!                                       [x,y,z] coordinates of each Voronoi nucleus
!               - number_of_cells       number of cells in the Voronoi tessellation
! OUTPUT:       - nucleus               number of Voronoi nucleus closest to point
!                                       (i.e., the point is located in this Voronoi
!                                       cell)

IMPLICIT NONE
DOUBLE PRECISION, DIMENSION(3), INTENT(IN)                      :: point
DOUBLE PRECISION, DIMENSION(number_of_cells,3), INTENT(IN)      :: voronoi_model
INTEGER, INTENT(IN)                                             :: number_of_cells
INTEGER, INTENT(OUT)                                            :: nucleus
DOUBLE PRECISION, DIMENSION(number_of_cells)                    :: distances

! Find distance from point to all Voronoi nuclei
distances=SQRT((voronoi_model(1:number_of_cells,1)-point(1))**2&
              +(voronoi_model(1:number_of_cells,2)-point(2))**2&
              +(voronoi_model(1:number_of_cells,3)-point(3))**2)

! Find nucleus at minimum distance
nucleus=MINLOC(distances,1)

END SUBROUTINE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE voronoi_to_2Dgrid(voronoi_model,number_of_cells,x_axis,y_axis,model2D)
! Creates a 2D grid of values from the Voronoi tessellation given as input.
!
! INPUT:        - voronoi model         2D array where each line contains the
!                                       [x,y] coordinates of each Voronoi nucleus
!                                       and their corresponding value
!               - number_of_cells       number of cells in the Voronoi tessellation
!               - x_axis                x-coordinates at which model values are to be evaluated
!               - y_axis                y-coordinates at which model values are to be evaluated
! OUTPUT:       - model2D               2D array of size [SIZE(x_axis) SIZE(y_axis)] where
!                                       each entry corresponds to the value of the Voronoi
!                                       cell which is closest to the [x,y] location

IMPLICIT NONE
DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN)                    :: voronoi_model
INTEGER, INTENT(IN)                                             :: number_of_cells
DOUBLE PRECISION, DIMENSION(:), INTENT(IN)                      :: x_axis,y_axis
DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE, INTENT(OUT)     :: model2D
INTEGER                                                         :: ii,jj
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE                   :: voronoi_temp
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE                     :: point
INTEGER                                                         :: nucleus

IF (ALLOCATED(model2D)) DEALLOCATE(model2D)
ALLOCATE(model2D(SIZE(y_axis),SIZE(x_axis)))

ALLOCATE(voronoi_temp(number_of_cells,4))
voronoi_temp(:,1)=voronoi_model(1:number_of_cells,1)
voronoi_temp(:,2)=0
voronoi_temp(:,3)=voronoi_model(1:number_of_cells,2)
voronoi_temp(:,4)=voronoi_model(1:number_of_cells,3)
ALLOCATE(point(3))

DO ii=1,SIZE(x_axis)
        DO jj=1,SIZE(y_axis)
                point=(/x_axis(ii),DBLE(0),y_axis(jj)/)
                CALL find_closest_nucleus(point,voronoi_temp(:,1:3),number_of_cells,nucleus)
                model2D(jj,ii)=voronoi_temp(nucleus,4)
        END DO
END DO

DEALLOCATE(voronoi_temp,point)

END SUBROUTINE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE voronoi_to_3Dgrid(voronoi_model,number_of_cells,x_axis,y_axis,z_axis,model3D)
! Creates a 3D grid of values from the Voronoi tessellation given as input.
!
! INPUT:        - voronoi model         2D array where each line contains the
!                                       [x,y,z] coordinates of each Voronoi nucleus
!                                       and their corresponding value
!               - number_of_cells       number of cells in the Voronoi tessellation
!               - x_axis                x-coordinates at which model values are to be evaluated
!               - y_axis                y-coordinates at which model values are to be evaluated
!               - z_axis                z-coordinates at which model values are to be evaluated
! OUTPUT:       - model3D               3D array of size [SIZE(x_axis) SIZE(y_axis) SIZE(z_axis)]
!                                       where each entry corresponds to the value of the Voronoi
!                                       cell which is closest to the [x,y,z] location

IMPLICIT NONE
DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN)                    :: voronoi_model
INTEGER, INTENT(IN)                                             :: number_of_cells
DOUBLE PRECISION, DIMENSION(:), INTENT(IN)                      :: x_axis,y_axis,z_axis
DOUBLE PRECISION, DIMENSION (:,:,:), ALLOCATABLE, INTENT(OUT)   :: model3D
INTEGER                                                         :: ii,jj,kk
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE                     :: point
INTEGER                                                         :: nucleus

IF (ALLOCATED(model3D)) DEALLOCATE(model3D)
ALLOCATE(model3D(SIZE(z_axis),SIZE(x_axis),SIZE(y_axis)))

ALLOCATE(point(3))

DO jj=1,SIZE(y_axis)
        DO ii=1,SIZE(x_axis)
                DO kk=1,SIZE(z_axis)
                        point=(/x_axis(ii),y_axis(jj),z_axis(kk)/)
                        CALL find_closest_nucleus(point,voronoi_model(1:number_of_cells,1:3),number_of_cells,nucleus)
                        model3D(kk,ii,jj)=voronoi_model(nucleus,4)
                END DO
        END DO
END DO

DEALLOCATE(point)

END SUBROUTINE

END MODULE
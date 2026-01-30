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


SUBROUTINE create_resistivity_grid2D(ncell,voronoi_model)

USE modelling_2_5D, ONLY: resistivity,x_centers,z_centers,x_number_of_cells,z_number_of_cells
USE voronoi_operations
IMPLICIT NONE
INTEGER, INTENT(IN)                                     :: ncell
DOUBLE PRECISION, DIMENSION(ncell,3), INTENT(IN)        :: voronoi_model
INTEGER                                                 :: ii,jj
DOUBLE PRECISION, DIMENSION(ncell,4)                    :: voronoi_temp
DOUBLE PRECISION, DIMENSION(3)                          :: point
INTEGER                                                 :: nucleus

voronoi_temp(:,1)=voronoi_model(1:ncell,1)
voronoi_temp(:,2)=0
voronoi_temp(:,3)=voronoi_model(1:ncell,2)
voronoi_temp(:,4)=voronoi_model(1:ncell,3)

DO jj=1,z_number_of_cells
        DO ii=1,x_number_of_cells
                point=(/x_centers(ii),DBLE(0),z_centers(jj)/)
                CALL find_closest_nucleus(point,voronoi_temp(:,1:3),ncell,nucleus)
                resistivity(ii,jj)=10**voronoi_temp(nucleus,4)
        END DO
END DO

END SUBROUTINE
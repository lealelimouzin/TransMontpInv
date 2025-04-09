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


SUBROUTINE interpolate_electrodes(electrode_list,Q_size,Q_indices,Q_values,Q_cpointer)

USE modelling_2_5D, ONLY: x_spacings,z_spacings,&
                          x_centers,z_centers,&
                          x_number_of_cells,z_number_of_cells
IMPLICIT NONE
DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN)                    :: electrode_list
INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT)                 :: Q_size
INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT)               :: Q_indices
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, INTENT(OUT)        :: Q_values
INTEGER, ALLOCATABLE, INTENT(OUT)                               :: Q_cpointer(:)
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE                   :: electrode_locations
INTEGER                                                         :: number_of_electrodes
INTEGER                                                         :: ii,jj,closest_index
INTEGER                                                         :: nn,row,column,entry_index
INTEGER, DIMENSION(2)                                           :: x_indices,z_indices
DOUBLE PRECISION                                                :: distance_left,distance_right
DOUBLE PRECISION                                                :: distance_top,distance_bottom
DOUBLE PRECISION                                                :: x_width,z_width
DOUBLE PRECISION                                                :: value
INTEGER, DIMENSION(:,:), ALLOCATABLE                            :: worki
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE                     :: workd
INTERFACE
        SUBROUTINE coordinate_to_sparse_column(M_size,M_indices,M_values,M_cpointer)
                INTEGER, DIMENSION(:), INTENT(IN)               :: M_size
                INTEGER, DIMENSION(:,:), INTENT(INOUT)          :: M_indices
                DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT)   :: M_values
                INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: M_cpointer
        END SUBROUTINE
END INTERFACE

IF (ALLOCATED(Q_size)) DEALLOCATE(Q_size)
IF (ALLOCATED(Q_indices)) DEALLOCATE(Q_indices)
IF (ALLOCATED(Q_values)) DEALLOCATE(Q_values)
IF (ALLOCATED(Q_cpointer)) DEALLOCATE(Q_cpointer)

number_of_electrodes=SIZE(electrode_list(:,1))
ALLOCATE(electrode_locations(number_of_electrodes,2))
electrode_locations=electrode_list

! Shift surface electrodes to the very top cell center
WHERE (electrode_locations(:,2).LT.MINVAL(z_centers))
        electrode_locations(:,2)=MINVAL(z_centers)
END WHERE

ALLOCATE(Q_size(2))
Q_size=(/ x_number_of_cells*z_number_of_cells,number_of_electrodes /)
ALLOCATE(Q_indices(4*number_of_electrodes,2))
ALLOCATE(Q_values(4*number_of_electrodes))

nn=0
DO ii=1,number_of_electrodes
        
        ! X AXIS
        ! Find cell center that is closest to the electrode
        closest_index=MINLOC(ABS(electrode_locations(ii,1)-x_centers),1)
        IF (electrode_locations(ii,1).GE.x_centers(closest_index)) THEN ! the closest cell center is on the left of or at the electrode location
                x_indices(1)=closest_index
                x_indices(2)=closest_index+1
        ELSE IF (electrode_locations(ii,1).LT.x_centers(closest_index)) THEN ! the closest cell center is on the right of the electrode location
                x_indices(1)=closest_index-1
                x_indices(2)=closest_index
        END IF
        ! Get distance from the electrode to the closest cell center left and right
        distance_left=electrode_locations(ii,1)-x_centers(x_indices(1))
        distance_right=x_centers(x_indices(2))-electrode_locations(ii,1)
        ! Get distance between the closest left and right cell centers
        x_width=x_centers(x_indices(2))-x_centers(x_indices(1))
        
        ! Z AXIS
        ! Find cell center that is closest to the electrode
        closest_index=MINLOC(ABS(electrode_locations(ii,2)-z_centers),1)
        IF (electrode_locations(ii,2).GE.z_centers(closest_index)) THEN ! the closest cell center is above or at the electrode location
                z_indices(1)=closest_index
                z_indices(2)=closest_index+1
        ELSE IF (electrode_locations(ii,2).LT.z_centers(closest_index)) THEN ! the closest cell center is below the electrode location
                z_indices(1)=closest_index-1
                z_indices(2)=closest_index
        END IF
        ! Get distance from the electrode to the closest cell center above and below
        distance_top=electrode_locations(ii,2)-z_centers(z_indices(1))
        distance_bottom=z_centers(z_indices(2))-electrode_locations(ii,2)
        ! Get distance between the closest top and bottom cell centers
        z_width=z_centers(z_indices(2))-z_centers(z_indices(1))
        
        ! Get nonzero elements
        ! Top-left
        value=(1-distance_left/x_width)*(1-distance_top/z_width)
        IF (value.NE.0) THEN
                nn=nn+1
                row=x_indices(1)
                column=z_indices(1)
                entry_index=x_number_of_cells*(column-1)+row
                Q_indices(nn,:)=(/ entry_index,ii /)
                Q_values(nn)=value
        END IF
        ! Top-right
        value=(1-distance_right/x_width)*(1-distance_top/z_width)
        IF (value.NE.0) THEN
                nn=nn+1
                row=x_indices(2)
                column=z_indices(1)
                entry_index=x_number_of_cells*(column-1)+row
                Q_indices(nn,:)=(/ entry_index,ii /)
                Q_values(nn)=value
        END IF
        ! Bottom-left
        value=(1-distance_left/x_width)*(1-distance_bottom/z_width)
        IF (value.NE.0) THEN
                nn=nn+1
                row=x_indices(1)
                column=z_indices(2)
                entry_index=x_number_of_cells*(column-1)+row
                Q_indices(nn,:)=(/ entry_index,ii /)
                Q_values(nn)=value
        END IF
        ! Bottom-right
        value=(1-distance_right/x_width)*(1-distance_bottom/z_width)
        IF (value.NE.0) THEN
                nn=nn+1
                row=x_indices(2)
                column=z_indices(2)
                entry_index=x_number_of_cells*(column-1)+row
                Q_indices(nn,:)=(/ entry_index,ii /)
                Q_values(nn)=value
        END IF
        
END DO

! Keep only nonzero elements
ALLOCATE(worki(nn,2),workd(nn))
worki=Q_indices(1:nn,:)
workd=Q_values(1:nn)
DEALLOCATE(Q_indices,Q_values)
ALLOCATE(Q_indices(nn,2))
ALLOCATE(Q_values(nn))
Q_indices=worki
Q_values=workd
DEALLOCATE(worki,workd)

CALL coordinate_to_sparse_column(Q_size,Q_indices,Q_values,Q_cpointer)

DEALLOCATE(electrode_locations)

END SUBROUTINE
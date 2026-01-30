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


SUBROUTINE get_data()

! This subroutine formats the acquisition data to make it compatible
! with the rest of the 2.5D modelling functions.
! From the list of quadrupoles, all unique pairs of current electrodes
! are extracted. Current and potential locations are then expressed in
! terms of their actual xz positions.
!
! Erica Galetti, November 2016
! erica.galetti@ed.ac.uk

USE m_mrgrnk
USE modelling_2_5D, ONLY: current_locations,potential_locations,current_electrode_numbers,number_of_sources
USE acquisition
IMPLICIT NONE
INTEGER, DIMENSION(:,:), ALLOCATABLE    :: current_electrodes
INTEGER, DIMENSION(:), ALLOCATABLE      :: source_numbers,source_order
INTEGER, ALLOCATABLE                    :: index_vector(:)
INTEGER                                 :: ii,jj
INTEGER, DIMENSION(1)                   :: nn
INTEGER                                 :: c1,c2,p1,p2
INTERFACE
        SUBROUTINE trim_sources(sources,index_vector,number_of_sources)
                INTEGER, DIMENSION(:,:), INTENT(IN)             :: sources
                INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: index_vector
                INTEGER, INTENT(OUT)                            :: number_of_sources
        END SUBROUTINE
END INTERFACE

IF (ALLOCATED(current_locations)) DEALLOCATE(current_locations)
IF (ALLOCATED(potential_locations)) DEALLOCATE(potential_locations)
IF (ALLOCATED(current_electrode_numbers)) DEALLOCATE(current_electrode_numbers)

! Keep only unique pairs of current electrodes
CALL trim_sources(quadrupoles(:,3:4),index_vector,number_of_sources)
ALLOCATE(current_electrodes(number_of_sources,2))
ALLOCATE(source_numbers(number_of_sources))
ALLOCATE(current_locations(number_of_sources,4))
ALLOCATE(potential_locations(number_of_measurements,4))
ALLOCATE(current_electrode_numbers(number_of_measurements))

 current_electrodes=quadrupoles(index_vector,3:4)
 source_numbers=(/ (ii,ii=1,number_of_sources) /)

! Loop over the unique current dipoles to get their locations
DO ii=1,number_of_sources
        c1=current_electrodes(ii,1)
        c2=current_electrodes(ii,2)
        current_locations(ii,1:2)=electrodes(c1,:)
        current_locations(ii,3:4)=electrodes(c2,:)
END DO

! Loop over the readings to get the locations of the 4 electrodes
! used for each reading
DO ii=1,number_of_measurements
        p1=quadrupoles(ii,1)
        p2=quadrupoles(ii,2)
        potential_locations(ii,1:2)=electrodes(p1,:)
        potential_locations(ii,3:4)=electrodes(p2,:)
        c1=quadrupoles(ii,3)
        c2=quadrupoles(ii,4)
        nn=PACK(source_numbers,current_electrodes(:,1).EQ.c1.AND.current_electrodes(:,2).EQ.c2)
        current_electrode_numbers(ii)=nn(1)
END DO

! Sort the data in increasing order of current electrode number
ALLOCATE(source_order(SIZE(current_electrode_numbers)))
CALL I_MRGRNK(current_electrode_numbers,source_order)
 current_electrode_numbers=current_electrode_numbers(source_order)
 potential_locations=potential_locations(source_order,:)
 observed_data=observed_data(source_order)
 data_noise=data_noise(source_order)
 dataset_associations=dataset_associations(source_order)
 
DEALLOCATE(index_vector,current_electrodes,source_numbers,source_order)

END SUBROUTINE


SUBROUTINE trim_sources(sources,index_vector,number_of_sources)

! This subroutine trims the list of current electrodes keeping only
! unique pairs. This ensures that the source term on the right-
! hand-side of the linear system of equations for the electrical
! potential has as few columns as possible.
! 
! The algorithm for this subroutine was inspired by 
! http://stackoverflow.com/questions/14137610/remove-repeated-elements-on-an-2d-array-in-fortran-90
!
! Erica Galetti, November 2016
! erica.galetti@ed.ac.uk

IMPLICIT NONE
INTEGER, DIMENSION(:,:), INTENT(IN)             :: sources
INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: index_vector
INTEGER, INTENT(OUT)                            :: number_of_sources
LOGICAL, DIMENSION(:), ALLOCATABLE              :: mask
INTEGER                                         :: ii,length_sources
INTEGER, DIMENSION(:), ALLOCATABLE              :: sources_index

IF (ALLOCATED(index_vector)) DEALLOCATE(index_vector)

length_sources=SIZE(sources(:,1))

ALLOCATE(mask(length_sources))
mask=.TRUE.
ALLOCATE(sources_index(length_sources))
sources_index=(/ (ii,ii=1,length_sources) /)

DO ii=length_sources,2,-1
        mask(ii)=.NOT.(ANY(sources(1:ii-1,1).EQ.sources(ii,1).AND.sources(1:ii-1,2).EQ.sources(ii,2)))
END DO

number_of_sources=SIZE(PACK(sources_index,mask))
ALLOCATE(index_vector(number_of_sources))
index_vector=PACK(sources_index,mask)

DEALLOCATE(mask,sources_index)

END SUBROUTINE
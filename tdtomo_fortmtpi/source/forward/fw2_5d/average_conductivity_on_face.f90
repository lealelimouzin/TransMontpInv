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


SUBROUTINE average_conductivity_on_face()

USE modelling_2_5D, ONLY: x_spacings,z_spacings,x_number_of_cells,z_number_of_cells,&
                          cell_areas,Nx,Nz,nax,naz,GRDax,GRDaz,resistivity,&
                          S_size,S_indices,S_values,resistivity_on_x_face,&
                          resistivity_on_z_face,x_face_size,z_face_size,&
                          x_rows,x_columns,z_rows,z_columns
IMPLICIT NONE
INTEGER, DIMENSION((x_number_of_cells-1)*z_number_of_cells)             :: lx
INTEGER, DIMENSION(x_number_of_cells*(z_number_of_cells-1))             :: lz
DOUBLE PRECISION, DIMENSION((x_number_of_cells-1)*z_number_of_cells)    :: kx
DOUBLE PRECISION, DIMENSION(x_number_of_cells*(z_number_of_cells-1))    :: kz
INTEGER                                                                 :: x_nonzero_length,z_nonzero_length
INTEGER                                                                 :: ini,fin

! Get average resistivity on x-faces
x_nonzero_length=x_number_of_cells-1
z_nonzero_length=z_number_of_cells

resistivity_on_x_face=((cell_areas(x_rows,x_columns)*resistivity(x_rows,x_columns))&
                      +(cell_areas(x_rows-1,x_columns)*resistivity(x_rows-1,x_columns)))/2

resistivity_on_x_face=resistivity_on_x_face/x_face_size

lx=RESHAPE(GRDax,(/x_nonzero_length*z_nonzero_length/))
kx=RESHAPE(resistivity_on_x_face,(/x_nonzero_length*z_nonzero_length/))

! Get average resistivity on z-faces
x_nonzero_length=x_number_of_cells
z_nonzero_length=z_number_of_cells-1

resistivity_on_z_face=((cell_areas(z_rows,z_columns-1)*resistivity(z_rows,z_columns-1))&
                      +(cell_areas(z_rows,z_columns)*resistivity(z_rows,z_columns)))/2

resistivity_on_z_face=resistivity_on_z_face/z_face_size

lz=RESHAPE(GRDaz,(/x_nonzero_length*z_nonzero_length/))
kz=RESHAPE(resistivity_on_z_face,(/x_nonzero_length*z_nonzero_length/))

! Create sparse matrix
ini=1
fin=nax
S_indices(ini:fin)=lx
S_values(ini:fin)=kx
ini=fin+1
fin=ini-1+naz
S_indices(ini:fin)=lz+nax
S_values(ini:fin)=kz
S_values=1/S_values

END SUBROUTINE
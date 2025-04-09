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


SUBROUTINE read_grid2D()

! Read the grid parameters from input file forward.in

USE modelling_2_5D
IMPLICIT NONE
INTEGER         :: ii,jj,kk

OPEN(UNIT=10,FILE=fw2_5d_file,STATUS='old')

READ(10,*)number_of_coefficients

READ(10,*)x_number_of_edges,z_number_of_edges
ALLOCATE(x_edges(x_number_of_edges),z_edges(z_number_of_edges))
DO ii=1,x_number_of_edges
        READ(10,*)x_edges(ii)
END DO
DO kk=1,z_number_of_edges
        READ(10,*)z_edges(kk)
END DO

CLOSE(10)

x_number_of_cells=x_number_of_edges-1
z_number_of_cells=z_number_of_edges-1

! Spacing between grid lines in the x and z direction
ALLOCATE(x_spacings(x_number_of_cells),z_spacings(z_number_of_cells))
x_spacings=x_edges(2:x_number_of_edges)-x_edges(1:x_number_of_edges-1)
z_spacings=z_edges(2:z_number_of_edges)-z_edges(1:z_number_of_edges-1)

! Location of cell centers in the x and z direction
ALLOCATE(x_centers(x_number_of_cells),z_centers(z_number_of_cells))
x_centers=x_edges(1:x_number_of_edges-1)+x_spacings/2
z_centers=z_edges(1:z_number_of_edges-1)+z_spacings/2

! Allocate arrays of cell areas and resistivity in each cell
ALLOCATE(cell_areas(x_number_of_cells,z_number_of_cells))
ALLOCATE(resistivity(x_number_of_cells,z_number_of_cells))
ALLOCATE(Dx(x_number_of_cells,z_number_of_cells),Dz(x_number_of_cells,z_number_of_cells))
DO ii=1,x_number_of_cells
        DO kk=1,z_number_of_cells
                cell_areas(ii,kk)=x_spacings(ii)*z_spacings(kk)
                Dx(ii,kk)=x_spacings(ii)
                Dz(ii,kk)=z_spacings(kk)
        END DO
END DO

! This part is taken from FW2_5D, and it sets up the parameters needed
! to calculate Divergence and Gradient operators.
Nx=x_number_of_cells-2
Nz=z_number_of_cells-2
np=(Nx+2)*(Nz+2)
nax=(Nx+1)*(Nz+2)
naz=(Nx+2)*(Nz+1)
ALLOCATE(GRDp(Nx+2,Nz+2))
ALLOCATE(GRDax(Nx+1,Nz+2))
ALLOCATE(GRDaz(Nx+2,Nz+1))
GRDp=RESHAPE((/ (ii,ii=1,np) /),(/Nx+2,Nz+2/))
GRDax=RESHAPE((/ (ii,ii=1,nax) /),(/(Nx+1),(Nz+2)/))
GRDaz=RESHAPE((/ (ii,ii=1,naz) /),(/(Nx+2),(Nz+1)/))


! This part allocates the arrays necessary to solve the forward problem.
! Allocating all of these arrays in advance ensures fast computation.

! S is a diagonal matrix of average face conductivities. It is used in
! subroutines average_conductivity_on_face and calculate_DS.
S_size=nax+naz
ALLOCATE(S_indices(nax+naz))
ALLOCATE(S_values(nax+naz))
! A is a 5-diagonal matrix that is obtained from the matrix product DSG having
! 'number_of_nonzeros' nonzero elements. It is calculated by subroutine
! calculate_A and subsequently used in subroutine calculate_data.
number_of_nonzeros=5*x_number_of_cells*z_number_of_cells-2*x_number_of_cells-2*z_number_of_cells
ALLOCATE(A_size(2))
A_size=(/x_number_of_cells*z_number_of_cells,x_number_of_cells*z_number_of_cells/)
ALLOCATE(A_indices(number_of_nonzeros,2))
ALLOCATE(A_values(number_of_nonzeros))
ALLOCATE(A_cpointer(x_number_of_cells*z_number_of_cells+1))
! C is a diagonal matrix of conductivities. It is used in subroutine calculate_data.
ALLOCATE(C_indices(x_number_of_cells*z_number_of_cells))
ALLOCATE(C_values(x_number_of_cells*z_number_of_cells))
 C_size=x_number_of_cells*z_number_of_cells
 C_indices=(/ (ii,ii=1,x_number_of_cells*z_number_of_cells) /)
! L is a sparse matrix of the same size as A which iteratively calculated in subroutine 
! calculate_data for each k and g parameter.
ALLOCATE(L_size(2))
L_size=(/x_number_of_cells*z_number_of_cells,x_number_of_cells*z_number_of_cells/)
ALLOCATE(L_indices(number_of_nonzeros,2))
ALLOCATE(L_values(number_of_nonzeros))
! These arrays are used by subroutine average_conductivity_on_face.
ALLOCATE(x_rows(x_number_of_cells-1),x_columns(z_number_of_cells))
x_rows=(/ (ii,ii=2,Nx+2) /)
x_columns=(/ (ii,ii=1,Nz+2) /)
ALLOCATE(z_rows(x_number_of_cells),z_columns(z_number_of_cells-1))
z_rows=(/ (ii,ii=1,Nx+2) /)
z_columns=(/ (ii,ii=2,Nz+2) /)
ALLOCATE(x_face_size(x_number_of_cells-1,z_number_of_cells))
x_face_size=(cell_areas(x_rows,x_columns)+cell_areas(x_rows-1,x_columns))/2
ALLOCATE(z_face_size(x_number_of_cells,z_number_of_cells-1))
z_face_size=(cell_areas(z_rows,z_columns-1)+cell_areas(z_rows,z_columns))/2
ALLOCATE(resistivity_on_x_face(x_number_of_cells-1,z_number_of_cells))
ALLOCATE(resistivity_on_z_face(x_number_of_cells,z_number_of_cells-1))


END SUBROUTINE
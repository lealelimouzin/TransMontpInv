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


SUBROUTINE calculate_A_slow(DS_size,DS_indices,DS_cpointer,DS_values,A_size,A_indices,A_cpointer,A_values)

USE modelling_2_5D, ONLY: G_size,G_indices,G_cpointer,G_values,x_number_of_cells,z_number_of_cells,nax,naz
IMPLICIT NONE
INTEGER, DIMENSION(:), INTENT(IN)                                :: DS_size,DS_cpointer
INTEGER, DIMENSION(:,:), INTENT(IN)                                :: DS_indices
DOUBLE PRECISION, DIMENSION(:), INTENT(IN)                        :: DS_values
INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT)                        :: A_size,A_cpointer
INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT)                :: A_indices
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, INTENT(OUT)        :: A_values
INTEGER                                                                :: ii,jj,nn,row,column,number_of_nonzeros,entry_number,nx,nz,cpoint
DOUBLE PRECISION, DIMENSION(1)                                        :: value_0,value_1u,value_1l,value_nxu,value_nxl

! To calculate A=-D*S*G we can exploit the fact that we know its structure: it is a sparse matrix with nonzero
! entries along diagonal 0 (the main diagonal), along diagonals +1 and -1 (the first upper and lower diagonals),
! and along diagonals +nx and -nx. Since we also know the structure of D*S and G, we can easily work out the value
! of each nonzero entry of A.

number_of_nonzeros=5*x_number_of_cells*z_number_of_cells-2*x_number_of_cells-2*z_number_of_cells
ALLOCATE(A_size(2))
ALLOCATE(A_indices(number_of_nonzeros,2))
ALLOCATE(A_values(number_of_nonzeros))
ALLOCATE(A_cpointer(x_number_of_cells*z_number_of_cells+1))

A_size=[x_number_of_cells*z_number_of_cells,x_number_of_cells*z_number_of_cells]

nx=x_number_of_cells
nz=z_number_of_cells

! Loop over the the diagonal entries
 nn=0
 entry_number=0
 cpoint=0
DO jj=1,z_number_of_cells
        DO ii=1,x_number_of_cells
                nn=nn+1
                value_0=0
                value_1u=0
                value_1l=0
                value_nxu=0
                value_nxl=0
                IF (jj.EQ.1) THEN
                        IF (ii.EQ.1) THEN
                                value_0=value_0+PACK(DS_values,DS_indices(:,1).EQ.nn.AND.DS_indices(:,2).EQ.nn)*PACK(G_values,G_indices(:,1).EQ.nn.AND.G_indices(:,2).EQ.nn)
                                value_0=value_0+PACK(DS_values,DS_indices(:,1).EQ.nn.AND.DS_indices(:,2).EQ.nn+nax)*PACK(G_values,G_indices(:,1).EQ.nn+nax.AND.G_indices(:,2).EQ.nn)
                                !value_1u=value_1u
                                value_1l=value_1l+PACK(DS_values,DS_indices(:,1).EQ.nn+1.AND.DS_indices(:,2).EQ.nn)*PACK(G_values,G_indices(:,1).EQ.nn.AND.G_indices(:,2).EQ.nn)
                                !value_nxu=value_nxu
                                value_nxl=value_nxl+PACK(DS_values,DS_indices(:,1).EQ.nn+nx.AND.DS_indices(:,2).EQ.nn+nax)*PACK(G_values,G_indices(:,1).EQ.nn+nax.AND.G_indices(:,2).EQ.nn)
                        ELSE IF (ii.GT.1.AND.ii.LT.x_number_of_cells) THEN
                                value_0=value_0+PACK(DS_values,DS_indices(:,1).EQ.nn.AND.DS_indices(:,2).EQ.nn-1)*PACK(G_values,G_indices(:,1).EQ.nn-1.AND.G_indices(:,2).EQ.nn)
                                value_0=value_0+PACK(DS_values,DS_indices(:,1).EQ.nn.AND.DS_indices(:,2).EQ.nn)*PACK(G_values,G_indices(:,1).EQ.nn.AND.G_indices(:,2).EQ.nn)
                                value_0=value_0+PACK(DS_values,DS_indices(:,1).EQ.nn.AND.DS_indices(:,2).EQ.nn+nax)*PACK(G_values,G_indices(:,1).EQ.nn+nax.AND.G_indices(:,2).EQ.nn)
                                value_1u=value_1u+PACK(DS_values,DS_indices(:,1).EQ.nn-1.AND.DS_indices(:,2).EQ.nn-1)*PACK(G_values,G_indices(:,1).EQ.nn-1.AND.G_indices(:,2).EQ.nn)
                                value_1l=value_1l+PACK(DS_values,DS_indices(:,1).EQ.nn+1.AND.DS_indices(:,2).EQ.nn)*PACK(G_values,G_indices(:,1).EQ.nn.AND.G_indices(:,2).EQ.nn)
                                !value_nxu=value_nxu
                                value_nxl=value_nxl+PACK(DS_values,DS_indices(:,1).EQ.nn+nx.AND.DS_indices(:,2).EQ.nn+nax)*PACK(G_values,G_indices(:,1).EQ.nn+nax.AND.G_indices(:,2).EQ.nn)
                        ELSE IF (ii.EQ.x_number_of_cells) THEN
                                value_0=value_0+PACK(DS_values,DS_indices(:,1).EQ.nn.AND.DS_indices(:,2).EQ.nn-1)*PACK(G_values,G_indices(:,1).EQ.nn-1.AND.G_indices(:,2).EQ.nn)
                                value_0=value_0+PACK(DS_values,DS_indices(:,1).EQ.nn.AND.DS_indices(:,2).EQ.nn+nax)*PACK(G_values,G_indices(:,1).EQ.nn+nax.AND.G_indices(:,2).EQ.nn)
                                value_1u=value_1u+PACK(DS_values,DS_indices(:,1).EQ.nn-1.AND.DS_indices(:,2).EQ.nn-1)*PACK(G_values,G_indices(:,1).EQ.nn-1.AND.G_indices(:,2).EQ.nn)
                                !value_1l=value_1l
                                !value_nxu=value_nxu
                                value_nxl=value_nxl+PACK(DS_values,DS_indices(:,1).EQ.nn+nx.AND.DS_indices(:,2).EQ.nn+nax)*PACK(G_values,G_indices(:,1).EQ.nn+nax.AND.G_indices(:,2).EQ.nn)
                        END IF
                ELSE IF (jj.GT.1.AND.jj.LT.z_number_of_cells) THEN
                        IF (ii.EQ.1) THEN
                                value_0=value_0+PACK(DS_values,DS_indices(:,1).EQ.nn.AND.DS_indices(:,2).EQ.nn-jj+1)*PACK(G_values,G_indices(:,1).EQ.nn-jj+1.AND.G_indices(:,2).EQ.nn)
                                value_0=value_0+PACK(DS_values,DS_indices(:,1).EQ.nn.AND.DS_indices(:,2).EQ.nn-nx+nax)*PACK(G_values,G_indices(:,1).EQ.nn-nx+nax.AND.G_indices(:,2).EQ.nn)
                                value_0=value_0+PACK(DS_values,DS_indices(:,1).EQ.nn.AND.DS_indices(:,2).EQ.nn+nax)*PACK(G_values,G_indices(:,1).EQ.nn+nax.AND.G_indices(:,2).EQ.nn)
                                !value_1u=value_1u
                                value_1l=value_1l+PACK(DS_values,DS_indices(:,1).EQ.nn+1.AND.DS_indices(:,2).EQ.nn-jj+1)*PACK(G_values,G_indices(:,1).EQ.nn-jj+1.AND.G_indices(:,2).EQ.nn)
                                value_nxu=value_nxu+PACK(DS_values,DS_indices(:,1).EQ.nn-nx.AND.DS_indices(:,2).EQ.nn-nx+nax)*PACK(G_values,G_indices(:,1).EQ.nn-nx+nax.AND.G_indices(:,2).EQ.nn)
                                value_nxl=value_nxl+PACK(DS_values,DS_indices(:,1).EQ.nn+nx.AND.DS_indices(:,2).EQ.nn+nax)*PACK(G_values,G_indices(:,1).EQ.nn+nax.AND.G_indices(:,2).EQ.nn)
                        ELSE IF (ii.GT.1.AND.ii.LT.x_number_of_cells) THEN
                                value_0=value_0+PACK(DS_values,DS_indices(:,1).EQ.nn.AND.DS_indices(:,2).EQ.nn-jj)*PACK(G_values,G_indices(:,1).EQ.nn-jj.AND.G_indices(:,2).EQ.nn)
                                value_0=value_0+PACK(DS_values,DS_indices(:,1).EQ.nn.AND.DS_indices(:,2).EQ.nn-jj+1)*PACK(G_values,G_indices(:,1).EQ.nn-jj+1.AND.G_indices(:,2).EQ.nn)
                                value_0=value_0+PACK(DS_values,DS_indices(:,1).EQ.nn.AND.DS_indices(:,2).EQ.nn-nx+nax)*PACK(G_values,G_indices(:,1).EQ.nn-nx+nax.AND.G_indices(:,2).EQ.nn)
                                value_0=value_0+PACK(DS_values,DS_indices(:,1).EQ.nn.AND.DS_indices(:,2).EQ.nn+nax)*PACK(G_values,G_indices(:,1).EQ.nn+nax.AND.G_indices(:,2).EQ.nn)
                                value_1u=value_1u+PACK(DS_values,DS_indices(:,1).EQ.nn-1.AND.DS_indices(:,2).EQ.nn-jj)*PACK(G_values,G_indices(:,1).EQ.nn-jj.AND.G_indices(:,2).EQ.nn)
                                value_1l=value_1l+PACK(DS_values,DS_indices(:,1).EQ.nn+1.AND.DS_indices(:,2).EQ.nn-jj+1)*PACK(G_values,G_indices(:,1).EQ.nn-jj+1.AND.G_indices(:,2).EQ.nn)
                                value_nxu=value_nxu+PACK(DS_values,DS_indices(:,1).EQ.nn-nx.AND.DS_indices(:,2).EQ.nn-nx+nax)*PACK(G_values,G_indices(:,1).EQ.nn-nx+nax.AND.G_indices(:,2).EQ.nn)
                                value_nxl=value_nxl+PACK(DS_values,DS_indices(:,1).EQ.nn+nx.AND.DS_indices(:,2).EQ.nn+nax)*PACK(G_values,G_indices(:,1).EQ.nn+nax.AND.G_indices(:,2).EQ.nn)
                        ELSE IF (ii.EQ.x_number_of_cells) THEN
                                value_0=value_0+PACK(DS_values,DS_indices(:,1).EQ.nn.AND.DS_indices(:,2).EQ.nn-jj)*PACK(G_values,G_indices(:,1).EQ.nn-jj.AND.G_indices(:,2).EQ.nn)
                                value_0=value_0+PACK(DS_values,DS_indices(:,1).EQ.nn.AND.DS_indices(:,2).EQ.nn-nx+nax)*PACK(G_values,G_indices(:,1).EQ.nn-nx+nax.AND.G_indices(:,2).EQ.nn)
                                value_0=value_0+PACK(DS_values,DS_indices(:,1).EQ.nn.AND.DS_indices(:,2).EQ.nn+nax)*PACK(G_values,G_indices(:,1).EQ.nn+nax.AND.G_indices(:,2).EQ.nn)
                                value_1u=value_1u+PACK(DS_values,DS_indices(:,1).EQ.nn-1.AND.DS_indices(:,2).EQ.nn-jj)*PACK(G_values,G_indices(:,1).EQ.nn-jj.AND.G_indices(:,2).EQ.nn)
                                !value_1l=value_1l
                                value_nxu=value_nxu+PACK(DS_values,DS_indices(:,1).EQ.nn-nx.AND.DS_indices(:,2).EQ.nn-nx+nax)*PACK(G_values,G_indices(:,1).EQ.nn-nx+nax.AND.G_indices(:,2).EQ.nn)
                                value_nxl=value_nxl+PACK(DS_values,DS_indices(:,1).EQ.nn+nx.AND.DS_indices(:,2).EQ.nn+nax)*PACK(G_values,G_indices(:,1).EQ.nn+nax.AND.G_indices(:,2).EQ.nn)
                        END IF
                ELSE IF (jj.EQ.z_number_of_cells) THEN
                        IF (ii.EQ.1) THEN
                                value_0=value_0+PACK(DS_values,DS_indices(:,1).EQ.nn.AND.DS_indices(:,2).EQ.nn-jj+1)*PACK(G_values,G_indices(:,1).EQ.nn-jj+1.AND.G_indices(:,2).EQ.nn)
                                value_0=value_0+PACK(DS_values,DS_indices(:,1).EQ.nn.AND.DS_indices(:,2).EQ.nn-nx+nax)*PACK(G_values,G_indices(:,1).EQ.nn-nx+nax.AND.G_indices(:,2).EQ.nn)
                                !value_1u=value_1u
                                value_1l=value_1l+PACK(DS_values,DS_indices(:,1).EQ.nn+1.AND.DS_indices(:,2).EQ.nn-jj+1)*PACK(G_values,G_indices(:,1).EQ.nn-jj+1.AND.G_indices(:,2).EQ.nn)
                                value_nxu=value_nxu+PACK(DS_values,DS_indices(:,1).EQ.nn-nx.AND.DS_indices(:,2).EQ.nn-nx+nax)*PACK(G_values,G_indices(:,1).EQ.nn-nx+nax.AND.G_indices(:,2).EQ.nn)
                                !value_nxl=value_nxl
                        ELSE IF (ii.GT.1.AND.ii.LT.x_number_of_cells) THEN
                                value_0=value_0+PACK(DS_values,DS_indices(:,1).EQ.nn.AND.DS_indices(:,2).EQ.nn-jj)*PACK(G_values,G_indices(:,1).EQ.nn-jj.AND.G_indices(:,2).EQ.nn)
                                value_0=value_0+PACK(DS_values,DS_indices(:,1).EQ.nn.AND.DS_indices(:,2).EQ.nn-jj+1)*PACK(G_values,G_indices(:,1).EQ.nn-jj+1.AND.G_indices(:,2).EQ.nn)
                                value_0=value_0+PACK(DS_values,DS_indices(:,1).EQ.nn.AND.DS_indices(:,2).EQ.nn-nx+nax)*PACK(G_values,G_indices(:,1).EQ.nn-nx+nax.AND.G_indices(:,2).EQ.nn)
                                value_1u=value_1u+PACK(DS_values,DS_indices(:,1).EQ.nn-1.AND.DS_indices(:,2).EQ.nn-jj)*PACK(G_values,G_indices(:,1).EQ.nn-jj.AND.G_indices(:,2).EQ.nn)
                                value_1l=value_1l+PACK(DS_values,DS_indices(:,1).EQ.nn+1.AND.DS_indices(:,2).EQ.nn-jj+1)*PACK(G_values,G_indices(:,1).EQ.nn-jj+1.AND.G_indices(:,2).EQ.nn)
                                value_nxu=value_nxu+PACK(DS_values,DS_indices(:,1).EQ.nn-nx.AND.DS_indices(:,2).EQ.nn-nx+nax)*PACK(G_values,G_indices(:,1).EQ.nn-nx+nax.AND.G_indices(:,2).EQ.nn)
                                !value_nxl=value_nxl
                        ELSE IF (ii.EQ.x_number_of_cells) THEN
                                value_0=value_0+PACK(DS_values,DS_indices(:,1).EQ.nn.AND.DS_indices(:,2).EQ.nn-jj)*PACK(G_values,G_indices(:,1).EQ.nn-jj.AND.G_indices(:,2).EQ.nn)
                                value_0=value_0+PACK(DS_values,DS_indices(:,1).EQ.nn.AND.DS_indices(:,2).EQ.nn-nx+nax)*PACK(G_values,G_indices(:,1).EQ.nn-nx+nax.AND.G_indices(:,2).EQ.nn)
                                value_1u=value_1u+PACK(DS_values,DS_indices(:,1).EQ.nn-1.AND.DS_indices(:,2).EQ.nn-jj)*PACK(G_values,G_indices(:,1).EQ.nn-jj.AND.G_indices(:,2).EQ.nn) 
                                !value_1l=value_1l
                                value_nxu=value_nxu+PACK(DS_values,DS_indices(:,1).EQ.nn-nx.AND.DS_indices(:,2).EQ.nn-nx+nax)*PACK(G_values,G_indices(:,1).EQ.nn-nx+nax.AND.G_indices(:,2).EQ.nn)
                                !value_nxl=value_nxl
                        END IF
                END IF
                A_cpointer(nn)=cpoint+1
                IF (value_nxu(1).NE.0) THEN
                        entry_number=entry_number+1
                        cpoint=cpoint+1
                        A_indices(entry_number,:)=[nn-nx,nn]
                        A_values(entry_number)=-value_nxu(1)
                END IF
                IF (value_1u(1).NE.0) THEN
                        entry_number=entry_number+1
                        cpoint=cpoint+1
                        A_indices(entry_number,:)=[nn-1,nn]
                        A_values(entry_number)=-value_1u(1)
                END IF
                IF (value_0(1).NE.0) THEN
                        entry_number=entry_number+1
                        cpoint=cpoint+1
                        A_indices(entry_number,:)=[nn,nn]
                        A_values(entry_number)=-value_0(1)
                END IF
                IF (value_1l(1).NE.0) THEN
                        entry_number=entry_number+1
                        cpoint=cpoint+1
                        A_indices(entry_number,:)=[nn+1,nn]
                        A_values(entry_number)=-value_1l(1)
                END IF
                IF (value_nxl(1).NE.0) THEN
                        entry_number=entry_number+1
                        cpoint=cpoint+1
                        A_indices(entry_number,:)=[nn+nx,nn]
                        A_values(entry_number)=-value_nxl(1)
                END IF
        END DO
END DO
A_cpointer(nn+1)=cpoint+1

END SUBROUTINE
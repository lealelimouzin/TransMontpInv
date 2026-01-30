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


SUBROUTINE calculate_gradient2D()

USE modelling_2_5D
IMPLICIT NONE
INTEGER, DIMENSION(:), ALLOCATABLE              :: lx,jx
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE     :: kx
INTEGER, DIMENSION(:), ALLOCATABLE              :: lz,jz
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE     :: kz
INTEGER                                         :: x_nonzero_length,z_nonzero_length
INTEGER                                         :: ini,fin,nn
INTERFACE
        SUBROUTINE coordinate_to_sparse_column(M_size,M_indices,M_values,M_cpointer)
                INTEGER, DIMENSION(:), INTENT(IN)               :: M_size
                INTEGER, DIMENSION(:,:), INTENT(INOUT)          :: M_indices
                DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT)   :: M_values
                INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: M_cpointer
        END SUBROUTINE
END INTERFACE

IF (ALLOCATED(G_size)) DEALLOCATE(G_size)
IF (ALLOCATED(G_indices)) DEALLOCATE(G_indices)
IF (ALLOCATED(G_values)) DEALLOCATE(G_values)

x_nonzero_length=2*nax
z_nonzero_length=2*naz

ALLOCATE(lx(x_nonzero_length),jx(x_nonzero_length))
ALLOCATE(kx(x_nonzero_length))
ALLOCATE(lz(z_nonzero_length),jz(z_nonzero_length))
ALLOCATE(kz(z_nonzero_length))
ALLOCATE(G_size(2))
ALLOCATE(G_indices(x_nonzero_length+z_nonzero_length,2))
ALLOCATE(G_values(x_nonzero_length+z_nonzero_length))

! Calculate d/dx
! Entries (l,j,k)
ini=1
fin=nax
lx(ini:fin)=RESHAPE(GRDax,(/nax/))
jx(ini:fin)=RESHAPE(GRDp(1:Nx+1,:),(/nax/))
kx(ini:fin)=RESHAPE(-2/(Dx(1:Nx+1,:)+Dx(2:Nx+2,:)),(/nax/))
! Entries (l+1,j,k)
ini=fin+1
fin=ini-1+nax
lx(ini:fin)=RESHAPE(GRDax,(/nax/))
jx(ini:fin)=RESHAPE(GRDp(2:Nx+2,:),(/nax/))
kx(ini:fin)=RESHAPE(2/(Dx(1:Nx+1,:)+Dx(2:Nx+2,:)),(/nax/))
! Calculate d/dz
! Entries (l,j,k)
ini=1
fin=naz
lz(ini:fin)=RESHAPE(GRDaz,(/naz/))
jz(ini:fin)=RESHAPE(GRDp(:,1:Nz+1),(/naz/))
kz(ini:fin)=RESHAPE(-2/(Dz(:,1:Nz+1)+Dz(:,2:Nz+2)),(/naz/))
! Entries (l,j,k+1)
ini=fin+1
fin=ini-1+naz
lz(ini:fin)=RESHAPE(GRDaz,(/naz/))
jz(ini:fin)=RESHAPE(GRDp(:,2:Nz+2),(/naz/))
kz(ini:fin)=RESHAPE(2/(Dz(:,1:Nz+1)+Dz(:,2:Nz+2)),(/naz/))

! Create gradient operator as a sparse array in coordinate format
G_size=(/nax+naz,np/)
ini=1
fin=x_nonzero_length
G_indices(ini:fin,1)=lx
G_indices(ini:fin,2)=jx
G_values(ini:fin)=kx
ini=fin+1
fin=fin+z_nonzero_length
G_indices(ini:fin,1)=lz+nax
G_indices(ini:fin,2)=jz
G_values(ini:fin)=kz

! Sort the nonzero entries in column order
CALL coordinate_to_sparse_column(G_size,G_indices,G_values,G_cpointer)

OPEN(UNIT=60,FILE='gradient.dat',STATUS='replace')
WRITE(60,*)G_size
DO nn=1,x_nonzero_length+z_nonzero_length
        WRITE(60,*)G_indices(nn,:),G_values(nn)
END DO
CLOSE(60)

DEALLOCATE(lx,jx,kx,lz,jz,kz)

END SUBROUTINE
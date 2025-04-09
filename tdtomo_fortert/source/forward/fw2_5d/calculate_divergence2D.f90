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


SUBROUTINE calculate_divergence2D()

USE modelling_2_5D
IMPLICIT NONE
INTEGER, DIMENSION(:), ALLOCATABLE              :: lx,jx
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE     :: kx
INTEGER, DIMENSION(:), ALLOCATABLE              :: lz,jz
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE     :: kz
INTEGER                                         :: x_nonzero_length,z_nonzero_length
INTEGER                                         :: ini,fin,nn
INTERFACE
        SUBROUTINE coordinate_to_sparse_row(M_size,M_indices,M_values,M_rpointer)
                INTEGER, DIMENSION(:), INTENT(IN)               :: M_size
                INTEGER, DIMENSION(:,:), INTENT(INOUT)          :: M_indices
                DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT)   :: M_values
                INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: M_rpointer
        END SUBROUTINE
END INTERFACE

IF (ALLOCATED(D_size)) DEALLOCATE(D_size)
IF (ALLOCATED(D_indices)) DEALLOCATE(D_indices)
IF (ALLOCATED(D_values)) DEALLOCATE(D_values)

x_nonzero_length=2*Nx*(Nz+2)+2*(Nz+2)
z_nonzero_length=2*Nz*(Nx+2)+2*(Nx+2)

ALLOCATE(lx(x_nonzero_length),jx(x_nonzero_length))
ALLOCATE(kx(x_nonzero_length))
ALLOCATE(lz(z_nonzero_length),jz(z_nonzero_length))
ALLOCATE(kz(z_nonzero_length))
ALLOCATE(D_size(2))
ALLOCATE(D_indices(x_nonzero_length+z_nonzero_length,2))
ALLOCATE(D_values(x_nonzero_length+z_nonzero_length))
ALLOCATE(DS_size(2))
ALLOCATE(DS_indices(x_nonzero_length+z_nonzero_length,2))
ALLOCATE(DS_values(x_nonzero_length+z_nonzero_length))

! Calculate d/dx
! Entries (l,j,k)
ini=1
fin=Nx*(Nz+2)
lx(ini:fin)=RESHAPE(GRDp(2:Nx+1,:),(/Nx*(Nz+2)/))
jx(ini:fin)=RESHAPE(GRDax(1:Nx,:),(/Nx*(Nz+2)/))
kx(ini:fin)=RESHAPE(-1/Dx(2:Nx+1,:),(/Nx*(Nz+2)/))
! Entries (l+1,j,k)
ini=fin+1
fin=ini-1+Nx*(Nz+2)
lx(ini:fin)=RESHAPE(GRDp(2:Nx+1,:),(/Nx*(Nz+2)/))
jx(ini:fin)=RESHAPE(GRDax(2:Nx+1,:),(/Nx*(Nz+2)/))
kx(ini:fin)=RESHAPE(1/Dx(2:Nx+1,:),(/Nx*(Nz+2)/))
! BC at x=0
ini=fin+1
fin=ini-1+Nz+2
lx(ini:fin)=RESHAPE(GRDp(1,:),(/Nz+2/))
jx(ini:fin)=RESHAPE(GRDax(1,:),(/Nz+2/))
kx(ini:fin)=RESHAPE(1/Dx(1,:),(/Nz+2/))
! BC at x=end
ini=fin+1
fin=ini-1+Nz+2
lx(ini:fin)=RESHAPE(GRDp(Nx+2,:),(/Nz+2/))
jx(ini:fin)=RESHAPE(GRDax(Nx+1,:),(/Nz+2/))
kx(ini:fin)=RESHAPE(-1/Dx(Nx+2,:),(/Nz+2/))

! Calculate d/dz
! Entries (l,j,k)
ini=1
fin=Nz*(Nx+2)
lz(ini:fin)=RESHAPE(GRDp(:,2:Nz+1),(/Nz*(Nx+2)/))
jz(ini:fin)=RESHAPE(GRDaz(:,1:Nz),(/Nz*(Nx+2)/))
kz(ini:fin)=RESHAPE(-1/Dz(:,2:Nz+1),(/Nz*(Nx+2)/))
! Entries (l+1,j,k)
ini=fin+1
fin=ini-1+Nz*(Nx+2)
lz(ini:fin)=RESHAPE(GRDp(:,2:Nz+1),(/Nz*(Nx+2)/))
jz(ini:fin)=RESHAPE(GRDaz(:,2:Nz+1),(/Nz*(Nx+2)/))
kz(ini:fin)=RESHAPE(1/Dz(:,2:Nz+1),(/Nz*(Nx+2)/))
! BC on z=0
ini=fin+1
fin=ini-1+Nx+2
lz(ini:fin)=RESHAPE(GRDp(:,1),(/Nx+2/))
jz(ini:fin)=RESHAPE(GRDaz(:,1),(/Nx+2/))
kz(ini:fin)=RESHAPE(1/Dz(:,1),(/Nx+2/))
! BC on z=end
ini=fin+1
fin=ini-1+Nx+2
lz(ini:fin)=RESHAPE(GRDp(:,Nz+2),(/Nx+2/))
jz(ini:fin)=RESHAPE(GRDaz(:,Nz+1),(/Nx+2/))
kz(ini:fin)=RESHAPE(-1/Dz(:,Nz+2),(/Nx+2/))

! Create divergence operator as a sparse array in coordinate format
D_size=(/np,nax+naz/)
ini=1
fin=x_nonzero_length
D_indices(ini:fin,1)=lx
D_indices(ini:fin,2)=jx
D_values(ini:fin)=kx
ini=fin+1
fin=fin+z_nonzero_length
D_indices(ini:fin,1)=lz
D_indices(ini:fin,2)=jz+nax
D_values(ini:fin)=kz

! Sort the nonzero entries in column order
CALL coordinate_to_sparse_row(D_size,D_indices,D_values,D_rpointer)

OPEN(UNIT=50,FILE='divergence.dat',STATUS='replace')
WRITE(50,*)D_size
DO nn=1,x_nonzero_length+z_nonzero_length
        WRITE(50,*)D_indices(nn,:),D_values(nn)
END DO
CLOSE(50)

DEALLOCATE(lx,jx,kx,lz,jz,kz)

END SUBROUTINE
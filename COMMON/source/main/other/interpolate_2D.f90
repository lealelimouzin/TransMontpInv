! ********************************************************************* !
! Program to create interpolated 2D models
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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This code was inspired by Nick Rawlinson's program tslicess.f90  !
! (part of the FMST package).                                      !
! It creates an interpolated 2D map from an input map file and     !
! outputs the map in a format compatible with GMT for plotting.    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

PROGRAM interpolate_2D

USE input_files, ONLY: x_number_of_pixels,y_number_of_pixels,z_number_of_pixels
USE processing, ONLY: x_pixels,y_pixels,z_pixels,x_pixel_width,y_pixel_width,z_pixel_width
USE voronoi_operations
USE chains, ONLY: x_min,x_max,y_min,y_max,z_min,z_max
USE acquisition
IMPLICIT NONE
CHARACTER(LEN=30) :: model_file_in
CHARACTER(LEN=30) :: model_file_out,bound_file,electrode_file_out
INTEGER :: vor_vtx2_vtx3,plot_plane,plot_electrodes,grid_dicing
DOUBLE PRECISION :: constant_value
INTEGER :: i,j,k,l,m,i1,j1,p,nn,ii,jj,kk,nucleus,ncell
INTEGER :: checkstat,conp,cont
INTEGER :: nnt,nnp
INTEGER :: nnx,nnz,stp,stt
INTEGER :: ddt,ddp
DOUBLE PRECISION :: sumi,sumj,u,rd1,rd2
DOUBLE PRECISION :: x_origin,y_origin,z_origin
DOUBLE PRECISION :: left,right,bottom,top
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: ui,vi
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: model_2d_rough,model_2d_fine
DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: model_3d
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: voronoi
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: nnx_pixels,nnz_pixels

OPEN(UNIT=10,FILE='interpolate_2D.in',STATUS='OLD')
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,*)vor_vtx2_vtx3
READ(10,1)model_file_in
READ(10,1)electrode_file
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,*)grid_dicing
READ(10,1)model_file_out
READ(10,1)bound_file
READ(10,*)plot_electrodes
READ(10,1)electrode_file_out
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,*)plot_plane
READ(10,*)constant_value
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,*)x_min,x_max
READ(10,*)y_min,y_max
READ(10,*)z_min,z_max
READ(10,*)x_number_of_pixels,y_number_of_pixels,z_number_of_pixels
1   FORMAT(a30)
CLOSE(10)

ddt=grid_dicing
ddp=grid_dicing

IF (vor_vtx2_vtx3.EQ.2) plot_plane=0

OPEN(UNIT=20,FILE=model_file_in,STATUS='OLD')

IF (vor_vtx2_vtx3.EQ.3) THEN

        ! Read in the 3D vtx file of the model
        ! [I read in the entire 3D model because in the future I might
        ! extend the code to plot more than one slice at a time]
        READ(20,*)y_number_of_pixels,x_number_of_pixels,z_number_of_pixels
        READ(20,*)y_origin,x_origin,z_origin
        READ(20,*)y_pixel_width,x_pixel_width,z_pixel_width
        
        ALLOCATE(model_3d(0:z_number_of_pixels+1,0:x_number_of_pixels+1,0:y_number_of_pixels+1))
        
        DO k=0,z_number_of_pixels+1
                DO i=0,x_number_of_pixels+1
                        DO j=y_number_of_pixels+1,0,-1
                                READ(20,*)model_3d(k,i,j)
                        END DO
                END DO
        END DO
        
        ! Calculate pixel arrays
        ALLOCATE(x_pixels(x_number_of_pixels))
        ALLOCATE(y_pixels(y_number_of_pixels))
        ALLOCATE(z_pixels(z_number_of_pixels))
        DO i=1,x_number_of_pixels
                x_pixels(i)=x_origin+x_pixel_width*(i-1)
        END DO
        DO j=1,y_number_of_pixels
                y_pixels(y_number_of_pixels+1-j)=y_origin-y_pixel_width*(j-1)
        END DO
        DO k=1,z_number_of_pixels
                z_pixels(k)=z_origin+z_pixel_width*(k-1)
        END DO
        
ELSE IF (vor_vtx2_vtx3.EQ.2) THEN
        
        ! Read in the 2D vtx file of the model
        READ(20,*)y_number_of_pixels,x_number_of_pixels
        READ(20,*)y_origin,x_origin
        READ(20,*)y_pixel_width,x_pixel_width
        
        ALLOCATE(model_2d_rough(0:y_number_of_pixels+1,0:x_number_of_pixels+1))
        
        DO i=0,x_number_of_pixels+1
                DO j=y_number_of_pixels+1,0,-1
                        READ(20,*)model_2d_rough(j,i)
                END DO
        END DO
        
        ! Calculate pixel arrays
        ALLOCATE(x_pixels(x_number_of_pixels))
        ALLOCATE(y_pixels(y_number_of_pixels))
        DO i=1,x_number_of_pixels
                x_pixels(i)=x_origin+x_pixel_width*(i-1)
        END DO
        DO j=1,y_number_of_pixels
                y_pixels(y_number_of_pixels+1-j)=y_origin-y_pixel_width*(j-1)
        END DO
        
ELSE IF (vor_vtx2_vtx3.EQ.1) THEN
        
        ! Read in the file of Voronoi cells and create 3D Voronoi model
        READ(20,*)ncell
        ALLOCATE(voronoi(ncell,4))
        DO nn=1,ncell
                READ(20,*)p,voronoi(nn,1),voronoi(nn,2),voronoi(nn,3),voronoi(nn,4)
        END DO
        
        ! Re-define the number of pixels on the refined (diced) grid
        x_number_of_pixels=grid_dicing*(x_number_of_pixels-1)+1
        y_number_of_pixels=grid_dicing*(y_number_of_pixels-1)+1
        z_number_of_pixels=grid_dicing*(z_number_of_pixels-1)+1
        
        ! Create pixel arrays on the refined (diced) grid
        ALLOCATE(x_pixels(x_number_of_pixels))
        ALLOCATE(y_pixels(y_number_of_pixels))
        ALLOCATE(z_pixels(z_number_of_pixels))
        CALL linspace(x_min,x_max,x_number_of_pixels,x_pixels)
        CALL linspace(y_min,y_max,y_number_of_pixels,y_pixels)
        CALL linspace(z_min,z_max,z_number_of_pixels,z_pixels)
        
        ! Get width of pixels on the refined (diced) grid
        x_pixel_width=ABS(x_max-x_min)/DBLE(x_number_of_pixels-1)
        y_pixel_width=ABS(y_max-y_min)/DBLE(y_number_of_pixels-1)
        z_pixel_width=ABS(z_max-z_min)/DBLE(z_number_of_pixels-1)
        
END IF
CLOSE(20)

! Set up the 2D array depending on which plane we want
IF (plot_plane.EQ.0) THEN ! xy plane
        
        nnt=y_number_of_pixels
        nnp=x_number_of_pixels
        
        IF (vor_vtx2_vtx3.EQ.3.OR.vor_vtx2_vtx3.EQ.2) THEN
                
                IF (vor_vtx2_vtx3.EQ.3) THEN
                        ALLOCATE(model_2d_rough(0:y_number_of_pixels+1,0:x_number_of_pixels+1))
                END IF
                
                nnx=(nnt-1)*ddt+1
                nnz=(nnp-1)*ddp+1
                
                rd1=x_pixel_width/DBLE(grid_dicing)
                rd2=y_pixel_width/DBLE(grid_dicing)
                
        ELSE IF (vor_vtx2_vtx3.EQ.1) THEN
                
                nnx=nnt
                nnz=nnp
                
                rd1=x_pixel_width
                rd2=y_pixel_width
                
        END IF
        
        IF (vor_vtx2_vtx3.NE.2) p=MINLOC(ABS(z_pixels-constant_value),1)
        
        left=x_pixels(1)
        right=x_pixels(x_number_of_pixels)
        bottom=y_pixels(1)
        top=y_pixels(y_number_of_pixels)
        
ELSE IF (plot_plane.EQ.1) THEN ! xz plane
        
        nnt=z_number_of_pixels
        nnp=x_number_of_pixels
        
        IF (vor_vtx2_vtx3.EQ.3) THEN
                
                ALLOCATE(model_2d_rough(0:z_number_of_pixels+1,0:x_number_of_pixels+1))
                
                nnx=(nnt-1)*ddt+1
                nnz=(nnp-1)*ddp+1
                
                rd1=x_pixel_width/DBLE(grid_dicing)
                rd2=z_pixel_width/DBLE(grid_dicing)
                
        ELSE IF (vor_vtx2_vtx3.EQ.1) THEN
                
                nnx=nnt
                nnz=nnp
                
                rd1=x_pixel_width
                rd2=z_pixel_width
                
        END IF
        
        p=MINLOC(ABS(y_pixels-constant_value),1)
        
        left=x_pixels(1)
        right=x_pixels(x_number_of_pixels)
        bottom=z_pixels(1)
        top=z_pixels(z_number_of_pixels)
        
ELSE IF (plot_plane.EQ.2) THEN ! yz plane
        
        nnt=z_number_of_pixels
        nnp=y_number_of_pixels
        
        IF (vor_vtx2_vtx3.EQ.3) THEN
                
                ALLOCATE(model_2d_rough(0:z_number_of_pixels+1,0:y_number_of_pixels+1))
                
                nnx=(nnt-1)*ddt+1
                nnz=(nnp-1)*ddp+1
                
                rd1=y_pixel_width/DBLE(grid_dicing)
                rd2=z_pixel_width/DBLE(grid_dicing)
                
        ELSE IF (vor_vtx2_vtx3.EQ.1) THEN
                
                nnx=nnt
                nnz=nnp
                
                rd1=y_pixel_width
                rd2=z_pixel_width
                
        END IF
        
        p=MINLOC(ABS(x_pixels-constant_value),1)
        
        left=y_pixels(1)
        right=y_pixels(y_number_of_pixels)
        bottom=z_pixels(1)
        top=z_pixels(z_number_of_pixels)
        
END IF

! Calculate refined grid
ALLOCATE(model_2d_fine(nnx,nnz))

IF (vor_vtx2_vtx3.EQ.3.OR.vor_vtx2_vtx3.EQ.2) THEN
        
        IF (vor_vtx2_vtx3.EQ.3) THEN
                
                IF (plot_plane.EQ.0) THEN ! xy plane
                        
                        model_2d_rough=model_3d(p,:,:)
                        
                ELSE IF (plot_plane.EQ.1) THEN ! xz plane
                        
                        model_2d_rough=model_3d(:,:,p)
                        
                ELSE IF (plot_plane.EQ.2) THEN ! yz plane
                
                        model_2d_rough=model_3d(:,p,:)
                        
                END IF
                
        END IF
        
        IF (grid_dicing.NE.1) THEN
                !!!!! --- This part has been taken and adapted from tslicess.f90 --- !!!!!
                ! Compute the values of the basis functions
                ALLOCATE(ui(ddt+1,4))
                DO i=1,ddt+1
                        u=ddt
                        u=(i-1)/u
                        ui(i,1)=(1.0-u)**3/6.0
                        ui(i,2)=(4.0-6.0*u**2+3.0*u**3)/6.0
                        ui(i,3)=(1.0+3.0*u+3.0*u**2-3.0*u**3)/6.0
                        ui(i,4)=u**3/6.0
                ENDDO
                ALLOCATE(vi(ddp+1,4))
                DO i=1,ddp+1
                        u=ddp
                        u=(i-1)/u
                        vi(i,1)=(1.0-u)**3/6.0
                        vi(i,2)=(4.0-6.0*u**2+3.0*u**3)/6.0
                        vi(i,3)=(1.0+3.0*u+3.0*u**2-3.0*u**3)/6.0
                        vi(i,4)=u**3/6.0
                ENDDO
                DO i=1,nnp-1
                        conp=ddp
                        IF(i==nnp-1)conp=ddp+1
                        DO j=1,nnt-1
                                cont=ddt
                                IF(j==nnt-1)cont=ddt+1
                                DO l=1,conp
                                        stp=ddp*(i-1)+l
                                        DO m=1,cont
                                                stt=ddt*(j-1)+m
                                                sumi=0.0
                                                DO i1=1,4
                                                        sumj=0.0
                                                        DO j1=1,4
                                                                sumj=sumj+ui(m,j1)*model_2d_rough(j-2+j1,i-2+i1)
                                                        ENDDO
                                                        sumi=sumi+vi(l,i1)*sumj
                                                ENDDO
                                                model_2d_fine(stt,stp)=sumi
                                        ENDDO
                                ENDDO
                        ENDDO
                ENDDO
                !!!!! --- End of part taken and adapted from tslicess.f90 --- !!!!!
        ELSE IF (grid_dicing.EQ.1) THEN
                
                model_2d_fine=model_2d_rough
                
        END IF

ELSE IF (vor_vtx2_vtx3.EQ.1) THEN

        IF (plot_plane.EQ.0) THEN ! xy plane
                
                ALLOCATE(nnz_pixels(x_number_of_pixels))
                ALLOCATE(nnx_pixels(y_number_of_pixels))
                nnz_pixels=x_pixels
                nnx_pixels=y_pixels
                
        ELSE IF (plot_plane.EQ.1) THEN ! xz plane
                
                ALLOCATE(nnz_pixels(x_number_of_pixels))
                ALLOCATE(nnx_pixels(z_number_of_pixels))
                nnz_pixels=x_pixels
                nnx_pixels=z_pixels
                
        ELSE IF (plot_plane.EQ.2) THEN ! yz plane
                
                ALLOCATE(nnz_pixels(y_number_of_pixels))
                ALLOCATE(nnx_pixels(z_number_of_pixels))
                nnz_pixels=y_pixels
                nnx_pixels=z_pixels
                
        END IF
        
        DO i=1,nnz
                DO j=1,nnx
                        
                        IF (plot_plane.EQ.0) THEN ! xy plane
                                
                                CALL find_closest_nucleus([nnz_pixels(i),nnx_pixels(j),z_pixels(p)],&
                                                           voronoi,ncell,nucleus)
                                
                        ELSE IF (plot_plane.EQ.1) THEN ! xz plane
                                
                                CALL find_closest_nucleus([nnz_pixels(i),y_pixels(p),nnx_pixels(j)],&
                                                           voronoi,ncell,nucleus)
                                
                        ELSE IF (plot_plane.EQ.2) THEN ! yz plane
                                
                                CALL find_closest_nucleus([x_pixels(p),nnz_pixels(i),nnx_pixels(j)],&
                                                           voronoi,ncell,nucleus)
                                
                        END IF
                         
                        model_2d_fine(j,i)=voronoi(nucleus,4)
                END DO
        END DO
        
END IF

! Write to file
OPEN(UNIT=50,FILE=bound_file,STATUS='REPLACE')
WRITE(50,'(f16.10)')left
WRITE(50,'(f16.10)')right
WRITE(50,'(f16.10)')bottom
WRITE(50,'(f16.10)')top
WRITE(50,'(f16.10)')rd1
WRITE(50,'(f16.10)')rd2
CLOSE(50)

OPEN(UNIT=50,FILE=model_file_out,STATUS='REPLACE')
DO i=1,nnz
        DO j=1,nnx
                WRITE(50,*)model_2d_fine(j,i)
        ENDDO
ENDDO
CLOSE(50)

IF (ALLOCATED(model_3d)) DEALLOCATE(model_3d)
IF (ALLOCATED(model_2d_rough)) DEALLOCATE(model_2d_rough)
IF (ALLOCATED(ui)) DEALLOCATE(ui)
IF (ALLOCATED(vi)) DEALLOCATE(vi)

IF (ALLOCATED(voronoi)) DEALLOCATE(voronoi)
IF (ALLOCATED(nnz_pixels)) DEALLOCATE(nnz_pixels)
IF (ALLOCATED(nnx_pixels)) DEALLOCATE(nnx_pixels)

DEALLOCATE(x_pixels,y_pixels)
IF (ALLOCATED(z_pixels)) DEALLOCATE(z_pixels)
DEALLOCATE(model_2d_fine)

IF (plot_electrodes.EQ.1) THEN

        CALL read_electrodes()
        
        OPEN(UNIT=60,FILE=electrode_file_out,STATUS='REPLACE')
        
        IF (SIZE(electrodes,2).EQ.2) THEN
        
                IF (plot_plane.EQ.0.AND.vor_vtx2_vtx3.NE.2) THEN ! xy plane
                
                        WRITE(*,*)'Plotting xy plane, but y-coordinates'
                        WRITE(*,*)'of electrodes missing in electrodes'
                        WRITE(*,*)'file ',electrode_file,'!'
                        WRITE(*,*)'TERMINATING PROGRAM!!'
                        STOP
                
                ELSE IF (plot_plane.EQ.1.OR.plot_plane.EQ.2.OR.vor_vtx2_vtx3.EQ.2) THEN
                        
                        DO nn=1,number_of_electrodes
                                
                                WRITE(60,*)electrodes(nn,1),electrodes(nn,2)
                                
                        END DO
                        
                END IF
                
        ELSE IF (SIZE(electrodes,2).EQ.2) THEN
                
                DO nn=1,number_of_electrodes
                        
                        IF (plot_plane.EQ.0) THEN ! xy plane
                                
                                IF (electrodes(nn,3).EQ.constant_value) WRITE(60,*)electrodes(nn,1),electrodes(nn,2)
                                
                        ELSE IF (plot_plane.EQ.1) THEN ! xz plane
                                
                                IF (electrodes(nn,2).EQ.constant_value) WRITE(60,*)electrodes(nn,1),electrodes(nn,3)
                                
                        ELSE IF (plot_plane.EQ.2) THEN ! yz plane
                                
                                IF (electrodes(nn,1).EQ.constant_value) WRITE(60,*)electrodes(nn,2),electrodes(nn,3)
                                
                        END IF
                        
                END DO
                
        END IF
        
        CLOSE(60)
        
END IF
        

END PROGRAM interpolate_2D
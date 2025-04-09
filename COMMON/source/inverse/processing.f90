! ********************************************************************* !
! Module for processing Markov chains in procsamples.f90
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


MODULE processing

USE chains
USE input_files
USE acquisition, ONLY: number_of_measurements,observed_data,data_noise
USE debug
IMPLICIT NONE
INTEGER                                                 :: iteration
DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE, SAVE   :: average,stdev,median,maximum,harmean,root_mean_square,&
                                                           entropy,skewness,kurtosis,node_density,pixel_size
INTEGER, DIMENSION(:,:,:), ALLOCATABLE, SAVE            :: node_count
DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE, SAVE   :: moment_1,moment_2,moment_3,moment_4
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, SAVE       :: ncells_iter,misfit_iter
INTEGER, DIMENSION(:), ALLOCATABLE, SAVE                :: posterior_ncells
INTEGER, DIMENSION(:,:), ALLOCATABLE, SAVE              :: posterior_noise
INTEGER, DIMENSION(:,:,:,:), ALLOCATABLE, SAVE          :: posterior_value
INTEGER, DIMENSION(:), ALLOCATABLE, SAVE                :: ncell_all_first,ncell_all_last
DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE, SAVE   :: voronoi_all_first,voronoi_all_last
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, SAVE     :: residuals_pixels,residuals_histo
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, SAVE       :: x_pixels,y_pixels,z_pixels
DOUBLE PRECISION, SAVE                                  :: x_pixel_width,y_pixel_width,z_pixel_width
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, SAVE     :: noise_pixels
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, SAVE       :: noise_pixel_width
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, SAVE       :: value_pixels
DOUBLE PRECISION, SAVE                                  :: value_pixel_width
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, SAVE     :: temperatures
INTEGER, SAVE                                           :: accept
DOUBLE PRECISION, SAVE                                  :: x_new,y_new,z_new,value_new
DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE, SAVE   :: xyz_model


CONTAINS


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE initialize_processing()
! Allocates all necessary arrays for processing
IMPLICIT NONE
INTEGER                 :: nn
INTEGER                 :: ii,jj,kk
DOUBLE PRECISION        :: xp_min,xp_max,yp_min,yp_max,zp_min,zp_max
DOUBLE PRECISION        :: xp_size,yp_size,zp_size
INTERFACE
        SUBROUTINE linspace(d1,d2,n,grid)
                INTEGER, INTENT(IN)                             :: n
                DOUBLE PRECISION, INTENT(IN)                    :: d1,d2
                DOUBLE PRECISION, DIMENSION(n), INTENT(OUT)     :: grid
        END SUBROUTINE
END INTERFACE

IF (debug_mode) WRITE(debug_unit,*)'--- Initializing arrays for processing ---'

! Models and maps
ALLOCATE(x_pixels(x_number_of_pixels))
IF (vtx_or_xyz.EQ.0) THEN
        CALL linspace(x_min,x_max,x_number_of_pixels,x_pixels)
        x_pixel_width=(x_max-x_min)/(x_number_of_pixels-1)
ELSE IF (vtx_or_xyz.EQ.1) THEN
        DO nn=1,x_number_of_pixels
                x_pixels(nn)=x_min+(DBLE(nn)-DBLE(0.5))*(x_max-x_min)/x_number_of_pixels
        END DO
        x_pixel_width=ABS(x_pixels(2)-x_pixels(1))
END IF
ALLOCATE(y_pixels(y_number_of_pixels))
IF (vtx_or_xyz.EQ.0) THEN
        CALL linspace(y_min,y_max,y_number_of_pixels,y_pixels)
        y_pixel_width=(y_max-y_min)/(y_number_of_pixels-1)
ELSE IF (vtx_or_xyz.EQ.1) THEN
        DO nn=1,y_number_of_pixels
                y_pixels(nn)=y_min+(DBLE(nn)-DBLE(0.5))*(y_max-y_min)/y_number_of_pixels
        END DO
        y_pixel_width=ABS(y_pixels(2)-y_pixels(1))
END IF
ALLOCATE(z_pixels(z_number_of_pixels))
IF (vtx_or_xyz.EQ.0) THEN
        CALL linspace(z_min,z_max,z_number_of_pixels,z_pixels)
        z_pixel_width=(z_max-z_min)/(z_number_of_pixels-1)
ELSE IF (vtx_or_xyz.EQ.1) THEN
        DO nn=1,z_number_of_pixels
                z_pixels(nn)=z_min+(DBLE(nn)-DBLE(0.5))*(z_max-z_min)/z_number_of_pixels
        END DO
        z_pixel_width=ABS(z_pixels(2)-z_pixels(1))
END IF
ALLOCATE(average(z_number_of_pixels,x_number_of_pixels,y_number_of_pixels))
average=0
ALLOCATE(stdev(z_number_of_pixels,x_number_of_pixels,y_number_of_pixels))
stdev=0
ALLOCATE(median(z_number_of_pixels,x_number_of_pixels,y_number_of_pixels))
median=0
ALLOCATE(maximum(z_number_of_pixels,x_number_of_pixels,y_number_of_pixels))
maximum=0
ALLOCATE(harmean(z_number_of_pixels,x_number_of_pixels,y_number_of_pixels))
harmean=0
ALLOCATE(root_mean_square(z_number_of_pixels,x_number_of_pixels,y_number_of_pixels))
root_mean_square=0
ALLOCATE(entropy(z_number_of_pixels,x_number_of_pixels,y_number_of_pixels))
entropy=0
ALLOCATE(skewness(z_number_of_pixels,x_number_of_pixels,y_number_of_pixels))
skewness=0
ALLOCATE(kurtosis(z_number_of_pixels,x_number_of_pixels,y_number_of_pixels))
kurtosis=0
ALLOCATE(node_count(z_number_of_pixels,x_number_of_pixels,y_number_of_pixels))
node_count=0
ALLOCATE(node_density(z_number_of_pixels,x_number_of_pixels,y_number_of_pixels))
node_density=0
ALLOCATE(moment_1(z_number_of_pixels,x_number_of_pixels,y_number_of_pixels))
moment_1=0
ALLOCATE(moment_2(z_number_of_pixels,x_number_of_pixels,y_number_of_pixels))
moment_2=0
ALLOCATE(moment_3(z_number_of_pixels,x_number_of_pixels,y_number_of_pixels))
moment_3=0
ALLOCATE(moment_4(z_number_of_pixels,x_number_of_pixels,y_number_of_pixels))
moment_4=0
ALLOCATE(pixel_size(z_number_of_pixels,x_number_of_pixels,y_number_of_pixels))
pixel_size=0

! Calculate size of pixels for node search
DO jj=1,y_number_of_pixels
        yp_min=y_pixels(jj)-y_pixel_size/2
        IF (yp_min.LT.y_min) yp_min=y_min
        yp_max=y_pixels(jj)+y_pixel_size/2
        IF (yp_max.GT.y_max) yp_max=y_max
        yp_size=ABS(yp_max-yp_min)
        IF (yp_size.EQ.0) yp_size=1
        DO ii=1,x_number_of_pixels
                xp_min=x_pixels(ii)-x_pixel_size/2
                IF (xp_min.LT.x_min) xp_min=x_min
                xp_max=x_pixels(ii)+x_pixel_size/2
                IF (xp_max.GT.x_max) xp_max=x_max
                xp_size=ABS(xp_max-xp_min)
                IF (xp_size.EQ.0) xp_size=1
                DO kk=1,z_number_of_pixels
                        zp_min=z_pixels(kk)-z_pixel_size/2
                        IF (zp_min.LT.z_min) zp_min=z_min
                        zp_max=z_pixels(kk)+z_pixel_size/2
                        IF (zp_max.GT.z_max) zp_max=z_max
                        zp_size=ABS(zp_max-zp_min)
                        IF (zp_size.EQ.0) zp_size=1
                        pixel_size(kk,ii,jj)=xp_size*yp_size*zp_size
                END DO
        END DO
END DO

! Chain convergence
ALLOCATE(ncells_iter(samples_so_far+1))
ALLOCATE(misfit_iter(samples_so_far+1))
ncells_iter=0
misfit_iter=0

! Posteriors
ALLOCATE(posterior_ncells(ncell_max))
posterior_ncells=0
ALLOCATE(value_pixels(value_number_of_pixels))
DO nn=1,value_number_of_pixels
        value_pixels(nn)=value_min+(DBLE(nn)-DBLE(0.5))*(value_max-value_min)/value_number_of_pixels
END DO
value_pixel_width=ABS(value_pixels(2)-value_pixels(1))
ALLOCATE(posterior_value(z_number_of_pixels,x_number_of_pixels,y_number_of_pixels,value_number_of_pixels))
posterior_value=0
IF (noise_type.EQ.0.OR.noise_type.EQ.3.OR.noise_type.EQ.6) THEN
        ALLOCATE(noise_pixels(noise_number_of_pixels,1))
        ALLOCATE(noise_pixel_width(1))
        DO nn=1,noise_number_of_pixels
                IF (noise_type.EQ.0) noise_pixels(nn,1)=sigma_min+(DBLE(nn)-DBLE(0.5))&
                                                        *(sigma_max-sigma_min)/noise_number_of_pixels
                IF (noise_type.EQ.3.OR.noise_type.EQ.6) noise_pixels(nn,1)=lambda_min+(DBLE(nn)-DBLE(0.5))&
                                                        *(lambda_max-lambda_min)/noise_number_of_pixels
        END DO
        noise_pixel_width(1)=ABS(noise_pixels(2,1)-noise_pixels(1,1))
        ALLOCATE(posterior_noise(noise_number_of_pixels,number_of_datasets))
ELSE IF (noise_type.EQ.1.OR.noise_type.EQ.2.OR.noise_type.EQ.5) THEN
        ALLOCATE(noise_pixels(noise_number_of_pixels,2))
        ALLOCATE(noise_pixel_width(2))
        DO nn=1,noise_number_of_pixels
                noise_pixels(nn,1)=a_min+(DBLE(nn)-DBLE(0.5))*(a_max-a_min)/noise_number_of_pixels
                noise_pixels(nn,2)=b_min+(DBLE(nn)-DBLE(0.5))*(b_max-b_min)/noise_number_of_pixels
        END DO
        noise_pixel_width(1)=ABS(noise_pixels(2,1)-noise_pixels(1,1))
        noise_pixel_width(2)=ABS(noise_pixels(2,2)-noise_pixels(1,2))
        ALLOCATE(posterior_noise(noise_number_of_pixels,2*number_of_datasets))
ELSE IF (noise_type.EQ.4) THEN
        ALLOCATE(noise_pixels(noise_number_of_pixels,1))
        ALLOCATE(noise_pixel_width(1))
        DO nn=1,noise_number_of_pixels
                noise_pixels(nn,1)=MINVAL(data_noise)+(DBLE(nn)-DBLE(0.5))&
                                   *(MAXVAL(data_noise)-MINVAL(data_noise))/noise_number_of_pixels
        END DO
        noise_pixel_width(1)=ABS(noise_pixels(2,1)-noise_pixels(1,1))
        ALLOCATE(posterior_noise(noise_number_of_pixels,1))
END IF
posterior_noise=0

! Voronoi models
ALLOCATE(voronoi_current(ncell_max,4))
ALLOCATE(ncell_all_first(number_of_chains))
ALLOCATE(ncell_all_last(number_of_chains))
ALLOCATE(voronoi_all_first(ncell_max,4,number_of_chains))
ALLOCATE(voronoi_all_last(ncell_max,4,number_of_chains))

! Noise parameters
IF (noise_type.EQ.1.OR.noise_type.EQ.2.OR.noise_type.EQ.5) THEN
        ALLOCATE(a_current(number_of_datasets))
        ALLOCATE(b_current(number_of_datasets))
ELSE IF (noise_type.EQ.3.OR.noise_type.EQ.6) THEN
        ALLOCATE(lambda_current(number_of_datasets))
ELSE IF (noise_type.EQ.0) THEN
        ALLOCATE(sigma_current(number_of_datasets))
END IF

! Residuals
ALLOCATE(residuals(number_of_measurements,2))
ALLOCATE(residuals_pixels(residuals_number_of_pixels,2))
CALL linspace(-max_abs_residual,max_abs_residual,residuals_number_of_pixels,residuals_pixels(:,1))
CALL linspace(-max_percentage_residual,max_percentage_residual,residuals_number_of_pixels,residuals_pixels(:,2))
ALLOCATE(residuals_histo(residuals_number_of_pixels,2))
residuals_histo=0

IF (debug_mode) WRITE(debug_unit,*)'--- Initialized arrays for processing ---'

END SUBROUTINE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE check_processing_parameters()
IMPLICIT NONE
CHARACTER(LEN=3)        :: section
input_error=.FALSE.
! Compare burn-in and number of samples
IF (burn_in.GE.samples_so_far) THEN
        input_error=.TRUE.
        section='(D)'
END IF
! Compare thinning and number of samples
IF (thinning.GE.samples_so_far) THEN
        input_error=.TRUE.
        section='(D)'
END IF
IF (input_error) THEN
        WRITE(*,*)
        WRITE(*,*)'********************************************'
        WRITE(*,*)'ERROR in procsamples.in'
        WRITE(*,*)'--------------------------------------------'
        WRITE(*,*)'Check section ',section,' of procsamples.in.'
        WRITE(*,*)'TERMINATING PROGRAM!!'
        WRITE(*,*)'********************************************'
        WRITE(*,*)
END IF
END SUBROUTINE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE open_files()
IMPLICIT NONE
! File of chain steps
OPEN(UNIT=samples_unit,FILE=TRIM(samples_file) // TRIM(chain_name),&
     STATUS='OLD',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
READ(samples_unit) ! ignore line of file type
READ(samples_unit) ! ignore line with number of cells
! File of residuals
OPEN(UNIT=residuals_unit,FILE=TRIM(residuals_file) // TRIM(chain_name),&
     STATUS='OLD',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
READ(residuals_unit) ! ignore line of file type
! File of number of cells
OPEN(UNIT=ncell_unit,FILE=TRIM(ncell_file) // TRIM(chain_name),&
     STATUS='OLD',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
READ(ncell_unit) ! ignore line of file type
! File of noise
IF (noise_type.NE.4)&
OPEN(UNIT=noise_unit,FILE=TRIM(noise_file) // TRIM(chain_name),&
     STATUS='OLD',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
IF (noise_type.EQ.4.AND.chain_number.EQ.1)&
OPEN(UNIT=noise_unit,FILE=TRIM(noise_file),&
     STATUS='OLD',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
READ(noise_unit) ! ignore line of file type
READ(noise_unit) ! ignore line of noise type
READ(noise_unit) ! ignore line with number of datasets/readings
! File of misfits
OPEN(UNIT=misfits_unit,FILE=TRIM(misfits_file) // TRIM(chain_name),&
     STATUS='OLD',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
READ(misfits_unit) ! ignore line of file type
! File of acceptance ratios
OPEN(UNIT=acceptance_unit,FILE=TRIM(acceptance_file) // TRIM(chain_name),&
     STATUS='OLD',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
READ(acceptance_unit) ! ignore line of file type
END SUBROUTINE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE close_files()
IMPLICIT NONE
LOGICAL         :: file_open
CLOSE(samples_unit)
CLOSE(residuals_unit)
CLOSE(ncell_unit)
INQUIRE(noise_unit,OPENED=file_open)
IF (file_open) CLOSE(noise_unit)
CLOSE(misfits_unit)
CLOSE(acceptance_unit)
END SUBROUTINE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE arrange_temperatures()

USE tempering, ONLY: T1_chain,T2_chain
IMPLICIT NONE
INTEGER                 :: nn,pp
INTEGER                 :: swap
DOUBLE PRECISION        :: dummy

OPEN(UNIT=temperatures_unit,FILE=TRIM(temperatures_file),&
     STATUS='OLD',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
READ(temperatures_unit)
READ(temperatures_unit)number_of_temperatures
ALLOCATE(temperatures(samples_so_far,number_of_temperatures))
temperatures=1
ALLOCATE(temperature_values0(number_of_temperatures))
DO nn=1,number_of_temperatures
        READ(temperatures_unit)temperature_values0(nn)
END DO
IF (tempering_start.EQ.0) tempering_start=1
IF (samples_so_far.GE.tempering_start) THEN
        temperatures(tempering_start,:)=temperature_values0
        DO nn=tempering_start,samples_so_far-1
                READ(temperatures_unit)pp,swap,T1_chain,T2_chain,dummy
                temperatures(nn+1,:)=temperatures(nn,:)
                IF (swap.EQ.1) temperatures(nn+1,(/T1_chain,T2_chain/))=temperatures(nn+1,(/T2_chain,T1_chain/))
        END DO
END IF
READ(temperatures_unit)pp,swap,T1_chain,T2_chain
CLOSE(temperatures_unit)

END SUBROUTINE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE update_model()
IMPLICIT NONE
IF (step.EQ.1.AND.accept.EQ.1) THEN
        voronoi_current(ncell_current,1)=x_new
        voronoi_current(ncell_current,2)=y_new
        voronoi_current(ncell_current,3)=z_new
        voronoi_current(ncell_current,4)=value_new
        voronoi_current(ncell_current+1:ncell_max,:)=0
ELSE IF (step.EQ.2.AND.accept.EQ.1) THEN
        voronoi_current(cell_number,:)=voronoi_current(ncell_current+1,:)
        voronoi_current(ncell_current+1:ncell_max,:)=0
ELSE IF (step.EQ.3.AND.accept.EQ.1) THEN
        voronoi_current(cell_number,1)=x_new
        voronoi_current(cell_number,2)=y_new
        voronoi_current(cell_number,3)=z_new
ELSE IF (step.EQ.4.AND.accept.EQ.1) THEN
        voronoi_current(cell_number,4)=value_new
ELSE IF (step.EQ.5.AND.accept.EQ.1) THEN
        voronoi_current=voronoi_current
ELSE
        voronoi_current=voronoi_current
END IF
END SUBROUTINE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE voronoi_to_xyz()
! Turns the Voronoi model into an xyz model and adds to the posterior on values
USE voronoi_operations
IMPLICIT NONE
INTEGER                         :: ii,jj,kk,ll
DOUBLE PRECISION, DIMENSION(3)  :: point
INTEGER                         :: nucleus
IF (.NOT.ALLOCATED(xyz_model)) ALLOCATE(xyz_model(z_number_of_pixels,&
                                                  x_number_of_pixels,&
                                                  y_number_of_pixels))
DO jj=1,y_number_of_pixels
        DO ii=1,x_number_of_pixels
                DO kk=1,z_number_of_pixels
                        ! Find which Voronoi cell the pixel falls into
                        point=(/x_pixels(ii),y_pixels(jj),z_pixels(kk)/)
                        CALL find_closest_nucleus(point,voronoi_current(1:ncell_current,1:3),ncell_current,nucleus)
                        ! Assign value of Voronoi cell to pixel
                        xyz_model(kk,ii,jj)=voronoi_current(nucleus,4)
                        ! Add value to posterior
                        ll=CEILING((voronoi_current(nucleus,4)-value_min)*value_number_of_pixels/(value_max-value_min))
                        IF (ll.EQ.0) ll=1
                        posterior_value(kk,ii,jj,ll)=posterior_value(kk,ii,jj,ll)+1
                END DO
        END DO
END DO
END SUBROUTINE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE evaluate_posteriors(pncell,pnoise)
! Adds the sample to the posteriors on number of cells and noise
IMPLICIT NONE
LOGICAL, INTENT(IN)     :: pncell,pnoise
INTEGER                 :: nn,ll
! Number of cells
IF (pncell) THEN
        posterior_ncells(ncell_current)=posterior_ncells(ncell_current)+1
END IF
! Noise
IF (pnoise) THEN
        IF (noise_type.EQ.0) THEN
                DO nn=1,number_of_datasets
                        ll=CEILING((sigma_current(nn)-sigma_min)*noise_number_of_pixels/(sigma_max-sigma_min))
                        IF (ll.EQ.0) ll=1
                        posterior_noise(ll,nn)=posterior_noise(ll,nn)+1
                END DO
        ELSE IF (noise_type.EQ.1.OR.noise_type.EQ.2.OR.noise_type.EQ.5) THEN
                DO nn=1,number_of_datasets
                        ll=CEILING((a_current(nn)-a_min)*noise_number_of_pixels/(a_max-a_min))
                        IF (ll.EQ.0) ll=1
                        posterior_noise(ll,nn) = posterior_noise(ll,nn)+1
                        ll=CEILING((b_current(nn)-b_min)*noise_number_of_pixels/(b_max-b_min))
                        IF (ll.EQ.0) ll=1
                        posterior_noise(ll,nn+number_of_datasets)=posterior_noise(ll,nn+number_of_datasets)+1
                END DO
        ELSE IF (noise_type.EQ.3.OR.noise_type.EQ.6) THEN
                DO nn=1,number_of_datasets
                        ll=CEILING((lambda_current(nn)-lambda_min)*noise_number_of_pixels/(lambda_max-lambda_min))
                        IF (ll.EQ.0) ll=1       
                        posterior_noise(ll,nn) = posterior_noise(ll,nn)+1
                END DO
        ELSE IF (noise_type.EQ.4.AND.chain_number.EQ.1.AND.iteration.EQ.0) THEN
                DO nn=1,number_of_measurements
                        ll=CEILING((data_noise(nn)-MINVAL(data_noise))&
                                   *noise_number_of_pixels/(MAXVAL(data_noise)-MINVAL(data_noise)))
                        IF (ll.EQ.0) ll=1
                        IF (ll.GT.noise_number_of_pixels) ll=noise_number_of_pixels
                        posterior_noise(ll,1) = posterior_noise(ll,1)+1
                END DO
        END IF
END IF
END SUBROUTINE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE evaluate_statistics()
! Adds the Voronoi model to the statistical moments
IMPLICIT NONE
moment_1=moment_1+xyz_model
moment_2=moment_2+xyz_model**2
moment_3=moment_3+xyz_model**3
moment_4=moment_4+xyz_model**4
WHERE (xyz_model.EQ.0)
        xyz_model=EPSILON(DBLE(0))
END WHERE
harmean=harmean+1/xyz_model
END SUBROUTINE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE calculate_node_density()
! Adds the model nodes to the node density count
IMPLICIT NONE
INTEGER                 :: nn
INTEGER                 :: ii,jj,kk
DOUBLE PRECISION        :: xp_min,xp_max,yp_min,yp_max,zp_min,zp_max
DOUBLE PRECISION        :: x_node,y_node,z_node
DO jj=1,y_number_of_pixels
        yp_min=y_pixels(jj)-y_pixel_size/2
        IF (yp_min.LT.y_min) yp_min=y_min
        yp_max=y_pixels(jj)+y_pixel_size/2
        IF (yp_max.GT.y_max) yp_max=y_max
        DO ii=1,x_number_of_pixels
                xp_min=x_pixels(ii)-x_pixel_size/2
                IF (xp_min.LT.x_min) xp_min=x_min
                xp_max=x_pixels(ii)+x_pixel_size/2
                IF (xp_max.GT.x_max) xp_max=x_max
                DO kk=1,z_number_of_pixels
                        zp_min=z_pixels(kk)-z_pixel_size/2
                        IF (zp_min.LT.z_min) zp_min=z_min
                        zp_max=z_pixels(kk)+z_pixel_size/2
                        IF (zp_max.GT.z_max) zp_max=z_max
                        DO nn=1,ncell_current
                                x_node=voronoi_current(nn,1)
                                y_node=voronoi_current(nn,2)
                                z_node=voronoi_current(nn,3)
                                IF (x_node.GE.xp_min.AND.x_node.LE.xp_max.AND.&
                                    y_node.GE.yp_min.AND.y_node.LE.yp_max.AND.&
                                    z_node.GE.zp_min.AND.z_node.LE.zp_max) THEN
                                        node_count(kk,ii,jj)=node_count(kk,ii,jj)+1
                                END IF
                        END DO
                END DO
        END DO
END DO
END SUBROUTINE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE calculate_entropy(samples_in_ensemble)
! Calculates entropy using eq. (2) by Wellmann & Regenauer-Lieb (2012)
IMPLICIT NONE
INTEGER, INTENT(IN)     :: samples_in_ensemble
INTEGER                 :: nn
INTEGER                 :: ii,jj,kk,ll
DOUBLE PRECISION        :: pm
DO jj=1,y_number_of_pixels
        DO ii=1,x_number_of_pixels
                DO kk=1,z_number_of_pixels
                        DO ll=1,value_number_of_pixels
                                pm=DBLE(posterior_value(kk,ii,jj,ll))/(DBLE(samples_in_ensemble))!*value_pixel_width)
                                IF (pm.EQ.0.OR.pm.EQ.1) THEN
                                        entropy(kk,ii,jj)=entropy(kk,ii,jj)-DBLE(0)
                                ELSE
                                        entropy(kk,ii,jj)=entropy(kk,ii,jj)-pm*LOG(pm)/LOG(DBLE(2))
                                END IF
                        END DO
                END DO
        END DO
END DO
END SUBROUTINE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE calculate_median_and_maximum(samples_in_ensemble)
! Calculates the median and maximum (mode) model from the posterior on values
IMPLICIT NONE
INTEGER, INTENT(IN)     :: samples_in_ensemble
INTEGER                 :: ii,jj,kk,ll,nn
DO jj=1,y_number_of_pixels
        DO ii=1,x_number_of_pixels
                DO kk=1,z_number_of_pixels
                        ll=MAXLOC(posterior_value(kk,ii,jj,:),1) ! find maximum
                        maximum(kk,ii,jj)=value_pixels(ll)
                        nn=0
                        DO ll=1,value_number_of_pixels ! find median
                                nn=nn+posterior_value(kk,ii,jj,ll)
                                IF (nn.GT.samples_in_ensemble/2) THEN
                                        median(kk,ii,jj)=value_pixels(ll)
                                        EXIT
                                END IF
                        END DO
                END DO
        END DO
END DO
END SUBROUTINE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE calculate_residuals()
! Adds the residuals to the residuals histogram array
IMPLICIT NONE
INTEGER         :: nn,ll
DO nn=1,number_of_measurements
        ll=NINT((residuals(nn,1)-(-max_abs_residual))*residuals_number_of_pixels/(2*max_abs_residual))
        IF (ll.LT.1) ll=1
        IF (ll.GT.residuals_number_of_pixels) ll=residuals_number_of_pixels
        residuals_histo(ll,1)=residuals_histo(ll,1)+1
        ll=NINT((residuals(nn,2)-(-max_percentage_residual))*residuals_number_of_pixels/(2*max_percentage_residual))
        IF (ll.LT.1) ll=1
        IF (ll.GT.residuals_number_of_pixels) ll=residuals_number_of_pixels
        residuals_histo(ll,2)=residuals_histo(ll,2)+1
END DO
END SUBROUTINE

END MODULE
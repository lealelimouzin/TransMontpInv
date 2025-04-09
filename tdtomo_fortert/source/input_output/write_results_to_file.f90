! ********************************************************************* !
! Subroutine to write processing results to file in procsamples.f90
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


SUBROUTINE write_results_to_file(samples_in_ensemble)

USE input_files, ONLY: average_file,stdev_file,median_file,maximum_file,harmean_file,&
                       rms_file,entropy_file,skewness_file,kurtosis_file,ncellsi_file,misfiti_file,&
                       postncells_file,postnoise_file,postvalue_file,node_file,first_file,last_file,&
                       residuals_histo_file
USE processing
USE chains
USE string_conversion
IMPLICIT NONE
INTEGER, INTENT(IN)             :: samples_in_ensemble
INTEGER                         :: nn,ll
INTEGER                         :: ii,jj,kk
CHARACTER(LEN=*), PARAMETER     :: outfmt0="(1F22.10)"
CHARACTER(LEN=*), PARAMETER     :: outfmt1="(4F22.10)"
CHARACTER(LEN=*), PARAMETER     :: outfmt2="(1I10,4F22.10)"
CHARACTER(LEN=*), PARAMETER     :: outfmt3="(1I12,1F32.10)"
CHARACTER(LEN=*), PARAMETER     :: outfmt4="(1F22.10,1F32.2,1F22.10,1F32.2)"
INTERFACE
        SUBROUTINE write_vtx(file_unit,map)
                INTEGER, INTENT(IN)                             :: file_unit
                DOUBLE PRECISION, DIMENSION(:,:,:), INTENT(IN)  :: map
        END SUBROUTINE
END INTERFACE

! Grid files
WRITE(*,*)
IF (vtx_or_xyz.EQ.0) WRITE(*,*)'Saving grid files in vtx format...'
IF (vtx_or_xyz.EQ.1) WRITE(*,*)'Saving grid files in xyz format...'
OPEN(100,FILE=average_file,STATUS='REPLACE')
OPEN(101,FILE=stdev_file,STATUS='REPLACE')
OPEN(102,FILE=median_file,STATUS='REPLACE')
OPEN(103,FILE=maximum_file,STATUS='REPLACE')
OPEN(104,FILE=harmean_file,STATUS='REPLACE')
OPEN(105,FILE=rms_file,STATUS='REPLACE')
OPEN(106,FILE=entropy_file,STATUS='REPLACE')
OPEN(107,FILE=skewness_file,STATUS='REPLACE')
OPEN(108,FILE=kurtosis_file,STATUS='REPLACE')
OPEN(109,FILE=node_file,STATUS='REPLACE')
IF (vtx_or_xyz.EQ.0) THEN
        CALL write_vtx(100,average)
        CALL write_vtx(101,stdev)
        CALL write_vtx(102,median)
        CALL write_vtx(103,maximum)
        CALL write_vtx(104,harmean)
        CALL write_vtx(105,root_mean_square)
        CALL write_vtx(106,entropy)
        CALL write_vtx(107,skewness)
        CALL write_vtx(108,kurtosis)
        CALL write_vtx(109,node_density/pixel_size)
ELSE IF (vtx_or_xyz.EQ.1) THEN
        DO jj=1,y_number_of_pixels
                DO ii=1,x_number_of_pixels
                        DO kk=1,z_number_of_pixels
                                WRITE(100,outfmt1)x_pixels(ii),y_pixels(jj),z_pixels(kk),average(kk,ii,jj)
                                WRITE(101,outfmt1)x_pixels(ii),y_pixels(jj),z_pixels(kk),stdev(kk,ii,jj)
                                WRITE(102,outfmt1)x_pixels(ii),y_pixels(jj),z_pixels(kk),median(kk,ii,jj)
                                WRITE(103,outfmt1)x_pixels(ii),y_pixels(jj),z_pixels(kk),maximum(kk,ii,jj)
                                WRITE(104,outfmt1)x_pixels(ii),y_pixels(jj),z_pixels(kk),harmean(kk,ii,jj)
                                WRITE(105,outfmt1)x_pixels(ii),y_pixels(jj),z_pixels(kk),root_mean_square(kk,ii,jj)
                                WRITE(106,outfmt1)x_pixels(ii),y_pixels(jj),z_pixels(kk),entropy(kk,ii,jj)
                                WRITE(107,outfmt1)x_pixels(ii),y_pixels(jj),z_pixels(kk),skewness(kk,ii,jj)
                                WRITE(108,outfmt1)x_pixels(ii),y_pixels(jj),z_pixels(kk),kurtosis(kk,ii,jj)
                                WRITE(109,outfmt1)x_pixels(ii),y_pixels(jj),z_pixels(kk),node_density(kk,ii,jj)/pixel_size(kk,ii,jj)
                        END DO
                END DO
        END DO
END IF
CLOSE(100)
CLOSE(101)
CLOSE(102)
CLOSE(103)
CLOSE(104)
CLOSE(105)
CLOSE(106)
CLOSE(107)
CLOSE(108)
CLOSE(109)
WRITE(*,*)'...done!'

! Posterior files
! - number of cells
WRITE(*,*)
WRITE(*,*)'Saving number of cells posterior...'
OPEN(110,FILE=postncells_file,STATUS='REPLACE')
DO nn=ncell_min,ncell_max
        WRITE(110,*)nn,DBLE(posterior_ncells(nn))/DBLE(samples_in_ensemble)
END DO
CLOSE(110)
WRITE(*,*)'...done!'
! - noise
WRITE(*,*)
WRITE(*,*)'Saving noise posterior...'
OPEN(111,FILE=postnoise_file,STATUS='REPLACE')
IF (noise_type.EQ.0.OR.noise_type.EQ.3.OR.noise_type.EQ.6) THEN
        DO ll=1,noise_number_of_pixels
                WRITE(111,outfmt0,ADVANCE='NO')noise_pixels(ll,1)
                DO nn=1,number_of_datasets
                        IF (nn.LT.number_of_datasets) WRITE(111,outfmt0,ADVANCE='NO')&
                                                DBLE(posterior_noise(ll,nn))&
                                              /(DBLE(samples_in_ensemble)*noise_pixel_width(1))
                        IF (nn.EQ.number_of_datasets) WRITE(111,outfmt0,ADVANCE='YES')&
                                                DBLE(posterior_noise(ll,nn))&
                                              /(DBLE(samples_in_ensemble)*noise_pixel_width(1))
                END DO
        END DO
ELSE IF (noise_type.EQ.1.OR.noise_type.EQ.2.OR.noise_type.EQ.5) THEN
        DO ll=1,noise_number_of_pixels
                WRITE(111,outfmt0,ADVANCE='NO')noise_pixels(ll,1)
                DO nn=1,number_of_datasets
                        WRITE(111,outfmt0,ADVANCE='NO')&
                                                DBLE(posterior_noise(ll,nn))&
                                              /(DBLE(samples_in_ensemble)*noise_pixel_width(1))
                END DO
                WRITE(111,outfmt0,ADVANCE='NO')noise_pixels(ll,2)
                DO nn=1,number_of_datasets
                        IF (nn.LT.number_of_datasets) WRITE(111,outfmt0,ADVANCE='NO')&
                                                DBLE(posterior_noise(ll,nn+number_of_datasets))&
                                              /(DBLE(samples_in_ensemble)*noise_pixel_width(2))
                        IF (nn.EQ.number_of_datasets) WRITE(111,outfmt0,ADVANCE='YES')&
                                                DBLE(posterior_noise(ll,nn+number_of_datasets))&
                                              /(DBLE(samples_in_ensemble)*noise_pixel_width(2))
                END DO
        END DO
ELSE IF (noise_type.EQ.4) THEN
        DO ll=1,noise_number_of_pixels
                WRITE(111,outfmt0,ADVANCE='NO')noise_pixels(ll,1)
                WRITE(111,outfmt0,ADVANCE='YES')&
                                                DBLE(posterior_noise(ll,1))&
                                              /(DBLE(samples_in_ensemble)*noise_pixel_width(1))
        END DO
END IF
CLOSE(111)
WRITE(*,*)'...done!'
! - values
WRITE(*,*)
WRITE(*,*)'Saving value posterior...'
OPEN(112,FILE=postvalue_file,STATUS='REPLACE')
WRITE(112,'(1I8)')x_number_of_pixels
WRITE(112,'(1I8)')y_number_of_pixels
WRITE(112,'(1I8)')z_number_of_pixels
WRITE(112,'(1I8)')value_number_of_pixels
DO ii=1,x_number_of_pixels
        WRITE(112,*)x_pixels(ii)
END DO
DO jj=1,y_number_of_pixels
        WRITE(112,*)y_pixels(jj)
END DO
DO kk=1,z_number_of_pixels
        WRITE(112,*)z_pixels(kk)
END DO
DO ll=1,value_number_of_pixels
        WRITE(112,*)value_pixels(ll)
END DO
DO jj=1,y_number_of_pixels
        DO ii=1,x_number_of_pixels
                DO kk=1,z_number_of_pixels
                        DO ll=1,value_number_of_pixels
                                WRITE(112,*)DBLE(posterior_value(kk,ii,jj,ll))/(DBLE(samples_in_ensemble)*value_pixel_width)
                        END DO
                END DO
        END DO
END DO
CLOSE(112)
WRITE(*,*)'...done!'

! Voronoi models
! - first
WRITE(*,*)
WRITE(*,*)'Saving initial models...'
DO ll=1,number_of_chains
        IF (number_of_bad_chains.GT.0) THEN
                IF (ANY(bad_chains.EQ.ll)) THEN
                        WRITE(*,*)'   chain ',TRIM(int2str(ll)),' ignored'
                        CYCLE
                END IF
        END IF
        OPEN(113,FILE=TRIM(first_file) // TRIM(int2str(ll)),STATUS='REPLACE')
        WRITE(113,*)ncell_all_first(ll)
        DO nn=1,ncell_all_first(ll)
                WRITE(113,outfmt2)nn,voronoi_all_first(nn,1,ll),&
                               voronoi_all_first(nn,2,ll),&
                               voronoi_all_first(nn,3,ll),&
                               voronoi_all_first(nn,4,ll)
        END DO
        CLOSE(113)
END DO
WRITE(*,*)'...done!'
! - last
WRITE(*,*)
WRITE(*,*)'Saving final models...'
DO ll=1,number_of_chains
        IF (number_of_bad_chains.GT.0) THEN
                IF (ANY(bad_chains.EQ.ll)) THEN
                        WRITE(*,*)'   chain ',TRIM(int2str(ll)),' ignored'
                        CYCLE
                END IF
        END IF
        OPEN(114,FILE=TRIM(last_file) // TRIM(int2str(ll)),STATUS='REPLACE')
        WRITE(114,*)ncell_all_last(ll)
        DO nn=1,ncell_all_last(ll)
                WRITE(114,outfmt2)nn,voronoi_all_last(nn,1,ll),&
                               voronoi_all_last(nn,2,ll),&
                               voronoi_all_last(nn,3,ll),&
                               voronoi_all_last(nn,4,ll)
        END DO
        CLOSE(114)
END DO
WRITE(*,*)'...done!'

! Average number of cells vs. iteration
WRITE(*,*)
WRITE(*,*)'Saving average number of cells vs. iteration...'
OPEN(115,FILE=ncellsi_file,STATUS='REPLACE')
DO nn=0,samples_so_far
        WRITE(115,outfmt3)nn,ncells_iter(nn+1)
END DO
CLOSE(115)
WRITE(*,*)'...done!'

! Average misfit vs. iteration
WRITE(*,*)
WRITE(*,*)'Saving average misfit vs. iteration...'
OPEN(116,FILE=misfiti_file,STATUS='REPLACE')
DO nn=0,samples_so_far
        WRITE(116,outfmt3)nn,misfit_iter(nn+1)
END DO
CLOSE(116)
WRITE(*,*)'...done!'

! Residuals
WRITE(*,*)
WRITE(*,*)'Saving residuals histogram...'
OPEN(117,FILE=residuals_histo_file,STATUS='REPLACE')
DO nn=1,residuals_number_of_pixels
        WRITE(117,outfmt4)residuals_pixels(nn,1),residuals_histo(nn,1),&
                          residuals_pixels(nn,2),residuals_histo(nn,2)
END DO
CLOSE(117)
WRITE(*,*)'...done!'

END SUBROUTINE
! ********************************************************************* !
! Subroutine to write inversion parameters to file in mksamples.f90
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

! ********************************************************************* !
! Scripts adapted to inversion of preferential pathways:
! Copyright (C) 2026 Léa Lelimouzin
! ********************************************************************* !


SUBROUTINE write_inversion_parameters()

USE chains, ONLY: number_of_chains,samples_so_far,ncell_min,ncell_max,value_min,value_max,&
                  x_min,x_max,y_min,y_max,z_min,z_max,noise_type,number_of_datasets,&
                  a_min,a_max,b_min,b_max,lambda_min,lambda_max,sigma_min,sigma_max,&
                  parallel_tempering,tempering_start,&
                  opt_min,opt_max,opt_parameters
USE input_files, ONLY: modelling_method,forward_file
USE acquisition, ONLY: number_of_measurements,observed_data,data_noise
IMPLICIT NONE
INTEGER         :: nn

! Open file containing inversion parameters
OPEN(UNIT=90,FILE='parameters.txt',STATUS='REPLACE') ! this is just for reference
OPEN(UNIT=80,FILE='parameters.dat',STATUS='REPLACE',FORM='UNFORMATTED') ! this will be read by procsamples

! Write number of chains
WRITE(90,*)'Number of chains'
WRITE(90,*)number_of_chains
WRITE(80)number_of_chains

! Write number of samples
WRITE(90,*)'Number of samples in chains'
WRITE(90,*)samples_so_far
WRITE(80)samples_so_far

! Write prior parameters
WRITE(90,*)'Prior'
WRITE(90,*)ncell_min,ncell_max
WRITE(90,*)value_min,value_max
WRITE(90,*)x_min,x_max
WRITE(90,*)y_min,y_max
WRITE(90,*)z_min,z_max
WRITE(90,*)noise_type,number_of_datasets
WRITE(90,*)a_min,a_max
WRITE(90,*)b_min,b_max
WRITE(90,*)lambda_min,lambda_max
WRITE(90,*)sigma_min,sigma_max
WRITE(90,*)'Parallel tempering'
WRITE(90,*)parallel_tempering,tempering_start
WRITE(90,*)'Forward modelling'
WRITE(90,*)modelling_method
WRITE(90,*)forward_file
IF (opt_parameters.GT.0) THEN 
        WRITE(90,*)'Option : inversion on angle/aperture allowed'
        WRITE(90,*)opt_min,opt_max
END IF
WRITE(80)ncell_min,ncell_max
WRITE(80)value_min,value_max
WRITE(80)x_min,x_max
WRITE(80)y_min,y_max
WRITE(80)z_min,z_max
WRITE(80)noise_type,number_of_datasets
WRITE(80)a_min,a_max
WRITE(80)b_min,b_max
WRITE(80)lambda_min,lambda_max
WRITE(80)sigma_min,sigma_max
WRITE(80)parallel_tempering,tempering_start
WRITE(80)modelling_method
WRITE(80)forward_file
IF (opt_parameters.GT.0) THEN 
        WRITE(80)opt_min,opt_max
END IF

! Write observed data
WRITE(90,*)'Data'
WRITE(90,*)number_of_measurements
WRITE(80)number_of_measurements
DO nn=1,number_of_measurements
        WRITE(90,*)observed_data(nn),data_noise(nn)
        WRITE(80)observed_data(nn),data_noise(nn)
END DO

CLOSE(90)
CLOSE(80)

END SUBROUTINE


SUBROUTINE read_inversion_parameters()

USE chains, ONLY: number_of_chains,samples_so_far,ncell_min,ncell_max,value_min,value_max,&
                  x_min,x_max,y_min,y_max,z_min,z_max,noise_type,number_of_datasets,&
                  a_min,a_max,b_min,b_max,lambda_min,lambda_max,sigma_min,sigma_max,&
                  parallel_tempering,tempering_start,&
                  opt_min,opt_max,opt_parameters
USE input_files, ONLY: modelling_method,forward_file,stop_sampling
USE acquisition, ONLY: number_of_measurements,observed_data,data_noise
IMPLICIT NONE
INTEGER         :: nn

! Open file containing inversion parameters
OPEN(UNIT=80,FILE='parameters.dat',STATUS='OLD',FORM='UNFORMATTED')

! Read number of chains
READ(80)number_of_chains

! Read number of samples
READ(80)samples_so_far
IF (stop_sampling.NE.-1) samples_so_far=stop_sampling

! Read prior parameters
READ(80)ncell_min,ncell_max
READ(80)value_min,value_max
READ(80)x_min,x_max
READ(80)y_min,y_max
READ(80)z_min,z_max
READ(80)noise_type,number_of_datasets
READ(80)a_min,a_max
READ(80)b_min,b_max
READ(80)lambda_min,lambda_max
READ(80)sigma_min,sigma_max
READ(80)parallel_tempering,tempering_start
READ(80)modelling_method
READ(80)forward_file
IF (opt_parameters.GT.0) THEN 
        READ(80)opt_min,opt_max
END IF

! Read observed data
READ(80)number_of_measurements
ALLOCATE(observed_data(number_of_measurements))
ALLOCATE(data_noise(number_of_measurements))
DO nn=1,number_of_measurements
        READ(80)observed_data(nn),data_noise(nn)
END DO

CLOSE(80)

END SUBROUTINE

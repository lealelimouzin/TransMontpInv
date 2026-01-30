! ********************************************************************* !
! Module to calculate acceptance parameter "alpha" in mksamples.f90
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
! Scripts adapted to inversion of preferential pathways in TransMontpInv
! package.
! Copyright (C) 2026 Léa Lelimouzin
! ********************************************************************* !


MODULE alpha_calculation

USE chains, ONLY: chain_name,l_change
USE debug
USE input_files, ONLY: modelling_method
IMPLICIT NONE


CONTAINS


SUBROUTINE calculate_alpha(likelihood_current,likelihood_proposed,&
                           likelihood_normalization,temper,alpha)

USE chains, ONLY: alpha_component,birth,death,move,value,noise,value_min,value_max,&
                  no_data,temperature,option
IMPLICIT NONE
DOUBLE PRECISION, INTENT(IN)            :: likelihood_current,likelihood_proposed,likelihood_normalization
LOGICAL, INTENT(IN)                     :: temper
DOUBLE PRECISION, INTENT(OUT)           :: alpha

IF (debug_mode) WRITE(debug_unit,*)'--- Chain ',TRIM(chain_name),': calculating alpha ---'


IF (modelling_method.EQ.1) THEN 
        IF (birth) THEN
                IF (no_data.EQ.0) THEN
                        IF (temper) THEN
                                alpha=MINVAL([DBLE(0),alpha_component&
                                                -LOG(value_max-value_min)&
                                                +(-likelihood_proposed+likelihood_current+likelihood_normalization)/temperature])
                        ELSE
                                alpha=MINVAL([DBLE(0),alpha_component&
                                                -LOG(value_max-value_min)&
                                                -likelihood_proposed+likelihood_current+likelihood_normalization])
                        END IF
                ELSE IF (no_data.EQ.1) THEN
                        alpha=MINVAL([DBLE(0),alpha_component-LOG(value_max-value_min)])
                END IF
        ELSE IF (death) THEN
                IF (no_data.EQ.0) THEN
                        IF (temper) THEN
                                alpha=MINVAL([DBLE(0),alpha_component&
                                                +LOG(value_max-value_min)&
                                                +(-likelihood_proposed+likelihood_current+likelihood_normalization)/temperature])
                        ELSE
                                alpha=MINVAL([DBLE(0),alpha_component&
                                                +LOG(value_max-value_min)&
                                                -likelihood_proposed+likelihood_current+likelihood_normalization])
                        END IF
                ELSE IF (no_data.EQ.1) THEN
                        alpha=MINVAL([DBLE(0),alpha_component-LOG(value_max-value_min)])
                END IF
        ELSE IF (move.OR.value.OR.noise) THEN
                IF (no_data.EQ.0) THEN
                        IF (temper) THEN
                                alpha=MINVAL([DBLE(0),(-likelihood_proposed+likelihood_current&
                                                        +likelihood_normalization)/temperature])
                        ELSE
                                alpha=MINVAL([DBLE(0),-likelihood_proposed+likelihood_current+likelihood_normalization])
                        END IF
                ELSE IF (no_data.EQ.1) THEN
                        alpha=DBLE(0)
                END IF
        END IF

        IF (debug_mode) WRITE(debug_unit,*)'--- Chain ',TRIM(chain_name),': calculated alpha ---'

ELSEIF (modelling_method.EQ.2) THEN
        IF (birth) THEN
                IF (no_data.EQ.0) THEN
                        IF (temper) THEN
                                alpha=MINVAL([DBLE(0),(-likelihood_proposed+likelihood_current&
                                                        +likelihood_normalization)/temperature])
                        ELSE
                                alpha=MINVAL([DBLE(0),-likelihood_proposed+likelihood_current+likelihood_normalization])
                        END IF
                ELSE IF (no_data.EQ.1) THEN
                        alpha=MINVAL([DBLE(0),-LOG(value_max-value_min)])
                END IF
                write(*,*) 'alpha (add) at chain ', TRIM(chain_name), ' is : ', alpha
        ELSE IF (death) THEN
                IF (no_data.EQ.0) THEN
                        IF (temper) THEN
                                alpha=MINVAL([DBLE(0),(-likelihood_proposed+likelihood_current&
                                                        +likelihood_normalization)/temperature])
                        ELSE
                                alpha=MINVAL([DBLE(0),-likelihood_proposed+likelihood_current+likelihood_normalization])
                        END IF
                ELSE IF (no_data.EQ.1) THEN
                        alpha=MINVAL([DBLE(0),-LOG(value_max-value_min)])
                END IF
        ELSE IF (move.OR.value.OR.noise.OR.option) THEN
                IF (no_data.EQ.0) THEN
                        IF (temper) THEN
                                alpha=MINVAL([DBLE(0),(-likelihood_proposed+likelihood_current+likelihood_normalization)&
                                                        /temperature])
                        ELSE
                                alpha=MINVAL([DBLE(0),-likelihood_proposed+likelihood_current+likelihood_normalization])
                        END IF
                ELSE IF (no_data.EQ.1) THEN
                        alpha=DBLE(0)
                END IF
        END IF
END IF

END SUBROUTINE


SUBROUTINE calculate_alpha_delayed(likelihood_current,likelihood_proposed,likelihood_delayed,&
                                   likelihood_normalization,temper,alpha_delayed)

USE chains, ONLY: move,value,inside_prior,inside_prior_delayed,no_data,&
                  voronoi_current,voronoi_proposed,voronoi_delayed,cell_number,&
                  value_change,x_change,y_change,z_change,alpha,temperature,&
                  option,opt_change
USE random_generator

IMPLICIT NONE
DOUBLE PRECISION, INTENT(IN)    :: likelihood_current,likelihood_proposed,likelihood_delayed,likelihood_normalization
LOGICAL, INTENT(IN)             :: temper
DOUBLE PRECISION, INTENT(OUT)   :: alpha_delayed
DOUBLE PRECISION                :: value_min,value_max
DOUBLE PRECISION                :: TT1,TT2,TT3,TT2x,TT2y,TT2z,TT3t,TT3b

IF (debug_mode) WRITE(debug_unit,*)'--- Chain ',TRIM(chain_name),': calculating alpha delayed ---'

IF (no_data.EQ.0) THEN
        IF (temper) THEN
                TT1=(-likelihood_delayed+likelihood_current+likelihood_normalization)/temperature
        ELSE
                TT1=-likelihood_delayed+likelihood_current+likelihood_normalization
        END IF
ELSE IF (no_data.EQ.1) THEN
        TT1=DBLE(0)
END IF

IF (move) THEN
        TT2x=((voronoi_proposed(cell_number,1)-voronoi_current(cell_number,1))**2&
                -(voronoi_proposed(cell_number,1)-voronoi_delayed(cell_number,1))**2)/(2*x_change**2)
        TT2y=((voronoi_proposed(cell_number,2)-voronoi_current(cell_number,2))**2&
                -(voronoi_proposed(cell_number,2)-voronoi_delayed(cell_number,2))**2)/(2*y_change**2)
        TT2z=((voronoi_proposed(cell_number,3)-voronoi_current(cell_number,3))**2&
                -(voronoi_proposed(cell_number,3)-voronoi_delayed(cell_number,3))**2)/(2*z_change**2)
        TT2=TT2x+TT2y+TT2z
ELSE IF (value) THEN
        TT2=(-(voronoi_proposed(cell_number,4)-voronoi_delayed(cell_number,4))**2&
                +(voronoi_proposed(cell_number,4)-voronoi_current(cell_number,4))**2)/(2*l_change**2)
ELSE IF (option) THEN
        TT2=(-(voronoi_proposed(cell_number,5)-voronoi_delayed(cell_number,5))**2&
                +(voronoi_proposed(cell_number,5)-voronoi_current(cell_number,5))**2)/(2*opt_change**2)
END IF



IF (.NOT.inside_prior.AND.inside_prior_delayed) THEN ! if the first proposal was rejected because out of bounds, then TT3=0
        TT3=0
        alpha_delayed=MINVAL([DBLE(0),TT1+TT2])
ELSE IF (inside_prior.AND.inside_prior_delayed) THEN ! if both the first and the delayed proposal are within the bounds
        IF (no_data.EQ.0) THEN
                IF (temper) THEN
                        TT3t=LOG(1-EXP(MINVAL([DBLE(0),&
                                                (-likelihood_proposed+likelihood_delayed+likelihood_normalization)/temperature])))
                ELSE
                        TT3t=LOG(1-EXP(MINVAL([DBLE(0),-likelihood_proposed+likelihood_delayed+likelihood_normalization])))
                END IF
        END IF
        TT3b=LOG(1-EXP(alpha))
        TT3=TT3t-TT3b
        alpha_delayed=MINVAL([DBLE(0),TT1+TT2+TT3])
END IF

IF (debug_mode) WRITE(debug_unit,*)'--- Chain ',TRIM(chain_name),': calculated alpha delayed ---'

END SUBROUTINE


SUBROUTINE calculate_likelihood(number_of_measurements,observed_data,modelled_data,uncertainty,likelihood_type,likelihood,misfit)

IMPLICIT NONE
INTEGER, INTENT(IN)                             :: number_of_measurements
DOUBLE PRECISION, DIMENSION(:), INTENT(IN)      :: observed_data,modelled_data,uncertainty
INTEGER, INTENT(IN)                             :: likelihood_type
DOUBLE PRECISION, INTENT(OUT)                   :: likelihood,misfit
INTEGER                                         :: ii

IF (debug_mode) WRITE(debug_unit,*)'--- Chain ',TRIM(chain_name),': calculating likelihood ---'

likelihood=0
misfit=0

DO ii=1,number_of_measurements

        IF (modelling_method.EQ.2) THEN
                IF (modelled_data(ii).GT.0) THEN
                        IF (likelihood_type.EQ.0) likelihood=likelihood+(observed_data(ii)-modelled_data(ii))**2&
                                                                /(2*(uncertainty(ii)**2))
                        IF (likelihood_type.EQ.1) likelihood=likelihood+ABS((observed_data(ii)-modelled_data(ii))&
                                                                                /uncertainty(ii))

                        misfit=misfit+((observed_data(ii)-modelled_data(ii))/uncertainty(ii))**2
                END IF
        ELSE
                IF (likelihood_type.EQ.0) likelihood=likelihood+(observed_data(ii)-modelled_data(ii))**2/(2*(uncertainty(ii)**2))
                IF (likelihood_type.EQ.1) likelihood=likelihood+ABS((observed_data(ii)-modelled_data(ii))/uncertainty(ii))
                
                misfit=misfit+((observed_data(ii)-modelled_data(ii))/uncertainty(ii))**2
        END IF
        
END DO

IF (debug_mode) WRITE(debug_unit,*)'--- Chain ',TRIM(chain_name),': calculated likelihood ---'

END SUBROUTINE


SUBROUTINE calculate_normalization(number_of_measurements,modelled_data,uncertainty_current,uncertainty_proposed,normalization)

IMPLICIT NONE
INTEGER, INTENT(IN)                             :: number_of_measurements
DOUBLE PRECISION, DIMENSION(:), INTENT(IN)      :: uncertainty_current,uncertainty_proposed
DOUBLE PRECISION, DIMENSION(:), INTENT(IN)      :: modelled_data
DOUBLE PRECISION, INTENT(OUT)                   :: normalization
INTEGER                                         :: ii

normalization=0

DO ii=1,number_of_measurements
        IF (modelled_data(ii).GT.0) THEN
                normalization=normalization+LOG(uncertainty_current(ii)/uncertainty_proposed(ii))
        END IF
END DO

END SUBROUTINE

END MODULE
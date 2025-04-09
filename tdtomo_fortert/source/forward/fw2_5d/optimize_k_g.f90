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


SUBROUTINE optimize_k_g()

USE modelling_2_5D, ONLY: k_parameters,g_parameters,number_of_coefficients,x_spacings,z_spacings,current_locations,number_of_sources
IMPLICIT NONE
INTEGER                                         :: max_iterations,max_radii,line_search_steps,number_of_radii
DOUBLE PRECISION                                :: line_search_low,line_search_high
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE     :: line_search,line_search_results
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE     :: x_radius,z_radius
INTEGER                                         :: ii,jj,nn,ini,fin
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE     :: x_potential_electrode,z_potential_electrode
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE   :: radii,workd4,dv_dk
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE     :: radii_inv,workd1
LOGICAL, DIMENSION(:), ALLOCATABLE              :: mask
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE     :: k_initial,k_updated,k_temp
DOUBLE PRECISION, ALLOCATABLE                   :: A(:,:)
DOUBLE PRECISION, ALLOCATABLE                   :: g_initial(:),g_updated(:),g_temp(:)
DOUBLE PRECISION, ALLOCATABLE                   :: v_initial(:),v_updated(:),v_temp(:)
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE     :: objective_function,dk
DOUBLE PRECISION                                :: obj,conditioning,reduction
INTEGER                                         :: iteration_counter,stopping
INTEGER                                         :: info
INTEGER, DIMENSION(:), ALLOCATABLE              :: ipiv
INTERFACE
        SUBROUTINE trim_radii(radii,number_of_radii)
                DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT)    :: radii
                INTEGER, INTENT(OUT)                                            :: number_of_radii
        END SUBROUTINE
END INTERFACE
INTERFACE
        SUBROUTINE linspace(d1,d2,n,grid)
                INTEGER, INTENT(IN)                             :: n
                DOUBLE PRECISION, INTENT(IN)                    :: d1, d2
                DOUBLE PRECISION, DIMENSION(n), INTENT(OUT)     :: grid
        END SUBROUTINE
END INTERFACE
INTERFACE
        SUBROUTINE create_A(number_of_coefficients,number_of_radii,radii,radii_inv,k_values,A)
                INTEGER, INTENT(IN)                                             :: number_of_coefficients,number_of_radii
                DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN)                    :: radii
                DOUBLE PRECISION, DIMENSION(:), INTENT(IN)                      :: radii_inv,k_values
                DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT)      :: A
        END SUBROUTINE
END INTERFACE
INTERFACE
        SUBROUTINE estimate_g(number_of_coefficients,A,g_values)
                INTEGER, INTENT(IN)                                             :: number_of_coefficients
                DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN)                    :: A
                DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, INTENT(OUT)        :: g_values
        END SUBROUTINE
END INTERFACE
INTERFACE
        SUBROUTINE calculate_objective_function(number_of_radii,A,g_values,v,obj)
                INTEGER, INTENT(IN)                                             :: number_of_radii
                DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN)                    :: A
                DOUBLE PRECISION, DIMENSION(:), INTENT(IN)                      :: g_values
                DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, INTENT(OUT)        :: v
                DOUBLE PRECISION, INTENT(OUT)                                   :: obj
        END SUBROUTINE
END INTERFACE
INTERFACE
        FUNCTION EYE(n,value)
                INTEGER, INTENT(IN)                     :: n
                DOUBLE PRECISION, INTENT(IN)            :: value
                DOUBLE PRECISION, DIMENSION(n,n)        :: EYE
        END FUNCTION
END INTERFACE
INTERFACE
        SUBROUTINE check_conditioning(M,conditioning)
                DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN)    :: M
                DOUBLE PRECISION, INTENT(OUT)                   :: conditioning
        END SUBROUTINE
END INTERFACE

ALLOCATE(k_parameters(number_of_coefficients),g_parameters(number_of_coefficients))

max_iterations=25 ! maximum number of iterations for the optimization routine
ALLOCATE(objective_function(max_iterations+1))
max_radii=2000 ! max number of radii to search over
line_search_steps=10 ! number of line search steps
line_search_low=0.01 ! lower bound for line search
line_search_high=1 ! upper bound for line search
ALLOCATE(line_search(line_search_steps))
CALL linspace(line_search_low,line_search_high,line_search_steps,line_search)

! Hard wired search radius for determining k and g
ALLOCATE(x_radius(14),z_radius(14))
x_radius(1:7)=0
x_radius(8:14)=(/ DBLE(0.1),DBLE(0.5),DBLE(1),DBLE(5),DBLE(10),DBLE(20),DBLE(30) /)
z_radius=x_radius(14:1:-1)

ALLOCATE(x_potential_electrode(SIZE(x_radius)),z_potential_electrode(SIZE(x_radius)))

ALLOCATE(radii(SIZE(x_radius)*number_of_sources,4))

fin=0
DO nn=1,number_of_sources
        
        ini=fin+1
        fin=ini-1+SIZE(x_radius)
        
        x_potential_electrode=x_radius+current_locations(nn,1)
        z_potential_electrode=z_radius+current_locations(nn,2)
        
        ! Norm of positive current electrode and 1st potential electrode
        radii(ini:fin,1)=SQRT((x_potential_electrode-current_locations(nn,1))**2&
                             +(z_potential_electrode-current_locations(nn,2))**2)
        
        ! Norm of negative current electrode and 1st potential electrode
        radii(ini:fin,2)=SQRT((x_potential_electrode-current_locations(nn,3))**2&
                             +(z_potential_electrode-current_locations(nn,4))**2)
        
        ! Norm of imaginary positive current electrode and 1st potential electrode
        radii(ini:fin,3)=SQRT((x_potential_electrode-current_locations(nn,1))**2&
                             +(z_potential_electrode+current_locations(nn,2))**2)
        
        ! Norm of imaginary negative current electrode and 1st potential electrode
        radii(ini:fin,4)=SQRT((x_potential_electrode-current_locations(nn,3))**2&
                             +(z_potential_electrode+current_locations(nn,4))**2)
        
END DO

! Keep only unique radii
CALL trim_radii(radii,number_of_radii)

! Initialize a starting guess for k
ALLOCATE(k_initial(number_of_coefficients))
CALL linspace(DBLE(-2),DBLE(0.5),number_of_coefficients,k_initial)
k_initial=10**k_initial

! Set up a matrix of radii
ALLOCATE(radii_inv(number_of_radii))
radii_inv=1/(1/radii(:,1)-1/radii(:,2)+1/radii(:,3)-1/radii(:,4))

! Remove any divide by zero
ALLOCATE(mask(number_of_radii))
mask=.TRUE.
DO nn=1,number_of_radii
        IF (SUM(1/radii(nn,:))+radii_inv(nn)-1.EQ.SUM(1/radii(nn,:))+radii_inv(nn)) THEN
                mask(nn)=.FALSE.
        END IF
END DO
number_of_radii=SIZE(PACK(radii_inv,mask))
ALLOCATE(workd4(number_of_radii,4),workd1(number_of_radii))
workd4(:,1)=PACK(radii(:,1),mask)
workd4(:,2)=PACK(radii(:,2),mask)
workd4(:,3)=PACK(radii(:,3),mask)
workd4(:,4)=PACK(radii(:,4),mask)
workd1=PACK(radii_inv,mask)
DEALLOCATE(mask,radii,radii_inv)
ALLOCATE(radii(number_of_radii,4),radii_inv(number_of_radii))
radii=workd4
radii_inv=workd1
DEALLOCATE(workd4,workd1)

! Calculate matrix A
CALL create_A(number_of_coefficients,number_of_radii,radii,radii_inv,k_initial,A)

! Estimate g for the initial k values
CALL estimate_g(number_of_coefficients,A,g_initial)

! Calculate objective function for the initial k and g estimates
CALL calculate_objective_function(number_of_radii,A,g_initial,v_initial,obj)
objective_function(1)=obj

! Start counter and initialize the optimization
iteration_counter=1
ALLOCATE(k_updated(number_of_coefficients))
k_updated=k_initial
stopping=0 ! stopping toggle in case A becomes ill-conditioned
reduction=1 ! variable to ensure sufficent decrease between iterations
ALLOCATE(dv_dk(number_of_radii,number_of_coefficients))
ALLOCATE(dk(number_of_coefficients))
ALLOCATE(ipiv(number_of_coefficients))

! Optimization terminates if the objective function is not reduced by at least 5% at each iteration
DO WHILE (objective_function(iteration_counter).GT.DBLE(1e-5)&
                .AND.iteration_counter.LT.max_iterations&
                .AND.stopping.EQ.0.AND.reduction.GT.DBLE(0.05))
        
        ! Create the derivative matrix
        ALLOCATE(k_temp(number_of_coefficients))
        DO nn=1,number_of_coefficients
                
                ! Calculate a new matrix A
                k_temp=k_initial
                k_temp(nn)=k_temp(nn)*DBLE(1.05)
                CALL create_A(number_of_coefficients,number_of_radii,radii,radii_inv,k_temp,A)
                
                ! Estimate g for the updated k values
                CALL estimate_g(number_of_coefficients,A,g_updated)
                
                ! Calculate objective function for the updated k and g estimates
                CALL calculate_objective_function(number_of_radii,A,g_updated,v_updated,obj)
                
                dv_dk(:,nn)=(v_updated-v_initial)/(k_temp(nn)-k_initial(nn))
                
        END DO
        DEALLOCATE(k_temp)
        
        ! Apply some smallness regularization
        dk=MATMUL(TRANSPOSE(dv_dk),1-v_initial)+DBLE(1e-8)*k_updated
        CALL DGESV(number_of_coefficients,1,MATMUL(TRANSPOSE(dv_dk),dv_dk)+EYE(number_of_coefficients,&
                   DBLE(1e-8)),number_of_coefficients,ipiv,dk,number_of_coefficients,info)
        
        ! Perform a line-search to maximize the descent
        ALLOCATE(line_search_results(line_search_steps))
        DO jj=1,line_search_steps
                
                ! Calculate a new matrix A
                ALLOCATE(k_temp(number_of_coefficients))
                k_temp=k_updated+line_search(jj)*dk
                CALL create_A(number_of_coefficients,number_of_radii,radii,radii_inv,k_temp,A)
                DEALLOCATE(k_temp)
                
                ! Estimate g for the temporary k values
                CALL estimate_g(number_of_coefficients,A,g_temp)
                
                ! Calculate objective function for the temporary k and g estimates
                CALL calculate_objective_function(number_of_radii,A,g_temp,v_temp,obj)
                line_search_results(jj)=obj
                
        END DO
        
        ! Find the smallest objective function from the line-search
        ii=MINLOC(line_search_results,1)
        DEALLOCATE(line_search_results)
        
        ! Create a new guess for k
        k_updated=k_updated+line_search(ii)*dk
        
        ! Calculate a new matrix A
        CALL create_A(number_of_coefficients,number_of_radii,radii,radii_inv,k_updated,A)
        
        ! Estimate g for the new k values
        CALL estimate_g(number_of_coefficients,A,g_updated)
        
        ! Calculate objective function for the new k and g estimates
        CALL calculate_objective_function(number_of_radii,A,g_updated,v_updated,obj)
        objective_function(iteration_counter+1)=obj
        
        reduction=objective_function(iteration_counter)/objective_function(iteration_counter+1)-1
        iteration_counter=iteration_counter+1
        
        k_initial=k_updated
        v_initial=v_updated
        
        ! Check the conditioning of the matrix
        CALL check_conditioning(MATMUL(TRANSPOSE(A),A),conditioning)
        
        IF (conditioning.LT.DBLE(1e-20)) THEN
                k_updated=k_updated-line_search(ii)*dk
                stopping=1
        END IF
        
END DO

! Get final k parameters
k_parameters=ABS(k_updated)

! Calculate a new matrix A
CALL create_A(number_of_coefficients,number_of_radii,radii,radii_inv,k_parameters,A)

! Get final g parameters
CALL estimate_g(number_of_coefficients,A,g_parameters)

DEALLOCATE(x_radius,z_radius)
DEALLOCATE(x_potential_electrode,z_potential_electrode)
DEALLOCATE(radii,radii_inv)
DEALLOCATE(k_initial,k_updated)
DEALLOCATE(A)
DEALLOCATE(g_initial,g_updated,g_temp)
DEALLOCATE(v_initial,v_updated,v_temp)
DEALLOCATE(objective_function)
DEALLOCATE(dv_dk,dk)
DEALLOCATE(line_search)
DEALLOCATE(ipiv)

END SUBROUTINE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE trim_radii(radii,number_of_radii)

! This subroutine trims the list of radii by keeping only unique
! combinations.
! 
! The algorithm for this subroutine was inspired by 
! http://stackoverflow.com/questions/14137610/remove-repeated-elements-on-an-2d-array-in-fortran-90
!
! Erica Galetti, November 2016
! erica.galetti@ed.ac.uk

IMPLICIT NONE
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT)    :: radii
INTEGER, INTENT(OUT)                                            :: number_of_radii
INTEGER, DIMENSION(:), ALLOCATABLE                              :: index_vector
LOGICAL, DIMENSION(:), ALLOCATABLE                              :: mask
INTEGER                                                         :: ii,length_radii
INTEGER, DIMENSION(:), ALLOCATABLE                              :: radii_index
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE                   :: workd

length_radii=SIZE(radii(:,1))

ALLOCATE(mask(length_radii))
mask=.TRUE.
ALLOCATE(radii_index(length_radii))
radii_index=(/ (ii,ii=1,length_radii) /)

DO ii=length_radii,2,-1
        mask(ii)=.NOT.(ANY(radii(1:ii-1,1).EQ.radii(ii,1).AND.radii(1:ii-1,2).EQ.radii(ii,2)&
                        .AND.radii(1:ii-1,3).EQ.radii(ii,3).AND.radii(1:ii-1,4).EQ.radii(ii,4)))
END DO

number_of_radii=SIZE(PACK(radii_index,mask))
ALLOCATE(index_vector(number_of_radii))
ALLOCATE(workd(number_of_radii,4))
index_vector=PACK(radii_index,mask)
workd=radii(index_vector,:)

DEALLOCATE(radii)
ALLOCATE(radii(number_of_radii,4))
radii=workd

DEALLOCATE(index_vector)
DEALLOCATE(mask,radii_index)
DEALLOCATE(workd)

END SUBROUTINE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE create_A(number_of_coefficients,number_of_radii,radii,radii_inv,k_values,A)

IMPLICIT NONE
INTEGER, INTENT(IN)                                             :: number_of_coefficients,number_of_radii
DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN)                    :: radii
DOUBLE PRECISION, DIMENSION(:), INTENT(IN)                      :: radii_inv,k_values
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT)      :: A
INTEGER                                                         :: ii,jj
DOUBLE PRECISION                                                :: BESSK

IF (ALLOCATED(A)) DEALLOCATE(A)

ALLOCATE(A(number_of_radii,number_of_coefficients))

DO jj=1,number_of_coefficients
        DO ii=1,number_of_radii
                A(ii,jj)=radii_inv(ii)*(BESSK(0,ABS(radii(ii,1)*k_values(jj)))&
                                       -BESSK(0,ABS(radii(ii,2)*k_values(jj)))&
                                       +BESSK(0,ABS(radii(ii,3)*k_values(jj)))&
                                       -BESSK(0,ABS(radii(ii,4)*k_values(jj))))
        END DO
END DO

END SUBROUTINE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE estimate_g(number_of_coefficients,A,g_values)

IMPLICIT NONE
INTEGER, INTENT(IN)                                             :: number_of_coefficients
DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN)                    :: A
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, INTENT(OUT)        :: g_values
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE                   :: L
INTEGER                                                         :: info
INTEGER, DIMENSION(:), ALLOCATABLE                              :: ipiv
INTEGER                                                         :: jj

IF (ALLOCATED(g_values)) DEALLOCATE(g_values)

ALLOCATE(g_values(number_of_coefficients))
DO jj=1,number_of_coefficients
        g_values(jj)=SUM(A(:,jj))
END DO

ALLOCATE(L(number_of_coefficients,number_of_coefficients))
L=MATMUL(TRANSPOSE(A),A)

ALLOCATE(ipiv(number_of_coefficients))

CALL DGESV(number_of_coefficients,1,L,number_of_coefficients,ipiv,g_values,number_of_coefficients,info)

DEALLOCATE(L,ipiv)

END SUBROUTINE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE calculate_objective_function(number_of_radii,A,g_values,v,obj)

IMPLICIT NONE
INTEGER, INTENT(IN)                                             :: number_of_radii
DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN)                    :: A
DOUBLE PRECISION, DIMENSION(:), INTENT(IN)                      :: g_values
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, INTENT(OUT)        :: v
DOUBLE PRECISION, INTENT(OUT)                                   :: obj

IF (ALLOCATED(v)) DEALLOCATE(v)

ALLOCATE(v(number_of_radii))

v=MATMUL(A,g_values)

obj=SUM((1-v)**2)

END SUBROUTINE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION EYE(n,value)

IMPLICIT NONE
INTEGER, INTENT(IN)                     :: n
DOUBLE PRECISION, INTENT(IN)            :: value
DOUBLE PRECISION, DIMENSION(n,n)        :: EYE
INTEGER                                 :: ii,jj

EYE=0

DO jj=1,n
        DO ii=1,n
                IF (ii.EQ.jj) EYE(ii,jj)=value
        END DO
END DO

RETURN

END FUNCTION


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE check_conditioning(M,conditioning)

IMPLICIT NONE
DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN)    :: M
DOUBLE PRECISION, INTENT(OUT)                   :: conditioning
DOUBLE PRECISION                                :: anorm
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE     :: work
INTEGER, DIMENSION(:), ALLOCATABLE              :: iwork
INTEGER                                         :: info

ALLOCATE(work(4*SIZE(M,2)))
ALLOCATE(iwork(SIZE(M,2)))

anorm=MAXVAL(SUM(ABS(M),1))

CALL DGECON('O',SIZE(M,2),M,SIZE(M,1),anorm,conditioning,work,iwork,info)

DEALLOCATE(work,iwork)

END SUBROUTINE


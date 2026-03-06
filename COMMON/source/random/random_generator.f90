! ********************************************************************* !
! Module for random number generation
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

MODULE random_generator

USE mt19937 
IMPLICIT NONE
INTEGER                 :: random_counter
PUBLIC                  :: unirand,gausrand,powerrand


CONTAINS


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Random counter initialisation routine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE initialize_random_generator(seed)

IMPLICIT NONE
INTEGER                 :: seed

CALL init_genrand(seed)
random_counter=0

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Generate a random number from a uniform distribution
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION unirand() RESULT(value)

IMPLICIT NONE
DOUBLE PRECISION        :: value

value=grnd()
random_counter=random_counter+1

RETURN
END FUNCTION

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Generate a random number from a Gaussian distribution
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION gausrand() RESULT(value)

IMPLICIT NONE
DOUBLE PRECISION        :: value

value=GASDEV()

RETURN
END FUNCTION

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Random number generator for a Gaussian distribution
! (adapted from Numerical Recipes)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION GASDEV()
      
DOUBLE PRECISION        :: GASDEV
DOUBLE PRECISION        :: v1,v2,r,fac
      
!IF (idum.lt.0) iset=0
10 v1=2*unirand()-1
v2=2*unirand()-1

r=v1**2+v2**2

IF(r.ge.1.or.r.eq.0) GOTO 10

fac=sqrt(-2*log(r)/r)

GASDEV=v2*fac

RETURN
END FUNCTION

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Generate a random number from a power law distribution (fracture length estimation)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION powerrand(value_min,value_max,a) RESULT(value)

IMPLICIT NONE
DOUBLE PRECISION        :: value_min, value_max, value, a 
!value=value_min*unirand()**(1/(1-a)) ! if truncated only by a l min
value=(unirand()*(value_max**(1-a) - value_min**(1-a)) + value_min**(1-a))**(1/(1-a))
write(*,*) value
!random_counter=random_counter+1

RETURN
END FUNCTION

FUNCTION powerrand_k(value_min,value_max,a) RESULT(powerlaw_k)

    IMPLICIT NONE
    DOUBLE PRECISION        :: value_min,value_max,a,powerlaw_k
    powerlaw_k=(1-a)/(value_max**(1-a) - value_min**(1-a))
    
RETURN
END FUNCTION

FUNCTION powerrand_std(value_min,value_max,a) RESULT(powerlaw_std)
    
    IMPLICIT NONE
    DOUBLE PRECISION        :: value_min,value_max,part1,part2,part3
    DOUBLE PRECISION        :: a, powerlaw_std

    if (a.EQ.2) then
        powerlaw_std = sqrt( value_min*value_max - ((log(value_max)-log(value_min))/(value_min**(-1)-value_max**(-1)))**2 )
    ELSEIF ((a.NE.2).AND.(a.NE.3)) then
        part1 = (1-a)/(value_max**(1-a)-value_min**(1-a))
        part2 = (value_max**(3-a)-value_min**(3-a))/(3-a)
        part3= (value_max**(2-a)-value_min**(2-a)/(2-a))**2*(1-a)/(value_max**(1-a)-value_min**(1-a))
        powerlaw_std= sqrt(part1*(part2-part3))
    ELSE
        write(*,*) 'When a>3, the standard deviation of a power-law truncated at both ends has not been implemented.'
    end if 
    
    RETURN
END FUNCTION

END MODULE
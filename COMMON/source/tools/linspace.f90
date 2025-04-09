!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Downloaded from https://github.com/lucaguerrieri/occbin_paper/tree/master/model_irrcap/Fortran/dp/sargent_fortran
! (direct link: https://github.com/lucaguerrieri/occbin_paper/blob/master/model_irrcap/Fortran/dp/sargent_fortran/linspace.f90)
! 13/11/2016
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE linspace(d1,d2,n,grid)

IMPLICIT NONE

INTEGER, INTENT(IN)                             :: n
DOUBLE PRECISION, INTENT(IN)                    :: d1,d2
DOUBLE PRECISION, DIMENSION(n), INTENT(OUT)     :: grid
INTEGER                                         :: indxi

grid(1) = d1
DO indxi= 0,n-2
        grid(indxi+1) = d1+(DBLE(indxi)*(d2-d1))/DBLE(n-1)
END DO
grid(n) = d2

END SUBROUTINE
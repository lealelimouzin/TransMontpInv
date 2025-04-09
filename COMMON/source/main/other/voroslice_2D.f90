! ********************************************************************* !
! Program to convert a Voronoi model to a GMT-friendly format
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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This program is designed to convert velocity from a file    !
! containing coordinates and velocities of Voronoi cells into !
! a form suitable for input to GMT.                           !
!                                                             !
! Adapted from tslicess.f90 in the FMST code package to plot  !
! Voronoi cells avoiding interpolation.                       !
!                                                             !
! Erica Galetti, March 2017                                   !
! erica.galetti@ed.ac.uk                                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

PROGRAM slice
USE voronoi_operations
IMPLICIT NONE
INTEGER :: i,j
INTEGER :: checkstat
INTEGER :: nnt,nnp
INTEGER :: nnx,nnz
INTEGER :: ddt,ddp
INTEGER :: prp,nr,nre
INTEGER, PARAMETER :: i5=SELECTED_REAL_KIND(5,10)
REAL(KIND=i5) :: rlat,rlon
DOUBLE PRECISION :: gotn,gots,gopw,gope,rgst,rgsp
DOUBLE PRECISION :: lft,rgt,btm,top
DOUBLE PRECISION :: rd1,rd2
CHARACTER (LEN=30) :: sep
CHARACTER (LEN=30) :: ifilev
CHARACTER (LEN=30) :: ofileb,ofilev
CHARACTER (LEN=30) :: irfile,ofiler
! Added by EG
DOUBLE PRECISION :: p1,p2,x_voro,y_voro,value_voro
INTEGER :: nn,ncell,nucleus
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: voronoi


OPEN(UNIT=10,FILE='voroslice.in',STATUS='old')
!
! Read in input file names
!
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,1)ifilev
READ(10,1)irfile
!
! Bounding box of plot
!
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,1)ofileb
!
! Now read in velocity grid parameters
!
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,*)nnt,nnp
READ(10,*)gotn,gots
READ(10,*)gopw,gope
rgst=(gotn-gots)/(nnt-1)
rgsp=(gope-gopw)/(nnp-1)
READ(10,*)ddt,ddp
READ(10,1)ofilev
!
! Now read in ray path parameters
!
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,*)prp
READ(10,1)ofiler
1   FORMAT(a30)
CLOSE(10)
!
! Read in the Voronoi tessellation file.
!
OPEN(UNIT=20,FILE=ifilev,status='old')
READ(20,*)ncell
ALLOCATE(voronoi(ncell,4), STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with ALLOCATE: PROGRAM voroslice: REAL voronoi'
ENDIF
DO i=1,ncell
        READ(20,*)nn,x_voro,y_voro,value_voro
        voronoi(i,1)=x_voro
        voronoi(i,2)=y_voro
        voronoi(i,4)=value_voro
END DO
CLOSE(20)
voronoi(:,3)=0
!
! Calculate GMT bounds file for depth slice if required.
!
lft=gopw
rgt=gope
btm=gots
top=gotn
rd1=rgsp/ddp
rd2=rgst/ddt
OPEN(UNIT=50,FILE=ofileb,STATUS='unknown')
WRITE(50,'(f16.10)')lft
WRITE(50,'(f16.10)')rgt
WRITE(50,'(f16.10)')btm
WRITE(50,'(f16.10)')top
WRITE(50,'(f16.10)')rd1
WRITE(50,'(f16.10)')rd2
CLOSE(50)
!
! Extract velocity slice 
!
nnx=(nnt-1)*ddt+1
nnz=(nnp-1)*ddp+1
OPEN(UNIT=30,FILE=ofilev,STATUS='unknown')
DO i=1,nnz
        DO j=nnx,1,-1
       
        p2=gotn-(j-1)*rd2
        p1=gopw+(i-1)*rd1

        CALL find_closest_nucleus([p1,p2,DBLE(0)],voronoi,ncell,nucleus)
        
        WRITE(30,*)voronoi(nucleus,4)
        
        ENDDO
ENDDO
CLOSE(30)
!
! Plot raypaths if required.
!
IF(prp.EQ.1)THEN
   sep='>'
   OPEN(UNIT=20,FILE=irfile,FORM='unformatted',STATUS='old')
   READ(20)nr
   OPEN(UNIT=30,FILE=ofiler,STATUS='unknown')
   DO i=1,nr
      READ(20)nre
      IF(nre.NE.0)THEN
         DO j=1,nre
            READ(20)rlat,rlon
            WRITE(30,'(2f10.4)')rlon,rlat
         ENDDO
         WRITE(30,'(a1)')sep
      ENDIF
   ENDDO
   CLOSE(20)
   CLOSE(30)
ENDIF

DEALLOCATE(voronoi, STAT=checkstat)
IF(checkstat > 0)THEN
        WRITE(6,*)'Error with DEALLOCATE: PROGRAM voroslice: voronoi'
ENDIF

STOP
END PROGRAM slice


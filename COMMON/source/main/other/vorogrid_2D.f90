! ********************************************************************* !
! Program to convert a Voronoi model to a vtx file
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
! a form suitable for input to fm2dss.                        !
!                                                             !
! Erica Galetti, March 2017                                   !
! erica.galetti@ed.ac.uk                                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM grid
USE voronoi_operations
USE string_conversion
IMPLICIT NONE
INTEGER :: i,j,ii,jj
INTEGER :: checkstat
INTEGER :: nnt,nnp
INTEGER :: nnx,nnz
DOUBLE PRECISION :: got,gop,rgst,rgsp
DOUBLE PRECISION :: lft,rgt,btm,top
DOUBLE PRECISION :: rd1,rd2
CHARACTER (LEN=30) :: ifilev
CHARACTER (LEN=30) :: ofilev
! Added by EG
DOUBLE PRECISION :: p1,p2,x_voro,y_voro,value_voro
INTEGER :: nn,ncell,nucleus
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: voronoi


OPEN(UNIT=10,FILE='vorogrid.in',STATUS='old')
!
! Read in input file names
!
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,1)ifilev
!
! Area (velocity grid parameters)
!
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,*)top,btm
READ(10,*)lft,rgt
READ(10,*)nnt,nnp
got=top
gop=lft
nnx=nnt
nnz=nnp
!
! Output file
!
READ(10,*)
READ(10,*)
READ(10,*)
READ(10,1)ofilev
1   FORMAT(a30)
CLOSE(10)
!
! Read in the Voronoi tessellation file.
!
OPEN(UNIT=20,FILE=ifilev,status='old')
READ(20,*)ncell
WRITE(*,*)'Model has ',TRIM(int2str(ncell)),' cells.'
ALLOCATE(voronoi(ncell,4), STAT=checkstat)
IF(checkstat > 0)THEN
   WRITE(6,*)'Error with ALLOCATE: PROGRAM vorogrid: REAL voronoi'
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
! Node spacing
!
rd1=(rgt-lft)/(nnp-1)
rd2=(top-btm)/(nnt-1)
!
! Write grid file
!
OPEN(UNIT=30,FILE=ofilev,STATUS='unknown')
WRITE(30,*)nnx,nnz
WRITE(30,*)got,gop
WRITE(30,*)rd2,rd1
WRITE(30,'(1X)')
DO i=0,nnz+1
        DO j=0,nnx+1 !nnx+1,0,-1
        
        ii=i
        jj=j
        IF  (i.EQ.0) ii=1
        IF  (j.EQ.0) jj=1
        IF  (i.EQ.nnz+1) ii=nnz
        IF  (j.EQ.nnx+1) jj=nnx
        
        p2=got-(jj-1)*rd2 ! theta (latitude)
        p1=gop+(ii-1)*rd1 ! phi (longitude)

        CALL find_closest_nucleus([p1,p2,DBLE(0)],voronoi,ncell,nucleus)
        
        WRITE(30,*)voronoi(nucleus,4)
        
        ENDDO
        WRITE(30,'(1X)')
ENDDO
CLOSE(30)

DEALLOCATE(voronoi, STAT=checkstat)
IF(checkstat > 0)THEN
	WRITE(6,*)'Error with DEALLOCATE: PROGRAM voroslice: voronoi'
ENDIF

STOP
END PROGRAM grid


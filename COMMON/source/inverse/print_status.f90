! ********************************************************************* !
! Subroutine to print chain status to screen in mksamples.f90
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


SUBROUTINE print_status(start_end,acceptances,iteration,number_of_chains)

USE string_conversion
USE chains, ONLY: noise_type,chain_name,samples_at_start,ncell_current,&
                  misfit_current,max_number_of_samples,number_of_datasets,&
                  noise_current,number_of_measurements,sigma_current,lambda_current,&
                  a_current,b_current,data_current,number_of_samples
IMPLICIT NONE
INTEGER, INTENT(IN)                                     :: start_end
DOUBLE PRECISION, DIMENSION(:), INTENT(IN), OPTIONAL    :: acceptances
INTEGER, INTENT(IN), OPTIONAL                           :: iteration,number_of_chains
INTEGER, DIMENSION(8)                                   :: values
CHARACTER(LEN=2)                                        :: day,month
CHARACTER(LEN=4)                                        :: year
CHARACTER(LEN=2)                                        :: hour,minute,second
INTEGER                                                 :: nn
CHARACTER(LEN=10)                                       :: tot_acc
CHARACTER(LEN=10)                                       :: b_acc,d_acc,m_acc,v_acc,n_acc
CHARACTER(LEN=10)                                       :: mm_acc,md_acc,vm_acc,vd_acc
CHARACTER(LEN=10)                                       :: o_acc,om_acc,od_acc
CHARACTER(LEN=13)                                       :: min_unc,max_unc
CHARACTER(LEN=10)                                       :: t_acc
CHARACTER(LEN=*), PARAMETER                             :: FMT1="(f9.4)"
CHARACTER(LEN=*), PARAMETER                             :: FMT2="(f12.4)"
INTEGER                                                 :: ii,k=1
DOUBLE PRECISION,DIMENSION(SIZE(data_current))                :: data_tempo,noise_tempo

IF (start_end.EQ.1.OR.start_end.EQ.2) THEN
        
        CALL date_and_time(VALUES=values)
        
        WRITE(day,'(i2.2)')values(3)
        WRITE(month,'(i2.2)')values(2)
        WRITE(year,'(i4.4)')values(1)
        WRITE(hour,'(i2.2)')values(5)
        WRITE(minute,'(i2.2)')values(6)
        WRITE(second,'(i2.2)')values(7)
        
        WRITE(*,*)
        IF (start_end.EQ.1) THEN
                WRITE(*,*)'***************************'
                WRITE(*,*)' Job started on ',day,'-',month,'-',year
                WRITE(*,*)' at ',hour,':',minute,':',second
                WRITE(*,*)'***************************'
        ELSE IF (start_end.EQ.2) THEN
                WRITE(*,*)'****************************'
                WRITE(*,*)' Job finished on ',day,'-',month,'-',year
                WRITE(*,*)' at ',hour,':',minute,':',second
                WRITE(*,*)'****************************'
        END IF
        WRITE(*,*)
        
ELSE IF (start_end.EQ.3) THEN
        
        WRITE(tot_acc,FMT1)acceptances(1)
        WRITE(b_acc,FMT1)acceptances(2)
        WRITE(d_acc,FMT1)acceptances(3)
        WRITE(m_acc,FMT1)acceptances(4)
        WRITE(v_acc,FMT1)acceptances(5)
        WRITE(n_acc,FMT1)acceptances(6)
        WRITE(mm_acc,FMT1)acceptances(7)
        WRITE(md_acc,FMT1)acceptances(8)
        WRITE(vm_acc,FMT1)acceptances(9)
        WRITE(vd_acc,FMT1)acceptances(10)
        !----Lea
        WRITE(o_acc,FMT1)acceptances(14)
        WRITE(om_acc,FMT1)acceptances(15)
        WRITE(od_acc,FMT1)acceptances(16)
        !-------
        WRITE(min_unc,FMT2)MINVAL(noise_current/ABS(data_current)*100)
        WRITE(max_unc,FMT2)MAXVAL(noise_current/ABS(data_current)*100)


        !--------------------------Modifs-------------------------------------------------
        !DO ii=1,number_of_measurements
        !        IF (data_current(ii).GT.0) THEN
        !                data_tempo(k) = data_current(ii)
        !                noise_tempo(k) = noise_current(ii)
        !                k=k+1
        !        END IF
        !END DO
        !DO ii=1,number_of_measurements
        !        IF (data_current(ii).LT.0) THEN
        !                data_current(ii)=0
        !                noise_current(ii)=0
        !        END IF
        !END DO

        !WRITE(min_unc,FMT2)MINVAL(noise_tempo/ABS(data_tempo)*100)
        !WRITE(max_unc,FMT2)MAXVAL(noise_tempo/ABS(data_tempo)*100)

        !WRITE(*,*)MINVAL(noise_current/ABS(data_current)*100),MINVAL(ABS(noise_current)),MAXVAL(data_current)
        !WRITE(*,*)MAXVAL(noise_current/ABS(data_current)*100),MAXVAL(noise_current),MINVAL(ABS(data_current))
        !WRITE(*,*)MINVAL(noise_tempo/ABS(data_tempo)*100),MINVAL(ABS(noise_tempo)),MAXVAL(data_tempo)
        !WRITE(*,*)MAXVAL(noise_tempo/ABS(data_tempo)*100),MAXVAL(noise_tempo),MINVAL(ABS(data_tempo))
        !----------------------------------------------------------------------

        WRITE(*,*)'--------------------------------------'
        WRITE(*,*)'Chain number: ',TRIM(chain_name),'/',TRIM(int2str(number_of_chains))
        WRITE(*,*)'Sample: ',TRIM(int2str(samples_at_start+iteration)),'/',&
                             TRIM(int2str(max_number_of_samples))
        WRITE(*,*)'Iteration in this run: ',TRIM(int2str(iteration)),'/',&
                                            TRIM(int2str(number_of_samples))
        WRITE(*,*)
        WRITE(*,*)'Number of cells: ',TRIM(int2str(ncell_current))
        WRITE(*,*)
        IF (noise_type.EQ.0) THEN
                DO nn=1,number_of_datasets
                        WRITE(*,*)'Dataset ',TRIM(int2str(nn)),', "sigma": ',&
                                        TRIM(dble2str(sigma_current(nn)))
                END DO
        ELSE IF (noise_type.EQ.3.OR.noise_type.EQ.6) THEN
                DO nn=1,number_of_datasets
                        WRITE(*,*)'Dataset ',TRIM(int2str(nn)),', "lambda": ',&
                                        TRIM(dble2str(lambda_current(nn)))
                END DO
        ELSE IF (noise_type.EQ.1.OR.noise_type.EQ.2.OR.noise_type.EQ.5) THEN
                DO nn=1,number_of_datasets
                        WRITE(*,*)'Dataset ',TRIM(int2str(nn)),', "a": ',&
                                        TRIM(dble2str(a_current(nn)))
                        WRITE(*,*)'Dataset ',TRIM(int2str(nn)),', "b": ',&
                                        TRIM(dble2str(b_current(nn)))
                END DO
        END IF
        WRITE(*,*)'Minimum uncertainty: ',min_unc,'%'
        WRITE(*,*)'Maximum uncertainty: ',max_unc,'%'
        WRITE(*,*)
        WRITE(*,*)'Misfit : ',TRIM(dble2str(misfit_current))
        WRITE(*,*)
        WRITE(*,*)'BIRTH acceptance: ',b_acc,'%'
        WRITE(*,*)'DEATH acceptance: ',d_acc,'%'
        WRITE(*,*)'MOVE  acceptance: ',m_acc,'%'
        WRITE(*,*)'      main: ',mm_acc,'%'
        WRITE(*,*)'      delayed: ',md_acc,'%'
        WRITE(*,*)'VALUE acceptance: ',v_acc,'%'
        WRITE(*,*)'      main: ',vm_acc,'%'
        WRITE(*,*)'      delayed: ',vd_acc,'%'
        WRITE(*,*)'NOISE acceptance: ',n_acc,'%'
        ! --- Lea
        WRITE(*,*)'OPTION acceptance: ',o_acc,'%'
        WRITE(*,*)'      main: ',om_acc,'%'
        WRITE(*,*)'      delayed: ',od_acc,'%'
        ! ---
        WRITE(*,*)
        WRITE(*,*)'TOTAL acceptance: ',tot_acc,'%'
        WRITE(*,*)'--------------------------------------'
        
ELSE IF (start_end.EQ.4) THEN
        
        WRITE(t_acc,FMT1)acceptances(13)
        WRITE(*,*)'--------------------------------------'
        WRITE(*,*)'PARALLEL TEMPERING acceptance: ',t_acc,'%'
        WRITE(*,*)'--------------------------------------'
        
END IF

END SUBROUTINE
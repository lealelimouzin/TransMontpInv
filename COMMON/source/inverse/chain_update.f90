! ********************************************************************* !
! Module to update Markov chain status after model acceptance/rejection
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


MODULE chain_update

USE chains
USE acquisition, ONLY: observed_data


CONTAINS


SUBROUTINE update_models(attempt)

IMPLICIT NONE
INTEGER, INTENT(IN)     :: attempt

IF (attempt.EQ.1) THEN ! we tried once
        
        likelihood_current=likelihood_proposed
        misfit_current=misfit_proposed
        noise_current=noise_proposed
        
        IF (.NOT.noise) THEN ! the Voronoi model was changed
                
                ncell_current=ncell_proposed
                voronoi_current=voronoi_proposed
                data_current=data_proposed
                IF (noise_type.EQ.1.OR.noise_type.EQ.2.OR.noise_type.EQ.5) noise_x_current=noise_x_proposed
                
        ELSE IF (noise) THEN ! a noise parameter was changed
                
                IF (noise_type.EQ.0) THEN
                        sigma_current=sigma_proposed
                ELSE IF (noise_type.EQ.1.OR.noise_type.EQ.2.OR.noise_type.EQ.5) THEN
                        a_current=a_proposed
                        b_current=b_proposed
                ELSE IF (noise_type.EQ.3.OR.noise_type.EQ.6) THEN
                        lambda_current=lambda_proposed
                END IF
                
        END IF
        
ELSE IF (attempt.EQ.2) THEN ! we used delayed rejection
        
        likelihood_current=likelihood_delayed
        misfit_current=misfit_delayed
        noise_current=noise_delayed
        ncell_current=ncell_proposed
        voronoi_current=voronoi_delayed
        data_current=data_delayed
        IF (noise_type.EQ.1.OR.noise_type.EQ.2.OR.noise_type.EQ.5) noise_x_current=noise_x_delayed
        
END IF

END SUBROUTINE


SUBROUTINE update_acceptance(accept,accept_main,accept_delayed)

IMPLICIT NONE
LOGICAL, INTENT(IN)     :: accept,accept_main,accept_delayed

IF (accept) THEN
        
        accepted_in_run=accepted_in_run+1
        
        acceptance_total(1)=acceptance_total(1)+1
        IF (birth) acceptance_birth(1)=acceptance_birth(1)+1
        IF (death) acceptance_death(1)=acceptance_death(1)+1
        IF (move) acceptance_move(1)=acceptance_move(1)+1
        IF (value) acceptance_value(1)=acceptance_value(1)+1
        IF (noise) acceptance_noise(1)=acceptance_noise(1)+1
        IF (noise_type.EQ.1.OR.noise_type.EQ.2.OR.noise_type.EQ.5) THEN
                IF (noise.AND.change_a) acceptance_a(1)=acceptance_a(1)+1
                IF (noise.AND.change_b) acceptance_b(1)=acceptance_b(1)+1
        END IF
        IF (option) acceptance_opt(1)=acceptance_opt(1)+1
ELSE
        
        acceptance_total(2)=acceptance_total(2)+1
        IF (birth) acceptance_birth(2)=acceptance_birth(2)+1
        IF (death) acceptance_death(2)=acceptance_death(2)+1
        IF (move) acceptance_move(2)=acceptance_move(2)+1
        IF (value) acceptance_value(2)=acceptance_value(2)+1
        IF (noise) acceptance_noise(2)=acceptance_noise(2)+1
        IF (noise_type.EQ.1.OR.noise_type.EQ.2.OR.noise_type.EQ.5) THEN
                IF (noise.AND.change_a) acceptance_a(2)=acceptance_a(2)+1
                IF (noise.AND.change_b) acceptance_b(2)=acceptance_b(2)+1
        END IF
        IF (option) acceptance_opt(2)=acceptance_opt(2)+1
END IF

IF (move.OR.value.OR.option) THEN
        IF (accept_main) THEN
                IF (move) acceptance_move_main(1)=acceptance_move_main(1)+1
                IF (value) acceptance_value_main(1)=acceptance_value_main(1)+1
                IF (option) acceptance_opt_main(1)=acceptance_opt_main(1)+1
        ELSE
                IF (move) acceptance_move_main(2)=acceptance_move_main(2)+1
                IF (value) acceptance_value_main(2)=acceptance_value_main(2)+1
                IF (option) acceptance_opt_main(2)=acceptance_opt_main(2)+1
        END IF
        
        IF (delayed_position.OR.delayed_value.OR.delayed_opt) THEN
                IF (accept_delayed) THEN
                        IF (move) acceptance_move_delayed(1)=acceptance_move_delayed(1)+1
                        IF (value) acceptance_value_delayed(1)=acceptance_value_delayed(1)+1
                        IF (option) acceptance_opt_delayed(1)=acceptance_opt_delayed(1)+1
                ELSE
                        IF (move) acceptance_move_delayed(2)=acceptance_move_delayed(2)+1
                        IF (value) acceptance_value_delayed(2)=acceptance_value_delayed(2)+1
                        IF (option) acceptance_opt_delayed(2)=acceptance_opt_delayed(2)+1
                END IF
        END IF
END IF

IF (SUM(acceptance_total(1:2)).GT.0) acceptance_total(3)=100*acceptance_total(1)/SUM(acceptance_total(1:2))
IF (birth.AND.SUM(acceptance_birth(1:2)).GT.0) acceptance_birth(3)=100*acceptance_birth(1)/SUM(acceptance_birth(1:2))
IF (death.AND.SUM(acceptance_death(1:2)).GT.0) acceptance_death(3)=100*acceptance_death(1)/SUM(acceptance_death(1:2))
IF (move.AND.SUM(acceptance_move(1:2)).GT.0) acceptance_move(3)=100*acceptance_move(1)/SUM(acceptance_move(1:2))
IF (move.AND.SUM(acceptance_move_main(1:2)).GT.0) acceptance_move_main(3)=100*acceptance_move_main(1)&
                                                                                /SUM(acceptance_move_main(1:2))
IF (move.AND.SUM(acceptance_move_delayed(1:2)).GT.0) acceptance_move_delayed(3)=100*acceptance_move_delayed(1)&
                                                                                /SUM(acceptance_move_delayed(1:2))
IF (value.AND.SUM(acceptance_value(1:2)).GT.0) acceptance_value(3)=100*acceptance_value(1)/SUM(acceptance_value(1:2))
IF (value.AND.SUM(acceptance_value_main(1:2)).GT.0) acceptance_value_main(3)=100*acceptance_value_main(1)&
                                                                                /SUM(acceptance_value_main(1:2))
IF (value.AND.SUM(acceptance_value_delayed(1:2)).GT.0) acceptance_value_delayed(3)=100*acceptance_value_delayed(1)&
                                                                                /SUM(acceptance_value_delayed(1:2))
IF (noise.AND.SUM(acceptance_noise(1:2)).GT.0) acceptance_noise(3)=100*acceptance_noise(1)/SUM(acceptance_noise(1:2))
IF (noise_type.EQ.1.OR.noise_type.EQ.2.OR.noise_type.EQ.5) THEN
        IF (noise.AND.SUM(acceptance_a(1:2)).GT.0) acceptance_a(3)=100*acceptance_a(1)/SUM(acceptance_a(1:2))
        IF (noise.AND.SUM(acceptance_b(1:2)).GT.0) acceptance_b(3)=100*acceptance_b(1)/SUM(acceptance_b(1:2))
END IF
! ----- Option -----
IF (option.AND.SUM(acceptance_opt(1:2)).GT.0) acceptance_opt(3)=100*acceptance_opt(1)/SUM(acceptance_opt(1:2))
IF (option.AND.SUM(acceptance_opt_main(1:2)).GT.0) acceptance_opt_main(3)=100*acceptance_opt_main(1)&
                                                                                /SUM(acceptance_opt_main(1:2))
IF (option.AND.SUM(acceptance_opt_delayed(1:2)).GT.0) acceptance_opt_delayed(3)=100*acceptance_opt_delayed(1)&
                                                                                /SUM(acceptance_opt_delayed(1:2))
acceptance_all(1)=acceptance_total(3)
acceptance_all(2)=acceptance_birth(3)
acceptance_all(3)=acceptance_death(3)
acceptance_all(4)=acceptance_move(3)
acceptance_all(5)=acceptance_value(3)
acceptance_all(6)=acceptance_noise(3)
acceptance_all(7)=acceptance_move_main(3)
acceptance_all(8)=acceptance_move_delayed(3)
acceptance_all(9)=acceptance_value_main(3)
acceptance_all(10)=acceptance_value_delayed(3)
acceptance_all(11)=acceptance_a(3)
acceptance_all(12)=acceptance_b(3)
IF (opt_parameters.GT.0) THEN
                acceptance_all(14)=acceptance_opt(3)
                acceptance_all(15)=acceptance_opt_main(3)
                acceptance_all(16)=acceptance_opt_delayed(3)
                ! acceptance_all(13) is for tempering
END IF

END SUBROUTINE


SUBROUTINE update_chain_arrays(accept,iteration)

IMPLICIT NONE
LOGICAL, INTENT(IN)     :: accept
INTEGER, INTENT(IN)     :: iteration
INTEGER                 :: sample_number,nn

IF (run_number.EQ.1) THEN
        sample_number=iteration+1
ELSE IF (run_number.GT.1) THEN
        sample_number=iteration
END IF

chain_steps(sample_number,1)=step
chain_ncell(sample_number)=ncell_current
chain_misfits(sample_number)=misfit_current
chain_acceptance(sample_number,1)=acceptance_total(3)
chain_acceptance(sample_number,2)=acceptance_birth(3)
chain_acceptance(sample_number,3)=acceptance_death(3)
chain_acceptance(sample_number,4)=acceptance_move(3)
chain_acceptance(sample_number,5)=acceptance_value(3)
chain_acceptance(sample_number,6)=acceptance_noise(3)
IF (opt_parameters.GT.0)  chain_acceptance(sample_number,8)=acceptance_opt(3)

IF (birth) THEN
        IF (accept) THEN
                chain_steps(sample_number,2)=1
                chain_indices(sample_number)=0
                chain_values(sample_number,1)=voronoi_current(ncell_current,1)
                chain_values(sample_number,2)=voronoi_current(ncell_current,2)
                chain_values(sample_number,3)=voronoi_current(ncell_current,3)
                chain_values(sample_number,4)=voronoi_current(ncell_current,4)
                IF (opt_parameters.GT.0)  chain_values(sample_number,5)=voronoi_current(ncell_current,5)
        END IF
ELSE IF (death) THEN
        IF (accept) THEN
                chain_steps(sample_number,2)=1
                chain_indices(sample_number)=cell_number
                chain_values(sample_number,1)=0
                chain_values(sample_number,2)=0
                chain_values(sample_number,3)=0
                chain_values(sample_number,4)=0
                IF (opt_parameters.GT.0)  chain_values(sample_number,5)=0
        END IF
ELSE IF (move) THEN
        IF (accept) THEN
                chain_steps(sample_number,2)=1
                chain_indices(sample_number)=cell_number
                chain_values(sample_number,1)=voronoi_current(cell_number,1)
                chain_values(sample_number,2)=voronoi_current(cell_number,2)
                chain_values(sample_number,3)=voronoi_current(cell_number,3)
                chain_values(sample_number,4)=voronoi_current(cell_number,4)
                IF (opt_parameters.GT.0)  chain_values(sample_number,5)=0
        END IF
ELSE IF (value) THEN
        IF (accept) THEN
                chain_steps(sample_number,2)=1
                chain_indices(sample_number)=cell_number
                chain_values(sample_number,1)=voronoi_current(cell_number,1) 
                chain_values(sample_number,2)=voronoi_current(cell_number,2) 
                chain_values(sample_number,3)=voronoi_current(cell_number,3) 
                chain_values(sample_number,4)=voronoi_current(cell_number,4)
                IF (opt_parameters.GT.0)  chain_values(sample_number,5)=0
        END IF
ELSE IF (option) THEN
        IF (accept) THEN
                chain_steps(sample_number,2)=1
                chain_indices(sample_number)=cell_number
                chain_values(sample_number,1)=voronoi_current(cell_number,1)
                chain_values(sample_number,2)=voronoi_current(cell_number,2)
                chain_values(sample_number,3)=voronoi_current(cell_number,3)
                chain_values(sample_number,4)=voronoi_current(cell_number,4)
                chain_values(sample_number,5)=voronoi_current(ncell_current,5)
        END IF
ELSE IF (noise) THEN
        IF (accept) THEN
                chain_steps(sample_number,2)=1
                chain_indices(sample_number)=cell_number
                chain_values(sample_number,1)=0
                chain_values(sample_number,2)=0
                chain_values(sample_number,3)=0
                chain_values(sample_number,4)=voronoi_current(cell_number,4)
                IF (opt_parameters.GT.0)  chain_values(sample_number,5)=0
        END IF
END IF
DO nn=1,number_of_datasets
        IF (noise_type.EQ.0) THEN
                chain_noise(sample_number,nn)=sigma_current(nn)
        ELSE IF (noise_type.EQ.1.OR.noise_type.EQ.2.OR.noise_type.EQ.5) THEN
                chain_noise(sample_number,nn)=a_current(nn)
                chain_noise(sample_number,nn+number_of_datasets)=b_current(nn)
        ELSE IF (noise_type.EQ.3.OR.noise_type.EQ.6) THEN
                chain_noise(sample_number,nn)=lambda_current(nn)
        END IF
END DO

IF (accept.AND..NOT.noise) residuals(:,sample_number)=observed_data-data_current

END SUBROUTINE

END MODULE
! --------------------------------------------------------------------
!
! Copyright (C) 2015 Rocco Meli
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
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
! ---------------------------------------------------------------------

MODULE CORE

    USE KINETIC, only: T_kinetic
    USE NUCLEAR, only: V_nuclear
    USE UTILS, only: print_real_matrix

    IMPLICIT NONE

    CONTAINS

        SUBROUTINE H_core(Kf,c,Nn,basis_D,basis_A,basis_L,basis_R,Rn,Zn,H,verbose)

            ! TODO Allow flexibility for basis sets other than STO-3G
            INTEGER, intent(in) :: c ! Number of contractions per basis function

            LOGICAL, intent(in):: verbose

            ! INPUT
            INTEGER, intent(in) :: Kf                       ! Number of basis functions
            INTEGER, intent(in) :: Nn                       ! Total number of nuclei
            REAL*8, dimension(Kf,3), intent(in) :: basis_R  ! Basis set niclear positions
            INTEGER, dimension(Kf,3), intent(in) :: basis_L ! Basis set angular momenta
            REAL*8, dimension(Kf,3), intent(in) :: basis_D  ! Basis set contraction coefficients
            REAL*8, dimension(Kf,3), intent(in) :: basis_A  ! Basis set exponential contraction coefficients
            REAL*8, dimension(Nn,3), intent(in) :: Rn       ! Nuclear positions
            INTEGER, dimension(Nn), intent(in) :: Zn        ! Nuclear charges (> 0)

            ! INTERMEDIATE VARIABLE
            REAL*8, dimension(Kf,Kf) :: V   ! Matrix storing total potential components
            REAL*8, dimension(Kf,Kf) :: M   ! Matrix storing one-atom potential conponents
            INTEGER :: i                    ! Loop index

            ! OUTPUT
            REAL*8, dimension(Kf,Kf), intent(out) :: H ! Core Hamiltonian

            H(:,:) = 0.0D0

            V(:,:) = 0.0D0

            CALL T_kinetic(Kf,c,basis_D,basis_A,basis_L,basis_R,H)

            IF (verbose) THEN
                WRITE(*,*) "Kinetic energy matrix T:"
                CALL print_real_matrix(Kf,Kf,H)
            END IF

            DO i = 1, Nn
                M(:,:) = 0.0D0

                CALL V_nuclear(Kf,c,basis_D,basis_A,basis_L,basis_R,M,Rn(i,1:3),Zn(i))

                IF (verbose) THEN
                    WRITE(*,*) "Nuclear attraction matrix V:", i
                    CALL print_real_matrix(Kf,Kf,M)
                END IF

                V = V + M
            END DO

            IF (verbose) THEN
                WRITE(*,*) "Total nuclear attraction matrix V:", i
                CALL print_real_matrix(Kf,Kf,V)
            END IF

            H = H + V

        END SUBROUTINE H_core

        SUBROUTINE density()

        END SUBROUTINE density

END MODULE CORE

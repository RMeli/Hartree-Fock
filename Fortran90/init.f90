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

MODULE INIT
    ! ---------------------------------
    ! INITIAL GUESS FOR THE FOCK MATRIX
    ! ---------------------------------

    IMPLICIT NONE

    CONTAINS

    ! ----------
    ! CORE GUESS
    ! ----------
    SUBROUTINE core_guess(Kf,Hc,F)
        ! ----------------------------
        ! Core Hamiltonian
        ! ----------------------------
        !
        ! Source:
        !   A. Szabo and N. S. Ostlund
        !   Modern Quantum Chemistry
        !   Dover
        !   1996
        !
        !-----------------------------

        IMPLICIT NONE

        ! INPUT
        INTEGER, intent(in) :: Kf                   ! Basis set size
        REAL*8, dimension(Kf,Kf), intent(in) :: Hc  ! Core Hamiltonian

        ! OUTPUT
        REAL*8, dimension(Kf,Kf), intent(out) :: F   ! Initial Guess for the Fock operator

        ! Initial guess: Core Hamiltonian
        F = Hc

    END SUBROUTINE core_guess



    ! ------------
    ! HÜCKEL GUESS
    ! ------------
    SUBROUTINE huckel_guess(Kf,Hc,S,F,cst)
        ! -----------------------------------------------------
        ! Extended Hückel Theory
        ! -----------------------------------------------------
        !
        ! Source:
        !   Extended Hückel and Slater’s rule initial guess for
        !       real space grid-based density functional theory
        !   M. Lee, K. Leiter, C. Eisner, J. Crone and J. Knap
        !   Computational and Theoretical Chemistry
        !   2015
        !
        !------------------------------------------------------

        IMPLICIT NONE

        ! INPUT
        INTEGER, intent(in) :: Kf                   ! Basis set size
        REAL*8, dimension(Kf,Kf),intent(in) :: S    ! Overlap matrix
        REAL*8, dimension(Kf,Kf), intent(in) :: Hc  ! Core Hamiltonian
        Real*8, intent(in) :: cst                   ! Multiplicative constant in Hückel model

        ! INTERMEDIATE VARIABLES
        INTEGER :: i, j                             ! Loop indices

        ! OUTPUT
        REAL*8, dimension(Kf,Kf), intent(out) :: F   ! Initial Guess for the Fock operator

        F(:,:) = 0.0D0

        DO i = 1,Kf
            DO j = 1,Kf
                F(i,j) = cst * S(i,j) * 0.5D0 * (Hc(i,i) + Hc(j,j))
            END DO
        END DO

    END SUBROUTINE huckel_guess


END MODULE INIT

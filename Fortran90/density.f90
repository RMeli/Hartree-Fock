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

MODULE DENSITY

    IMPLICIT NONE

    CONTAINS

        ! --------------
        ! DENSITY MATRIX
        ! --------------
        SUBROUTINE P_density(Kf,Ne,C,P)
            ! ----------------------------------------------------------
            ! Compute the density matrix from the optimized coefficients
            ! ----------------------------------------------------------
            !
            ! Source:
            !   Szabo and Ostlund
            !   Modern Quantum Chemistry
            !   Doever
            !   1989
            !
            ! -----------------------------------------------------------

            IMPLICIT NONE

            ! INPUT
            INTEGER, intent(in) :: Kf                           ! Number of basis functions
            INTEGER, intent(in) :: Ne                           ! Number of electrons
            REAL*8, dimension(Kf,Kf), intent(in) :: C           ! Coefficients matrix

            ! INTERMEDIATE VARIABLES
            INTEGER :: i, j, k                                  ! Loop coefficients

            ! OUTPUT
            REAL*8, dimension(Kf,Kf), intent(out) :: P          ! Density matrix

            P(:,:) = 0.0D0

            DO i = 1, Kf
                DO j = 1, Kf
                    DO k = 1, Ne / 2 ! TODO Only in RHF
                        P(i,j) = P(i,j) + 2.0D0 * C(i,k) * C(j,k)
                    END DO
                END DO
            END DO

        END SUBROUTINE P_density


        ! -------------------------
        ! DENSITY MATRIX DIFFERENCE
        ! -------------------------
        FUNCTION delta_P(Kf,Pold,Pnew) result(delta)
            ! ----------------------------------------------------------
            ! Compute density matrix difference as convergence criterion
            ! ----------------------------------------------------------
            !
            ! Source:
            !   A. Szabo and N. S. Ostlund
            !   Modern Quantum Chemistry
            !   Dover
            !   1996
            !
            ! -----------------------------------------------------------

            IMPLICIT NONE

            ! INPUT
            INTEGER, intent(in) :: Kf                       ! Basis set size
            REAL*8, dimension(Kf,Kf), intent(in) :: Pold    ! Old density matrix
            REAL*8, dimension(Kf,Kf), intent(in) :: Pnew    ! New density matix

            ! INTERMEDIATE VARIABLES
            INTEGER :: i, j

            ! OUTPUT
            REAL*8 :: delta                                 ! Sum of matrix elements square differences

            delta = 0


            DO i = 1, Kf
                DO j = 1, Kf
                    delta = delta + (Pold(i,j)-Pnew(i,j))**2
                END DO
            END DO

            delta = DSQRT(delta / Kf**2)

        END FUNCTION delta_P

END MODULE DENSITY

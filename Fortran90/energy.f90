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

MODULE ENERGY

    IMPLICIT NONE

    CONTAINS

        FUNCTION E_electronic(Kf,P,F,H) result(E_el)
            ! -
            !
            ! -
            !
            ! Source:
            !   A. Szabo and N. Ostlund
            !
            ! -

            IMPLICIT NONE

            ! INPUT
            INTEGER, intent(in) :: Kf
            REAL*8, dimension(Kf,Kf), intent(in) :: P   ! Density matrix
            REAL*8, dimension(Kf,Kf), intent(in) :: F   ! Fock operator
            REAL*8, dimension(Kf,Kf), intent(in) :: H   ! Core Hamiltonian

            ! INTERMEDIATE VARIABLES
            INTEGER :: i, j                             ! Loop indices

            ! OUTPUT
            REAL*8 :: E_el                              ! Electronic energy

            E_el = 0.0D0

            DO i = 1, Kf
                DO j = 1, Kf
                    E_el = E_el + 0.5D0 * P(j,i) * (H(i,j) + F(i,j))
                END DO
            END DO

        END FUNCTION E_electronic

        FUNCTION E_electronic_spin(Kf,Pa,Fa,Pb,Fb,H) result(E_el)

            IMPLICIT NONE

            ! INPUT
            INTEGER, intent(in) :: Kf
            REAL*8, dimension(Kf,Kf), intent(in) :: Pa      ! Density matrix (alpha spin)
            REAL*8, dimension(Kf,Kf), intent(in) :: Pb      ! Density matrix (beta spin)
            REAL*8, dimension(Kf,Kf), intent(in) :: Fa      ! Fock operator (alpha spin)
            REAL*8, dimension(Kf,Kf), intent(in) :: Fb      ! Fock operator (beta spin)
            REAL*8, dimension(Kf,Kf), intent(in) :: H       ! Core Hamiltonian

            ! INTERMEDIATE VARIABLES
            INTEGER :: i, j                                 ! Loop indices
            REAL*8, dimension(Kf,Kf) :: Ptot                ! Total density matrix

            ! OUTPUT
            REAL*8 :: E_el                              ! Electronic energy

            Ptot = Pa + Pb

            E_el = 0.0D0

            DO i = 1, Kf
                DO j = 1, Kf
                    E_el = E_el + Ptot(j,i) * H(i,j) + Pa(j,i) * Fa(i,j) + Pb(j,i) * Fb(i,j)
                END DO
            END DO

            E_el = 0.5D0 * E_el

        END FUNCTION E_electronic_spin


        FUNCTION E_nuclear(Nn,Rn,Zn) result(E_nu)

            IMPLICIT NONE

            ! INPUT
            INTEGER, intent(in) :: Nn                       ! Total number of nuclei
            INTEGER, dimension(Nn), intent(in) :: Zn        ! Nuclear charges
            REAL*8, dimension(Nn,3), intent(in) :: Rn       ! Nuclear positions

            ! INTERMEDIATE VARIABLES
            INTEGER :: i, j                 ! Loop variables

            ! OUTPUT
            REAl*8 :: E_nu                  ! Nuclear energy

            E_nu = 0.0D0

            DO i = 1, Nn
                DO j = i + 1, Nn
                    E_nu = E_nu + Zn(i) * Zn(j) / NORM2(Rn(i,1:3) - Rn(j,1:3))
                END DO
            END DO

        END FUNCTION E_nuclear

        FUNCTION E_tot(Kf,Nn,Rn,Zn,P,F,H)
            IMPLICIT NONE

            ! INPUT
            INTEGER, intent(in) :: Kf                   ! Basis set size
            INTEGER, intent(in) :: Nn                   ! Total number of nuclei
            INTEGER, dimension(Nn), intent(in) :: Zn    ! Nuclear charges
            REAL*8, dimension(Nn,3), intent(in) :: Rn   ! Nuclear positions
            REAL*8, dimension(Kf,Kf), intent(in) :: P   ! Density matrix
            REAL*8, dimension(Kf,Kf), intent(in) :: F   ! Fock operator
            REAL*8, dimension(Kf,Kf), intent(in) :: H   ! Core Hamiltonian

            ! OUTPUT
            REAL*8 :: E_tot

            E_tot = E_nuclear(Nn,Rn,Zn) + E_electronic(Kf,P,F,H)

        END FUNCTION E_tot

        FUNCTION E_tot_spin(Kf,Nn,Rn,Zn,Pa,Pb,Fa,Fb,H) result(E_tot)
            IMPLICIT NONE

            ! INPUT
            INTEGER, intent(in) :: Kf                   ! Basis set size
            INTEGER, intent(in) :: Nn                   ! Total number of nuclei
            INTEGER, dimension(Nn), intent(in) :: Zn    ! Nuclear charges
            REAL*8, dimension(Nn,3), intent(in) :: Rn   ! Nuclear positions
            REAL*8, dimension(Kf,Kf), intent(in) :: Pa   ! Density matrix (alpha spin)
            REAL*8, dimension(Kf,Kf), intent(in) :: Fa   ! Fock operator (alpha spin)
            REAL*8, dimension(Kf,Kf), intent(in) :: Pb   ! Density matrix (beta spin)
            REAL*8, dimension(Kf,Kf), intent(in) :: Fb   ! Fock operator (beta spin)
            REAL*8, dimension(Kf,Kf), intent(in) :: H   ! Core Hamiltonian

            ! OUTPUT
            REAL*8 :: E_tot

            E_tot = E_nuclear(Nn,Rn,Zn) + E_electronic_spin(Kf,Pa,Fa,Pb,Fb,H)

        END FUNCTION E_tot_spin



END MODULE ENERGY

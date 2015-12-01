MODULE ENERGY

    IMPLICIT NONE

    CONTAINS

        FUNCTION E_electronic(Kf,P,F,H) result(E_el)

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
                    E_el = E_el + 0.5D0 * P(i,j) * (H(i,j) + F(i,j))
                END DO
            END DO

        END FUNCTION E_electronic


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



END MODULE ENERGY

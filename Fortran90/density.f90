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

            delta = DSQRT(delta)

        END FUNCTION delta_P

END MODULE DENSITY

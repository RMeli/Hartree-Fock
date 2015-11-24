MODULE UTILS

    IMPLICIT NONE

    CONTAINS

    ! ------------
    ! PRINT MATRIX
    ! ------------
    SUBROUTINE print_real_matrix(r,c,M)

        IMPLICIT NONE

        ! INPUT
        INTEGER, intent(in) :: r, c
        REAL*8, dimension(r,c), intent(in) :: M

        ! VARIABLES
        INTEGER :: i

        DO i = 1, r
            WRITE(*,'(20G16.6)') M(i,1:c)
        END DO

    END SUBROUTINE print_real_matrix

    ! --------------------------------------
    ! PRINT ELECTRON-ELECTRON INTEGRALS LIST
    ! --------------------------------------
    SUBROUTINE print_ee_list(Kf,ee)

        IMPLICIT NONE

        ! INPUT
        INTEGER, intent(in) :: Kf               ! Basis set size
        REAL*8, dimension(Kf,Kf,Kf,Kf) :: ee    ! List of electron-electron integrals

        ! INTERMEDIATE VARIABLES
        INTEGER :: i, j, k, l   ! Loop indices

        DO i = 1,Kf
            DO j = 1,Kf
                DO k = 1,Kf
                    DO l = 1,Kf

                        WRITE(*,*) "(", i, j, k, l, ") ", ee(i,j,k,l)

                    END DO ! l
                END DO ! k
            END DO ! j
        END DO ! i

    END SUBROUTINE print_ee_list

    ! ------------------------
    ! SOLVE EIGENVALUE PROBLEM
    ! ------------------------
    SUBROUTINE EIGS(d,M,V,lambda)

        IMPLICIT NONE

        ! INPUT
        INTEGER, intent(in) :: d                        ! Dimension of M (d x d)
        REAL*8, dimension(d,d),intent(in) :: M          ! Matrix M

        ! INTERMEDIATE VARIABLES
        INTEGER, PARAMETER :: LWMAX = 1000             ! Maximal workspace size
        INTEGER :: LWORK                           ! Query optimal workspace size ()
        INTEGER :: INFO                                 ! Information flag for DSYEV
        REAL*8, dimension(LWMAX) :: WORK

        ! OUTPUT
        REAL*8, dimension(d,d), intent(out) :: V        ! Eigenvectors
        REAL*8, dimension(d), intent(out) :: lambda     ! Eigenvalues

        V = M

        LWORK = -1
        CALL  DSYEV('V','U', d, V, d, lambda, WORK, LWORK, INFO )       ! Query optimal workspace size
        LWORK = MIN(LWMAX, INT(WORK(1))) * 2 * d                        ! Optimal workspace size
                                                                        ! See (http://www.netlib.org/lapack/lug/node120.html)

        CALL  DSYEV('V','U', d, V, d, lambda, WORK, LWORK, INFO )       ! Solve the eigenvalue problem

        IF (INFO .NE. 0) THEN
            WRITE(*,*) "ERROR: IMPOSSIBLE TO SOLVE THE EIGENVALUE PROBLEM!"
        END IF

    END SUBROUTINE EIGS


END MODULE UTILS

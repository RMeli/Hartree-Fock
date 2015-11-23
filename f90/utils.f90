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

END MODULE UTILS
